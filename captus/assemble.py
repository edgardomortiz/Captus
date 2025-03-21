#!/usr/bin/env python3
"""
Copyright 2020-2025 Edgardo M. Ortiz (e.ortiz.v@gmail.com)
https://github.com/edgardomortiz/Captus

This file is part of Captus. Captus is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Captus is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Captus. If
not, see <http://www.gnu.org/licenses/>.
"""

import math
import platform
import shutil
import statistics
import subprocess
import time
from pathlib import Path

from . import log, settings
from .bioformats import dict_to_fasta, fasta_headers_to_spades, fasta_to_dict, get_read_stats
from .misc import (
    bbtools_path_version,
    bold,
    dim,
    elapsed_time,
    find_and_match_fastqs,
    format_dep_msg,
    make_output_dir,
    make_tmp_dir_within,
    megahit_path_version,
    megahit_tk_path_version,
    python_library_check,
    quit_with_error,
    red,
    salmon_path_version,
    set_ram,
    set_threads,
    successful_exit,
    tqdm_parallel_async_run,
    tqdm_serial_run,
)
from .version import __version__


def assemble(full_command, args):
    captus_start = time.time()
    out_dir, out_dir_msg = make_output_dir(args.out)
    log.logger = log.Log(Path(args.out, "captus-assemble.log"), stdout_verbosity_level=1)

    mar = 21  # Margin for aligning parameters and values

    ################################################################################################
    ############################################################################### STARTING SECTION
    log.log_section_header("Starting Captus-assembly: ASSEMBLE", single_newline=False)
    log.log_explanation(
        "Welcome to the de novo assembly step of Captus-assembly. In this step, Captus will use"
        " MEGAHIT to assemble your input reads. It is also possible to subsample a number of reads"
        " using reformat.sh from BBTools prior to assembly, this is useful while performing tests or"
        " when including samples with considerably higher sequencing depth in a dataset.",
        extra_empty_lines_after=0,
    )
    if len(args.reads) == 1 and Path(args.reads[0]).is_dir():
        intro_msg = (
            "Since you provided a directory name, Captus will look in that location for all the"
            " FASTQ files that contain the string '_R1' in their names and match them with their"
            " respective '_R2' pairs. If the '_R2' can not be found, the sample is treated as"
            " single-end. Sample names are derived from the text found before the '_R1' string."
        )
    else:
        intro_msg = (
            "Since you provided a file name or a list of file names, Captus will only process the"
            " FASTQ files that contain the string '_R1' in their names and match them with their"
            " respective '_R2' pairs. If the '_R2' can not be found, the sample is treated as"
            " single-end. Sample names are derived from the text found before the '_R1' string."
        )
    skipped_subsample = []
    skipped_assemble = []
    if args.sample_reads_target > 0:
        fastqs_to_subsample, skipped_subsample = find_and_match_fastqs(args.reads)
        skip_subsampling = False
        _, reformat_version, reformat_status = bbtools_path_version(args.reformat_path)
        intro_msg += (
            f" MEGAHIT de novo assemblies will start after subsampling {args.sample_reads_target}"
            " read pairs (or single-end reads)."
        )
    else:
        fastqs_to_assemble, skipped_assemble = find_and_match_fastqs(args.reads, recursive=False)
        skip_subsampling = True
        reformat_version, reformat_status = "", "not used"
        intro_msg += (
            "The full set of reads per sample will be assembled, no subsampling will be performed."
        )
    log.log_explanation(intro_msg, extra_empty_lines_after=0)
    log.log_explanation("For more information, please see https://github.com/edgardomortiz/Captus")

    log.log(f"{'Captus version':>{mar}}: {bold(f'v{__version__}')}")
    log.log(f"{'Command':>{mar}}: {bold(full_command)}")
    log.log(f"{'OS':>{mar}}: {bold(platform.platform())}")
    log.log(f"{'Host':>{mar}}: {bold(platform.node())}")
    tsv_comment = f"#Captus v{__version__}\n#Command: {full_command}\n"
    ram_B, ram_MB, ram_GB, ram_GB_total = set_ram(args.ram)
    log.log(f"{'Max. RAM':>{mar}}: {bold(f'{ram_GB:.1f}GB')} {dim(f'(out of {ram_GB_total:.1f}GB)')}")
    threads_max, threads_total = set_threads(args.threads)
    log.log(f"{'Max. Threads':>{mar}}: {bold(threads_max)} {dim(f'(out of {threads_total})')}")
    log.log("")

    log.log(f"{'Dependencies':>{mar}}:")
    _, megahit_version, megahit_status = megahit_path_version(args.megahit_path)
    log.log(format_dep_msg(f"{'MEGAHIT':>{mar}}: ", megahit_version, megahit_status))
    _, megahit_tk_version, megahit_tk_status = megahit_tk_path_version(args.megahit_toolkit_path)
    log.log(format_dep_msg(f"{'megahit_toolkit':>{mar}}: ", megahit_tk_version, megahit_tk_status))
    log.log(format_dep_msg(f"{'BBTools':>{mar}}: ", reformat_version, reformat_status))
    if args.disable_mapping:
        salmon_version, salmon_status = "", "not used"
    else:
        _, salmon_version, salmon_status = salmon_path_version(args.salmon_path)
    log.log(format_dep_msg(f"{'Salmon':>{mar}}: ", salmon_version, salmon_status))
    log.log("")

    log.log(f"{'Python libraries':>{mar}}:")
    numpy_found, numpy_version, numpy_status = python_library_check("numpy")
    pandas_found, pandas_version, pandas_status = python_library_check("pandas")
    plotly_found, plotly_version, plotly_status = python_library_check("plotly")
    log.log(format_dep_msg(f"{'numpy':>{mar}}: ", numpy_version, numpy_status))
    log.log(format_dep_msg(f"{'pandas':>{mar}}: ", pandas_version, pandas_status))
    log.log(format_dep_msg(f"{'plotly':>{mar}}: ", plotly_version, plotly_status))
    log.log("")

    log.log(f"{'Output directory':>{mar}}: {bold(out_dir)}")
    log.log(f"{'':>{mar}}  {dim(out_dir_msg)}")
    log.log("")

    if reformat_status == "not found":
        skip_subsampling = True
        log.log(
            f"{bold('WARNING:')} reformat.sh from BBTools could not be found, the reads will not be"
            " subsampled. Please verify you have it installed or provide the full path to the"
            " program with '--reformat_path'"
        )
    if megahit_status == "not found":
        quit_with_error(
            "MEGAHIT could not be found, please verify you have it installed or provide the full"
            " path to the progam with '--megahit_path'"
        )
    if salmon_status == "not found":
        args.disable_mapping = True
        log.log(
            f"{bold('WARNING:')} Salmon could not be found, the reads will not be mapped back to the"
            " contigs for accurate depth of coverage estimation. Please verify you have it installed"
            " or provide the full path to the program with '--salmon_path'"
        )

    ################################################################################################
    ############################################################################ SUBSAMPLING SECTION
    log.log_section_header("Subsampling Reads with reformat.sh")
    log.log_explanation(
        f"Now Captus will randomly subsample {args.sample_reads_target} read pairs (or single-end"
        " reads) from each sample prior to de novo assembly with MEGAHIT."
    )
    if skip_subsampling:
        log.log(
            red(
                "Skipping read subsampling step... (to enable provide a number of reads to subsample"
                " with '--sample_reads_target')"
            )
        )
        log.log("")
    else:
        if not fastqs_to_subsample:
            quit_with_error(
                "Captus could not find FASTQs to subsample, please verify your '--reads' argument"
            )

        log.log(f"{'Reads to subsample':>{mar}}: {bold(args.sample_reads_target)}")
        log.log(f"{'Overwrite files':>{mar}}: {bold(args.overwrite)}")
        log.log(f"{'Keep all files':>{mar}}: {bold(args.keep_all)}")
        log.log(f"{'Samples to process':>{mar}}: {bold(len(fastqs_to_subsample))}")
        log.log("")
        log.log(
            f"{'Output directories':>{mar}}:"
            f" {bold(f'{out_dir}/[Sample_name]__captus-asm/00_subsampled_reads')}"
        )
        log.log(f"{'':>{mar}}  {dim('A directory will be created for every sample')}")
        log.log("")

        if skipped_subsample:
            log.log(f"{bold('WARNING:')} {len(skipped_subsample)} sample(s) will be skipped")
            for msg in skipped_subsample:
                log.log(msg)
            log.log("")

        ram_B, ram_MB, ram_GB, ram_GB_total = set_ram(args.ram, java=True)
        subsample_reads_params = []
        for fastq_r1 in sorted(fastqs_to_subsample):
            in_fastq = fastq_r1
            if fastqs_to_subsample[fastq_r1]["fastq_r2"] is not None:
                if "_R1." in in_fastq:
                    in_fastq = fastq_r1.replace("_R1.", "_R#.")
                elif "_R1_" in in_fastq:
                    in_fastq = fastq_r1.replace("_R1_", "_R#_")
            subsample_reads_params.append(
                (
                    args.reformat_path,
                    ram_MB,
                    threads_max,
                    fastqs_to_subsample[fastq_r1]["fastq_dir"],
                    in_fastq,
                    out_dir,
                    args.sample_reads_target,
                    args.overwrite,
                )
            )
        tqdm_serial_run(
            subsample_reads,
            subsample_reads_params,
            "Subsampling reads",
            "Reads subsampling completed",
            "sample",
            args.show_less,
        )
        log.log("")

        fastqs_to_assemble = find_and_match_subsampled_fastqs(fastqs_to_subsample, out_dir)

    ################################################################################################
    ####################################################################### MEGAHIT ASSEMBLY SECTION
    log.log_section_header("De Novo Assembly with MEGAHIT")
    assembly_msg = (
        "Now Captus will perform de novo assembly with your input reads using MEGAHIT. Both"
        " '--min_contig_len' (when set to 'auto') and '--k_list' will be adjusted sample-wise"
        " according to the mean read length of the sample's FASTQ files."
    )
    if args.sample_reads_target > 0:
        assembly_msg += (
            "Since you specified a value for '--sample_reads_target', MEGAHIT will use the"
            " subsampled reads."
        )
    log.log_explanation(assembly_msg)

    if not fastqs_to_assemble:
        quit_with_error(
            "Captus could not find FASTQs to assemble, please verify your '--reads' argument or"
            " verify that the reads were correctly subsampled"
        )

    ram_B, ram_MB, ram_GB, ram_GB_total = set_ram(args.ram)
    concurrent, threads_per_assembly, ram_B_per_assembly = adjust_megahit_concurrency(
        args.concurrent, threads_max, ram_B, len(fastqs_to_assemble), args.preset
    )
    log.log(f"{'Concurrent assemblies':>{mar}}: {bold(concurrent)}")
    log.log(f"{'RAM per assembly':>{mar}}: {bold(f'{ram_B_per_assembly / 1024**3:.1f}GB')}")
    log.log(f"{'Threads per assembly':>{mar}}: {bold(threads_per_assembly)}")
    log.log("")

    if not args.preset:
        args.preset = "CAPSKIM"
    if args.preset.upper() not in settings.MEGAHIT_PRESETS:
        invalid_preset = args.preset
        args.preset = "CAPSKIM"
        log.log(f"{'preset':>{mar}}: {bold(args.preset)} ({invalid_preset} is not a valid preset)")
    else:
        args.preset = args.preset.upper()
        log.log(f"{'preset':>{mar}}: {bold(args.preset.upper())}")
    if not args.k_list:
        args.k_list = settings.MEGAHIT_PRESETS[args.preset]["k_list"]
    log.log(f"{'k_list':>{mar}}: {bold(args.k_list)}")
    if not args.min_count:
        args.min_count = settings.MEGAHIT_PRESETS[args.preset]["min_count"]
    log.log(f"{'min_count':>{mar}}: {bold(args.min_count)}")
    if not args.prune_level:
        args.prune_level = settings.MEGAHIT_PRESETS[args.preset]["prune_level"]
    log.log(f"{'prune_level':>{mar}}: {bold(args.prune_level)}")
    log.log(f"{'merge_level':>{mar}}: {bold(args.merge_level)}")
    log.log(f"{'min_contig_len':>{mar}}: {bold(args.min_contig_len)}")
    extra_options = settings.MEGAHIT_PRESETS[args.preset]["extra_options"]
    log.log(f"{'extra_options':>{mar}}: {bold(extra_options)}")
    tmp_dir = make_tmp_dir_within(args.tmp_dir, "captus_megahit_tmp")
    log.log(f"{'tmp_dir':>{mar}}: {bold(tmp_dir)}")
    log.log(f"{'max_contig_gc':>{mar}}: {bold(args.max_contig_gc)}%")
    log.log(f"{'disable_mapping':>{mar}}: {bold(args.disable_mapping)}")
    if args.min_contig_depth == "auto":
        msg_min_contig_depth = "1.0x (auto)" if args.disable_mapping is True else "1.5 (auto)"
    else:
        try:
            msg_min_contig_depth = f"{float(args.min_contig_depth)}x"
        except ValueError:
            quit_with_error(
                f"'--min_contig_depth should be a decimal, you provided '{args.min_contig_depth}'"
            )
    log.log(f"{'min_contig_depth':>{mar}}: {bold(msg_min_contig_depth)}")
    salmon_ran = bool(list(Path(out_dir).resolve().rglob("quant.sf")))
    if args.disable_mapping is False and salmon_ran is False:
        args.redo_filtering = True
    log.log(f"{'redo_filtering':>{mar}}: {bold(args.redo_filtering)}")
    log.log("")
    log.log(f"{'Overwrite files':>{mar}}: {bold(args.overwrite)}")
    log.log(f"{'Keep all files':>{mar}}: {bold(args.keep_all)}")
    log.log(f"{'Samples to assemble':>{mar}}: {bold(len(fastqs_to_assemble))}")
    log.log("")
    log.log(f"{'Output directories':>{mar}}: {bold(f'{out_dir}/[Sample_name]__captus-asm/01_assembly')}")
    log.log(f"{'':>{mar}}  {dim('A directory will be created for each sample')}")
    log.log("")

    if skipped_assemble:
        log.log(f"{bold('WARNING:')} {len(skipped_assemble)} sample(s) will be skipped")
        for msg in skipped_assemble:
            log.log(msg)
        log.log("")

    megahit_params = []
    for fastq_r1 in sorted(fastqs_to_assemble):
        sample_name = "_R1".join(fastq_r1.split("_R1")[:-1])
        if "_R1." in fastq_r1:
            sample_name = "_R1.".join(fastq_r1.split("_R1.")[:-1])
        elif "_R1_" in fastq_r1:
            sample_name = "_R1_".join(fastq_r1.split("_R1_")[:-1])
        fastqs_to_assemble[fastq_r1]["sample_name"] = sample_name
        read_stats = get_read_stats(
            Path(fastqs_to_assemble[fastq_r1]["fastq_dir"], fastq_r1),
            settings.NUM_READS_TO_CALCULATE_STATS,
        )
        fastqs_to_assemble[fastq_r1]["min_read_length"] = read_stats["min_read_length"]
        fastqs_to_assemble[fastq_r1]["max_read_length"] = read_stats["max_read_length"]
        fastqs_to_assemble[fastq_r1]["mean_read_length"] = read_stats["mean_read_length"]
        megahit_params.append(
            (
                args.megahit_path,
                args.megahit_toolkit_path,
                fastqs_to_assemble[fastq_r1]["fastq_dir"],
                fastq_r1,
                fastqs_to_assemble[fastq_r1]["fastq_r2"],
                sample_name,
                read_stats["mean_read_length"],
                args.k_list,
                args.min_count,
                args.prune_level,
                args.merge_level,
                ram_B_per_assembly,
                threads_per_assembly,
                out_dir,
                args.min_contig_len,
                extra_options,
                tmp_dir,
                tsv_comment,
                args.keep_all,
                args.overwrite,
            )
        )

    if args.debug:
        tqdm_serial_run(
            megahit,
            megahit_params,
            "De novo assembling with MEGAHIT",
            "De novo assembly completed",
            "sample",
            args.show_less,
        )
    else:
        tqdm_parallel_async_run(
            megahit,
            megahit_params,
            "De novo assembling with MEGAHIT",
            "De novo assembly completed",
            "sample",
            concurrent,
            args.show_less,
        )
    log.log("")

    if args.disable_mapping is False:
        salmon_params = []
        for fastq_r1 in sorted(fastqs_to_assemble):
            salmon_params.append(
                (
                    args.salmon_path,
                    fastqs_to_assemble[fastq_r1]["fastq_dir"],
                    fastq_r1,
                    fastqs_to_assemble[fastq_r1]["fastq_r2"],
                    fastqs_to_assemble[fastq_r1]["sample_name"],
                    args.k_list,
                    fastqs_to_assemble[fastq_r1]["max_read_length"],
                    threads_max,
                    out_dir,
                    tsv_comment,
                    args.keep_all,
                    args.overwrite,
                )
            )

        tqdm_serial_run(
            salmon,
            salmon_params,
            "Depth of coverage estimation with Salmon",
            "Depth of coverage estimation completed",
            "sample",
            args.show_less,
        )
        log.log("")

    filter_assembly_params = []
    for fastq_r1 in sorted(fastqs_to_assemble):
        filter_assembly_params.append(
            (
                fastqs_to_assemble[fastq_r1]["sample_name"],
                args.max_contig_gc,
                args.disable_mapping,
                args.min_contig_depth,
                args.redo_filtering,
                out_dir,
                tsv_comment,
                args.overwrite,
            )
        )
    symbol = "<"
    if args.disable_mapping is True and args.min_contig_depth == "auto":
        symbol = "<="
    filtering_msg = (
        f"Filtering contigs >{args.max_contig_gc}% GC and {symbol}{msg_min_contig_depth} depth"
    )
    if args.debug:
        tqdm_serial_run(
            filter_assembly,
            filter_assembly_params,
            filtering_msg,
            "Filtering completed",
            "sample",
            args.show_less,
        )
    else:
        tqdm_parallel_async_run(
            filter_assembly,
            filter_assembly_params,
            filtering_msg,
            "Filtering completed",
            "sample",
            concurrent,
            args.show_less,
        )
    log.log("")

    asm_stats_tsv, dep_stats_tsv, len_stats_tsv = collect_asm_stats(out_dir, tsv_comment)
    if asm_stats_tsv:
        log.log(f"{'Assembly statistics':>{mar}}: {bold(asm_stats_tsv)}")
        log.log(f"{'Depth statistics':>{mar}}: {bold(dep_stats_tsv)}")
        log.log(f"{'Length statistics':>{mar}}: {bold(len_stats_tsv)}")
        log.log("")
        if all([numpy_found, pandas_found, plotly_found]):
            from .report import build_assembly_report

            log.log_explanation("Generating Assembly report...")
            asm_html_report, asm_html_msg = build_assembly_report(
                out_dir, asm_stats_tsv, len_stats_tsv, dep_stats_tsv
            )
            log.log(f"{'Assembly report':>{mar}}: {bold(asm_html_report)}")
            log.log(f"{'':>{mar}}  {dim(asm_html_msg)}")
            log.log("")
        else:
            log.log(
                f"{bold('WARNING:')} Captus uses 'numpy', 'pandas', and 'plotly' to generate  an HTML"
                " report based on the assembly statistics. At least one of these libraries could not be"
                " found, please verify these libraries are installed and available."
            )
            log.log("")
    else:
        log.log(red("Skipping summarization step... (no assembly statistics files were produced)"))
        log.log("")

    shutil.rmtree(tmp_dir, ignore_errors=True)
    log.log(f"MEGAHIT temporary directory '{tmp_dir}' deleted")
    log.log("")

    ################################################################################################
    ################################################################################# ENDING SECTION
    successful_exit(
        "Captus-assembly: ASSEMBLE -> successfully completed"
        f" [{elapsed_time(time.time() - captus_start)}]"
    )


def subsample_reads(
    reformat_path, ram_MB, threads, in_dir, in_fastq, out_dir, sample_reads_target, overwrite
):
    """
    Subsample reads using reformat.sh from BBTools
    """
    start = time.time()

    sample_name = "_R".join(in_fastq.split("_R")[:-1])
    if "_R#" in in_fastq:
        sample_name = "_R#".join(in_fastq.split("_R#")[:-1])
    elif "_R1." in in_fastq:
        sample_name = "_R1.".join(in_fastq.split("_R1.")[:-1])
    elif "_R1_" in in_fastq:
        sample_name = "_R1_".join(in_fastq.split("_R1_")[:-1])

    sample_subsampled_out_dir, _ = make_output_dir(
        Path(out_dir, f"{sample_name}__captus-asm", "00_subsampled_reads")
    )
    reformat_cmd = [
        reformat_path,
        f"-Xmx{ram_MB}m",
        f"threads={threads}",
        f"in={Path(in_dir, in_fastq)}",
        f"out={Path(sample_subsampled_out_dir, in_fastq)}",
        f"samplereadstarget={sample_reads_target}",
        f"ziplevel={settings.REFORMAT_ZIPLEVEL}",
    ]
    reformat_log_file = Path(sample_subsampled_out_dir, f"{sample_name}_subsampling.log")
    if overwrite is True or not reformat_log_file.exists():
        reformat_cmd += ["overwrite=t"]
        with open(reformat_log_file, "w") as reformat_log:
            subprocess.run(reformat_cmd, stderr=reformat_log)
        message = f"'{sample_name}': reads subsampled [{elapsed_time(time.time() - start)}]"
    else:
        message = dim(f"'{sample_name}': SKIPPED (output files already exist)")

    return message


def find_and_match_subsampled_fastqs(fastqs_to_subsample, out_dir):
    """
    Find subsampled reads to build the dictionary 'fastqs_to_assemble'
    """
    fastqs_to_assemble = {}
    for fastq_r1 in fastqs_to_subsample:
        fastq_r2 = fastqs_to_subsample[fastq_r1]["fastq_r2"]
        sample_name = "_R1".join(fastq_r1.split("_R1")[:-1])
        if "_R1." in fastq_r1:
            sample_name = "_R1.".join(fastq_r1.split("_R1.")[:-1])
        elif "_R1_" in fastq_r1:
            sample_name = "_R1_".join(fastq_r1.split("_R1_")[:-1])
        sample_subsampled_out_dir = Path(out_dir, f"{sample_name}__captus-asm", "00_subsampled_reads")
        if Path(sample_subsampled_out_dir, fastq_r1).exists():
            fastqs_to_assemble[fastq_r1] = {
                "fastq_dir": f"{sample_subsampled_out_dir}",
                "fastq_r2": None,
            }
        if fastq_r2 is not None and Path(sample_subsampled_out_dir, fastq_r2).exists():
            fastqs_to_assemble[fastq_r1]["fastq_r2"] = fastq_r2

    return fastqs_to_assemble


def adjust_megahit_concurrency(concurrent, threads_max, ram_B, num_samples, preset):
    """
    Adjust the proposed number of 'concurrent' MEGAHIT processes so 'threads_per_assembly' are never
    fewer than 'settings.MEGAHIT_MIN_THREADS' or 'ram_B_per_assembly' is smaller than
    'settings.MEGAHIT_MIN_RAM_B'. Once the right 'concurrent' number has been found
    'threads_per_assembly' and 'ram_per_assembly' are readjusted. If the computer has fewer than 4
    CPUs or less than 4GB of RAM 'concurrent' defaults to 1 and MEGAHIT uses whatever resources are
    available. If assembly presets are used, then adjust minimum RAM per assembly according to the
    preset (e.g. min 8GB for RNAseq or WGS)
    """

    if concurrent == "auto":
        concurrent = min(threads_max, num_samples)
    else:
        try:
            concurrent = int(concurrent)
        except ValueError:
            quit_with_error("Invalid value for '--concurrent', set it to 'auto' or use a number")

    min_ram_B = settings.MEGAHIT_MIN_RAM_B
    min_threads = settings.MEGAHIT_MIN_THREADS
    if preset:
        if preset.upper() in settings.MEGAHIT_PRESETS:
            min_ram_B = settings.MEGAHIT_PRESETS[preset.upper()]["min_ram_B"]
            min_threads = settings.MEGAHIT_PRESETS[preset.upper()]["min_threads"]

    if threads_max < min_threads or ram_B < min_ram_B:
        return 1, threads_max, ram_B

    if concurrent > 1:
        while threads_max // concurrent < min_threads:
            concurrent -= 1
    ram_B_per_assembly = ram_B // concurrent

    if ram_B_per_assembly < min_ram_B:
        while ram_B // concurrent < min_ram_B:
            concurrent -= 1

    threads_per_assembly = threads_max // concurrent
    ram_B_per_assembly = ram_B // concurrent

    return concurrent, threads_per_assembly, ram_B_per_assembly


def megahit(
    megahit_path,
    megahit_toolkit_path,
    fastq_dir,
    fastq_r1,
    fastq_r2,
    sample_name,
    mean_read_length,
    k_list,
    min_count,
    prune_level,
    merge_level,
    ram_B,
    threads,
    out_dir,
    min_contig_len,
    extra_options,
    tmp_dir,
    tsv_comment,
    keep_all,
    overwrite,
):
    """
    De novo assembly with MEGAHIT >= v1.2.9
    """
    start = time.time()
    sample_out_dir, _ = make_output_dir(Path(out_dir, f"{sample_name}__captus-asm"))

    # MEGAHIT won't run if 'out_dir' already exists, it needs to create it by itself, we can only
    # create or verify that 'sample_out_dir' exists
    incomplete_log = Path(sample_out_dir, "megahit_brief.log")
    sample_megahit_out_dir = Path(sample_out_dir, "01_assembly")
    if overwrite is True or incomplete_log.exists() or not sample_megahit_out_dir.exists():
        if sample_megahit_out_dir.exists():
            shutil.rmtree(sample_megahit_out_dir, ignore_errors=True)
        if mean_read_length is False:
            message = red(
                f"'{sample_name}': SKIPPED (FASTQ files have .gz"
                " extension but are not compressed, please verify)"
            )
            return message
        adjusted_k_list = adjust_megahit_k_list(
            k_list, mean_read_length, settings.DELTA_MEAN_READ_LENGTH_TO_MAX_KMER_SIZE
        )
        adjusted_min_contig_len = adjust_megahit_min_contig_len(
            min_contig_len, mean_read_length, adjusted_k_list
        )
        megahit_cmd = [megahit_path]
        if fastq_r2:
            megahit_cmd += [
                "-1",
                f"{Path(fastq_dir, fastq_r1)}",
                "-2",
                f"{Path(fastq_dir, fastq_r2)}",
            ]
        else:
            megahit_cmd += ["--read", f"{Path(fastq_dir, fastq_r1)}"]
        megahit_cmd += [
            "--min-count",
            f"{min_count}",
            "--k-list",
            adjusted_k_list,
            "--merge-level",
            merge_level,
            "--prune-level",
            f"{prune_level}",
            "--memory",
            f"{ram_B}",
            "--num-cpu-threads",
            f"{threads}",
            "--out-dir",
            f"{sample_megahit_out_dir}",
            "--min-contig-len",
            f"{adjusted_min_contig_len}",
            "--tmp-dir",
            f"{tmp_dir}",
        ]
        if extra_options:
            megahit_cmd += [extra_options]
        megahit_log_file = Path(sample_out_dir, "megahit_brief.log")
        with open(megahit_log_file, "w") as megahit_log:
            megahit_log.write(f"Captus' MEGAHIT Command:\n  {' '.join(megahit_cmd)}\n\n\n")
        with open(megahit_log_file, "a") as megahit_log:
            subprocess.run(megahit_cmd, stdout=megahit_log, stderr=megahit_log)
        if cleanup_megahit_out_dir(
            sample_megahit_out_dir, megahit_log_file, megahit_toolkit_path, keep_all
        ):
            write_depth_coverage_tsv(sample_megahit_out_dir, "megahit", tsv_comment)
            message = f"'{sample_name}': assembled [{elapsed_time(time.time() - start)}]"
        else:
            message = red(f"'{sample_name}': FAILED assembly (no contigs were produced)")
    else:
        message = dim(f"'{sample_name}': SKIPPED (output files already exist)")

    return message


def adjust_megahit_k_list(k_list, mean_read_length, delta):
    """
    Remove largest kmers sizes from 'k_list' if they exceed 'mean_read_length' by at least 'delta'
    """
    return ",".join([k for k in k_list.split(",") if int(k) - mean_read_length <= delta])


def adjust_megahit_min_contig_len(min_contig_len, mean_read_length, k_list):
    """
    If 'min_contig_len' is set to 'auto' then adjust it to 'min_read_length' + smallest kmer size in
    'k_list'
    """
    if min_contig_len == "auto":
        return min(200, mean_read_length + int(k_list.split(",")[0]))
    elif int(min_contig_len) > 0:
        return min_contig_len
    else:
        return 200


def cleanup_megahit_out_dir(sample_megahit_out_dir, megahit_log_file, megahit_toolkit_path, keep_all):
    """
    Tidy up the output folder, MEGAHIT leaves many unnecesary trace files, especially in the
    'intermediate_contigs' directory
    """

    # Move/rename logs into place
    megahit_log_file.replace(Path(sample_megahit_out_dir, "megahit_brief.log"))
    Path(sample_megahit_out_dir, "log").replace(Path(sample_megahit_out_dir, "megahit_full.log"))
    intermediate_contigs_dir = Path(sample_megahit_out_dir, "intermediate_contigs")

    # Generate 'fastg' assembly graph for 'final.contigs.fa' and reformat FASTA headers
    assembly_file = Path(sample_megahit_out_dir, "final.contigs.fa")
    assembly_size = assembly_file.stat().st_size
    if assembly_size > 0:
        if megahit_toolkit_path:
            megahit_contig2fastg(
                megahit_toolkit_path,
                assembly_file,
                Path(sample_megahit_out_dir, "assembly_graph.fastg"),
            )
        fasta_reheaded, _ = fasta_headers_to_spades(fasta_to_dict(assembly_file))
        dict_to_fasta(fasta_reheaded, Path(sample_megahit_out_dir, "assembly.fasta"), wrap=80)

    # Delete unimportant/already-processed files
    assembly_file.unlink()
    if not keep_all:
        shutil.rmtree(intermediate_contigs_dir, ignore_errors=True)
        Path(sample_megahit_out_dir, "checkpoints.txt").unlink()
        Path(sample_megahit_out_dir, "done").unlink()
        Path(sample_megahit_out_dir, "options.json").unlink()

    return bool(assembly_size)


def salmon(
    salmon_path,
    fastq_dir,
    fastq_r1,
    fastq_r2,
    sample_name,
    k_list,
    max_read_length,
    threads,
    out_dir,
    tsv_comment,
    keep_all,
    overwrite,
):
    start = time.time()
    sample_megahit_out_dir = Path(out_dir, f"{sample_name}__captus-asm", "01_assembly")
    assembly_fasta = Path(sample_megahit_out_dir, "assembly.fasta")
    removed_fasta = Path(sample_megahit_out_dir, "removed_contigs.fasta")
    sample_index_dir = Path(sample_megahit_out_dir, settings.SALMON_INDEX_DIR)
    sample_quant_dir = Path(sample_megahit_out_dir, settings.SALMON_QUANT_DIR)
    k_min = k_list.split(",")[0]

    if not assembly_fasta.exists() or assembly_fasta.stat().st_size == 0:
        message = dim(f"'{sample_name}': SKIPPED (assembly file empty or not found)")
        return message

    if overwrite is True or not sample_quant_dir.exists():
        if removed_fasta.is_file():
            assembly = fasta_to_dict(assembly_fasta)
            removed = fasta_to_dict(removed_fasta)
            for seq_name in removed:
                assembly[seq_name] = removed[seq_name]
            dict_to_fasta(assembly, assembly_fasta, wrap=80)
            removed_fasta.unlink()
        if sample_index_dir.exists():
            shutil.rmtree(sample_index_dir, ignore_errors=True)
        if sample_quant_dir.exists():
            shutil.rmtree(sample_quant_dir, ignore_errors=True)
        salmon_index_cmd = [
            salmon_path,
            "index",
            "--transcripts",
            f"{assembly_fasta}",
            "--index",
            f"{sample_index_dir}",
            "--kmerLen",
            f"{k_min}",
        ]
        salmon_quant_cmd = [
            salmon_path,
            "quant",
            "--index",
            f"{sample_index_dir}",
            "--output",
            f"{sample_quant_dir}",
            "--validateMappings",
            "--recoverOrphans",
            "--threads",
            f"{threads}",
        ]
        if fastq_r2 is None:
            salmon_quant_cmd += [
                "--libType",
                "SF",
                "--unmatedReads",
                f"{Path(fastq_dir, fastq_r1)}",
            ]
        else:
            max_read_length *= 2
            salmon_quant_cmd += [
                "--libType",
                "IU",
                "--mates1",
                f"{Path(fastq_dir, fastq_r1)}",
                "--mates2",
                f"{Path(fastq_dir, fastq_r2)}",
            ]
        salmon_log_file = Path(sample_megahit_out_dir, "salmon.log")
        with open(salmon_log_file, "w") as salmon_log:
            salmon_log.write(f"Captus' Salmon Index Command:\n  {' '.join(salmon_index_cmd)}\n\n\n")
        with open(salmon_log_file, "a") as salmon_log:
            subprocess.run(salmon_index_cmd, stdout=salmon_log, stderr=salmon_log)
        with open(salmon_log_file, "a") as salmon_log:
            salmon_log.write(
                f"\n\n\nCaptus' Salmon Quant Command:\n  {' '.join(salmon_quant_cmd)}\n\n\n"
            )
        with open(salmon_log_file, "a") as salmon_log:
            subprocess.run(salmon_quant_cmd, stdout=salmon_log, stderr=salmon_log)
        if cleanup_salmon_out_dirs(sample_index_dir, sample_quant_dir, keep_all):
            write_depth_coverage_tsv(sample_megahit_out_dir, "salmon", tsv_comment, max_read_length)
            message = (
                f"'{sample_name}': depth of coverage estimated [{elapsed_time(time.time() - start)}]"
            )
        else:
            message = red(f"'{sample_name}': FAILED depth estimation (Salmon output not found)")
    else:
        message = dim(f"'{sample_name}': SKIPPED (output files already exist)")

    return message


def cleanup_salmon_out_dirs(sample_index_dir, sample_quant_dir, keep_all):
    """
    Remove Salmon index and other Salmon quant unnecessary files
    """
    # If the file quant.sf is not produced we should skip the cleanup and flag the task as failed
    quant_sf = Path(sample_quant_dir, "quant.sf")
    if not quant_sf.exists():
        return False
    if not keep_all:
        shutil.rmtree(sample_index_dir, ignore_errors=True)
        shutil.rmtree(Path(sample_quant_dir, "aux_info"), ignore_errors=True)
        shutil.rmtree(Path(sample_quant_dir, "libParams"), ignore_errors=True)
        shutil.rmtree(Path(sample_quant_dir, "logs"), ignore_errors=True)
        Path(sample_quant_dir, "cmd_info.json").unlink()
        Path(sample_quant_dir, "lib_format_counts.json").unlink()
    return True


def write_depth_coverage_tsv(sample_megahit_out_dir, depth_estimator, tsv_comment, max_read_length=150):
    fasta_asm_file = Path(sample_megahit_out_dir, "assembly.fasta")
    fasta_asm = fasta_to_dict(fasta_asm_file)

    tsv_header = "\t".join(
        [
            "megahit_contig_name",
            "megahit_depth",
            "length",
            "salmon_contig_name",
            "salmon_num_reads",
            "salmon_depth",
        ]
    )
    if depth_estimator == "megahit":
        with open(Path(sample_megahit_out_dir, settings.CONTIGS_DEPTH), "w") as tsv_out:
            tsv_out.write(tsv_comment)
            tsv_out.write(f"{tsv_header}\n")
            for ctg_name in fasta_asm:
                name_parts = ctg_name.split("_")
                tsv_record = "\t".join(
                    [
                        ctg_name,
                        f"{float(name_parts[5])}",
                        f"{int(name_parts[3])}",
                        "NA",
                        "NA",
                        "NA",
                    ]
                )
                tsv_out.write(f"{tsv_record}\n")
    elif depth_estimator == "salmon":
        reheaded_fasta_asm = {}
        salmon_quant_file = Path(sample_megahit_out_dir, settings.SALMON_QUANT_DIR, "quant.sf")
        with open(Path(sample_megahit_out_dir, settings.CONTIGS_DEPTH), "w") as tsv_out:
            tsv_out.write(tsv_comment)
            tsv_out.write(f"{tsv_header}\n")
            with open(salmon_quant_file, "r") as quant:
                for line in quant:
                    if not line.startswith("Name"):
                        record = line.strip().split()
                        name_parts = record[0].split("_")
                        salmon_num_reads = float(record[4])
                        length = int(record[1])
                        salmon_depth = round((salmon_num_reads * max_read_length / length), 4)
                        salmon_contig_name = (
                            f"{'_'.join(name_parts[0:5])}_{salmon_depth}_{'_'.join(name_parts[6:])}"
                        )
                        tsv_record = "\t".join(
                            [
                                record[0],
                                f"{float(name_parts[5])}",
                                f"{length}",
                                f"{salmon_contig_name}",
                                f"{salmon_num_reads}",
                                f"{salmon_depth}",
                            ]
                        )
                        tsv_out.write(f"{tsv_record}\n")
                        reheaded_fasta_asm[salmon_contig_name] = fasta_asm[record[0]]
        dict_to_fasta(reheaded_fasta_asm, fasta_asm_file, wrap=80)

    return


def filter_assembly(
    sample_name,
    max_contig_gc,
    disable_mapping,
    min_contig_depth,
    redo_filtering,
    out_dir,
    tsv_comment,
    overwrite,
):
    start = time.time()
    sample_megahit_out_dir = Path(out_dir, f"{sample_name}__captus-asm", "01_assembly")
    assembly_fasta = Path(sample_megahit_out_dir, "assembly.fasta")
    removed_fasta = Path(sample_megahit_out_dir, "removed_contigs.fasta")
    assembly_stats_tsv = Path(sample_megahit_out_dir, "assembly_stats.tsv")
    depth_stats_tsv = Path(sample_megahit_out_dir, "depth_stats.tsv")
    length_stats_tsv = Path(sample_megahit_out_dir, "length_stats.tsv")

    if not assembly_fasta.exists() or assembly_fasta.stat().st_size == 0:
        message = dim(f"'{sample_name}': SKIPPED (assembly file empty or not found)")
        return message

    skip_gc = False
    if max_contig_gc >= 100 or max_contig_gc <= 0:
        skip_gc = True
    if min_contig_depth == "auto":
        min_contig_depth = 1.0001 if disable_mapping is True else 1.5
    else:
        min_contig_depth = float(min_contig_depth)
    skip_depth = False
    if min_contig_depth <= 0:
        skip_depth = True

    if overwrite is True or redo_filtering is True or not removed_fasta.is_file():
        if removed_fasta.is_file():
            assembly = fasta_to_dict(assembly_fasta)
            removed = fasta_to_dict(removed_fasta)
            for seq_name in removed:
                assembly[seq_name] = removed[seq_name]
            dict_to_fasta(assembly, assembly_fasta, wrap=80)
            removed_fasta.unlink()
        if assembly_stats_tsv.is_file():
            assembly_stats_tsv.unlink()
        if depth_stats_tsv.is_file():
            depth_stats_tsv.unlink()
        if length_stats_tsv.is_file():
            length_stats_tsv.unlink()
        unfiltered = fasta_to_dict(assembly_fasta)
        accepted = {}
        rejected = {}
        for seq_name in unfiltered:
            if skip_gc is True:
                gc = 100
            else:
                seq = unfiltered[seq_name]["sequence"].upper().replace("-", "").replace("N", "")
                gc = round(seq.replace("C", "G").count("G") / len(seq) * 100, 5)
            if skip_depth is True:
                depth = float("inf")
            else:
                depth = float(seq_name.split("_")[5])
            if gc > max_contig_gc or depth < min_contig_depth:
                rejected[seq_name] = unfiltered[seq_name]
            else:
                accepted[seq_name] = unfiltered[seq_name]
        dict_to_fasta(accepted, assembly_fasta, wrap=80, write_if_empty=True)
        dict_to_fasta(rejected, removed_fasta, wrap=80, write_if_empty=True)
        message = get_asm_stats(
            sample_name,
            unfiltered,
            accepted,
            assembly_stats_tsv,
            depth_stats_tsv,
            length_stats_tsv,
            tsv_comment,
        )
        message = f"{message}\n'{sample_name}': assembly filtered [{elapsed_time(time.time() - start)}]"
    else:
        message = dim(f"'{sample_name}': SKIPPED (output files already exist)")

    return message


def get_asm_stats(
    sample_name,
    asm_before,
    asm_after,
    asm_stats_tsv_path,
    depth_stats_tsv_path,
    length_stats_tsv_path,
    tsv_comment,
):
    def calc_asm_stats(
        sample_name, asm, stage, asm_stats_tsv_path, depth_stats_tsv_path, length_stats_tsv_path
    ):
        """
        4847 contigs, total 2895219 bp, min 183 bp, max 4250 bp, avg 597 bp, N50 673 bp
        """
        lengths = []
        depths = []
        depths_hist = {}
        for k in range(5, 105, 5):
            depths_hist[str(k / 10)] = {"num_contigs": 0, "length": 0, "fraction": 0}
        for k in range(15, 105, 5):
            depths_hist[str(k)] = {"num_contigs": 0, "length": 0, "fraction": 0}
        for k in range(150, 1050, 50):
            depths_hist[str(k)] = {"num_contigs": 0, "length": 0, "fraction": 0}
        for k in range(1500, 10500, 500):
            depths_hist[str(k)] = {"num_contigs": 0, "length": 0, "fraction": 0}
        depths_hist["inf"] = {"num_contigs": 0, "length": 0, "fraction": 0}
        lengths_hist = {}
        for k in range(50, 1050, 50):
            lengths_hist[str(k)] = {"num_contigs": 0, "length": 0, "fraction": 0}
        for k in range(500, 10500, 500):
            lengths_hist[str(k)] = {"num_contigs": 0, "length": 0, "fraction": 0}
        for k in range(5000, 105000, 5000):
            lengths_hist[str(k)] = {"num_contigs": 0, "length": 0, "fraction": 0}
        lengths_hist["inf"] = {"num_contigs": 0, "length": 0, "fraction": 0}
        gc_count = 0

        for seq in asm:
            length = len(asm[seq]["sequence"])
            depth = float(seq.split("_")[5])
            lengths.append(length)
            depths.append(depth)
            gc_count += asm[seq]["sequence"].count("G")
            gc_count += asm[seq]["sequence"].count("C")
            if depth == 0:
                depths_hist["0.5"]["num_contigs"] += 1
                depths_hist["0.5"]["length"] += length
            elif depth <= 10:
                key = str(float(math.ceil(depth / 0.5) * 0.5))
                depths_hist[key]["num_contigs"] += 1
                depths_hist[key]["length"] += length
            elif depth <= 100:
                key = str(math.ceil(depth / 5) * 5)
                depths_hist[key]["num_contigs"] += 1
                depths_hist[key]["length"] += length
            elif depth <= 1000:
                key = str(math.ceil(depth / 50) * 50)
                depths_hist[key]["num_contigs"] += 1
                depths_hist[key]["length"] += length
            elif depth <= 10000:
                key = str(math.ceil(depth / 500) * 500)
                depths_hist[key]["num_contigs"] += 1
                depths_hist[key]["length"] += length
            else:
                depths_hist["inf"]["num_contigs"] += 1
                depths_hist["inf"]["length"] += length
            if length <= 1000:
                key = str(math.ceil(length / 50) * 50)
                lengths_hist[key]["num_contigs"] += 1
                lengths_hist[key]["length"] += length
            elif length <= 10000:
                key = str(math.ceil(length / 500) * 500)
                lengths_hist[key]["num_contigs"] += 1
                lengths_hist[key]["length"] += length
            elif length <= 100000:
                key = str(math.ceil(length / 5000) * 5000)
                lengths_hist[key]["num_contigs"] += 1
                lengths_hist[key]["length"] += length
            else:
                lengths_hist["inf"]["num_contigs"] += 1
                lengths_hist["inf"]["length"] += length

        num_contigs = len(asm)
        tot_length = sum(lengths)
        min_length = min(lengths)
        max_length = max(lengths)
        avg_length = math.ceil(statistics.mean(lengths))
        med_length = math.ceil(statistics.median(lengths))
        avg_depth = round(statistics.mean(depths), 4)
        med_depth = round(statistics.median(depths), 4)
        gc = round(gc_count / tot_length * 100, 3)
        for D in depths_hist:
            depths_hist[D]["fraction"] = round(depths_hist[D]["length"] / tot_length * 100, 3)
        for L in lengths_hist:
            lengths_hist[L]["fraction"] = round(lengths_hist[L]["length"] / tot_length * 100, 3)
        p50_length, cum_length, n50, l50 = tot_length * 0.50, 0, 0, 0
        for L in sorted(lengths, reverse=True):
            cum_length += L
            l50 += 1
            if cum_length >= p50_length:
                n50 = L
                break
        p75_length, cum_length, n75, l75 = tot_length * 0.75, 0, 0, 0
        for L in sorted(lengths, reverse=True):
            cum_length += L
            l75 += 1
            if cum_length >= p75_length:
                n75 = L
                break

        num_1kbp = round(len(list(filter(lambda L: L >= 1000, lengths))) / num_contigs * 100, 3)
        num_2kbp = round(len(list(filter(lambda L: L >= 2000, lengths))) / num_contigs * 100, 3)
        num_5kbp = round(len(list(filter(lambda L: L >= 5000, lengths))) / num_contigs * 100, 3)
        num_10kbp = round(len(list(filter(lambda L: L >= 10000, lengths))) / num_contigs * 100, 3)
        num_20kbp = round(len(list(filter(lambda L: L >= 20000, lengths))) / num_contigs * 100, 3)
        num_50kbp = round(len(list(filter(lambda L: L >= 50000, lengths))) / num_contigs * 100, 3)

        len_1kbp = round(sum(list(filter(lambda L: L >= 1000, lengths))) / tot_length * 100, 3)
        len_2kbp = round(sum(list(filter(lambda L: L >= 2000, lengths))) / tot_length * 100, 3)
        len_5kbp = round(sum(list(filter(lambda L: L >= 5000, lengths))) / tot_length * 100, 3)
        len_10kbp = round(sum(list(filter(lambda L: L >= 10000, lengths))) / tot_length * 100, 3)
        len_20kbp = round(sum(list(filter(lambda L: L >= 20000, lengths))) / tot_length * 100, 3)
        len_50kbp = round(sum(list(filter(lambda L: L >= 50000, lengths))) / tot_length * 100, 3)

        asm_stats_record = "\t".join(
            [
                sample_name,
                stage,
                f"{num_contigs}",
                f"{num_1kbp}",
                f"{num_2kbp}",
                f"{num_5kbp}",
                f"{num_10kbp}",
                f"{num_20kbp}",
                f"{num_50kbp}",
                f"{tot_length}",
                f"{len_1kbp}",
                f"{len_2kbp}",
                f"{len_5kbp}",
                f"{len_10kbp}",
                f"{len_20kbp}",
                f"{len_50kbp}",
                f"{min_length}",
                f"{max_length}",
                f"{avg_length}",
                f"{med_length}",
                f"{avg_depth}",
                f"{med_depth}",
                f"{gc}",
                f"{n50}",
                f"{n75}",
                f"{l50}",
                f"{l75}",
            ]
        )
        with open(asm_stats_tsv_path, "a") as stats_tsv:
            stats_tsv.write(f"{asm_stats_record}\n")

        with open(depth_stats_tsv_path, "a") as stats_tsv:
            for D in depths_hist:
                depth_stats_record = "\t".join(
                    [
                        sample_name,
                        stage,
                        f"{D}",
                        f"{depths_hist[D]['length']}",
                        f"{depths_hist[D]['fraction']}",
                        f"{depths_hist[D]['num_contigs']}",
                    ]
                )
                stats_tsv.write(f"{depth_stats_record}\n")

        with open(length_stats_tsv_path, "a") as stats_tsv:
            for L in lengths_hist:
                length_stats_record = "\t".join(
                    [
                        sample_name,
                        stage,
                        f"{L}",
                        f"{lengths_hist[L]['length']}",
                        f"{lengths_hist[L]['fraction']}",
                        f"{lengths_hist[L]['num_contigs']}",
                    ]
                )
                stats_tsv.write(f"{length_stats_record}\n")

        msg = (
            f"'{sample_name}': {stage.upper()} {num_contigs:,} contigs, total {tot_length:,} bp,"
            f" min {min_length:,} bp, max {max_length:,} bp, avg {avg_length:,} bp, N50 {n50:,} bp"
        )
        return msg

    asm_stats_header = "\t".join(
        [
            "sample_name",
            "stage",
            "num_contigs",
            "pct_contigs_1kbp",
            "pct_contigs_2kbp",
            "pct_contigs_5kbp",
            "pct_contigs_10kbp",
            "pct_contigs_20kbp",
            "pct_contigs_50kbp",
            "total_length",
            "pct_length_1kbp",
            "pct_length_2kbp",
            "pct_length_5kbp",
            "pct_length_10kbp",
            "pct_length_20kbp",
            "pct_length_50kbp",
            "shortest_contig",
            "longest_contig",
            "avg_length",
            "median_length",
            "avg_depth",
            "median_depth",
            "gc",
            "N50",
            "N75",
            "L50",
            "L75",
        ]
    )
    depth_stats_header = "\t".join(
        [
            "sample_name",
            "stage",
            "depth_bin",
            "length",
            "fraction",
            "num_contigs",
        ]
    )
    length_stats_header = "\t".join(
        [
            "sample_name",
            "stage",
            "length_bin",
            "length",
            "fraction",
            "num_contigs",
        ]
    )
    with open(asm_stats_tsv_path, "w") as stats_tsv:
        stats_tsv.write(tsv_comment)
        stats_tsv.write(f"{asm_stats_header}\n")
    with open(depth_stats_tsv_path, "w") as stats_tsv:
        stats_tsv.write(tsv_comment)
        stats_tsv.write(f"{depth_stats_header}\n")
    with open(length_stats_tsv_path, "w") as stats_tsv:
        stats_tsv.write(tsv_comment)
        stats_tsv.write(f"{length_stats_header}\n")

    msg_before = calc_asm_stats(
        sample_name,
        asm_before,
        "before",
        asm_stats_tsv_path,
        depth_stats_tsv_path,
        length_stats_tsv_path,
    )
    msg_after = calc_asm_stats(
        sample_name,
        asm_after,
        "after",
        asm_stats_tsv_path,
        depth_stats_tsv_path,
        length_stats_tsv_path,
    )

    return f"{msg_before}\n{msg_after}"


def collect_asm_stats(out_dir, tsv_comment):
    assembly_tsv_files = sorted(list(Path(out_dir).resolve().rglob("assembly_stats.tsv")))
    depth_tsv_files = sorted(list(Path(out_dir).resolve().rglob("depth_stats.tsv")))
    length_tsv_files = sorted(list(Path(out_dir).resolve().rglob("length_stats.tsv")))

    assembly_stats_tsv = Path(out_dir, "captus-assemble_assembly_stats.tsv")
    depth_stats_tsv = Path(out_dir, "captus-assemble_depth_stats.tsv")
    length_stats_tsv = Path(out_dir, "captus-assemble_length_stats.tsv")

    if not assembly_tsv_files or not depth_tsv_files or not length_tsv_files:
        return None, None, None
    else:
        assembly_stats_header = "\t".join(
            [
                "sample_name",
                "stage",
                "num_contigs",
                "pct_contigs_1kbp",
                "pct_contigs_2kbp",
                "pct_contigs_5kbp",
                "pct_contigs_10kbp",
                "pct_contigs_20kbp",
                "pct_contigs_50kbp",
                "total_length",
                "pct_length_1kbp",
                "pct_length_2kbp",
                "pct_length_5kbp",
                "pct_length_10kbp",
                "pct_length_20kbp",
                "pct_length_50kbp",
                "shortest_contig",
                "longest_contig",
                "avg_length",
                "median_length",
                "avg_depth",
                "median_depth",
                "gc",
                "N50",
                "N75",
                "L50",
                "L75",
            ]
        )
        depth_stats_header = "\t".join(
            [
                "sample_name",
                "stage",
                "depth_bin",
                "length",
                "fraction",
                "num_contigs",
            ]
        )
        length_stats_header = "\t".join(
            [
                "sample_name",
                "stage",
                "length_bin",
                "length",
                "fraction",
                "num_contigs",
            ]
        )
        assembly_stats_header = f"{assembly_stats_header}\n"
        depth_stats_header = f"{depth_stats_header}\n"
        length_stats_header = f"{length_stats_header}\n"

        with open(assembly_stats_tsv, "w") as tsv_out:
            tsv_out.write(tsv_comment)
            tsv_out.write(assembly_stats_header)
            for tsv in assembly_tsv_files:
                with open(tsv, "rt") as tsv_in:
                    for line in tsv_in:
                        if not line.startswith("#") and line != assembly_stats_header:
                            tsv_out.write(line)

        with open(depth_stats_tsv, "w") as tsv_out:
            tsv_out.write(tsv_comment)
            tsv_out.write(depth_stats_header)
            for tsv in depth_tsv_files:
                with open(tsv, "rt") as tsv_in:
                    for line in tsv_in:
                        if not line.startswith("#") and line != depth_stats_header:
                            tsv_out.write(line)

        with open(length_stats_tsv, "w") as tsv_out:
            tsv_out.write(tsv_comment)
            tsv_out.write(length_stats_header)
            for tsv in length_tsv_files:
                with open(tsv, "rt") as tsv_in:
                    for line in tsv_in:
                        if not line.startswith("#") and line != length_stats_header:
                            tsv_out.write(line)

        return assembly_stats_tsv, depth_stats_tsv, length_stats_tsv


def megahit_contig2fastg(megahit_toolkit_path, in_fasta_path, out_fastg_path):
    """
    Given the 'in_fasta_path' to a 'contigs.fa' produced by MEGAHIT, find out kmer size from header
    and then run 'contig2fastg' to create a '.fastg.' into 'out_fastg_path'
    """
    kmer = ""
    with open(in_fasta_path, "r") as fasta_in:
        for line in fasta_in:
            if line.startswith(">"):
                kmer = line[2:].split("_")[0]
                break
    if int(kmer) >= settings.MIN_KMER_SIZE_FOR_FASTG:
        contig2fastg_cmd = [megahit_toolkit_path, "contig2fastg", kmer, f"{in_fasta_path}"]
        with open(out_fastg_path, "w") as fastg_out:
            subprocess.run(contig2fastg_cmd, stdout=fastg_out, stderr=None)
