#!/usr/bin/env python3
"""
Copyright 2020-2023 Edgardo M. Ortiz (e.ortiz.v@gmail.com)
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
import shutil
import statistics
import subprocess
import time
from pathlib import Path

from . import log, settings
from .bioformats import dict_to_fasta, fasta_headers_to_spades, fasta_to_dict, get_mean_read_length
from .misc import (bbtools_path_version, bold, dim, elapsed_time, find_and_match_fastqs,
                   format_dep_msg, make_output_dir, make_tmp_dir_within, megahit_path_version,
                   megahit_tk_path_version, python_library_check, quit_with_error, red, set_ram,
                   set_threads, successful_exit, tqdm_parallel_async_run, tqdm_serial_run)
from .version import __version__


def assemble(full_command, args):

    captus_start = time.time()
    out_dir, out_dir_msg = make_output_dir(args.out)
    log.logger = log.Log(Path(args.out, "captus-assembly_assemble.log"), stdout_verbosity_level=1)

    mar = 21  # Margin for aligning parameters and values

    ################################################################################################
    ############################################################################### STARTING SECTION
    log.log_section_header("Starting Captus-assembly: ASSEMBLE", single_newline=False)
    log.log_explanation(
        "Welcome to the de novo assembly step of Captus-assembly. In this step, Captus will use"
        " MEGAHIT to assemble your input reads. It is also possible to subsample a number of reads"
        " using reformat.sh from BBTools prior to assembly, this is useful while performing tests or"
        " when including samples with considerably higher sequencing depth in a dataset.",
        extra_empty_lines_after=0
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
            f"MEGAHIT de novo assemblies will start after subsampling {args.sample_reads_target}"
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

    log.log(f'{"Captus version":>{mar}}: {bold(f"v{__version__}")}')
    log.log(f'{"Command":>{mar}}: {bold(full_command)}')
    ram_B, ram_MB, ram_GB, ram_GB_total = set_ram(args.ram)
    log.log(f'{"Max. RAM":>{mar}}: {bold(f"{ram_GB:.1f}GB")} {dim(f"(out of {ram_GB_total:.1f}GB)")}')
    threads_max, threads_total = set_threads(args.threads)
    log.log(f'{"Max. Threads":>{mar}}: {bold(threads_max)} {dim(f"(out of {threads_total})")}')
    log.log("")

    log.log(f'{"Dependencies":>{mar}}:')
    _, megahit_version, megahit_status = megahit_path_version(args.megahit_path)
    log.log(format_dep_msg(f'{"MEGAHIT":>{mar}}: ', megahit_version, megahit_status))
    _, megahit_tk_version, megahit_tk_status = megahit_tk_path_version(args.megahit_toolkit_path)
    log.log(format_dep_msg(f'{"megahit_toolkit":>{mar}}: ', megahit_tk_version, megahit_tk_status))
    log.log(format_dep_msg(f'{"BBTools":>{mar}}: ', reformat_version, reformat_status))
    log.log("")

    log.log(f'{"Python libraries":>{mar}}:')
    numpy_found, numpy_version, numpy_status = python_library_check("numpy")
    pandas_found, pandas_version, pandas_status = python_library_check("pandas")
    plotly_found, plotly_version, plotly_status = python_library_check("plotly")
    log.log(format_dep_msg(f'{"numpy":>{mar}}: ', numpy_version, numpy_status))
    log.log(format_dep_msg(f'{"pandas":>{mar}}: ', pandas_version, pandas_status))
    log.log(format_dep_msg(f'{"plotly":>{mar}}: ', plotly_version, plotly_status))
    log.log("")

    log.log(f'{"Output directory":>{mar}}: {bold(out_dir)}')
    log.log(f'{"":>{mar}}  {dim(out_dir_msg)}')
    log.log("")

    if megahit_status == "not found":
        quit_with_error(
            "MEGAHIT could not be found, please verify you have it installed or provide the full"
            " path to the progam with '--megahit_path'"
        )
    if reformat_status == "not found":
        skip_subsampling = True
        log.log(
            f"{bold('WARNING:')} reformat.sh from BBTools could not be found, the reads will not be"
            " subsampled. Please verify you have it installed or provide the full path to the"
            " program with '--reformat_path'"
        )


    ################################################################################################
    ############################################################################ SUBSAMPLING SECTION
    log.log_section_header("Subsampling Reads with reformat.sh")
    log.log_explanation(
        f"Now Captus will randomly subsample {args.sample_reads_target} read pairs (or single-end"
        " reads) from each sample prior to de novo assembly with MEGAHIT."
    )
    if skip_subsampling:
        log.log(red(
            "Skipping read subsampling step... (to enable provide a number of reads to subsample"
            " with '--sample_reads_target')"
        ))
        log.log("")
    else:
        if not fastqs_to_subsample:
            quit_with_error(
                "Captus could not find FASTQs to subsample, please verify your '--reads' argument"
            )

        log.log(f'{"Reads to subsample":>{mar}}: {bold(args.sample_reads_target)}')
        log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
        log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
        log.log(f'{"Samples to process":>{mar}}: {bold(len(fastqs_to_subsample))}')
        log.log("")
        log.log(f'{"Output directories":>{mar}}:'
                f' {bold(f"{out_dir}/[Sample_name]__captus-asm/00_subsampled_reads")}')
        log.log(f'{"":>{mar}}  {dim("A directory will be created for every sample")}')
        log.log("")

        if skipped_subsample:
            log.log(f'{bold("WARNING:")} {len(skipped_subsample)} sample(s) will be skipped')
            for msg in skipped_subsample:
                log.log(msg)
            log.log("")

        subsample_reads_params = []
        for fastq_r1 in sorted(fastqs_to_subsample):
            in_fastq = fastq_r1
            if fastqs_to_subsample[fastq_r1]["fastq_r2"] is not None:
                if "_R1." in in_fastq:
                    in_fastq = fastq_r1.replace("_R1.", "_R#.")
                elif "_R1_" in in_fastq:
                    in_fastq = fastq_r1.replace("_R1_", "_R#_")
            subsample_reads_params.append((
                args.reformat_path,
                ram_MB,
                threads_max,
                fastqs_to_subsample[fastq_r1]["fastq_dir"],
                in_fastq,
                out_dir,
                args.sample_reads_target,
                args.overwrite
            ))
        tqdm_serial_run(subsample_reads, subsample_reads_params, "Subsampling reads",
                        "Reads subsampling completed", "sample", args.show_less)
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

    concurrent, threads_per_assembly, ram_B_per_assembly = adjust_megahit_concurrency(
        args.concurrent, threads_max, ram_B, len(fastqs_to_assemble), args.preset
    )
    log.log(f'{"Concurrent assemblies":>{mar}}: {bold(concurrent)}')
    log.log(f'{"RAM per assembly":>{mar}}: {bold(f"{ram_B_per_assembly / 1024 ** 3:.1f}GB")}')
    log.log(f'{"Threads per assembly":>{mar}}: {bold(threads_per_assembly)}')
    log.log("")

    if not args.preset:
        args.preset = "CAPSKIM"
    if args.preset.upper() not in settings.MEGAHIT_PRESETS:
        invalid_preset = args.preset
        args.preset = "CAPSKIM"
        log.log(f'{"preset":>{mar}}: {bold(args.preset)} ({invalid_preset} is not a valid preset)')
    else:
        args.preset = args.preset.upper()
        log.log(f'{"preset":>{mar}}: {bold(args.preset.upper())}')
    if not args.k_list:
        args.k_list = settings.MEGAHIT_PRESETS[args.preset]["k_list"]
    log.log(f'{"k_list":>{mar}}: {bold(args.k_list)}')
    if not args.min_count:
        args.min_count = settings.MEGAHIT_PRESETS[args.preset]["min_count"]
    log.log(f'{"min_count":>{mar}}: {bold(args.min_count)}')
    if not args.prune_level:
        args.prune_level = settings.MEGAHIT_PRESETS[args.preset]["prune_level"]
    log.log(f'{"prune_level":>{mar}}: {bold(args.prune_level)}')
    log.log(f'{"merge_level":>{mar}}: {bold(args.merge_level)}')
    log.log(f'{"min_contig_len":>{mar}}: {bold(args.min_contig_len)}')
    log.log(f'{"max_contig_gc":>{mar}}: {bold(args.max_contig_gc)}%')
    extra_options = settings.MEGAHIT_PRESETS[args.preset]["extra_options"]
    log.log(f'{"extra_options":>{mar}}: {bold(extra_options)}')
    tmp_dir = make_tmp_dir_within(args.tmp_dir, "captus_megahit_tmp")
    log.log(f'{"tmp_dir":>{mar}}: {bold(tmp_dir)}')
    log.log("")
    log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
    log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
    log.log(f'{"Samples to assemble":>{mar}}: {bold(len(fastqs_to_assemble))}')
    log.log("")
    log.log(f'{"Output directories":>{mar}}: {bold(f"{out_dir}/[Sample_name]__captus-asm/01_assembly")}')
    log.log(f'{"":>{mar}}  {dim("A directory will be created for each sample")}')
    log.log("")

    if skipped_assemble:
        log.log(f'{bold("WARNING:")} {len(skipped_assemble)} sample(s) will be skipped')
        for msg in skipped_assemble:
            log.log(msg)
        log.log("")

    megahit_params = []
    for fastq_r1 in sorted(fastqs_to_assemble):
        fastq_r2 = fastqs_to_assemble[fastq_r1]["fastq_r2"]
        fastq_dir = fastqs_to_assemble[fastq_r1]["fastq_dir"]
        megahit_params.append((
            args.megahit_path,
            args.megahit_toolkit_path,
            fastq_dir,
            fastq_r1,
            fastq_r2,
            args.k_list,
            args.min_count,
            args.prune_level,
            args.merge_level,
            ram_B_per_assembly,
            threads_per_assembly,
            out_dir,
            args.min_contig_len,
            args.max_contig_gc,
            extra_options,
            tmp_dir,
            args.keep_all,
            args.overwrite
        ))

    if args.debug:
        tqdm_serial_run(megahit, megahit_params,
                        "De novo assembling with MEGAHIT",
                        "De novo assembly completed",
                        "sample", args.show_less)
    else:
        tqdm_parallel_async_run(megahit, megahit_params,
                                "De novo assembling with MEGAHIT",
                                "De novo assembly completed",
                                "sample", concurrent, args.show_less)
    log.log("")

    asm_stats_tsv = collect_asm_stats(out_dir)
    if asm_stats_tsv:
        log.log(f'{"Assembly statistics":>{mar}}: {bold(asm_stats_tsv)}')
        log.log("")
        if all([numpy_found, pandas_found, plotly_found]):

            from .report import build_assembly_report

            log.log_explanation(
                "Generating Assembly report..."
            )
            asm_html_report, asm_html_msg = build_assembly_report(out_dir, asm_stats_tsv)
            log.log(f'{"Assembly report":>{mar}}: {bold(asm_html_report)}')
            log.log(f'{"":>{mar}}  {dim(asm_html_msg)}')
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
        f"ziplevel={settings.REFORMAT_ZIPLEVEL}"
    ]
    reformat_log_file = Path(sample_subsampled_out_dir, f"{sample_name}.subsampling.log")
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
                "fastq_r2": None
            }
        if Path(sample_subsampled_out_dir, fastq_r2).exists():
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
    if preset:
        if preset.upper() in settings.MEGAHIT_PRESETS:
            min_ram_B = settings.MEGAHIT_PRESETS[preset.upper()]["min_ram_B"]

    if threads_max < settings.MEGAHIT_MIN_THREADS or ram_B < min_ram_B:
        return 1, threads_max, ram_B

    if concurrent > 1:
        while threads_max // concurrent < settings.MEGAHIT_MIN_THREADS:
            concurrent -= 1
    ram_B_per_assembly = ram_B // concurrent

    if ram_B_per_assembly < min_ram_B:
        while ram_B // concurrent < min_ram_B:
            concurrent -= 1

    threads_per_assembly = threads_max // concurrent
    ram_B_per_assembly = ram_B // concurrent

    return concurrent, threads_per_assembly, ram_B_per_assembly


def megahit(
    megahit_path, megahit_toolkit_path, fastq_dir, fastq_r1, fastq_r2, k_list, min_count,
    prune_level, merge_level, ram_B, threads, out_dir, min_contig_len, max_contig_gc, extra_options,
    tmp_dir, keep_all, overwrite
):
    """
    De novo assembly with MEGAHIT >= v1.2.9
    """
    start = time.time()
    sample_name = "_R1".join(fastq_r1.split("_R1")[:-1])
    if "_R1." in fastq_r1:
        sample_name = "_R1.".join(fastq_r1.split("_R1.")[:-1])
    elif "_R1_" in fastq_r1:
        sample_name = "_R1_".join(fastq_r1.split("_R1_")[:-1])
    sample_out_dir, _ = make_output_dir(Path(out_dir, f"{sample_name}__captus-asm"))

    # MEGAHIT won't run if 'out_dir' already exists, it needs to create it by itself, we can only
    # create or verify that 'sample_out_dir' exists
    sample_megahit_out_dir = Path(sample_out_dir, "01_assembly")
    if overwrite is True or not sample_megahit_out_dir.exists():
        if sample_megahit_out_dir.exists():
            shutil.rmtree(sample_megahit_out_dir, ignore_errors=True)
        mean_read_length = get_mean_read_length(Path(fastq_dir, fastq_r1),
                                                settings.NUM_READS_TO_CALCULATE_MEAN_READ_LENGTH)
        if mean_read_length is False:
            message = red(f"'{sample_name}': SKIPPED (FASTQ files have .gz"
                           " extension but are not compressed, please verify)")
            return message
        adjusted_k_list = adjust_megahit_k_list(k_list, mean_read_length,
                                                settings.DELTA_MEAN_READ_LENGTH_TO_MAX_KMER_SIZE)
        adjusted_min_contig_len = adjust_megahit_min_contig_len(min_contig_len, mean_read_length,
                                                                adjusted_k_list)
        megahit_command = [megahit_path]
        if fastq_r2:
            megahit_command += [
                "-1", f"{Path(fastq_dir, fastq_r1)}",
                "-2", f"{Path(fastq_dir, fastq_r2)}"
            ]
        else:
            megahit_command += ["--read", f"{Path(fastq_dir, fastq_r1)}"]
        megahit_command += [
            "--min-count", f"{min_count}",
            "--k-list", adjusted_k_list,
            "--merge-level", merge_level,
            "--prune-level", f"{prune_level}",
            "--memory", f"{ram_B}",
            "--num-cpu-threads", f"{threads}",
            "--out-dir", f"{sample_megahit_out_dir}",
            "--min-contig-len", f"{adjusted_min_contig_len}",
            "--tmp-dir", f"{tmp_dir}",
        ]
        if extra_options:
            megahit_command += [extra_options]
        megahit_log_file = Path(sample_out_dir, "megahit.brief.log")
        with open(megahit_log_file, "w") as megahit_log:
            megahit_log.write(f"Captus' MEGAHIT Command:\n  {' '.join(megahit_command)}\n\n\n")
        with open(megahit_log_file, "a") as megahit_log:
            subprocess.run(megahit_command, stdout=megahit_log, stderr=megahit_log)
        cleanup_megahit_out_dir(sample_megahit_out_dir, megahit_log_file,
                                megahit_toolkit_path, keep_all)
        filter_assembly_by_gc(sample_megahit_out_dir, max_contig_gc)
        asm_stats = get_asm_stats(sample_megahit_out_dir)
        message = f"'{sample_name}': {asm_stats} [{elapsed_time(time.time() - start)}]"
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
    megahit_log_file.replace(Path(sample_megahit_out_dir, "megahit.brief.log"))
    Path(sample_megahit_out_dir, "log").replace(Path(sample_megahit_out_dir, "megahit.full.log"))
    intermediate_contigs_dir = Path(sample_megahit_out_dir, "intermediate_contigs")

    # Now that assembly has been optimized, the intermedite contigs don't need to be saved or
    # processed further

    # Generate 'fastg' asembly graphs for 'intermediate_contigs' and reformat FASTA headers
    # intermediate_contigs_graphs_dir = Path(sample_megahit_out_dir, "01_intermediate_contigs_graphs")
    # intermediate_contigs_graphs_dir.mkdir(parents=True)
    # fastas_to_process = list(intermediate_contigs_dir.glob("*[!final].contigs.fa"))
    # for fasta_path in fastas_to_process:
    #     if megahit_toolkit_path:
    #         megahit_contig2fastg(megahit_toolkit_path,
    #                              fasta_path,
    #                              Path(intermediate_contigs_graphs_dir, f"{fasta_path.stem}.fastg"))
    #     fasta_reheaded, _ = fasta_headers_to_spades(fasta_to_dict(fasta_path))
    #     dict_to_fasta(fasta_reheaded,
    #                   Path(intermediate_contigs_graphs_dir, fasta_path.name),
    #                   wrap=80)

    # Generate 'fastg' assembly graph for 'final.contigs.fa' and reformat FASTA headers
    if megahit_toolkit_path:
        megahit_contig2fastg(megahit_toolkit_path,
                             Path(sample_megahit_out_dir, "final.contigs.fa"),
                             Path(sample_megahit_out_dir, "assembly_graph.fastg"))
    fasta_reheaded, _ = fasta_headers_to_spades(fasta_to_dict(Path(sample_megahit_out_dir,
                                                                   "final.contigs.fa")))
    dict_to_fasta(fasta_reheaded, Path(sample_megahit_out_dir, "assembly.fasta"), wrap=80)

    # Delete unimportant/already-processed files
    Path(sample_megahit_out_dir, "final.contigs.fa").unlink()
    if not keep_all:
        shutil.rmtree(intermediate_contigs_dir, ignore_errors=True)
        Path(sample_megahit_out_dir, "checkpoints.txt").unlink()
        Path(sample_megahit_out_dir, "done").unlink()
        Path(sample_megahit_out_dir, "options.json").unlink()


def filter_assembly_by_gc(sample_megahit_out_dir, max_contig_gc):
    if max_contig_gc >= 100 or max_contig_gc <= 0:
        return
    assembly_path = Path(sample_megahit_out_dir, "assembly.fasta")
    unfiltered = fasta_to_dict(assembly_path)
    accepted = {}
    rejected = {}
    for seq_name in unfiltered:
        seq = unfiltered[seq_name]["sequence"].upper().replace("-", "").replace("N", "")
        gc = round(seq.replace("C", "G").count("G") / len(seq) * 100, 5)
        if gc > max_contig_gc:
            rejected[seq_name] = unfiltered[seq_name]
        else:
            accepted[seq_name] = unfiltered[seq_name]
    dict_to_fasta(accepted, assembly_path, wrap=80)
    dict_to_fasta(rejected, Path(sample_megahit_out_dir, "filtered_contigs.fasta"), wrap=80)
    return


def get_asm_stats(sample_megahit_out_dir):
    asm = fasta_to_dict(Path(sample_megahit_out_dir, "assembly.fasta"))

    sample = Path(sample_megahit_out_dir).parts[-2].replace("__captus-asm", "")

    lengths = []
    depths = []
    gc_count = 0
    for seq in asm:
        lengths.append(len(asm[seq]["sequence"]))
        gc_count += asm[seq]["sequence"].count("G")
        gc_count += asm[seq]["sequence"].count("C")
        depths.append(float(seq.split("_cov_")[1].split("_")[0]))

    n_least_0bp = len(lengths)
    n_least_1kbp = round(len(list(filter(lambda L: L >= 1000, lengths))) / n_least_0bp * 100, 3)
    n_least_2kbp = round(len(list(filter(lambda L: L >= 2000, lengths))) / n_least_0bp * 100, 3)
    n_least_5kbp = round(len(list(filter(lambda L: L >= 5000, lengths))) / n_least_0bp * 100, 3)
    n_least_10kbp = round(len(list(filter(lambda L: L >= 10000, lengths))) / n_least_0bp * 100, 3)

    longest = max(lengths)
    shortest = min(lengths)
    s_least_0bp = sum(lengths)
    s_least_1kbp = round(sum(list(filter(lambda L: L >= 1000, lengths))) / s_least_0bp * 100, 3)
    s_least_2kbp = round(sum(list(filter(lambda L: L >= 2000, lengths))) / s_least_0bp * 100, 3)
    s_least_5kbp = round(sum(list(filter(lambda L: L >= 5000, lengths))) / s_least_0bp * 100, 3)
    s_least_10kbp = round(sum(list(filter(lambda L: L >= 10000, lengths))) / s_least_0bp * 100, 3)
    avg_length = math.ceil(statistics.mean(lengths))
    median_length = math.ceil(statistics.median(lengths))
    asm_half, size, n50 = s_least_0bp / 2, 0, 0
    for L in sorted(lengths, reverse=True):
        size += L
        if size >= asm_half:
            n50 = L
            break

    gc = round(gc_count / s_least_0bp * 100, 3)

    avg_depth = round(statistics.mean(depths), 2)
    d_least_1x = round(len(list(filter(lambda d: d >= 1, depths))) / n_least_0bp * 100, 3)
    d_least_2x = round(len(list(filter(lambda d: d >= 2, depths))) / n_least_0bp * 100, 3)
    d_least_5x = round(len(list(filter(lambda d: d >= 5, depths))) / n_least_0bp * 100, 3)
    d_least_10x = round(len(list(filter(lambda d: d >= 10, depths))) / n_least_0bp * 100, 3)

    filt_fasta = Path(sample_megahit_out_dir, "filtered_contigs.fasta")
    filt_n_contigs, filt_total_length, filt_avg_length, filt_gc_content = 0, 0, 0, 0
    if filt_fasta.exists():
        filt_asm = fasta_to_dict(Path(sample_megahit_out_dir, "filtered_contigs.fasta"))
        filt_lengths = []
        filt_gc_count = 0
        for seq in filt_asm:
            filt_lengths.append(len(filt_asm[seq]["sequence"]))
            filt_gc_count += filt_asm[seq]["sequence"].count("G")
            filt_gc_count += filt_asm[seq]["sequence"].count("C")
        filt_n_contigs = len(filt_lengths)
        filt_total_length = sum(filt_lengths)
        filt_avg_length = math.ceil(statistics.mean(filt_lengths))
        filt_gc_content = round(filt_gc_count / filt_total_length * 100, 3)

    stats_tsv = [
        "sample", sample,
        "n_contigs", n_least_0bp,
        "pct_contigs_>=_1kbp", n_least_1kbp,
        "pct_contigs_>=_2kbp", n_least_2kbp,
        "pct_contigs_>=_5kbp", n_least_5kbp,
        "pct_contigs_>=_10kbp", n_least_10kbp,
        "longest_contig", longest,
        "shortest_contig", shortest,
        "total_length", s_least_0bp,
        "pct_length_>=_1kbp", s_least_1kbp,
        "pct_length_>=_2kbp", s_least_2kbp,
        "pct_length_>=_5kbp", s_least_5kbp,
        "pct_length_>=_10kbp", s_least_10kbp,
        "avg_length", avg_length,
        "median_length", median_length,
        "N50", n50,
        "GC_content",gc,
        "avg_depth", avg_depth,
        "pct_contigs_>=_1x", d_least_1x,
        "pct_contigs_>=_2x", d_least_2x,
        "pct_contigs_>=_5x", d_least_5x,
        "pct_contigs_>=_10x", d_least_10x,
        "filtered_n_contigs", filt_n_contigs,
        "filtered_total_length", filt_total_length,
        "filtered_avg_length", filt_avg_length,
        "filtered_gc_content", filt_gc_content,
    ]

    with open(Path(sample_megahit_out_dir, "assembly.stats.tsv"), "wt") as tsv_out:
        for i in range(0, len(stats_tsv), 2):
            tsv_out.write(f"{stats_tsv[i].rjust(21)} : {stats_tsv[i+1]}\n")
    with open(Path(sample_megahit_out_dir, "assembly.stats.t.tsv"), "wt") as tsv_out:
        tsv_out.write("\t".join(f"{stats_tsv[i+1]}" for i in range(0, len(stats_tsv), 2)) + "\n")

    asm_msg = (
        f"{s_least_0bp:,} bp in {n_least_0bp:,} contigs, from {shortest:,} to {longest:,} bp,"
        f" avg {avg_length:,} bp, median {median_length:,} bp, N50 {n50:,} bp"
    )

    return asm_msg


def collect_asm_stats(out_dir):
    tsv_files = sorted(list(Path(out_dir).resolve().rglob("*assembly.stats.t.tsv")))
    if not tsv_files:
        return None
    else:
        tsv = ["\t".join([
            "sample",
            "n_contigs",
            "pct_contigs_>=_1kbp",
            "pct_contigs_>=_2kbp",
            "pct_contigs_>=_5kbp",
            "pct_contigs_>=_10kbp",
            "longest_contig",
            "shortest_contig",
            "total_length",
            "pct_length_>=_1kbp",
            "pct_length_>=_2kbp",
            "pct_length_>=_5kbp",
            "pct_length_>=_10kbp",
            "avg_length",
            "median_length",
            "N50",
            "GC_content",
            "avg_depth",
            "pct_contigs_>=_1x",
            "pct_contigs_>=_2x",
            "pct_contigs_>=_5x",
            "pct_contigs_>=_10x",
            "filtered_n_contigs",
            "filtered_total_length",
            "filtered_avg_length",
            "filtered_gc_content",
        ])]
        for file in tsv_files:
            with open(file, "rt") as tsv_in:
                for line in tsv_in:
                    tsv.append(line.strip("\n"))
        stats_tsv_file = Path(out_dir, "captus-assembly_assemble.stats.tsv")
        with open(stats_tsv_file, "wt") as tsv_out:
            tsv_out.write("\n".join(tsv))
        return stats_tsv_file


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

