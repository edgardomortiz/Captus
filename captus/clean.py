#!/usr/bin/env python3
"""
Copyright 2020 Edgardo M. Ortiz (e.ortiz.v@gmail.com)
https://github.com/edgardomortiz/Captus

This file is part of Captus. Captus is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Captus is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Captus. If
not, see <http://www.gnu.org/licenses/>.
"""


import multiprocessing
import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path

from . import log
from . import settings_assembly as settings
from .misc import (bbtools_path_version, bold, dim, elapsed_time, execute_jupyter_report,
                   fastqc_path_version, find_and_match_fastqs, format_dep_msg, has_valid_ext,
                   make_output_dir, quit_with_error, red, set_ram, set_threads,
                   tqdm_parallel_async_run, tqdm_serial_run)
from .version import __version__


def clean(full_command, args):

    captus_start = time.time()
    out_dir, out_dir_msg = make_output_dir(args.out)
    log.logger = log.Log(Path(out_dir, "captus-assembly_clean.log"), stdout_verbosity_level=1)

    mar = 21  # Margin for aligning parameters and values

    ################################################################################################
    ############################################################################### STARTING SECTION
    log.log_section_header("Starting Captus-assembly: Clean", single_newline=False)
    log.log_explanation(
        "Welcome to the read cleaning step of Captus-assembly. In this step, Captus will perform"
        " adaptor trimming and quality filtering/trimming on your input reads using bbduk.sh from"
        " BBTools.", extra_empty_lines_after=0
    )
    if len(args.reads) == 1 and Path(args.reads[0]).is_dir():
        intro_msg = (
            "Since you provided a directory name, Captus will look in that location for all the"
            " FASTQ files that contain the string _R1 in their names and match them with their"
            " respective _R2 pairs. If the _R2 can not be found, the sample is treated as single-end."
            " Sample names are derived from the text found before the _R1 string."
        )
    else:
        intro_msg = (
            "Since you provided a file name or a list of file names, Captus will only process the"
            " FASTQ files that contain the string _R1 in their names and match them with their"
            " respective _R2 pairs. If the _R2 can not be found, the sample is treated as single-end."
            " Sample names are derived from the text found before the _R1 string."
        )
    log.log_explanation(intro_msg, extra_empty_lines_after=0)
    if args.skip_fastqc:
        fastqc_version, fastqc_status = "", "not used"
        intro_msg = ("FastQC will not be run after the cleaning has been completed.")
    else:
        _, fastqc_version, fastqc_status = fastqc_path_version(args.fastqc_path)
        intro_msg = ("FastQC will run after the cleaning has been completed.")
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
    _, bbduk_version, bbduk_status = bbtools_path_version(args.bbduk_path)
    log.log(format_dep_msg(f'{"BBTools":>{mar}}: ', bbduk_version, bbduk_status))
    log.log(format_dep_msg(f'{"FastQC":>{mar}}: ', fastqc_version, fastqc_status))
    log.log("")
    log.log(f'{"Output directory":>{mar}}: {bold(out_dir)}')
    log.log(f'{"":>{mar}}  {dim(out_dir_msg)}')
    log.log("")

    if bbduk_status == "not found":
        quit_with_error(
            "bbduk.sh could not be found, please verify you have it installed or provide the full"
            " path to the progam with '--bbduk_path'"
        )
    if fastqc_status == "not found":
        args.skip_fastqc = True
        log.log(
            f"{bold('WARNING:')} FastQC could not be found, the program will be skipped after"
            " cleaning the reads. Please verify you have it installed or provide the full path to"
            " the program with '--fastqc_path'\n"
        )


    ################################################################################################
    ####################################################################### ADAPTOR TRIMMING SECTION
    log.log_section_header("Adaptor Trimming with bbduk.sh")
    adaptor_msg = (
        "Now Captus will perform two simultaneous rounds of adaptor removal using the specified set"
        " of adaptors. Most of the parameters are hard-coded but can be modified directly in the"
        f" file: '{settings.SETTINGS_ASSEMBLY_PATH}'"
    )
    if args.rna:
        adaptor_msg += ". Since you enabled the '--rna' option, poly-A tails will also be trimmed."
    log.log_explanation(adaptor_msg)

    if args.adaptor_set.lower() == "illumina":
        adaptor_set = settings.ILLUMINA_ADAPTORS
        adaptor_set_msg = f'{bold("Illumina")} {dim(adaptor_set)}'
    elif args.adaptor_set.lower() in ["bgi", "bgiseq", "dnbseq", "mgiseq"]:
        adaptor_set = settings.BGISEQ_ADAPTORS
        adaptor_set_msg = f'{bold("BGI")} {dim(adaptor_set)}'
    elif Path(args.adaptor_set).is_file():
        adaptor_set = Path(args.adaptor_set).resolve()
        adaptor_set_msg = bold(adaptor_set)
    else:
        quit_with_error(f"{args.adaptor_set} not found. Please check the path")

    log.log(f'{"Adaptor set":>{mar}}: {adaptor_set_msg}')
    log.log(f'{"Trim poly-A tails":>{mar}}: {bold(args.rna)}')
    log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
    log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
    fastqs_raw = find_and_match_fastqs(args.reads)
    log.log("")
    adaptors_trimmed_dir, adaptors_trimmed_msg = make_output_dir(Path(out_dir, "00_adaptors_trimmed"))
    log.log(f'{"Output directory":>{mar}}: {bold(adaptors_trimmed_dir)}')
    log.log(f'{"":>{mar}}  {dim(adaptors_trimmed_msg)}')
    log.log("")

    if not fastqs_raw:
        quit_with_error("No FASTQ files to process. Please verify your '--reads' argument.")

    bbduk_trim_adaptors_params = []
    for fastq_r1 in sorted(fastqs_raw):
        if fastqs_raw[fastq_r1]["fastq_r2"] is not None:
            in_fastq = fastq_r1.replace("_R1", "_R#")
        else:
            in_fastq = fastq_r1
        bbduk_trim_adaptors_params.append((
            args.bbduk_path,
            ram_MB,
            threads_max,
            fastqs_raw[fastq_r1]["fastq_dir"],
            in_fastq,
            adaptors_trimmed_dir,
            adaptor_set,
            args.rna,
            args.overwrite
        ))
    tqdm_serial_run(bbduk_trim_adaptors, bbduk_trim_adaptors_params,
                    "Trimming adaptors", "Adaptor trimming completed", "sample", args.show_less)
    log.log("")


    ################################################################################################
    ######################################################### QUALITY TRIMMING AND FILTERING SECTION
    log.log_section_header("Quality Trimming and Filtering with bbduk.sh")
    log.log_explanation(
        "Now Captus will trim the regions of reads with PHRED quality score lower than '--trimq' and"
        " afterwards will remove reads with average PHRED quality lower than '--maq'. Common"
        " sequencing artifacts and reads belonging to the viral genome phiX174 will also be"
        " filtered out. Other parameters are hard-coded but can be modified directly in the file:"
        f" '{settings.SETTINGS_ASSEMBLY_PATH}'"
    )

    log.log(f'{"trimq":>{mar}}: {bold(args.trimq)}')
    log.log(f'{"maq":>{mar}}: {bold(args.maq)}')
    log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
    log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
    fastqs_no_adaptors = find_and_match_fastqs(adaptors_trimmed_dir)
    log.log(f'{"Samples to clean":>{mar}}: {bold(len(fastqs_no_adaptors))}')
    log.log("")
    log.log(f'{"Output directory":>{mar}}: {bold(out_dir)}')
    log.log("")

    if fastqs_no_adaptors:
        bbduk_filter_quality_params = []
        for fastq_r1 in sorted(fastqs_no_adaptors):
            if fastqs_no_adaptors[fastq_r1]["fastq_r2"] is not None:
                in_fastq = fastq_r1.replace("_R1", "_R#")
            else:
                in_fastq = fastq_r1
            bbduk_filter_quality_params.append((
                args.bbduk_path,
                ram_MB,
                threads_max,
                fastqs_no_adaptors[fastq_r1]["fastq_dir"],
                in_fastq,
                out_dir,
                args.trimq,
                args.maq,
                args.overwrite
            ))
        tqdm_serial_run(bbduk_filter_quality, bbduk_filter_quality_params,
                        "Trimming/filtering quality", "Quality trimming/filtering completed",
                        "sample", args.show_less)
        log.log("")
    else:
        log.log(red(
            "Quality Trimming and Filtering will be skipped. The reads must have been cleaned"
            " already. If you want to clean the reads again enable '--overwrite'"
        ))
        log.log("")


    ################################################################################################
    ################################################################################# FASTQC SECTION
    log.log_section_header("FastQC Analysis: Raw and Clean Reads")
    log.log_explanation("Now Captus will run FastQC on the raw reads as well as on the clean reads.")
    if  args.skip_fastqc:
        log.log(red("Skipping FastQC analyses... (to enable omit '--skip_fastqc')"))
        log.log("")
    else:
        fastqc_before_dir, fastqc_before_msg = make_output_dir(Path(out_dir, "01_fastqc_before"))
        fastqc_after_dir, fastqc_after_msg = make_output_dir(Path(out_dir, "02_fastqc_after"))
        fastqc_params = []
        for fastq_r1 in sorted(fastqs_raw):
            fastqc_params.append((
                args.fastqc_path,
                Path(fastqs_raw[fastq_r1]["fastq_dir"], fastq_r1),
                fastqc_before_dir,
                args.overwrite,
                "BEFORE"
            ))
            if fastqs_raw[fastq_r1]["fastq_r2"] is not None:
                fastqc_params.append((
                    args.fastqc_path,
                    Path(fastqs_raw[fastq_r1]["fastq_dir"], fastqs_raw[fastq_r1]["fastq_r2"]),
                    fastqc_before_dir,
                    args.overwrite,
                    "BEFORE"
                ))
        clean_fastqs = find_and_match_fastqs(out_dir)
        for fastq_r1 in sorted(clean_fastqs):
            fastqc_params.append((
                args.fastqc_path,
                Path(clean_fastqs[fastq_r1]["fastq_dir"], fastq_r1),
                fastqc_after_dir,
                args.overwrite,
                "AFTER"
            ))
            if clean_fastqs[fastq_r1]["fastq_r2"] is not None:
                fastqc_params.append((
                    args.fastqc_path,
                    Path(clean_fastqs[fastq_r1]["fastq_dir"], clean_fastqs[fastq_r1]["fastq_r2"]),
                    fastqc_after_dir,
                    args.overwrite,
                    "AFTER"
                ))

        log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
        log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
        log.log(f'{"Files to analyze":>{mar}}: {bold(len(fastqc_params))}')
        log.log("")
        log.log(f'{"BEFORE cleaning":>{mar}}: {bold(fastqc_before_dir)}')
        log.log(f'{"":>{mar}}  {dim(fastqc_before_msg)}')
        log.log(f'{"AFTER cleaning":>{mar}}: {bold(fastqc_after_dir)}')
        log.log(f'{"":>{mar}}  {dim(fastqc_after_msg)}')
        log.log("")
        tqdm_parallel_async_run(fastqc, fastqc_params,
                                "Running FastQC", "FastQC analysis completed", "file",
                                min(threads_max, settings.FASTQC_MAX_INSTANCES), args.show_less)
        log.log("")


    ################################################################################################
    ################################################################################ CLEANUP SECTION
    log.log_section_header("Logs Summarization and File Cleanup")
    log.log_explanation(
        "Now Captus will summarize the logs from 'bbduk.sh' and 'FastQC' to produce a final report"
        " in HTML format. The data tables and individual graphs associated with the HTML report will"
        " be placed in a subdirectory  called '03_qc_extras'."
    )
    qc_extras_dir, qc_extras_msg = make_output_dir(Path(out_dir, "03_qc_extras"))
    log.log("")
    log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
    log.log(f'{"QC extras directory":>{mar}}: {bold(qc_extras_dir)}')
    log.log(f'{"":>{mar}}  {dim(qc_extras_msg)}')
    log.log("")
    bbduk_summarize_logs(out_dir, qc_extras_dir, mar)
    fastqc_summarize_logs(out_dir, qc_extras_dir, fastqc_before_dir, fastqc_after_dir, mar)
    log.log("")
    log.log_explanation(
        "Generating Jupyter Notebook report..."
    )
    qc_html_report, qc_html_msg = execute_jupyter_report(out_dir,
                                                         settings.CAPTUS_QC_REPORT,
                                                         "Captus_Cleaning_Report.ipynb",
                                                         "captus-assembly_clean")
    log.log(f'{"Final HTML QC report":>{mar}}: {bold(qc_html_report)}')
    log.log(f'{"":>{mar}}  {dim(qc_html_msg)}')
    log.log("")

    if not args.keep_all:
        start = time.time()
        log.log("")
        log.log_explanation(
            "Deleting intermediate FASTQ files and other unnecessary files..."
        )
        reclaimed_bytes = 0
        files_to_delete = [fq for fq in adaptors_trimmed_dir.resolve().glob("*")
                           if has_valid_ext(fq, settings.FASTQ_VALID_EXTENSIONS)]
        if not args.skip_fastqc:
            files_to_delete += [html for html in fastqc_before_dir.resolve().glob("*")
                                if has_valid_ext(html, [".html"])]
            files_to_delete += [html for html in fastqc_after_dir.resolve().glob("*")
                                if has_valid_ext(html, [".html"])]
        for del_file in files_to_delete:
            reclaimed_bytes += del_file.stat().st_size
            del_file.unlink()
        log.log(
            f'A total of {len(files_to_delete)} files'
            f' amounting to {reclaimed_bytes / 1024 ** 3:.2f}GB'
            f' were deleted in {elapsed_time(time.time() - start)}'
        )
    else:
        log.log(bold("No files were removed, '--keep_all' was enabled"))
    log.log("")


    ################################################################################################
    ################################################################################# ENDING SECTION
    log.log_section_header(
        "Captus-assembly: Clean -> successfully completed"
        f" [{elapsed_time(time.time() - captus_start)}]"
    )
    log.log("")


def bbduk_trim_adaptors(
        bbduk_path, ram_MB, threads, in_dir, in_fastq, out_dir, adaptor_set, rna, overwrite
):
    """
    Trims adaptors from FASTQ reads using bbduk.sh from BBTools. Runs two simultaneous rounds of
    cleaning, each with slightly different options. Accepts single- or paired-end reads. Also trims
    poly-A tails if the `--rna` flag is enabled.

    Args:
        bbduk_path (str): Path to bbduk.sh
        ram_MB (int): Amount of RAM to use in MB
        threads (int): Number of threads to use
        in_dir (Path): Input files directory
        in_fastq (str): FASTQ input filename
        out_dir (Path): Output files directory
        adaptor_set (str): Path to FASTA file with adaptor sequences to trim from reads
        rna (bool): Trims poly-A tails when enabled
        overwrite (bool): Enable overwriting of previous results
    """

    start = time.time()

    sample_name = "_R".join(in_fastq.split("_R")[:-1])

    fixed = [
        bbduk_path,
        f"-Xmx{ram_MB}m",
        f"threads={threads}",
        "ktrim=r",
        f"minlength={settings.BBDUK_MIN_LENGTH}",
        "interleaved=f"
    ]

    ref = [f"ref={adaptor_set}"]

    if rna is True:
        ref[0] += f",{settings.POLYA_ADAPTORS}"
        fixed += ["trimpolya=4"]

    if "#" in in_fastq:
        fixed += ["trimpairsevenly=t", "trimbyoverlap=t"]

    round_1_stdout_file = Path(out_dir, f"{sample_name}.stdout1.log")
    round_2_stdout_file = Path(out_dir, f"{sample_name}.stdout2.log")

    round_1 = [
        f"in={Path(in_dir, in_fastq)}",
        "out=stdout.fq",
        f"k={settings.BBDUK_ADAPTOR_ROUND1_K}",
        f"mink={settings.BBDUK_ADAPTOR_ROUND1_MINK}",
        f"hdist={settings.BBDUK_ADAPTOR_ROUND1_HDIST}",
        f"stats={Path(out_dir, f'{sample_name}.round1.stats.txt')}",
        f"2>{round_1_stdout_file}",
    ]

    round_2 = [
        "in=stdin.fq",
        f"out={Path(out_dir, in_fastq)}",
        f"k={settings.BBDUK_ADAPTOR_ROUND2_K}",
        f"mink={settings.BBDUK_ADAPTOR_ROUND2_MINK}",
        f"hdist={settings.BBDUK_ADAPTOR_ROUND2_HDIST}",
        f"stats={Path(out_dir, f'{sample_name}.round2.stats.txt')}",
        f"2>{round_2_stdout_file}",
    ]

    bbduk_stderr_file = Path(out_dir, f"{sample_name}.stderr.log")
    bbduk_log1_file = Path(out_dir, f"{sample_name}.round1.log")
    bbduk_log2_file = Path(out_dir, f"{sample_name}.round2.log")
    if overwrite is True or not bbduk_log2_file.exists():
        fixed += ["overwrite=t"]
        bbduk_cmd = fixed + ref + round_1 + ["|"] + fixed + ref + round_2
        with open(bbduk_stderr_file, "w") as bbduk_cmd_stderr:
            subprocess.run(bbduk_cmd, stderr=bbduk_cmd_stderr)
        header_cmd = f"Captus' BBDuk Command for BOTH rounds:\n  {' '.join(bbduk_cmd)}\n\n\n"
        with open(bbduk_log1_file, "w") as bbduk_log1:
            bbduk_log1.write(header_cmd)
            bbduk_log1.write("ROUND 1 LOG:\n")
            with open(round_1_stdout_file, "r") as real_log1_data:
                for line in real_log1_data:
                    if "calcmem.sh" not in line and not line.startswith("java -ea "):
                        bbduk_log1.write(line)
        with open(bbduk_log2_file, "w") as bbduk_log2:
            bbduk_log2.write(header_cmd)
            bbduk_log2.write("ROUND 2 LOG:\n")
            with open(round_2_stdout_file, "r") as real_log2_data:
                for line in real_log2_data:
                    if "calcmem.sh" not in line and not line.startswith("java -ea "):
                        bbduk_log2.write(line)
        round_1_stdout_file.unlink()
        round_2_stdout_file.unlink()
        bbduk_stderr_file.unlink()

        message = f"'{sample_name}': adaptors trimmed [{elapsed_time(time.time() - start)}]"
    else:
        message = dim(f"'{sample_name}': skipped (output files already exist)")

    return message


def bbduk_filter_quality(
        bbduk_path, ram_MB, threads, in_dir, in_fastq, out_dir, trimq, maq, overwrite
):
    """
    Use bbduk.sh from BBTools to quality-filter reads and filter sequencing artifacts and reads
    matching the phiX viral genome
    """

    start = time.time()

    sample_name = "_R".join(in_fastq.split("_R")[:-1])

    bbduk_cmd = [
        bbduk_path,
        f"-Xmx{ram_MB}m",
        f"threads={threads}",
        f"in={Path(in_dir, in_fastq)}",
        f"out={Path(out_dir, in_fastq)}",
        f"ref={settings.PHIX_REF_GENOME},{settings.SEQUENCING_ARTIFACTS}",
        f"k={settings.BBDUK_QUALITY_K}",
        f"hdist={settings.BBDUK_QUALITY_HDIST}",
        f"qtrim={settings.BBDUK_QUALITY_QTRIM}",
        f"trimq={trimq}",
        f"maq={maq}",
        f"minlength={settings.BBDUK_MIN_LENGTH}",
        f"maxns={settings.BBDUK_QUALITY_MAXNS}",
        f"ziplevel={settings.BBDUK_QUALITY_ZIPLEVEL}"
    ]

    stats_log = [
        f"stats={Path(out_dir, f'{sample_name}.cleaning.stats.txt')}",
        f"2>{Path(out_dir, f'{sample_name}.stdout.log')}"
    ]

    bbduk_stderr_file = Path(out_dir, f"{sample_name}.stderr.log")
    bbduk_log_file = Path(out_dir, f"{sample_name}.cleaning.log")
    if overwrite is True or not bbduk_log_file.exists():
        bbduk_cmd += ["overwrite=t"]
        bbduk_cmd += stats_log
        with open(bbduk_stderr_file, "w") as bbduk_cmd_stderr:
            subprocess.run(bbduk_cmd, stderr=bbduk_cmd_stderr)
        with open(bbduk_log_file, "w") as bbduk_log:
            bbduk_log.write(f"Captus' BBDuk Command:\n  {' '.join(bbduk_cmd)}\n\n\n")
            with open(Path(out_dir, f"{sample_name}.stdout.log"), "r") as real_log_data:
                for line in real_log_data:
                    bbduk_log.write(line)
        bbduk_stderr_file.unlink()
        Path(out_dir, f"{sample_name}.stdout.log").unlink()
        message = f"'{sample_name}': quality filtered [{elapsed_time(time.time() - start)}]"
    else:
        message = dim(f"'{sample_name}': skipped (output files already exist)")

    return message


def fastqc(fastqc_path, in_fastq, fastqc_out_dir, overwrite, stage):
    start = time.time()

    fastqc_cmd = [
        fastqc_path,
        "--outdir", f"{fastqc_out_dir}",
        "--nogroup",
        "--threads", "1",
        f"{in_fastq}"
    ]

    if "_R1" in in_fastq.name:
        idx = in_fastq.name.find("_R1") + 3
    elif "_R2" in in_fastq.name:
        idx = in_fastq.name.find("_R2") + 3
    fastqc_log_file = Path(fastqc_out_dir, f"{in_fastq.name[:idx]}.fastqc.log")
    file_name_stage = f"'{in_fastq.name}' ({stage} cleaning)"

    if overwrite is True or not fastqc_log_file.exists():
        with open(fastqc_log_file, "w") as fastqc_log:
            subprocess.call(fastqc_cmd, stdout=fastqc_log, stderr=fastqc_log)
        message = f"{file_name_stage}: FastQC finished [{elapsed_time(time.time() - start)}]"
    else:
        message = dim(f"{file_name_stage}: skipped (FastQC report already exists)")

    return message


def bbduk_summarize_logs(out_dir, qc_extras_dir, margin):
    round1_logs = list(out_dir.rglob("*.round1.log"))
    round2_logs = list(out_dir.rglob("*.round2.log"))
    cleaning_logs = list(out_dir.rglob("*.cleaning.log"))
    round1_stats = list(out_dir.rglob("*.round1.stats.txt"))
    round2_stats = list(out_dir.rglob("*.round2.stats.txt"))
    cleaning_stats = list(out_dir.rglob("*.cleaning.stats.txt"))

    sample_names = [file.name.replace(".round1.log", "") for file in round1_logs]
    sample_names += [file.name.replace(".round2.log", "") for file in round2_logs]
    sample_names += [file.name.replace(".cleaning.log", "") for file in cleaning_logs]
    sample_names += [file.name.replace(".round1.stats.txt", "") for file in round1_stats]
    sample_names += [file.name.replace(".round2.stats.txt", "") for file in round2_stats]
    sample_names += [file.name.replace(".cleaning.stats.txt", "") for file in cleaning_stats]
    sample_names = sorted(list(set(sample_names)))

    stats = {}
    for sample_name in sample_names:
        stats[sample_name] = {
            "reads_input": "NA",
            "bases_input": "NA",
            "reads_trimmed_round1": "NA",
            "reads_passed_round1": "NA",
            "bases_passed_round1": "NA",
            "reads_trimmed_round2": "NA",
            "reads_passed_round2": "NA",
            "bases_passed_round2": "NA",
            "reads_filtered_cleaning": "NA",
            "reads_passed_cleaning": "NA",
            "bases_passed_cleaning": "NA",
        }
    adaptors_round1 = {}
    for file in round1_stats:
        with open(file, "rt") as stat_file:
            for line in stat_file:
                if not line.startswith("#") and bool(line):
                    adaptor_name = " ".join(line.split()[:-2])
                    adaptors_round1[adaptor_name] = {sample: "NA" for sample in sample_names}
                    adaptors_round1[adaptor_name]["total"] = 0
    adaptors_round2 = {}
    for file in round2_stats:
        with open(file, "rt") as stat_file:
            for line in stat_file:
                if not line.startswith("#") and bool(line):
                    adaptor_name = " ".join(line.split()[:-2])
                    adaptors_round2[adaptor_name] = {sample: "NA" for sample in sample_names}
                    adaptors_round2[adaptor_name]["total"] = 0
    contaminants = {}
    for file in cleaning_stats:
        with open(file, "rt") as stat_file:
            for line in stat_file:
                if not line.startswith("#") and bool(line):
                    contaminant_name = " ".join(line.split()[:-2])
                    contaminants[contaminant_name] = {sample: "NA" for sample in sample_names}
                    contaminants[contaminant_name]["total"] = 0

    for file in round1_logs:
        sample_name = file.name.replace(".round1.log", "")
        with open(file, "rt") as log_file:
            for line in log_file:
                fields = line.split()
                if line.startswith("Input:"):
                    stats[sample_name]["reads_input"] = int(fields[1])
                    stats[sample_name]["bases_input"] = int(fields[3])
                if line.startswith("Result:"):
                    stats[sample_name]["reads_passed_round1"] = int(fields[1])
                    stats[sample_name]["bases_passed_round1"] = int(fields[4])
    for file in round2_logs:
        sample_name = file.name.replace(".round2.log", "")
        with open(file, "rt") as log_file:
            for line in log_file:
                fields = line.split()
                if line.startswith("Result:"):
                    stats[sample_name]["reads_passed_round2"] = int(fields[1])
                    stats[sample_name]["bases_passed_round2"] = int(fields[4])
    for file in cleaning_logs:
        sample_name = file.name.replace(".cleaning.log", "")
        with open(file, "rt") as log_file:
            for line in log_file:
                fields = line.split()
                if line.startswith("Result:"):
                    stats[sample_name]["reads_passed_cleaning"] = int(fields[1])
                    stats[sample_name]["bases_passed_cleaning"] = int(fields[4])
    for file in round1_stats:
        sample_name = file.name.replace(".round1.stats.txt", "")
        with open(file, "rt") as stat_file:
            for line in stat_file:
                fields = line.split()
                if line.startswith("#Matched"):
                    stats[sample_name]["reads_trimmed_round1"] = int(fields[-2])
                if not line.startswith("#") and bool(line):
                    adaptor_name = " ".join(line.split()[:-2])
                    adaptors_round1[adaptor_name][sample_name] = int(fields[-2])
                    adaptors_round1[adaptor_name]["total"] += int(fields[-2])
    for file in round2_stats:
        sample_name = file.name.replace(".round2.stats.txt", "")
        with open(file, "rt") as stat_file:
            for line in stat_file:
                fields = line.split()
                if line.startswith("#Matched"):
                    stats[sample_name]["reads_trimmed_round2"] = int(fields[-2])
                if not line.startswith("#") and bool(line):
                    adaptor_name = " ".join(line.split()[:-2])
                    adaptors_round2[adaptor_name][sample_name] = int(fields[-2])
                    adaptors_round2[adaptor_name]["total"] += int(fields[-2])
    for file in cleaning_stats:
        sample_name = file.name.replace(".cleaning.stats.txt", "")
        with open(file, "rt") as stat_file:
            for line in stat_file:
                fields = line.split()
                if line.startswith("#Matched"):
                    stats[sample_name]["reads_filtered_cleaning"] = int(fields[-2])
                if not line.startswith("#") and bool(line):
                    contaminant_name = " ".join(line.split()[:-2])
                    contaminants[contaminant_name][sample_name] = int(fields[-2])
                    contaminants[contaminant_name]["total"] += int(fields[-2])

    reads_bases_summary = Path(qc_extras_dir, "reads_bases.tsv")
    adaptors_round1_summary = Path(qc_extras_dir, "adaptors_round1.tsv")
    adaptors_round2_summary = Path(qc_extras_dir, "adaptors_round2.tsv")
    contaminants_summary = Path(qc_extras_dir, "contaminants.tsv")
    with open(reads_bases_summary, "wt") as summary_out:
        summary_out.write(
            f"sample\treads_input\tbases_input\t"
            f"reads_trimmed_round1\treads_passed_round1\tbases_passed_round1\t"
            f"reads_trimmed_round2\treads_passed_round2\tbases_passed_round2\t"
            f"reads_filtered_cleaning\treads_passed_cleaning\tbases_passed_cleaning\n"
        )
        for sample_name in sorted(stats):
            summary_out.write(
                f'{sample_name}\t'
                f'{stats[sample_name]["reads_input"]}\t'
                f'{stats[sample_name]["bases_input"]}\t'
                f'{stats[sample_name]["reads_trimmed_round1"]}\t'
                f'{stats[sample_name]["reads_passed_round1"]}\t'
                f'{stats[sample_name]["bases_passed_round1"]}\t'
                f'{stats[sample_name]["reads_trimmed_round2"]}\t'
                f'{stats[sample_name]["reads_passed_round2"]}\t'
                f'{stats[sample_name]["bases_passed_round2"]}\t'
                f'{stats[sample_name]["reads_filtered_cleaning"]}\t'
                f'{stats[sample_name]["reads_passed_cleaning"]}\t'
                f'{stats[sample_name]["bases_passed_cleaning"]}\n'
            )
    log.log(f'{"Reads/Bases stats":>{margin}}: {bold(reads_bases_summary)}')
    with open(adaptors_round1_summary, "wt") as summary_out:
        samples_header = "\t".join(sample_names)
        summary_out.write(
            f"adaptor\ttotal_reads\t{samples_header}\n"
        )
        adaptor_totals = {k: adaptors_round1[k]["total"] for k in adaptors_round1}
        for adaptor in {k: v for k, v in sorted(adaptor_totals.items(),
                                                key=lambda item: item[1], reverse=True)}:
            samples_counts = "\t".join([f"{adaptors_round1[adaptor][sample_name]}" for
                                        sample_name in sample_names])
            summary_out.write(
                f"{adaptor}\t{adaptors_round1[adaptor]['total']}\t{samples_counts}\n"
            )
    log.log(f'{"Adaptors 1st round":>{margin}}: {bold(adaptors_round1_summary)}')
    with open(adaptors_round2_summary, "wt") as summary_out:
        samples_header = "\t".join(sample_names)
        summary_out.write(
            f"adaptor\ttotal_reads\t{samples_header}\n"
        )
        adaptor_totals = {k: adaptors_round2[k]["total"] for k in adaptors_round2}
        for adaptor in {k: v for k, v in sorted(adaptor_totals.items(),
                                                key=lambda item: item[1], reverse=True)}:
            samples_counts = "\t".join([f"{adaptors_round2[adaptor][sample_name]}" for
                                        sample_name in sample_names])
            summary_out.write(
                f"{adaptor}\t{adaptors_round2[adaptor]['total']}\t{samples_counts}\n"
            )
    log.log(f'{"Adaptors 2nd round":>{margin}}: {bold(adaptors_round2_summary)}')
    with open(contaminants_summary, "wt") as summary_out:
        samples_header = "\t".join(sample_names)
        summary_out.write(
            f"contaminant\ttotal_reads\t{samples_header}\n"
        )
        contaminant_totals = {k: contaminants[k]["total"] for k in contaminants}
        for contaminant in {k: v for k, v in sorted(contaminant_totals.items(),
                                                    key=lambda item: item[1], reverse=True)}:
            samples_counts = "\t".join([f"{contaminants[contaminant][sample_name]}" for
                                        sample_name in sample_names])
            summary_out.write(
                f"{contaminant}\t{contaminants[contaminant]['total']}\t{samples_counts}\n"
            )
    log.log(f'{"Contaminants":>{margin}}: {bold(contaminants_summary)}')


def fastqc_summarize_logs(out_dir, qc_extras_dir, fastqc_before_dir, fastqc_after_dir, margin):
    """
    Summarize FastQC reports into a nice multisample .html, replacement for NGSReports which is
    cumbersome to install and very slow
    """
    reports = list(Path(out_dir, fastqc_before_dir).rglob("*_fastqc.zip"))
    reports += list(Path(out_dir, fastqc_after_dir).rglob("*_fastqc.zip"))

    fastqc_stats = {
        "per_base_seq_qual_data":    ["\t".join(["sample_name", "read", "stage",
                                                 "base", "mean", "median",
                                                 "upper_quartile", "lower_quartile",
                                                 "percentile_10", "percentile_90"])],
        "per_seq_qual_scores_data":  ["\t".join(["sample_name", "read", "stage",
                                                 "quality", "count"])],
        "per_base_seq_content_data": ["\t".join(["sample_name", "read", "stage",
                                                 "base", "G", "A", "T", "C"])],
        "per_seq_gc_content_data":   ["\t".join(["sample_name", "read", "stage",
                                                 "gc_content", "count"])],
        "seq_len_dist_data":         ["\t".join(["sample_name", "read", "stage",
                                                 "length", "count"])],
        "seq_dup_levels_data":       ["\t".join(["sample_name", "read", "stage",
                                                 "duplication_level",
                                                 "percentage_of_deduplicated",
                                                 "percentage_of_total"])],
        "adaptor_content_data":      ["\t".join(["sample_name", "read", "stage",
                                                 "position", "Illumina_universal_adaptor",
                                                 "Illumina_small_RNA_3'_adaptor",
                                                 "Illumina_small_RNA_5'_adaptor",
                                                 "Nextera_transposase_sequence",
                                                 "SOLID_small_RNA_adaptor"])],
    }

    for file in sorted(reports):
        if "_R1_fastqc.zip" in file.name:
            sample_name = file.name.replace("_R1_fastqc.zip", "")
            read = "R1"
        elif "_R2_fastqc.zip" in file.name:
            sample_name = file.name.replace("_R2_fastqc.zip", "")
            read = "R2"
        else:
            sample_name = file.name.replace("_fastqc.zip", "")
            read = "UNK"
        if "fastqc_before" in str(file):
            stage = "before"
        elif "fastqc_after" in str(file):
            stage = "after"
        else:
            stage = "UNK"
        prepend = [sample_name, read, stage]

        with zipfile.ZipFile(file) as z:
            add_to = ""
            with z.open(str(Path(file.stem, "fastqc_data.txt"))) as data_txt:
                for line in data_txt:
                    line = line.decode()
                    if not line.startswith("#") and not line.startswith(">>END_MODULE"):
                        if line.startswith(">>Basic Statistics"):
                            add_to = ""
                        if line.startswith(">>Per base sequence quality"):
                            add_to = "per_base_seq_qual_data"
                            continue
                        if line.startswith(">>Per tile sequence quality"):
                            add_to = ""
                        if line.startswith(">>Per sequence quality scores"):
                            add_to = "per_seq_qual_scores_data"
                            continue
                        if line.startswith(">>Per base sequence content"):
                            add_to = "per_base_seq_content_data"
                            continue
                        if line.startswith(">>Per sequence GC content"):
                            add_to = "per_seq_gc_content_data"
                            continue
                        if line.startswith(">>Per base N content"):
                            add_to = ""
                        if line.startswith(">>Sequence Length Distribution"):
                            add_to = "seq_len_dist_data"
                            continue
                        if line.startswith(">>Sequence Duplication Levels"):
                            add_to = "seq_dup_levels_data"
                            continue
                        if line.startswith(">>Overrepresented sequences"):
                            add_to = ""
                        if line.startswith(">>Adapter Content"):
                            add_to = "adaptor_content_data"
                            continue
                        if add_to == "":
                            continue
                        else:
                            record = "\t".join(prepend + line.strip("\n").split())
                            fastqc_stats[add_to].append(record)

    per_base_seq_qual = Path(qc_extras_dir, "per_base_seq_qual.tsv")
    with open(per_base_seq_qual, "wt") as stats_out:
        stats_out.write("\n".join(fastqc_stats["per_base_seq_qual_data"]))
    log.log(f'{"Per base seq. qual.":>{margin}}: {bold(per_base_seq_qual)}')

    per_seq_qual_scores = Path(qc_extras_dir, "per_seq_qual_scores.tsv")
    with open(per_seq_qual_scores, "wt") as stats_out:
        stats_out.write("\n".join(fastqc_stats["per_seq_qual_scores_data"]))
    log.log(f'{"Per seq. qual. scores":>{margin}}: {bold(per_seq_qual_scores)}')

    per_base_seq_content = Path(qc_extras_dir, "per_base_seq_content.tsv")
    with open(per_base_seq_content, "wt") as stats_out:
        stats_out.write("\n".join(fastqc_stats["per_base_seq_content_data"]))
    log.log(f'{"Per base seq. content":>{margin}}: {bold(per_base_seq_content)}')

    per_seq_gc_content = Path(qc_extras_dir, "per_seq_gc_content.tsv")
    with open(per_seq_gc_content, "wt") as stats_out:
        stats_out.write("\n".join(fastqc_stats["per_seq_gc_content_data"]))
    log.log(f'{"Per seq. GC content":>{margin}}: {bold(per_seq_gc_content)}')

    seq_len_dist = Path(qc_extras_dir, "seq_len_dist.tsv")
    with open(seq_len_dist, "wt") as stats_out:
        stats_out.write("\n".join(fastqc_stats["seq_len_dist_data"]))
    log.log(f'{"Seq. length distr.":>{margin}}: {bold(seq_len_dist)}')

    seq_dup_levels = Path(qc_extras_dir, "seq_dup_levels.tsv")
    with open(seq_dup_levels, "wt") as stats_out:
        stats_out.write("\n".join(fastqc_stats["seq_dup_levels_data"]))
    log.log(f'{"Seq. duplication":>{margin}}: {bold(seq_dup_levels)}')

    adaptor_content = Path(qc_extras_dir, "adaptor_content.tsv")
    with open(adaptor_content, "wt") as stats_out:
        stats_out.write("\n".join(fastqc_stats["adaptor_content_data"]))
    log.log(f'{"Adaptor content":>{margin}}: {bold(adaptor_content)}')

