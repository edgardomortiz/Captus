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

import gzip
import math
import shutil
import subprocess
import time
from pathlib import Path

from tqdm import tqdm

from . import log, settings
from .bioformats import bait_stats, dict_to_fasta, fasta_to_dict, resolve_iupac
from .misc import (bbtools_path_version, bold, compress_list_files, dim, dir_is_empty, elapsed_time,
                   file_is_empty, format_dep_msg, gzip_compress, has_valid_ext, make_output_dir,
                   pigz_compress, python_library_check, quit_with_error, red, set_ram, set_threads,
                   successful_exit, tqdm_parallel_async_run, tqdm_serial_run, vsearch_path_version)
from .version import __version__


def bait(full_command, args):

    captus_start = time.time()
    out_dir, out_dir_msg = make_output_dir(args.out)
    log.logger = log.Log(Path(out_dir, "captus-design_bait.log"), stdout_verbosity_level=1)

    mar = 24  # Margin for aligning parameters and values

    ################################################################################################
    ############################################################################### STARTING SECTION
    log.log_section_header("Starting Captus-design: BAIT", single_newline=False)
    log.log_explanation(
        "Welcome to the bait creation step of Captus-design. In this step, Captus will process"
        " the clusters automatically selected during the previous step as well as any other"
        " manually selected alignments and derive baits from them. Potential baits will be filtered"
        " and reduced in number through clustering. If the markers were imported from annotated"
        " genomes, the exon data will be used to create baits that do not cross exon boundaries",
        extra_empty_lines_after=0
    )
    log.log_explanation("For more information, please see https://github.com/edgardomortiz/Captus")

    log.log(f'{"Captus version":>{mar}}: {bold(f"v{__version__}")}')
    log.log(f'{"Command":>{mar}}: {bold(full_command)}')
    _, ram_MB, ram_GB, ram_GB_total = set_ram(args.ram)
    log.log(f'{"Max. RAM":>{mar}}: {bold(f"{ram_GB:.1f}GB")} {dim(f"(out of {ram_GB_total:.1f}GB)")}')
    threads_max, threads_total = set_threads(args.threads)
    log.log(f'{"Max. Threads":>{mar}}: {bold(threads_max)} {dim(f"(out of {threads_total})")}')
    log.log("")

    log.log(f'{"Dependencies":>{mar}}:')
    _, bbmap_version, bbmap_status = bbtools_path_version(args.bbmap_path)
    log.log(format_dep_msg(f'{"BBTools":>{mar}}: ', bbmap_version, bbmap_status))
    _, vsearch_version, vsearch_status = vsearch_path_version(args.vsearch_path)
    log.log(format_dep_msg(f'{"VSEARCH":>{mar}}: ', vsearch_version, vsearch_status))
    log.log("")

    log.log(f'{"Python libraries":>{mar}}:')
    _, numpy_version, numpy_status = python_library_check("numpy")
    _, pandas_version, pandas_status = python_library_check("pandas")
    _, plotly_version, plotly_status = python_library_check("plotly")
    log.log(format_dep_msg(f'{"numpy":>{mar}}: ', numpy_version, numpy_status))
    log.log(format_dep_msg(f'{"pandas":>{mar}}: ', pandas_version, pandas_status))
    log.log(format_dep_msg(f'{"plotly":>{mar}}: ', plotly_version, plotly_status))
    log.log("")

    log.log(f'{"Output directory":>{mar}}: {bold(out_dir)}')
    log.log(f'{"":>{mar}}  {dim(out_dir_msg)}')
    log.log("")

    if bbmap_status == "not_found":
        quit_with_error(
            "Captus could not find BBTools' bbmap.sh, please provide a valid path with '--bbmap_path'"
        )
    if vsearch_status == "not_found":
        quit_with_error(
            "Captus could not find VSEARCH, please provide a valid path with '--vsearch_path'"
        )


    ################################################################################################
    ########################################################################## BAIT CREATION SECTION
    log.log_section_header("Creating baits from selected loci alignments")
    log.log_explanation(
        "Now Captus will attempt to import the alignments from the selected clusters directory"
        " (called './02_selected_markers' if default options were used). Both automatically"
        " and manually selected clusters will be processed and all potential baits will be"
        " generated in this step. Finally the baits will be dereplicated."
    )

    log.log(f'{"Bait length":>{mar}}: {bold(args.bait_length)}')
    log.log("")
    log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
    log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
    fastas_auto = find_fastas(Path(args.captus_selected_dir, settings.DES_DIRS["AUT"]))
    log.log(f'{"FASTAs automatic":>{mar}}: {bold(len(fastas_auto))}')
    fastas_manual = find_fastas(Path(args.captus_selected_dir, settings.DES_DIRS["MAN"]))
    log.log(f'{"FASTAs manual":>{mar}}: {bold(len(fastas_manual))}')
    log.log(f'{"FASTAs total":>{mar}}: {bold(len(fastas_auto) + len(fastas_manual))}')
    log.log("")

    raw_baits_dir_path, raw_baits_dir_msg = make_output_dir(Path(out_dir, settings.DES_DIRS["RAW"]))
    log.log(f'{"Raw baits directory":>{mar}}: {bold(raw_baits_dir_path)}')
    log.log(f'{"":>{mar}}  {dim(raw_baits_dir_msg)}')
    log.log("")

    baits_exons_gz_path = Path(raw_baits_dir_path, f'{settings.DES_FILES["BDEX"]}.gz')
    baits_no_exons_gz_path = Path(raw_baits_dir_path, f'{settings.DES_FILES["BDNE"]}.gz')
    long_exons_path = Path(raw_baits_dir_path, settings.DES_FILES["LONG"])

    if baits_exons_gz_path.exists() or baits_no_exons_gz_path.exists():
        log.log(bold("The following output files were already found, nothing to process:"))
        for output_file in [baits_exons_gz_path, baits_no_exons_gz_path, long_exons_path]:
            if output_file.exists():
                log.log(f"'{output_file}'")
            else:
                output_file = None
    else:
        baits_full_exons_path, baits_full_no_exons_path, long_exons_path = create_baits(
            raw_baits_dir_path, Path(args.captus_clusters_dir), args.bait_length,
            fastas_auto, fastas_manual, args.overwrite, args.show_more
        )
        baits_exons_gz_path, baits_no_exons_gz_path = dereplicate_compress_baits(
            raw_baits_dir_path, baits_full_exons_path, baits_full_no_exons_path,
            args.bait_length, args.overwrite, args.vsearch_path, threads_max
        )

    log.log("")


    ################################################################################################
    ######################################################################### BAIT FILTERING SECTION
    log.log_section_header("Filtering baits")
    log.log_explanation(
        "Now Captus will filter the potential bait sequences according to several chemistry"
        " parameters. If you want to exclude potential bait sequences that map to a specific"
        " genomic sequence (e.g. if you want to avoid potential capture of regions of chloroplasts"
        " or mitochondria) you can provide a single FASTA file containing the genomic sequence(s)"
        " to exclude with '--exclude_reference'"
    )

    show_less =  not args.show_more
    log.log(f'{"Concurrent tasks":>{mar}}: {bold(threads_max)}')
    log.log("")
    if args.exclude_reference and not Path(args.exclude_reference).exists():
        quit_with_error(f"The provided FASTA file '{args.exclude_reference}' does not exist")
    log.log(f'{"Exclude if bait maps to":>{mar}}: {bold(args.exclude_reference)}')
    log.log(f'{"Allow Ns in baits":>{mar}}: {bold(args.include_n)}')
    log.log(f'{"Keep IUPAC ambiguities":>{mar}}: {bold(args.keep_iupac)}')
    log.log(f'{"GC content range":>{mar}}: {bold(args.gc_content)}%')
    log.log(f'{"Melting temp. range":>{mar}}: {bold(args.melting_temperature)}˚C')
    log.log(f'{"Hybridization chemistry":>{mar}}: {bold(args.hybridization_chemistry)}')
    log.log(f'{"Sodium concentration":>{mar}}: {bold(args.sodium)}')
    log.log(f'{"Formamide concentration":>{mar}}: {bold(args.formamide)}')
    log.log(f'{"Max. low complexity":>{mar}}: {bold(args.max_masked_percentage)}%')
    log.log(f'{"Max. homopolymer length":>{mar}}: {bold(args.max_homopolymer_length)} bp')
    log.log("")
    log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
    log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
    log.log("")

    fil_baits_dir_path, fil_baits_dir_msg = make_output_dir(Path(out_dir, settings.DES_DIRS["FIL"]))
    log.log(f'{"Filtered baits directory":>{mar}}: {bold(fil_baits_dir_path)}')
    log.log(f'{"":>{mar}}  {dim(fil_baits_dir_msg)}')
    log.log("")

    baits_exons_mapped_gz_path = map_exonic_baits(
        fil_baits_dir_path, baits_exons_gz_path, long_exons_path, args.bbmap_path,
        ram_MB, threads_max, args.bait_length, args.keep_all, args.overwrite
    )

    if baits_exons_mapped_gz_path: baits_exons_gz_path = baits_exons_mapped_gz_path

    baits_concat_mask_gz_path = concat_refex_mask_baits(
        fil_baits_dir_path, baits_exons_gz_path, baits_no_exons_gz_path, args.exclude_reference,
        args.bbmap_path, args.vsearch_path, ram_MB, threads_max, args.bait_length, args.overwrite
    )

    baits_filtered_gz_path = filter_baits(
        fil_baits_dir_path, baits_concat_mask_gz_path, args.include_n, args.keep_iupac,
        args.gc_content, args.melting_temperature, args.hybridization_chemistry, args.sodium,
        args.formamide, args.max_masked_percentage, args.max_homopolymer_length,
        args.bait_length, args.keep_all, args.overwrite, threads_max, args.debug, show_less
    )


    ################################################################################################
    ######################################################################## BAIT CLUSTERING SECTION
    log.log_section_header("Clustering and tiling baits")
    log.log_explanation(
        "Now Captus will cluster the filtered baits at '--clustering_threshold' percent identity"
        " and tiled at '--tiling_percentage_overlap' to create the final set of baits ready to be"
        " sent for synthesis. Increase the clustering threshold to increase the final number of"
        " or reduce it to achieve the opposite result."
    )
    log.log(f'{"Tiling overlap":>{mar}}: {bold(args.tiling_percentage_overlap)}%')
    log.log(f'{"Clustering id.":>{mar}}: {bold(args.clustering_threshold)}%')
    log.log("")
    log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
    log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
    log.log("")

    clu_baits_dir_path, clu_baits_dir_msg = make_output_dir(Path(out_dir, settings.DES_DIRS["BAI"]))
    log.log(f'{"Final baits directory":>{mar}}: {bold(clu_baits_dir_path)}')
    log.log(f'{"":>{mar}}  {dim(clu_baits_dir_msg)}')
    log.log("")

    cluster_tile_baits(
        clu_baits_dir_path, baits_filtered_gz_path, args.bait_length, args.vsearch_path,
        args.tiling_percentage_overlap, args.clustering_threshold, args.overwrite
    )


    ################################################################################################
    ################################################################################# ENDING SECTION
    successful_exit(
        "Captus-design: BAIT -> successfully completed"
        f" [{elapsed_time(time.time() - captus_start)}]"
    )


def find_fastas(fastas_dir: Path):
    valid_exts = settings.FASTA_VALID_EXTENSIONS
    if fastas_dir.exists():
        return [file for file in fastas_dir.resolve().rglob("*")
                if has_valid_ext(file, valid_exts)]
    else:
        return []


def find_fastas_in_dir(cluster_dir_path: Path, extension: str):
    if cluster_dir_path.exists() and not dir_is_empty(cluster_dir_path):
        fastas_to_process = list(Path(cluster_dir_path).rglob(f"*{extension}"))
        return fastas_to_process
    else:
        quit_with_error(
            f"The directory {cluster_dir_path} does not exist or it is empty, no files to process"
        )


def find_and_merge_exon_data(cluster_dir_path: Path):
    tsvs = []
    if cluster_dir_path.exists():
        sample_dirs = list(Path(cluster_dir_path).rglob("*__captus-clr"))
        for sample_dir in sample_dirs:
            tsvs += list(Path(sample_dir).rglob(f'*{settings.DES_SUFFIXES["DATA"]}'))
    exons_data = {}
    for tsv in tsvs:
        with open(tsv, "rt") as tsv_in:
            for line in tsv_in:
                if not line.startswith("cds_id"):
                    record = line.strip().split()
                    try:
                        ple = int(record[5]) / int(record[3])
                    except ZeroDivisionError:
                        ple = 0
                    try:
                        pse = int(record[7]) / int(record[3])
                    except ZeroDivisionError:
                        pse = 0
                    exons_data[record[0]] = {
                        "position": record[1],
                        "exons": int(record[2]),
                        "exons_len": int(record[3]),
                        "long_exons": int(record[4]),
                        "long_exons_len": int(record[5]),
                        "short_exons": int(record[6]),
                        "short_exons_len": int(record[7]),
                        "introns_len": int(record[8]),
                        "gene_len": int(record[9]),
                        "prop_long_exons": ple,
                        "prop_short_exons": pse,
                    }
    return exons_data


def create_baits(
    raw_baits_dir_path: Path, cluster_dir_path: Path, bait_length: int,
    fastas_auto: list, fastas_manual: list, overwrite: bool, show_more: bool
):
    baits_full_exons_path = Path(raw_baits_dir_path, settings.DES_FILES["BFEX"])
    baits_full_no_exons_path = Path(raw_baits_dir_path, settings.DES_FILES["BFNE"])
    long_exons_path = Path(raw_baits_dir_path, settings.DES_FILES["LONG"])
    fastas = fastas_auto + fastas_manual

    if overwrite or not baits_full_exons_path.exists() or not baits_full_no_exons_path.exists():
        exons_data = find_and_merge_exon_data(cluster_dir_path)
        all_cds_ids = []
        log.log(bold(f"Segmenting selected loci alignments into {bait_length} bp baits:"))
        tqdm_cols = min(shutil.get_terminal_size().columns, 120)
        with tqdm(total=len(fastas), ncols=tqdm_cols, unit="loci") as pbar:
            for fasta_path in fastas:
                start = time.time()
                fasta_in = fasta_to_dict(fasta_path)
                cds_ids = []
                bait_count = 1
                baits = {}
                locus, ext = "", ""
                for e in settings.FASTA_VALID_EXTENSIONS:
                    if f"{fasta_path}".lower().endswith(e.lower()):
                        locus = f"{fasta_path.name}".rstrip(e)
                        ext = e
                        break
                for seq_name in fasta_in:
                    sample = seq_name.split(settings.SEQ_NAME_SEP)[0]
                    seq_id = seq_name.split(settings.SEQ_NAME_SEP)[1]
                    if seq_id in exons_data: cds_ids.append(seq_id)
                    seq = fasta_in[seq_name]["sequence"].replace("-","").upper()
                    des = (f'[sample={sample}] [seq_id={seq_id}]'
                        f' {fasta_in[seq_name]["description"]}').strip()
                    for p in range(0, len(seq) - (bait_length - 1)):
                        bait_name = f"{locus}{settings.SEQ_NAME_SEP}b{bait_count:05}"
                        baits[bait_name] = {
                            "sequence": seq[p:p+bait_length],
                            "description": des,
                        }
                        bait_count += 1
                if cds_ids:
                    all_cds_ids += cds_ids
                    dict_to_fasta(baits, baits_full_exons_path, append=True)
                else:
                    dict_to_fasta(baits, baits_full_no_exons_path, append=True)
                msg = f"'{fasta_path.name}': processed [{elapsed_time(time.time() - start)}]"
                if show_more: tqdm.write(msg)
                log.log(msg, print_to_screen=False)
                pbar.update()
        log.log(bold(
            f" \u2514\u2500\u2192 Successfully processed {len(fastas)} loci"
            f" [{elapsed_time(time.time() - start)}]"
        ))
        log.log("")

        if long_exons_path.exists(): long_exons_path.unlink()
        long_exons_fastas = find_fastas_in_dir(cluster_dir_path, settings.DES_SUFFIXES["LONG"])
        if all_cds_ids and long_exons_fastas:
            start = time.time()
            log.log(bold(f"Saving FASTA file including long exons found in selected loci:"))
            tqdm_cols = min(shutil.get_terminal_size().columns, 120)
            with tqdm(total=len(long_exons_fastas), ncols=tqdm_cols, unit="loci") as pbar:
                inner_start = time.time()
                for fasta_path in long_exons_fastas:
                    long_exons = fasta_to_dict(fasta_path)
                    lef = {}
                    for seq_name in long_exons:
                        cds_id = seq_name.split("_exon")[0]
                        if cds_id in all_cds_ids:
                            lef[seq_name] = long_exons[seq_name]
                    dict_to_fasta(lef, long_exons_path, append=True)
                msg = f"'{long_exons_path}' saved in {elapsed_time(time.time() - inner_start)}"
                if show_more: tqdm.write(msg)
                log.log(msg, print_to_screen=False)
                pbar.update()
        log.log(bold(
            f" \u2514\u2500\u2192 FASTA file '{long_exons_path.name}' succesfully created"
            f" [{elapsed_time(time.time() - start)}]"
        ))
        log.log("")
        if not long_exons_path.exists() or file_is_empty(long_exons_path):
            long_exons_path = None

    else:
        log.log(bold("The following output files were already found:"))
        for output_file in [baits_full_exons_path, baits_full_no_exons_path, long_exons_path]:
            if output_file.exists(): log.log(f"'{output_file}'")
        log.log("")

    if not baits_full_exons_path.exists() or file_is_empty(baits_full_exons_path):
        baits_full_exons_path = None
    if not baits_full_no_exons_path.exists() or file_is_empty(baits_full_no_exons_path):
        baits_full_no_exons_path = None

    return baits_full_exons_path, baits_full_no_exons_path, long_exons_path


def dereplicate_compress_baits(
    raw_baits_dir_path: Path, baits_full_exons_path: Path, baits_full_no_exons_path: Path,
    bait_length: int, overwrite: bool, vsearch_path: str, threads_max: int
):
    start = time.time()

    baits_exons_path = Path(raw_baits_dir_path, settings.DES_FILES["BDEX"])
    baits_no_exons_path = Path(raw_baits_dir_path, settings.DES_FILES["BDNE"])
    baits_exons_gz_path = Path(raw_baits_dir_path, f'{settings.DES_FILES["BDEX"]}.gz')
    baits_no_exons_gz_path = Path(raw_baits_dir_path, f'{settings.DES_FILES["BDNE"]}.gz')
    baits_to_compress = []

    if overwrite or (not baits_exons_gz_path.exists() and not baits_no_exons_gz_path.exists()):
        log.log(bold(f"Dereplicating bait sequences:"))
        if baits_exons_path.exists(): baits_exons_path.unlink()
        if baits_no_exons_path.exists(): baits_no_exons_path.unlink()
        tqdm_cols = min(shutil.get_terminal_size().columns, 120)
        with tqdm(total=2, ncols=tqdm_cols, unit="file") as pbar:
            for bait_file in [baits_full_exons_path, baits_full_no_exons_path]:
                if bait_file is None:
                    pbar.update()
                else:
                    inner_start = time.time()
                    bait_derep_file = Path(f"{bait_file}".replace("_full", ""))
                    derep_log_file = Path(f"{bait_file}".replace(".fasta", ".derep.log"))
                    derep_cmd = [
                        vsearch_path,
                        "--derep_fulllength", f"{bait_file}",
                        "--output", f"{bait_derep_file}",
                        "--fasta_width", f"{bait_length}",
                        "--strand", "both",
                        "--notrunclabels",
                    ]
                    with open(derep_log_file, "w") as derep_log:
                        derep_log.write(f"Captus' Dereplication Command:\n  {' '.join(derep_cmd)}\n\n\n")
                    with open(derep_log_file, "a") as derep_log:
                        subprocess.run(derep_cmd, stdout=derep_log, stdin=derep_log, stderr=derep_log)
                    if not file_is_empty(bait_derep_file):
                        baits_to_compress.append(bait_derep_file)
                        bait_file.unlink()
                        msg = (f"'{bait_file.name}': dereplicated, saved as "
                               f"'{bait_derep_file.name}' [{elapsed_time(time.time() - inner_start)}]")
                    else:
                        bait_derep_file.unlink()
                        msg = red(f"'{bait_file.name}': FAILED dereplication")
                    tqdm.write(msg)
                    log.log(msg, print_to_screen=False)
                    pbar.update()
        log.log(bold(
            f" \u2514\u2500\u2192 Bait files successfully dereplicated"
            f" [{elapsed_time(time.time() - start)}]"
        ))
        log.log("")
    if baits_to_compress:
        compress_list_files(baits_to_compress, threads_max)
        for bait_gz_file in [baits_exons_gz_path, baits_no_exons_gz_path]:
            if not bait_gz_file.exists():
                bait_gz_file = None
    else:
        quit_with_error("The dereplication of the bait files failed, double check your command!")

    return baits_exons_gz_path, baits_no_exons_gz_path


def map_exonic_baits(
    filtered_baits_dir_path: Path, baits_exons_gz_path: Path, long_exons_path: Path, bbmap_path: str,
    ram_MB: int, threads_max: int, bait_length: int, keep_all: bool, overwrite: bool
):
    start = time.time()

    if not baits_exons_gz_path or not long_exons_path: return None
    baits_exons_mapped_path = Path(filtered_baits_dir_path, f'{settings.DES_FILES["BEXM"]}')
    baits_exons_mapped_gz_path = Path(f"{baits_exons_mapped_path}.gz")

    if overwrite or not baits_exons_mapped_gz_path.exists():
        if baits_exons_mapped_path.exists(): baits_exons_mapped_path.unlink()
        if baits_exons_mapped_gz_path.exists(): baits_exons_mapped_gz_path.unlink()
        baits_sam_gz_path = Path(f"{baits_exons_mapped_path}".replace(".fasta", ".sam.gz"))
        baits_sam_log_path = Path(f"{baits_exons_mapped_path}".replace(".fasta", ".bbmap.log"))
        baits_sam_stdout_path = Path(f"{baits_exons_mapped_path}".replace(".fasta", ".stdout.log"))
        baits_sam_stderr_path = Path(f"{baits_exons_mapped_path}".replace(".fasta", ".stderr.log"))
        log.log(bold(f"Mapping baits to long exon sequences:"))
        tqdm_cols = min(shutil.get_terminal_size().columns, 120)
        with tqdm(total=2, ncols=tqdm_cols, unit="task") as pbar:

            inner_start = time.time()
            bbmap_cmd = [
                bbmap_path,
                f"-Xmx{ram_MB}m",
                "nodisk=t",
                f"ref={long_exons_path}",
                f"in={baits_exons_gz_path}",
                "local=t",
                f"threads={threads_max}",
                f"outm={baits_sam_gz_path}",
                "sam=1.3",
                f"fastawrap={bait_length}",
                f"2>{baits_sam_stdout_path}"
            ]
            with open(baits_sam_stderr_path, "w") as bbmap_err:
                subprocess.run(bbmap_cmd, stderr=bbmap_err)
            with open(baits_sam_log_path, "w") as bbmap_log:
                bbmap_log.write(f"Captus' BBMap Command:\n  {' '.join(bbmap_cmd)}\n\n\n")
                with open(baits_sam_stdout_path, "r") as real_log_data:
                    for line in real_log_data:
                        bbmap_log.write(line)
            baits_sam_stderr_path.unlink()
            baits_sam_stdout_path.unlink()
            msg = (f"'{baits_sam_gz_path.name}': baits mapped to exons"
                   f" [{elapsed_time(time.time() - inner_start)}]")
            tqdm.write(msg)
            log.log(msg, print_to_screen=False)
            pbar.update()

            inner_start = time.time()
            total_baits = 0
            long_exon_baits = 0
            with gzip.open(baits_sam_gz_path, "rt") as sam:
                while 1:
                    sam_chunk = sam.readlines(10000)
                    if not sam_chunk: break
                    baits_fasta = {}
                    for line in sam_chunk:
                        if not line.startswith("@"):
                            total_baits += 1
                            record = line.split("\t")
                            if not "S" in record[5]:
                                long_exon_baits += 1
                                baits_fasta[record[0].split(" ")[0]] = {
                                    "sequence": record[9],
                                    "description": " ".join(record[0].split(" ")[1:]),
                                }
                    dict_to_fasta(baits_fasta, baits_exons_mapped_path, append=True)
            if shutil.which("pigz"):
                pigz_compress(baits_exons_mapped_path, threads_max)
            elif shutil.which("gzip"):
                gzip_compress(baits_exons_mapped_path)
            if not keep_all: baits_sam_gz_path.unlink()
            if baits_exons_mapped_gz_path.exists():
                pct_retained = long_exon_baits / total_baits * 100
                msg = (f"'{baits_exons_mapped_gz_path.name}': {long_exon_baits:n} of {total_baits:n}"
                       f" baits ({pct_retained:.2f}%) fully contained in long exons"
                       f" [{elapsed_time(time.time() - inner_start)}]")
            else:
                msg = red(f"'{baits_exons_mapped_gz_path.name}': FAILED")
            tqdm.write(msg)
            log.log(msg, print_to_screen=False)
            pbar.update()

        log.log(bold(
            f" \u2514\u2500\u2192 Exonic baits processed, only baits fully contained"
            f" in long exons retained [{elapsed_time(time.time() - start)}]"
        ))
        log.log("")
    else:
        log.log(f"'{bold(baits_exons_mapped_gz_path.name)}': output already"
                f" exists, SKIPPED mapping baits to long exons sequences")
        log.log("")

    if baits_exons_mapped_gz_path.exists():
        return baits_exons_mapped_gz_path
    else:
        return None


def concat_refex_mask_baits(
    filtered_baits_dir_path: Path, baits_exons_gz_path: Path, baits_no_exons_gz_path: Path,
    exclude_reference: str, bbmap_path: str, vsearch_path: str, ram_MB: int, threads_max: int,
    bait_length: int, overwrite: bool
):
    start = time.time()

    if not baits_exons_gz_path.exists() and not baits_no_exons_gz_path.exists():
        quit_with_error("No bait files were found, nothing else to process...")
    baits_concat_path = Path(filtered_baits_dir_path, f'{settings.DES_FILES["BCAT"]}')
    baits_mask_path = Path(f"{baits_concat_path}".replace(".fasta", "_mask.fasta"))
    baits_mask_gz_path = Path(f"{baits_concat_path}".replace(".fasta", "_mask.fasta.gz"))

    if overwrite or not baits_mask_gz_path.exists():
        if baits_mask_path.exists(): baits_mask_path.unlink()
        if baits_mask_gz_path.exists(): baits_mask_gz_path.unlink()
        baits_tomap_path = Path(f"{baits_concat_path}".replace(".fasta", "_tomap.fasta.gz"))
        baits_unmap_path = Path(f"{baits_concat_path}".replace(".fasta", "_unmap.fasta.gz"))
        baits_unmap_log_path = Path(f"{baits_concat_path}".replace(".fasta", "_unmap.bbmap.log"))
        baits_unmap_stdout_path = Path(f"{baits_concat_path}".replace(".fasta", "_unmap.stdout.log"))
        baits_unmap_stderr_path = Path(f"{baits_concat_path}".replace(".fasta", "_unmap.stderr.log"))
        baits_mask_log_path = Path(f"{baits_concat_path}".replace(".fasta", "_mask.vsearch.log"))

        log.log(bold(f"Excluding baits matching a given reference and masking low-complexity:"))
        tqdm_cols = min(shutil.get_terminal_size().columns, 120)
        with tqdm(total=3, ncols=tqdm_cols, unit="task") as pbar:

            inner_start = time.time()
            cat_cmd = ["cat"]
            if baits_exons_gz_path.exists(): cat_cmd += [f"{baits_exons_gz_path}"]
            if baits_no_exons_gz_path.exists(): cat_cmd += [f"{baits_no_exons_gz_path}"]
            with open(baits_tomap_path, "w") as baits_tomap:
                subprocess.run(cat_cmd, stdout=baits_tomap)
            if baits_tomap_path.exists():
                msg = (f"'{baits_tomap_path}': baits concatenated"
                       f" [{elapsed_time(time.time() - inner_start)}]")
            else:
                quit_with_error("No bait files were found, nothing else to process...")
            tqdm.write(msg)
            log.log(msg, print_to_screen=False)
            pbar.update()

            inner_start = time.time()
            if not exclude_reference or not Path(exclude_reference).exists():
                msg = (f"'{baits_tomap_path}': SKIPPED mapping against excluding"
                       f" reference [{elapsed_time(time.time() - inner_start)}]")
                baits_tomap_path.replace(baits_unmap_path)
            else:
                bbmap_cmd = [
                    bbmap_path,
                    f"-Xmx{ram_MB}m",
                    "nodisk=t",
                    f"ref={exclude_reference}",
                    f"in={baits_tomap_path}",
                    "local=t",
                    f"threads={threads_max}",
                    f"outu={baits_unmap_path}",
                    f"fastawrap={bait_length}",
                    f"2>{baits_unmap_stdout_path}"
                ]
                with open(baits_unmap_stderr_path, "w") as bbmap_err:
                    subprocess.run(bbmap_cmd, stderr=bbmap_err)
                with open(baits_unmap_log_path, "w") as bbmap_log:
                    bbmap_log.write(f"Captus' BBMap Command:\n  {' '.join(bbmap_cmd)}\n\n\n")
                    with open(baits_unmap_stdout_path, "r") as real_log_data:
                        for line in real_log_data:
                            bbmap_log.write(line)
                baits_unmap_stderr_path.unlink()
                baits_unmap_stdout_path.unlink()
                msg = (f"'{baits_unmap_path.name}': retained baits that do not match the"
                       f" excluding reference [{elapsed_time(time.time() - inner_start)}]")
            tqdm.write(msg)
            log.log(msg, print_to_screen=False)
            pbar.update()

            inner_start = time.time()
            mask_cmd = [
                vsearch_path,
                "--maskfasta", f"{baits_unmap_path}",
                "--output", f"{baits_mask_path}",
                "--qmask", "dust",
                "--fasta_width", f"{bait_length}",
                "--notrunclabels",
            ]
            with open(baits_mask_log_path, "w") as mask_log:
                mask_log.write(f"Captus' Low Complexity Masking Command:\n  {' '.join(mask_cmd)}\n\n\n")
            with open(baits_mask_log_path, "a") as mask_log:
                subprocess.run(mask_cmd, stdout=mask_log, stdin=mask_log, stderr=mask_log)
            if not file_is_empty(baits_mask_path):
                if shutil.which("pigz"):
                    pigz_compress(baits_mask_path, threads_max)
                elif shutil.which("gzip"):
                    gzip_compress(baits_mask_path)
            if baits_mask_gz_path.exists():
                msg = (f"'{baits_mask_gz_path.name}': masked low-complexity regions"
                       f" in baits [{elapsed_time(time.time() - inner_start)}]")
            else:
                quit_with_error("No bait files were found, nothing else to process...")
            tqdm.write(msg)
            log.log(msg, print_to_screen=False)
            pbar.update()

        log.log(bold(
            f" \u2514\u2500\u2192 Baits concatenated, filtered against excluding"
            f" reference and masked [{elapsed_time(time.time() - start)}]"
        ))
        log.log("")
    else:
        log.log(f"'{bold(baits_mask_gz_path.name)}': output already exists, SKIPPED"
                f" concatenation, filtering by reference, and masking of baits")
        log.log("")

    if baits_mask_gz_path.exists():
        return baits_mask_gz_path
    else:
        return None


def split_baits_file(filtered_baits_dir_path: Path, baits_path: Path, threads_max: int):
    baits = fasta_to_dict(baits_path)
    num_seqs = len(baits)
    chunks_floor = max(threads_max, math.floor(num_seqs / settings.BAITS_SPLIT_SIZE))
    chunks_ceiling = max(threads_max, math.ceil(num_seqs / settings.BAITS_SPLIT_SIZE))
    try:
        if (abs(settings.BAITS_SPLIT_SIZE - (num_seqs / chunks_ceiling))
            <= abs(settings.BAITS_SPLIT_SIZE - (num_seqs / chunks_floor))):
            seqs_per_chunk = math.ceil(num_seqs / chunks_ceiling)
        else:
            seqs_per_chunk = math.ceil(num_seqs / chunks_floor)
    except ZeroDivisionError:
        seqs_per_chunk = num_seqs
    bait_parts_paths = []
    part = 1
    split_bait = {}
    seq_count = 0
    for seq_name in baits:
        split_bait[seq_name] = baits[seq_name]
        seq_count += 1
        if seq_count < seqs_per_chunk:
            continue
        else:
            split_bait_path = Path(filtered_baits_dir_path, f"{baits_path.stem}_part{part}.fasta.gz")
            dict_to_fasta(split_bait, split_bait_path)
            bait_parts_paths.append(split_bait_path)
            part += 1
            split_bait = {}
            seq_count = 0
    return bait_parts_paths


def filter_baits_chunk(
    bait_part_path: Path, include_n: bool, keep_iupac: bool, gc_content: str,
    melting_temperature: str, hybridization_chemistry: str, sodium: float, formamide: float,
    max_masked_percentage: float, max_homopolymer_length: int, bait_length: int
):
    start = time.time()

    bait_part_accepted = {}
    bait_part_accepted_path = Path(f"{bait_part_path}".replace(".fasta.gz", "_accepted.fasta.gz"))
    bait_part_rejected = {}
    bait_part_rejected_path = Path(f"{bait_part_path}".replace(".fasta.gz", "_rejected.fasta.gz"))

    max_Ns = bait_length if include_n else 0
    min_gc, max_gc = (float(gc) for gc in gc_content.split(","))
    min_melt_temp, max_melt_temp = (float(mt) for mt in melting_temperature.split(","))
    sodium = max(0.00001, sodium)

    bait_part = fasta_to_dict(bait_part_path)
    for bait_name in bait_part:
        bs = bait_stats(bait_part[bait_name]["sequence"], hybridization_chemistry, sodium, formamide)
        bait_desc = (f'{bait_part[bait_name]["description"]} [gc={bs["gc"]:.2f}%]'
                     f' [melt={bs["melt_temp"]:.2f}C] [low_complex={bs["low_complexity"]:.2f}%]'
                     f' [homopol={bs["max_homopolymer_len"]}bp] [N={bs["Ns"]}]')
        if (min_gc <= bs["gc"] <= max_gc
            and min_melt_temp <= bs["melt_temp"] <= max_melt_temp
            and bs["low_complexity"] <= max_masked_percentage
            and bs["max_homopolymer_len"] <= max_homopolymer_length
            and bs["Ns"] <= max_Ns):
            bait_seq = bait_part[bait_name]["sequence"]
            if keep_iupac: bait_seq = resolve_iupac(bait_seq)
            bait_part_accepted[bait_name] = {
                "sequence": bait_seq,
                "description": bait_desc,
            }
        else:
            bait_part_rejected[bait_name] = {
                "sequence": bait_part[bait_name]["sequence"],
                "description": bait_desc,
            }
    dict_to_fasta(bait_part_accepted, bait_part_accepted_path)
    dict_to_fasta(bait_part_rejected, bait_part_rejected_path)

    return f"'{bait_part_accepted_path.name}': filtered [{elapsed_time(time.time() - start)}]"


def filter_baits(
    filtered_baits_dir_path: Path, baits_concat_mask_gz_path: Path, include_n: bool,
    keep_iupac: bool, gc_content: str, melting_temperature: str, hybridization_chemistry: str,
    sodium: float, formamide: float, max_masked_percentage: float, max_homopolymer_length: int,
    bait_length: int, keep_all: bool, overwrite: bool, threads_max: int, debug: bool,
    show_less: bool
):

    if not baits_concat_mask_gz_path.exists():
        quit_with_error("No bait files were found, nothing else to process...")
    baits_accepted_path = Path(filtered_baits_dir_path, f'{settings.DES_FILES["BACC"]}')
    baits_accepted_gz_path = Path(f"{baits_accepted_path}.gz")
    baits_rejected_path = Path(filtered_baits_dir_path, f'{settings.DES_FILES["BREJ"]}')
    baits_rejected_gz_path = Path(f"{baits_rejected_path}.gz")

    if overwrite or not baits_accepted_gz_path.exists():
        if baits_accepted_path.exists(): baits_accepted_path.unlink()
        if baits_accepted_gz_path.exists(): baits_accepted_gz_path.unlink()
        baits_chunks_paths = split_baits_file(filtered_baits_dir_path,
                                              baits_concat_mask_gz_path,
                                              threads_max)
        filter_params = []
        for baits_chunk_path in baits_chunks_paths:
            filter_params.append((
                baits_chunk_path,
                include_n,
                keep_iupac,
                gc_content,
                melting_temperature,
                hybridization_chemistry,
                sodium,
                formamide,
                max_masked_percentage,
                max_homopolymer_length,
                bait_length
            ))
        if debug:
            tqdm_serial_run(filter_baits_chunk, filter_params,
                            f"Filtering baits", f"Baits filtering completed",
                            "file", show_less)
        else:
            tqdm_parallel_async_run(filter_baits_chunk, filter_params,
                                    f"Filtering baits", f"Baits filtering completed",
                                    "file", threads_max, show_less)
        log.log("")

        if not keep_all:
            for baits_chunk_path in baits_chunks_paths:
                if baits_chunk_path.exists(): baits_chunk_path.unlink()

        start = time.time()
        rejected_bait_paths = list(filtered_baits_dir_path.rglob("*_rejected.fasta.gz"))
        for baits_chunk_path in rejected_bait_paths:
            rejected_baits = fasta_to_dict(baits_chunk_path)
            baits_chunk_path.unlink()
            dict_to_fasta(rejected_baits, baits_rejected_path, append=True)
            if shutil.which("pigz"):
                pigz_compress(baits_rejected_path, threads_max)
            elif shutil.which("gzip"):
                gzip_compress(baits_rejected_path)
        log.log(f"'{bold(baits_rejected_gz_path.name)}': concatenated [{elapsed_time(time.time() - start)}]")
        log.log("")

        accepted_bait_paths = list(filtered_baits_dir_path.rglob("*_accepted.fasta.gz"))
        baits_accepted_unsort_path = Path(f"{baits_accepted_path}".replace(".fasta", "_unsort.fasta"))
        for baits_chunk_path in accepted_bait_paths:
            accepted_baits = fasta_to_dict(baits_chunk_path)
            baits_chunk_path.unlink()
            dict_to_fasta(accepted_baits, baits_accepted_unsort_path, append=True)
        baits_accepted = fasta_to_dict(baits_accepted_unsort_path)
        baits_accepted_unsort_path.unlink()
        dict_to_fasta(baits_accepted, baits_accepted_path, sort=True)
        if shutil.which("pigz"):
            pigz_compress(baits_accepted_path, threads_max)
        elif shutil.which("gzip"):
            gzip_compress(baits_accepted_path)
        log.log(f"'{bold(baits_accepted_gz_path.name)}': concatenated [{elapsed_time(time.time() - start)}]")
        log.log("")
    else:
        log.log(f"'{bold(baits_accepted_gz_path.name)}': output"
                f" already exists, SKIPPED bait filtering")
        log.log("")

    if baits_accepted_gz_path.exists():
        return baits_accepted_gz_path
    else:
        return None


def cluster_tile_baits(
    clust_baits_dir_path: Path, baits_filtered_gz_path: Path, bait_length: int, vsearch_path: str,
    tiling_percentage_overlap:float, clustering_threshold: float, overwrite: bool
):

    mincols = math.ceil(tiling_percentage_overlap / 100 * bait_length)
    clust_baits_unsorted_file = f"baits_clust{clustering_threshold:.2f}_mincols{mincols}_unsorted.fasta"
    clust_baits_unsorted_path = Path(clust_baits_dir_path, clust_baits_unsorted_file)
    clust_baits_final_file = f"baits_clust{clustering_threshold:.2f}_mincols{mincols}.fasta"
    clust_baits_final_path = Path(clust_baits_dir_path, clust_baits_final_file)
    clust_baits_log_path = Path(f"{clust_baits_final_path}".replace(".fasta", ".log"))
    clust_baits_tsv_path = Path(f"{clust_baits_final_path}".replace(".fasta", ".tsv"))
    if clustering_threshold > 1.0: clustering_threshold /= 100

    if overwrite or not clust_baits_final_path.exists():
        if clust_baits_final_path.exists(): clust_baits_final_path.unlink()
        if clust_baits_log_path.exists(): clust_baits_log_path.unlink()
        if clust_baits_tsv_path.exists(): clust_baits_tsv_path.unlink()
        log.log(bold(f"Clustering and tiling baits:"))
        tqdm_cols = min(shutil.get_terminal_size().columns, 120)
        with tqdm(total=1, ncols=tqdm_cols, unit="task") as pbar:
            inner_start = time.time()
            clust_cmd = [
                vsearch_path,
                "--cluster_smallmem", f"{baits_filtered_gz_path}",
                "--usersort",
                "--strand", "both",
                "--id", f"{clustering_threshold}",
                "--mincols", f"{mincols}",
                "--fasta_width", f"{bait_length}",
                "--notrunclabels",
                "--centroids", f"{clust_baits_unsorted_path}"
            ]
            with open(clust_baits_log_path, "w") as clust_log:
                clust_log.write(f"Captus' Bait Clustering Command:\n  {' '.join(clust_cmd)}\n\n\n")
            with open(clust_baits_log_path, "a") as clust_log:
                subprocess.run(clust_cmd, stdout=clust_log, stdin=clust_log, stderr=clust_log)
            if clust_baits_unsorted_path.exists() and not file_is_empty(clust_baits_unsorted_path):
                msg = (f"'{clust_baits_final_path.name}': baits clustered"
                       f" [{elapsed_time(time.time() - inner_start)}]")
            else:
                quit_with_error("No bait files were found, nothing else to process...")
            tqdm.write(msg)
            log.log(msg, print_to_screen=False)
            pbar.update()
        log.log("")
        loci = {}
        baits = fasta_to_dict(clust_baits_unsorted_path)
        num_baits = len(baits)
        clust_baits_unsorted_path.unlink()
        dict_to_fasta(baits, clust_baits_final_path, sort=True)
        for bait_name in baits:
            locus = bait_name.split(settings.SEQ_NAME_SEP)[0]
            if locus in loci:
                loci[locus] += 1
            else:
                loci[locus] = 1
        log.log(f"'{bold(clust_baits_final_path.name)}': {num_baits} baits covering {len(loci)} loci")
        log.log("")
        with open(clust_baits_tsv_path, "wt") as tsv_out:
            tsv_out.write(f"locus\tnum_baits\n")
            for locus in sorted(loci):
                tsv_out.write(f"{locus}\t{loci[locus]}\n")
        log.log(f"'{bold(clust_baits_tsv_path.name)}': saved table with number of baits per locus")
        log.log("")
    else:
        log.log(f"'{bold(clust_baits_final_path.name)}': output"
                f" already exists, SKIPPED bait clustering")
        log.log("")

    return
