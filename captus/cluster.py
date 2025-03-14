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
import shutil
import subprocess
import time
from multiprocessing import Manager
from pathlib import Path

from tqdm import tqdm

from . import log, settings
from .bioformats import (
    alignment_stats,
    cds_from_gff,
    dict_to_fasta,
    fasta_to_dict,
    mmseqs_cluster,
    vsearch_cluster,
)
from .misc import (
    bold,
    dim,
    dir_is_empty,
    elapsed_time,
    file_is_empty,
    find_and_match_fastas_gffs,
    format_dep_msg,
    mafft_path_version,
    make_output_dir,
    mmseqs_path_version,
    python_library_check,
    quit_with_error,
    red,
    set_ram,
    set_threads,
    successful_exit,
    tqdm_parallel_async_run,
    tqdm_serial_run,
    vsearch_path_version,
)
from .version import __version__


def cluster(full_command, args):
    captus_start = time.time()
    out_dir, out_dir_msg = make_output_dir(args.out)
    log.logger = log.Log(Path(out_dir, "captus-cluster.log"), stdout_verbosity_level=1)

    mar = 23  # Margin for aligning parameters and values

    ################################################################################################
    ############################################################################### STARTING SECTION
    log.log_section_header("Starting Captus-design: CLUSTER", single_newline=False)
    log.log_explanation(
        "Welcome to the marker clustering step of Captus-design. In this step, Captus will process"
        " the files provided with the '--markers_to_cluster' option. If a FASTA file is found"
        " together with a matching GFF annotation file, Captus will extract the CDS for clustering,"
        " otherwise every sequence in the FASTA file that is shorter than '--max_seq_len' will be"
        " clustered",
        extra_empty_lines_after=0,
    )
    if len(args.markers_to_cluster) == 1 and Path(args.markers_to_cluster[0]).is_dir():
        intro_msg = (
            "Since you provided a directory name, Captus will look in that location for all the"
            " FASTA files and their matching GFF annotation files whenever posssible. To match"
            " a FASTA file with its GFF annotation file, the names should be identical except for"
            " their extension. Sample names are derived from the text before the extension. Use '_'"
            " to separate parts of the sample name, first part will be taken as Genus and second as"
            " Species, if multiple samples of the same species are used, add a identifier as third"
            " part to distinguish them."
        )
    else:
        intro_msg = (
            "Since you provided a file name or a list of file names, Captus will only process those"
            " FASTA files and their matching GFF annotation files whenever posssible. To match"
            " a FASTA file with its GFF annotation file, the names should be identical except for"
            " their extension. Sample names are derived from the text before the extension. Use '_'"
            " to separate parts of the sample name, first part will be taken as Genus and second as"
            " Species, if multiple samples of the same species are used, add a identifier as third"
            " part to distinguish them."
        )
    log.log_explanation(intro_msg, extra_empty_lines_after=0)

    intro_msg = ""
    no_clust_program = False
    _, mmseqs_version, mmseqs_status = mmseqs_path_version(args.mmseqs_path)
    _, vsearch_version, vsearch_status = vsearch_path_version(args.vsearch_path)
    if args.clust_program.lower() == "mmseqs" and mmseqs_status == "not found":
        if vsearch_status == "OK":
            args.clust_program = "vsearch"
            intro_msg = (
                "MMseqs2 was not found, but VSEARCH will be used for dereplicating and"
                " clustering instead."
            )
        else:
            no_clust_program = True
    elif args.clust_program.lower() == "vsearch" and vsearch_status == "not found":
        if mmseqs_status == "OK":
            args.clust_program = "mmseqs"
            intro_msg = (
                "VSEARCH was not found, but MMseqs2 will be used for dereplicating and"
                " clustering instead."
            )
        else:
            no_clust_program = True
    show_less = not args.show_more

    log.log_explanation(intro_msg, extra_empty_lines_after=0)
    log.log_explanation("For more information, please see https://github.com/edgardomortiz/Captus")

    max_seq_len = adjust_max_seq_len(args.clust_program, args.max_seq_len)

    log.log(f"{'Captus version':>{mar}}: {bold(f'v{__version__}')}")
    log.log(f"{'Command':>{mar}}: {bold(full_command)}")
    _, ram_MB, ram_GB, ram_GB_total = set_ram(args.ram)
    log.log(f"{'Max. RAM':>{mar}}: {bold(f'{ram_GB:.1f}GB')} {dim(f'(out of {ram_GB_total:.1f}GB)')}")
    threads_max, threads_total = set_threads(args.threads)
    log.log(f"{'Max. Threads':>{mar}}: {bold(threads_max)} {dim(f'(out of {threads_total})')}")
    log.log("")

    log.log(f"{'Dependencies':>{mar}}:")
    if args.clust_program.lower() == "mmseqs":
        log.log(format_dep_msg(f"{'MMseqs2':>{mar}}: ", mmseqs_version, mmseqs_status))
        log.log(format_dep_msg(f"{'VSEARCH':>{mar}}: ", "", "not used"))
    if args.clust_program.lower() == "vsearch":
        log.log(format_dep_msg(f"{'MMseqs2':>{mar}}: ", "", "not used"))
        log.log(format_dep_msg(f"{'VSEARCH':>{mar}}: ", vsearch_version, vsearch_status))
    _, mafft_version, mafft_status = mafft_path_version(args.mafft_path)
    log.log(format_dep_msg(f"{'MAFFT':>{mar}}: ", mafft_version, mafft_status))
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

    if args.redo_from:
        log.log(f"{'Redo from':>{mar}}: {bold(args.redo_from)}")
        prepare_redo(out_dir, args.redo_from)
        log.log("")

    if no_clust_program:
        quit_with_error(
            "Captus could not find a program for clustering, please provide a valid path with either"
            " '--mmseqs_path' or '--vsearch_path'"
        )
    if mafft_status == "not found":
        quit_with_error("Captus could not find MAFFT, please provide a valid path with '--mafft_path'")

    ################################################################################################
    ########################################################################## MARKER IMPORT SECTION
    log.log_section_header("Processing and import of markers for clustering")
    import_msg = (
        "Now Captus will import the markers to be clustered. It will look for valid FASTA files in"
        " the location provided with '--markers_to_cluster'. Additionally it will try to find"
        " matching GFF annotation files, names must be identical except for the extension. For"
        " example 'Poa_annua.fasta' and 'Poa_annua.gff.gz' will be taken as a valid FASTA file with"
        " a matching GFF file."
    )
    log.log_explanation(import_msg)

    log.log(f"{'Concurrent imports':>{mar}}: {bold(threads_max)}")
    log.log("")
    log.log(f"{'Bait length':>{mar}}: {bold(args.bait_length)}")
    log.log("")
    log.log(f"{'Overwrite files':>{mar}}: {bold(args.overwrite)}")
    log.log(f"{'Keep all files':>{mar}}: {bold(args.keep_all)}")
    fastas_to_import = find_and_match_fastas_gffs(args.markers_to_cluster)
    log.log(f"{'FASTAs to import':>{mar}}: {bold(len(fastas_to_import))}")
    log.log("")

    import_params = []
    for fasta in fastas_to_import:
        import_params.append(
            (
                fasta,
                fastas_to_import[fasta]["fasta_dir"],
                fastas_to_import[fasta]["gff_path"],
                args.bait_length,
                out_dir,
                args.overwrite,
            )
        )

    if args.debug:
        tqdm_serial_run(
            import_fasta,
            import_params,
            "Importing and preprocessing FASTAs and GFFs",
            "Data import completed",
            "sample",
            show_less,
        )
    else:
        tqdm_parallel_async_run(
            import_fasta,
            import_params,
            "Importing and preprocessing FASTAs and GFFs",
            "Data import completed",
            "sample",
            threads_max,
            show_less,
        )
    log.log("")

    ################################################################################################
    #################################################### MARKER DEDUPLICATION AND CLUSTERING SECTION
    log.log_section_header("Deduplication and clustering of markers")
    clust_msg = (
        "Now Captus will exclude sequences longer than '--max_seq_len', then the sequences will be"
        " deduplicated sample-wise at the identity threshold given by '--dedup_threshold'. Finally,"
        " markers will be clustered across samples using '--clust_threshold' as identity threshold."
    )
    log.log_explanation(clust_msg)

    log.log(f"{'Max. threads':>{mar}}: {bold(threads_max)}")
    log.log("")
    log.log(f"{'Clustering program':>{mar}}: {bold(args.clust_program)}")
    if args.clust_program == "mmseqs":
        log.log(f"{'MMseqs2 sensitivity':>{mar}}: {bold(args.mmseqs_sensitivity)}")
        log.log(f"{'MMseqs2 cluster mode':>{mar}}: {bold(args.mmseqs_cluster_mode)}")
    log.log(f"{'Max. sequence length':>{mar}}: {bold(max_seq_len)}")
    log.log(f"{'Min. sequence length':>{mar}}: {bold(args.bait_length)} (= bait length)")
    log.log(f"{'Strand':>{mar}}: {bold(args.strand)}")
    log.log(f"{'Deduplication id.':>{mar}}: {bold(args.dedup_threshold)}%")
    log.log(f"{'Clustering id.':>{mar}}: {bold(args.clust_threshold)}%")
    log.log(f"{'Align singletons':>{mar}}: {bold(args.align_singletons)}")
    log.log("")
    log.log(f"{'Overwrite files':>{mar}}: {bold(args.overwrite)}")
    log.log(f"{'Keep all files':>{mar}}: {bold(args.keep_all)}")
    samples_to_filter = find_fastas_in_sample_dirs(out_dir, settings.DES_SUFFIXES["MARKERS"])
    log.log(f"{'FASTAs to process':>{mar}}: {bold(len(samples_to_filter))}")
    log.log("")

    filter_params = []
    for sample_name in samples_to_filter:
        filter_params.append(
            (
                sample_name,
                samples_to_filter[sample_name]["sample_dir"],
                samples_to_filter[sample_name]["fasta_path"],
                out_dir,
                max_seq_len,
                args.bait_length,
                args.overwrite,
            )
        )
    task_title = (
        f"Removing sequences shorter than {args.bait_length}"
        f" bp or longer than {max_seq_len} bp from FASTAs"
    )
    if args.debug:
        tqdm_serial_run(
            filter_fasta,
            filter_params,
            task_title,
            "Length filtering completed",
            "sample",
            show_less,
        )
    else:
        tqdm_parallel_async_run(
            filter_fasta,
            filter_params,
            task_title,
            "Length filtering completed",
            "sample",
            threads_max,
            show_less,
        )

    log.log("")

    samples_to_dedup = find_fastas_in_sample_dirs(out_dir, settings.DES_SUFFIXES["FILTER"])
    dedup_params = []
    for sample_name in samples_to_dedup:
        dedup_params.append(
            (
                sample_name,
                samples_to_dedup[sample_name]["sample_dir"],
                samples_to_dedup[sample_name]["fasta_path"],
                out_dir,
                args.clust_program,
                args.mmseqs_path,
                args.mmseqs_sensitivity,
                args.mmseqs_cluster_mode,
                args.vsearch_path,
                args.dedup_threshold,
                args.strand,
                threads_max,
                args.overwrite,
                args.keep_all,
            )
        )
    tqdm_serial_run(
        dedup_fasta,
        dedup_params,
        f"Deduplicating sequences from each sample at {args.dedup_threshold}% identity",
        "Deduplication completed",
        "sample",
        show_less,
    )
    log.log("")

    concat_dir_path, concat_dir_msg = make_output_dir(Path(out_dir, settings.DES_DIRS["CAT"]))
    log.log(f"{'Concatenation directory':>{mar}}: {bold(concat_dir_path)}")
    log.log(f"{'':>{mar}}  {dim(concat_dir_msg)}")
    log.log("")
    samples_to_concat = find_fastas_in_sample_dirs(out_dir, settings.DES_SUFFIXES["DEDUPED"])
    fasta_concat_path = rehead_and_concatenate_fastas(
        samples_to_concat, concat_dir_path, args.overwrite, show_less
    )
    log.log("")

    cluster_dir_path, cluster_dir_msg = make_output_dir(Path(out_dir, settings.DES_DIRS["CLR"]))
    log.log(f"{'Clustering directory':>{mar}}: {bold(cluster_dir_path)}")
    log.log(f"{'':>{mar}}  {dim(cluster_dir_msg)}")
    log.log("")
    cluster_markers(
        fasta_concat_path,
        cluster_dir_path,
        args.clust_program,
        args.mmseqs_path,
        args.mmseqs_sensitivity,
        args.mmseqs_cluster_mode,
        args.vsearch_path,
        args.strand,
        args.clust_threshold,
        args.align_singletons,
        threads_max,
        args.overwrite,
    )
    log.log("")

    ################################################################################################
    ######################################################################### LOCI ALIGNMENT SECTION
    log.log_section_header("Alignment of clustered markers")
    align_msg = (
        "Now Captus will align the clusters using MAFFT, the algorithm used is 'genafpair' a.k.a"
        " 'E-INS-i' but it can be changed in 'settings.py' inside Captus' installation directory."
    )
    log.log_explanation(align_msg)

    fastas_to_align = find_fastas_in_dir(cluster_dir_path, ".fasta")
    concurrent, threads_per_alignment = adjust_align_concurrency(args.concurrent, threads_max)
    align_dir_path, align_dir_msg = make_output_dir(Path(out_dir, settings.DES_DIRS["ALN"]))
    align_params = []
    for fasta_path in fastas_to_align:
        align_params.append(
            (
                args.mafft_path,
                settings.DESIGN_ALIGN_ALGORITHM,
                threads_per_alignment,
                args.timeout,
                fasta_path,
                align_dir_path,
                args.overwrite,
            )
        )

    log.log(f"{'Concurrent alignments':>{mar}}: {bold(concurrent)}")
    log.log(f"{'Threads per alignment':>{mar}}: {bold(threads_per_alignment)}")
    log.log("")
    log.log(
        f"{'Algorithm':>{mar}}: {bold(settings.DESIGN_ALIGN_ALGORITHM)}"
        f" {dim(settings.ALIGN_ALGORITHMS[settings.DESIGN_ALIGN_ALGORITHM]['aka'])}"
    )
    log.log(f"{'Timeout':>{mar}}: {bold(args.timeout)} {dim(f'[{elapsed_time(args.timeout)}]')}")
    log.log("")
    log.log(f"{'Overwrite files':>{mar}}: {bold(args.overwrite)}")
    log.log(f"{'Keep all files':>{mar}}: {bold(args.keep_all)}")
    log.log(f"{'FASTAs to align':>{mar}}: {bold(len(fastas_to_align))}")
    log.log("")
    log.log(f"{'Alignment directory':>{mar}}: {bold(align_dir_path)}")
    log.log(f"{'':>{mar}}  {dim(align_dir_msg)}")
    log.log("")

    if args.debug:
        tqdm_serial_run(
            mafft_assembly,
            align_params,
            "Aligning clusters",
            "Cluster alignment completed",
            "alignment",
            show_less,
        )
    else:
        tqdm_parallel_async_run(
            mafft_assembly,
            align_params,
            "Aligning clusters",
            "Cluster alignment completed",
            "alignment",
            concurrent,
            show_less,
        )
    log.log("")

    ################################################################################################
    ##################################################################### ALIGNMENT CURATION SECTION
    log.log_section_header("Curation of aligned clusters")
    curate_msg = (
        "Now Captus will curate the alignments of the marker clusters. First, it will try to stitch"
        " together partial sequences from a single sample when these sequences do not overlap or"
        " the overlapping section is identical. Then it will remove leading and trailing columns"
        " with a proportion of missing data greater than '--gaps', and annotate the presence of"
        " '--focal_species', '--outgroup_species', and '--addon_samples' in each alignment."
    )
    log.log_explanation(curate_msg)

    fastas_to_curate = find_fastas_in_dir(align_dir_path, ".fasta")
    args.focal_species = format_sample_list(args.focal_species)
    args.outgroup_species = format_sample_list(args.outgroup_species)
    args.addon_samples = format_sample_list(args.addon_samples, list_format="samples")
    exons_data = find_and_merge_exon_data(out_dir)
    curate_dir_path, curate_dir_msg = make_output_dir(Path(out_dir, settings.DES_DIRS["CUR"]))
    manager = Manager()
    shared_aln_stats = manager.list()
    curate_params = []
    for fasta_path in fastas_to_curate:
        curate_params.append(
            (
                args.gaps,
                args.bait_length,
                args.focal_species,
                args.outgroup_species,
                args.addon_samples,
                exons_data,
                fasta_path,
                curate_dir_path,
                shared_aln_stats,
                args.overwrite,
            )
        )

    log.log(f"{'Concurrent curations':>{mar}}: {bold(threads_max)}")
    log.log("")
    log.log(f"{'Max. gap proportion':>{mar}}: {bold(args.gaps)}")
    log.log(f"{'Bait length':>{mar}}: {bold(args.bait_length)}")
    log.log(f"{'Focal species':>{mar}}: {bold(args.focal_species)}")
    log.log(f"{'Outgroup species':>{mar}}: {bold(args.outgroup_species)}")
    log.log(f"{'Add-on samples':>{mar}}: {bold(args.addon_samples)}")
    log.log("")
    log.log(f"{'Overwrite files':>{mar}}: {bold(args.overwrite)}")
    log.log(f"{'Keep all files':>{mar}}: {bold(args.keep_all)}")
    log.log(f"{'FASTAs to curate':>{mar}}: {bold(len(fastas_to_curate))}")
    log.log("")
    log.log(f"{'Curation directory':>{mar}}: {bold(curate_dir_path)}")
    log.log(f"{'':>{mar}}  {dim(curate_dir_msg)}")
    log.log("")

    tqdm_serial_run(
        curate,
        curate_params,
        "Curating cluster alignments",
        "Curation of alignments completed",
        "alignment",
        show_less,
    )

    # manager.list() slows down the curate() function when used in parallel, until better solution
    # we will run it serially
    # if args.debug:
    #     tqdm_serial_run(curate, curate_params,
    #                     "Curating cluster alignments", "Curation of alignments completed",
    #                     "alignment", show_less)
    # else:
    #     tqdm_parallel_async_run(curate, curate_params,
    #                             "Curating cluster alignments", "Curation of alignments completed",
    #                             "alignment", threads_max, show_less)

    log.log("")

    ################################################################################################
    ################################################################################ CLEANUP SECTION
    log.log_section_header("Statistics Summarization and File Cleanup")
    log.log_explanation(
        "Captus will calculate and then summarize alignment statistics into an HTML report."
        " Additionally, the sequence-to-sample equivalence file needed for ASTRAL-Pro will be"
        " generated. Finally, empty directories will be removed as well as MAFFT and ClipKIT logs"
        " unless the flag '--keep_all' is enabled."
    )

    start = time.time()
    log.log_explanation("Writing curated alignments statistics...")
    aln_stats_tsv = write_aln_stats(out_dir, shared_aln_stats)
    if aln_stats_tsv:
        log.log(f"{'Alignment statistics':>{mar}}: {bold(aln_stats_tsv)}")
        log.log(f"{'':>{mar}}  {dim(f'File saved in {elapsed_time(time.time() - start)}')}")
        log.log("")
        if all([numpy_found, pandas_found, plotly_found]):
            from .report import build_design_report

            log.log("")
            log.log_explanation("Generating Alignment Statistics report...")
            aln_html_report, aln_html_msg = build_design_report(out_dir, aln_stats_tsv, "cluster")
            log.log(f"{'Alignment report':>{mar}}: {bold(aln_html_report)}")
            log.log(f"{'':>{mar}}  {dim(aln_html_msg)}")
        else:
            log.log(
                f"{bold('WARNING:')} Captus uses 'numpy', 'pandas', and 'plotly' to generate an HTML"
                " report based on the marker recovery statistics. At least one of these libraries"
                " could not be found, please verify these libraries are installed and available."
            )
    else:
        log.log(red("Skipping summarization step... (no alignment statistics files were produced)"))
    log.log("")

    if not args.keep_all:
        start = time.time()
        log.log("")
        log.log_explanation("Deleting MAFFT's logs and other unnecessary files...")
        reclaimed_bytes = 0
        files_to_delete = list(out_dir.resolve().rglob("*.mafft.log"))
        for del_file in files_to_delete:
            reclaimed_bytes += del_file.stat().st_size
            del_file.unlink()
        log.log(
            f"    A total of {len(files_to_delete)} files"
            f" amounting to {reclaimed_bytes / 1024**2:.2f}MB"
            f" were deleted in {elapsed_time(time.time() - start)}"
        )
    else:
        log.log(bold("No files were removed, '--keep_all' was enabled"))

    log.log("")

    ################################################################################################
    ################################################################################# ENDING SECTION
    successful_exit(
        f"Captus-design: CLUSTER -> successfully completed [{elapsed_time(time.time() - captus_start)}]"
    )


def prepare_redo(out_dir: Path, redo_from: str):
    importation = list(Path(out_dir).rglob("*__captus-clr"))
    dereplication = list(Path(out_dir).rglob(f"*{settings.DES_SUFFIXES['DEDUPED']}"))
    concatenation = [Path(out_dir, settings.DES_DIRS["CAT"])]
    clustering = [Path(out_dir, settings.DES_DIRS["CLR"])]
    alignment = [Path(out_dir, settings.DES_DIRS["ALN"])]
    curation = [Path(out_dir, settings.DES_DIRS["CUR"])]
    dirs_to_delete = []
    if redo_from == "importation":
        dirs_to_delete = importation + concatenation + clustering + alignment + curation
    elif redo_from == "dereplication":
        dirs_to_delete = concatenation + clustering + alignment + curation
    elif redo_from == "clustering":
        dirs_to_delete = concatenation + clustering + alignment + curation
    elif redo_from == "alignment":
        dirs_to_delete = alignment + curation
    elif redo_from == "curation":
        dirs_to_delete = curation

    start = time.time()
    log.log("")
    log.log(bold(f"Deleting directories and files to redo from the '{redo_from}' stage:"))
    tqdm_cols = min(shutil.get_terminal_size().columns, 120)
    items_to_delete = len(dirs_to_delete)
    if redo_from == "dereplication":
        items_to_delete += len(dereplication)
    with tqdm(total=items_to_delete, ncols=tqdm_cols, unit="items") as pbar:
        if redo_from == "dereplication":
            for del_file in dereplication:
                del_file.unlink()
                tqdm.write(f"'{del_file}': deleted")
                pbar.update()
        for del_dir in dirs_to_delete:
            shutil.rmtree(del_dir, ignore_errors=True)
            tqdm.write(f"'{del_dir}': deleted")
            pbar.update()
    log.log(
        bold(
            f" \u2514\u2500\u2192 Ready for redo from the '{redo_from}' stage, {items_to_delete}"
            f" directories and files deleted [{elapsed_time(time.time() - start)}]"
        )
    )


def import_fasta(
    fasta_name: str,
    fasta_parent: Path,
    gff_path: Path,
    bait_length: int,
    out_dir: Path,
    overwrite: bool,
):
    start = time.time()

    fasta_path = Path(fasta_parent, fasta_name)
    sample_name = fasta_path.name.replace("".join(fasta_path.suffixes), "")
    sample_out_dir = Path(out_dir, f"{sample_name}__captus-clr")

    if overwrite is True or not sample_out_dir.exists():
        sample_out_dir, _ = make_output_dir(Path(out_dir, f"{sample_name}__captus-clr"))
        if gff_path:
            cds, long, short, data = cds_from_gff(gff_path, fasta_path, bait_length)
            cds_path = Path(sample_out_dir, f"{sample_name}{settings.DES_SUFFIXES['CDS']}")
            long_path = Path(sample_out_dir, f"{sample_name}{settings.DES_SUFFIXES['LONG']}")
            short_path = Path(sample_out_dir, f"{sample_name}{settings.DES_SUFFIXES['SHORT']}")
            data_path = Path(sample_out_dir, f"{sample_name}{settings.DES_SUFFIXES['DATA']}")
            dict_to_fasta(cds, cds_path, write_if_empty=True)
            dict_to_fasta(long, long_path, write_if_empty=True)
            dict_to_fasta(short, short_path, write_if_empty=True)
            with open(data_path, "wt") as data_out:
                for line in data:
                    data_out.write("\t".join(line) + "\n")
            message = (
                f"'{sample_name}': CDS and exon data imported [{elapsed_time(time.time() - start)}]"
            )
        else:
            markers = fasta_to_dict(fasta_path)
            markers_path = Path(sample_out_dir, f"{sample_name}{settings.DES_SUFFIXES['MARKERS']}")
            dict_to_fasta(markers, markers_path, write_if_empty=True)
            message = f"'{sample_name}': FASTA imported [{elapsed_time(time.time() - start)}]"
    else:
        message = dim(f"'{sample_name}': SKIPPED import (output files already exist)")

    return message


def find_fastas_in_sample_dirs(out_dir: Path, suffix: str):
    fastas_to_process = {}
    if out_dir.exists():
        sample_dirs = list(Path(out_dir).rglob("*__captus-clr"))
        for sample_dir in sample_dirs:
            fastas = list(Path(sample_dir).rglob(f"*{suffix}"))
            for fasta in fastas:
                if f"{fasta.parent}".endswith("__captus-clr"):
                    sample_dir = fasta.parent.parts[-1]
                    sample_name = sample_dir.replace("__captus-clr", "")
                    fastas_to_process[sample_name] = {
                        "sample_dir": sample_dir,
                        "fasta_path": fasta,
                    }

    return fastas_to_process


def adjust_max_seq_len(clust_program: str, max_seq_len: str):
    if max_seq_len == "auto":
        if clust_program == "vsearch":
            return settings.MAX_SEQ_LEN["vsearch"]
        elif clust_program == "mmseqs":
            return settings.MAX_SEQ_LEN["mmseqs"]
    else:
        try:
            return int(max_seq_len)
        except ValueError:
            quit_with_error("Please enter a valid number for '--max_seq_len'")


def filter_fasta(
    sample_name: str,
    sample_dir: str,
    fasta_path: Path,
    out_dir: Path,
    max_seq_len: int,
    min_seq_len: int,
    overwrite: bool,
):
    start = time.time()
    fasta_out_path = Path(out_dir, sample_dir, f"{sample_name}{settings.DES_SUFFIXES['FILTER']}")

    if overwrite is True or not fasta_out_path.exists():
        fasta_in = fasta_to_dict(fasta_path)
        fasta_out = {}
        for seq_name in fasta_in:
            if min_seq_len <= len(fasta_in[seq_name]["sequence"]) <= max_seq_len:
                fasta_out[seq_name] = fasta_in[seq_name]
        dict_to_fasta(fasta_out, fasta_out_path)
        if fasta_out_path.exists():
            message = f"'{sample_name}': FASTA filtered by length [{elapsed_time(time.time() - start)}]"
        else:
            message = red(f"'{sample_name}': FAILED filtering by length")
    else:
        message = dim(f"'{sample_name}': SKIPPED filtering by length (output file already exist)")

    return message


def dedup_fasta(
    sample_name: str,
    sample_dir: str,
    fasta_path: Path,
    out_dir: Path,
    clust_program: str,
    mmseqs_path: str,
    mmseqs_sensitivity: float,
    mmseqs_cluster_mode: int,
    vsearch_path: str,
    dedup_threshold: float,
    strand: str,
    threads_max: int,
    overwrite: bool,
    keep_all: bool,
):
    start = time.time()
    fasta_out_path = Path(out_dir, sample_dir, f"{sample_name}{settings.DES_SUFFIXES['DEDUPED']}")
    if dedup_threshold > 1:
        dedup_threshold = dedup_threshold / 100

    if overwrite is True or not fasta_out_path.exists():
        result_prefix = f"{fasta_out_path}".replace(".fasta", "")
        tmp_dir = Path(out_dir, sample_dir, "mmseqs_tmp")
        if clust_program == "vsearch":
            dedup_cmd = [
                vsearch_path,
                "--cluster_fast",
                f"{fasta_path}",
                "--strand",
                f"{strand}",
                "--iddef",
                "0",
                "--id",
                f"{dedup_threshold}",
                "--query_cov",
                f"{min(0.99, dedup_threshold)}",
                "--notrunclabels",
                "--centroids",
                f"{fasta_out_path}",
                "--threads",
                f"{threads_max}",
            ]
        elif clust_program == "mmseqs":
            dedup_cmd = [
                mmseqs_path,
                "easy-cluster",
                f"{fasta_path}",
                f"{result_prefix}",
                f"{tmp_dir}",
                "-s",
                f"{mmseqs_sensitivity}",
                "--min-seq-id",
                f"{dedup_threshold}",
                "--seq-id-mode",
                "1",
                "-c",
                f"{min(0.99, dedup_threshold)}",
                "--cov-mode",
                "1",
                "--cluster-mode",
                f"{mmseqs_cluster_mode}",
                "--gap-open",
                f"{max(1, settings.MMSEQS_GAP_OPEN)}",
                "--gap-extend",
                f"{max(1, settings.MMSEQS_GAP_EXTEND)}",
                "--kmer-per-seq-scale",
                f"{settings.MMSEQS_KMER_PER_SEQ_SCALE}",
                "--threads",
                f"{threads_max}",
            ]
            if mmseqs_cluster_mode != 0:
                dedup_cmd += ["--cluster-reassign"]
        dedup_log_file = Path(f"{fasta_out_path}".replace(".fasta", ".log"))
        with open(dedup_log_file, "w") as dedup_log:
            dedup_log.write(f"Captus' Deduplication Command:\n  {' '.join(dedup_cmd)}\n\n\n")
        with open(dedup_log_file, "a") as dedup_log:
            subprocess.run(dedup_cmd, stdout=dedup_log, stdin=dedup_log, stderr=dedup_log)
        rep_seq_path = Path(f"{result_prefix}_rep_seq.fasta")
        if not keep_all:
            to_delete = [
                Path(f"{result_prefix}_all_seqs.fasta"),
                Path(f"{result_prefix}_cluster.tsv"),
            ]
            for file in to_delete:
                if file.exists():
                    file.unlink()
            if tmp_dir.exists():
                shutil.rmtree(tmp_dir, ignore_errors=True)
        if rep_seq_path.exists():
            rep_seq_path.replace(fasta_out_path)
        if file_is_empty(fasta_out_path):
            message = red(f"'{sample_name}': FAILED deduplication")
        else:
            message = f"'{sample_name}': FASTA deduplicated [{elapsed_time(time.time() - start)}]"
    else:
        message = dim(f"'{sample_name}': SKIPPED deduplication (output file already exist)")

    return message


def rehead_and_concatenate_fastas(
    samples_to_concat: dict, concat_dir: Path, overwrite: bool, show_less: bool
):
    """
    Since all the FASTAs coming from all the samples will be joined in a single file for clustering,
    we have to include the sample names in the headers to be able to distinguish them later.
    The sample's name + '__' + original contig name will be used as new headers, so try to avoid
    the use of '__' to name the samples. The descriptions will be lost (MMseqs ignores them afaik).
    """
    start = time.time()
    log.log(bold(f"Concatenating FASTAs from {len(samples_to_concat)} samples for clustering:"))
    fasta_concat_path = Path(concat_dir, "concatenated.fasta")
    if overwrite is True or not fasta_concat_path.exists():
        tqdm_cols = min(shutil.get_terminal_size().columns, 120)
        with tqdm(total=len(samples_to_concat), ncols=tqdm_cols, unit="file") as pbar:
            for sample in samples_to_concat:
                fasta_rehead = {}
                fasta_sample = fasta_to_dict(samples_to_concat[sample]["fasta_path"])
                for seq_name in fasta_sample:
                    fasta_rehead[f"{sample}{settings.SEQ_NAME_SEP}{seq_name}"] = fasta_sample[seq_name]
                dict_to_fasta(fasta_rehead, fasta_concat_path, append=True)
                message = f"'{sample}': FASTA reheaded [{elapsed_time(time.time() - start)}]"
                log.log(message, print_to_screen=False)
                if not show_less:
                    tqdm.write(message)
                pbar.update()
        log.log(
            bold(
                f" \u2514\u2500\u2192 File '{fasta_concat_path.name}'"
                f" prepared in {elapsed_time(time.time() - start)}(s)"
            )
        )
    else:
        log.log(dim(f"SKIPPED concatenation, '{fasta_concat_path.name}' already exists"))
    return fasta_concat_path


def cluster_markers(
    fasta_concat_path: Path,
    cluster_dir_path: Path,
    clust_program: str,
    mmseqs_path: str,
    mmseqs_sensitivity: float,
    mmseqs_cluster_mode: int,
    vsearch_path: str,
    strand: str,
    clust_threshold: float,
    align_singletons: bool,
    threads: int,
    overwrite: bool,
):
    log.log(bold(f"Clustering '{fasta_concat_path.name}' at {clust_threshold}% identity:"))
    if clust_threshold > 1:
        clust_threshold = clust_threshold / 100
    cluster_prefix = "concatenated"
    cluster_tsv_path = Path(fasta_concat_path.parent, f"{cluster_prefix}_cluster.tsv")
    if overwrite is True or not cluster_tsv_path.exists() or file_is_empty(cluster_tsv_path):
        if clust_program == "mmseqs":
            mmseqs_tmp_dir = Path(fasta_concat_path.parent, "mmseqs_tmp")
            message = mmseqs_cluster(
                mmseqs_path,
                "easy-cluster",
                Path(fasta_concat_path.parent),
                fasta_concat_path,
                cluster_prefix,
                mmseqs_tmp_dir,
                mmseqs_sensitivity,
                clust_threshold,
                1,
                clust_threshold,
                1,
                mmseqs_cluster_mode,
                threads,
            )
        elif clust_program == "vsearch":
            message = vsearch_cluster(
                vsearch_path,
                "--cluster_fast",
                Path(fasta_concat_path.parent),
                fasta_concat_path,
                cluster_prefix,
                clust_threshold,
                0,
                clust_threshold,
                strand,
                threads,
            )
        log.log(message)
    else:
        log.log(dim(f"SKIPPED clustering, '{cluster_tsv_path.name}' already exists and is not empty"))

    start = time.time()
    log.log("")
    log.log(bold("Writing clusters in FASTA format:"))
    clusters_all = {}
    clusters_pass = {}
    if file_is_empty(cluster_tsv_path) or not cluster_tsv_path.exists():
        quit_with_error(f"Clustering FAILED, the file {cluster_tsv_path} was empty or not found")
    else:
        with open(cluster_tsv_path) as cluster_tsv:
            for line in cluster_tsv:
                record = line.strip().split()
                centroid = record[0]
                member = record[1]
                if centroid not in clusters_all:
                    clusters_all[centroid] = [centroid]
                    if member not in clusters_all[centroid]:
                        clusters_all[centroid] += [member]
                else:
                    clusters_all[centroid] += [member]
        single_sample_clusters = 0
        single_species_clusters = 0
        for centroid in clusters_all:
            samples = []
            species = []
            for member in clusters_all[centroid]:
                sample_name = member.split(settings.SEQ_NAME_SEP)[0]
                if sample_name not in samples:
                    samples.append(sample_name)
                species_name = "_".join(sample_name.split("_")[0:2])
                if species_name not in species:
                    species.append(species_name)
            if len(samples) > 1:
                if len(species) > 1:
                    clusters_pass[centroid] = clusters_all[centroid]
                else:
                    single_species_clusters += 1
                    if align_singletons:
                        clusters_pass[centroid] = clusters_all[centroid]
            else:
                single_sample_clusters += 1
    if overwrite is True or dir_is_empty(cluster_dir_path):
        fasta_concat = fasta_to_dict(fasta_concat_path)
        tqdm_cols = min(shutil.get_terminal_size().columns, 120)
        with tqdm(total=len(clusters_pass), ncols=tqdm_cols, unit="cluster") as pbar:
            cluster_count = 1
            num_digits = len(str(len(clusters_pass)))
            for centroid in sorted(clusters_pass):
                fasta_out = {}
                for member in clusters_pass[centroid]:
                    fasta_out[member] = fasta_concat[member]
                fasta_out_path = Path(cluster_dir_path, f"locus{cluster_count:0{num_digits}}.fasta")
                dict_to_fasta(fasta_out, fasta_out_path)
                cluster_count += 1
                pbar.update()
        if align_singletons:
            message = (
                f"Clusters written, {single_species_clusters} clusters contained a single"
                f" species out of a total of {len(clusters_pass)} valid clusters"
                f" [{elapsed_time(time.time() - start)}]"
            )
        else:
            message = (
                f"Clusters written, {len(clusters_pass)} total valid clusters"
                f" [{elapsed_time(time.time() - start)}]"
            )
        log.log(bold(f" \u2514\u2500\u2192 {message}"))
    else:
        log.log(dim(f"SKIPPED writing of FASTA clusters, '{cluster_dir_path}' is not empty"))

    return


def find_fastas_in_dir(cluster_dir_path: Path, extension: str):
    if cluster_dir_path.exists() and not dir_is_empty(cluster_dir_path):
        fastas_to_process = list(Path(cluster_dir_path).rglob(f"*{extension}"))
        return fastas_to_process
    else:
        quit_with_error(
            f"The directory {cluster_dir_path} does not exist or it is empty, no files to process"
        )


def adjust_align_concurrency(concurrent, threads_max):
    if concurrent == "auto":
        concurrent = max(threads_max // 2, 1)
    else:
        try:
            concurrent = int(concurrent)
        except ValueError:
            quit_with_error("Invalid value for '--concurrent', set it to 'auto' or use a number")
    if concurrent > threads_max:
        concurrent = min(threads_max, concurrent)
    threads_per_alignment = threads_max // concurrent

    return concurrent, threads_per_alignment


def mafft_assembly(
    mafft_path: Path,
    mafft_algorithm: str,
    threads_per_alignment: int,
    timeout: int,
    input_fasta_path: Path,
    align_dir_path: Path,
    overwrite: bool,
):
    def cleanup_fastas(all_fastas: list):
        for file in all_fastas:
            if file.exists():
                file.unlink()
        return

    start = time.time()

    output_fasta_path = Path(align_dir_path, input_fasta_path.name)

    if overwrite is True or not output_fasta_path.exists():
        long_fasta_path = Path(align_dir_path, f"{input_fasta_path.stem}_L.fasta")
        short_fasta_path = Path(align_dir_path, f"{input_fasta_path.stem}_S.fasta")
        long_aligned_fasta_path = Path(align_dir_path, f"{input_fasta_path.stem}_A.fasta")
        intermediate_fasta_path = Path(align_dir_path, f"{input_fasta_path.stem}_I.fasta")
        all_fastas = [
            long_fasta_path,
            short_fasta_path,
            long_aligned_fasta_path,
            intermediate_fasta_path,
        ]
        mafft_log_file = Path(align_dir_path, f"{input_fasta_path.stem}.mafft.log")
        input_fasta = fasta_to_dict(input_fasta_path)
        long_fasta = {}
        short_fasta = {}

        max_seq_length = 0
        for seq_name in input_fasta:
            seq_length = len(input_fasta[seq_name]["sequence"])
            if seq_length > max_seq_length:
                max_seq_length = seq_length
        for seq_name in input_fasta:
            seq_length = len(input_fasta[seq_name]["sequence"])
            if seq_length > max_seq_length * settings.MIN_SEQ_LEN_PROP:
                long_fasta[seq_name] = input_fasta[seq_name]
            else:
                short_fasta[seq_name] = input_fasta[seq_name]
        dict_to_fasta(long_fasta, long_fasta_path)
        if len(short_fasta) > 0:
            dict_to_fasta(short_fasta, short_fasta_path)

        mafft_long_cmd = [
            mafft_path,
            settings.ALIGN_ALGORITHMS[mafft_algorithm]["arg"],
            "--nuc",
            "--maxiterate",
            "1000",
            "--adjustdirection",
            "--reorder",
            "--thread",
            f"{threads_per_alignment}",
            f"{long_fasta_path}",
        ]

        mafft_add_cmd = [
            mafft_path,
            "--auto",
            "--nuc",
            "--maxiterate",
            "1000",
            "--adjustdirection",
            "--reorder",
            "--addfragments",
            f"{short_fasta_path}",
            "--thread",
            f"{threads_per_alignment}",
            f"{long_aligned_fasta_path}",
        ]

        mid_fasta_path = long_aligned_fasta_path if len(short_fasta) > 0 else intermediate_fasta_path

        if len(long_fasta) > 1:
            with open(mid_fasta_path, "w") as mafft_mid_out:
                with open(mafft_log_file, "w") as mafft_log:
                    mafft_log.write(
                        f"Captus' MAFFT Command:\n  {' '.join(mafft_long_cmd)} > {mid_fasta_path}\n\n\n"
                    )
                with open(mafft_log_file, "a") as mafft_log:
                    try:
                        subprocess.run(
                            mafft_long_cmd, stdout=mafft_mid_out, stderr=mafft_log, timeout=timeout
                        )
                        if file_is_empty(mid_fasta_path):
                            cleanup_fastas(all_fastas)
                            message = red(
                                f"'{input_fasta_path.name}': FAILED alignment, empty output file"
                            )
                            return message
                    except subprocess.TimeoutExpired:
                        cleanup_fastas(all_fastas)
                        message = red(
                            f"'{input_fasta_path.name}': FAILED alignment, timeout"
                            f" exceeded [{elapsed_time(time.time() - start)}]"
                        )
                        return message
        else:
            long_fasta_path.replace(long_aligned_fasta_path)
        if len(short_fasta) > 0:
            with open(intermediate_fasta_path, "w") as mafft_intermediate_out:
                with open(mafft_log_file, "a") as mafft_log:
                    mafft_log.write(
                        f"\n\n\nCaptus' MAFFT Command:\n  {' '.join(mafft_add_cmd)}"
                        f" > {intermediate_fasta_path}\n\n\n"
                    )
                with open(mafft_log_file, "a") as mafft_log:
                    try:
                        subprocess.run(
                            mafft_add_cmd,
                            stdout=mafft_intermediate_out,
                            stderr=mafft_log,
                            timeout=timeout,
                        )
                        if file_is_empty(intermediate_fasta_path):
                            cleanup_fastas(all_fastas)
                            message = red(
                                f"'{input_fasta_path.name}': FAILED alignment, empty output file"
                            )
                            return message
                    except subprocess.TimeoutExpired:
                        cleanup_fastas(all_fastas)
                        message = red(
                            f"'{input_fasta_path.name}': FAILED alignment, timeout"
                            f" exceeded [{elapsed_time(time.time() - start)}]"
                        )
                        return message
        intermediate_fasta = fasta_to_dict(intermediate_fasta_path)
        output_fasta = {}
        for seq_name in intermediate_fasta:
            if seq_name.startswith("_R_"):
                new_seq_name = seq_name[3:]
            else:
                new_seq_name = seq_name
            output_fasta[new_seq_name] = {
                "sequence": intermediate_fasta[seq_name]["sequence"],
                "description": input_fasta[new_seq_name]["description"],
            }
        dict_to_fasta(output_fasta, output_fasta_path)
        cleanup_fastas(all_fastas)
        message = f"'{input_fasta_path.name}': aligned [{elapsed_time(time.time() - start)}]"
    else:
        message = dim(f"'{input_fasta_path.name}': SKIPPED (output FASTA already exists)")

    return message


def format_sample_list(sample_list: list, list_format="species"):
    output_list = []
    if sample_list is None:
        return output_list
    elif list_format == "species":
        for sample in sample_list.split(","):
            output_list.append("_".join(sample.split("_")[0:2]))
    elif list_format == "samples":
        output_list = sample_list.split(",")
    return list(set(output_list))


def find_and_merge_exon_data(out_dir: Path):
    tsvs = []
    if out_dir.exists():
        sample_dirs = list(Path(out_dir).rglob("*__captus-clr"))
        for sample_dir in sample_dirs:
            tsvs += list(Path(sample_dir).rglob(f"*{settings.DES_SUFFIXES['DATA']}"))
    exons_data = {}
    for tsv in tsvs:
        with open(tsv, "rt") as tsv_in:
            for line in tsv_in:
                if not line.startswith("cds_id"):
                    record = line.strip().split()
                    try:
                        prop_long_exons = int(record[5]) / int(record[3])
                    except ZeroDivisionError:
                        prop_long_exons = 0
                    try:
                        prop_short_exons = int(record[7]) / int(record[3])
                    except ZeroDivisionError:
                        prop_short_exons = 0
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
                        "prop_long_exons": prop_long_exons,
                        "prop_short_exons": prop_short_exons,
                    }
    return exons_data


def curate(
    gaps: float,
    bait_length: int,
    focal_species: list,
    outgroup_species: list,
    addon_samples: list,
    exons_data: dict,
    input_fasta_path: Path,
    curate_dir_path: Path,
    shared_aln_stats: list,
    overwrite: bool,
):
    def can_be_merged(hap: str, seq: str):
        """
        A sequence can be merged into the haplotype if they don't overlap or they don't differ in
        nucleotides, ignoring gaps

        Parameters
        ----------
        hap : str
            The haplotype
        seq : str
            The sequence to be merged in the haplotype
        """
        if len(hap.strip("-")) == 0:
            return True
        hap = hap.upper()
        seq = seq.upper()
        for p in range(len(hap)):
            if hap[p] == "-" or seq[p] == "-":
                continue
            elif hap[p] != seq[p]:
                return False
        return True

    def min_copies(aln_trimmed: dict, aln_width: int):
        haplos = {}
        for seq_name in aln_trimmed:
            sample_name = aln_trimmed[seq_name]["sample_name"]
            if sample_name in haplos:
                haplos[sample_name].append("-" * aln_width)
            else:
                haplos[sample_name] = ["-" * aln_width]
        for seq_name in aln_trimmed:
            sample_name = aln_trimmed[seq_name]["sample_name"]
            seq = aln_trimmed[seq_name]["sequence"]
            if len(haplos[sample_name]) == 1:
                haplos[sample_name][0] = seq
            else:
                for hap in range(len(haplos[sample_name])):
                    if can_be_merged(haplos[sample_name][hap], seq):
                        merged_hap = ""
                        for pos in range(len(seq)):
                            if haplos[sample_name][hap][pos] == "-":
                                merged_hap += seq[pos]
                            else:
                                merged_hap += haplos[sample_name][hap][pos]
                        haplos[sample_name][hap] = merged_hap
                        break
                    else:
                        continue
        min_copies = 0
        for sample_name in haplos:
            for haplo in haplos[sample_name]:
                if len(haplo.strip("-")) > 0:
                    min_copies += 1
        return min_copies

    start = time.time()
    curated_fasta_path = Path(curate_dir_path, input_fasta_path.name)

    if overwrite is True or not curated_fasta_path.exists():
        aln = fasta_to_dict(input_fasta_path)
        aln_height = len(aln)
        aln_width = None
        num_addons = 0
        col_depth = {}
        for seq_name in aln:
            if aln_width is None:
                aln_width = len(aln[seq_name]["sequence"])
                for p in range(aln_width):
                    col_depth[p] = 0
            else:
                if len(aln[seq_name]["sequence"]) != aln_width:
                    message = red(f"'{input_fasta_path.name}': FAILED (sequences were not aligned)")
                    return message
            aln[seq_name]["sample_name"] = seq_name.split(settings.SEQ_NAME_SEP)[0]
            aln[seq_name]["cds_id"] = seq_name.split(settings.SEQ_NAME_SEP)[1]
            aln[seq_name]["genus_name"] = aln[seq_name]["sample_name"].split("_")[0]
            aln[seq_name]["species_name"] = "_".join(aln[seq_name]["sample_name"].split("_")[0:2])
            aln[seq_name]["start"] = aln_width - len(aln[seq_name]["sequence"].lstrip("-"))
            aln[seq_name]["end"] = len(aln[seq_name]["sequence"].rstrip("-"))
            if aln[seq_name]["sample_name"] not in addon_samples:
                for p in range(aln[seq_name]["start"], aln[seq_name]["end"]):
                    col_depth[p] += 1
            else:
                num_addons += 1  # needed for calculating missing data excluding the add-on samples

        min_data = math.floor((aln_height - num_addons) * (1 - gaps))
        trim_start = 0
        trim_end = aln_width
        for p in sorted(col_depth):
            if col_depth[p] < min_data:
                trim_start = p + 1
            else:
                break
        for p in sorted(col_depth, reverse=True):
            if col_depth[p] < min_data:
                trim_end = p
            else:
                break

        aln_width = trim_end - trim_start
        aln_trimmed = {}
        aln_samples = []
        aln_genera = []
        aln_species = []
        aln_focal_species = []
        aln_outgroup_species = []
        aln_addon_samples = []
        cds_ids = {}
        for seq_name in aln:
            aln[seq_name]["sequence"] = aln[seq_name]["sequence"][trim_start:trim_end]
            seq_len_ungapped = len(aln[seq_name]["sequence"].replace("-", ""))
            if seq_len_ungapped >= bait_length:
                aln_trimmed[seq_name] = aln[seq_name]
                aln_samples.append(aln[seq_name]["sample_name"])
                aln_genera.append(aln[seq_name]["genus_name"])
                aln_species.append(aln[seq_name]["species_name"])
                if aln[seq_name]["species_name"] in focal_species:
                    aln_focal_species.append(aln[seq_name]["species_name"])
                if aln[seq_name]["species_name"] in outgroup_species:
                    aln_outgroup_species.append(aln[seq_name]["species_name"])
                if aln[seq_name]["sample_name"] in addon_samples:
                    aln_addon_samples.append(aln[seq_name]["sample_name"])
                if aln[seq_name]["cds_id"] in exons_data:
                    cds_ids[aln[seq_name]["cds_id"]] = seq_len_ungapped

        if len(aln_trimmed) > 1:
            ast = alignment_stats(aln_trimmed, "NT", coding=False)
            ast["genera"] = len(set(aln_genera))
            ast["species"] = len(set(aln_species))
            ast["focal"] = len(set(aln_focal_species))
            ast["outgroup"] = len(set(aln_outgroup_species))
            ast["samples"] = len(set(aln_samples))
            ast["addons"] = len(set(aln_addon_samples))
            ast["copies"] = min_copies(aln_trimmed, aln_width)
            ast["avg_copies"] = ast["copies"] / ast["samples"]
            ast["cds_id"] = "NA"
            ast["exons"] = "NA"
            ast["exons_len"] = "NA"
            ast["long_exons"] = "NA"
            ast["long_exons_len"] = "NA"
            ast["short_exons"] = "NA"
            ast["short_exons_len"] = "NA"
            ast["introns_len"] = "NA"
            ast["gene_len"] = "NA"
            ast["prop_exons_retained"] = "NA"
            ast["prop_long_exons"] = "NA"
            ast["prop_long_exons_retained"] = "NA"
            ast["len_long_exons_retained"] = "NA"
            ast["prop_short_exons"] = "NA"
            ast["prop_short_exons_retained"] = "NA"
            ast["len_short_exons_retained"] = "NA"
            ast["perc_exons_retained"] = "NA"
            ast["perc_long_exons_retained"] = "NA"
            ast["perc_short_exons_retained"] = "NA"

            if cds_ids:
                if len(cds_ids) > 1:
                    median_len = list(sorted(cds_ids.values()))[max(0, (len(cds_ids) // 2) - 1)]
                else:
                    median_len = list(sorted(cds_ids.values()))[0]
                cds_id = ""
                for cds in cds_ids:
                    if cds_ids[cds] == median_len:
                        cds_id = cds

                ast["cds_id"] = cds_id
                ast["exons"] = exons_data[cds_id]["exons"]
                ast["exons_len"] = exons_data[cds_id]["exons_len"]
                ast["long_exons"] = exons_data[cds_id]["long_exons"]
                ast["long_exons_len"] = exons_data[cds_id]["long_exons_len"]
                ast["short_exons"] = exons_data[cds_id]["short_exons"]
                ast["short_exons_len"] = exons_data[cds_id]["short_exons_len"]
                ast["introns_len"] = exons_data[cds_id]["introns_len"]
                ast["gene_len"] = exons_data[cds_id]["gene_len"]
                ast["prop_exons_retained"] = median_len / exons_data[cds_id]["exons_len"]
                ast["prop_long_exons"] = exons_data[cds_id]["prop_long_exons"]
                ast["prop_long_exons_retained"] = ast["prop_long_exons"] * ast["prop_exons_retained"]
                ast["len_long_exons_retained"] = ast["prop_long_exons_retained"] * ast["exons_len"]
                ast["prop_short_exons"] = exons_data[cds_id]["prop_short_exons"]
                ast["prop_short_exons_retained"] = ast["prop_short_exons"] * ast["prop_exons_retained"]
                ast["len_short_exons_retained"] = ast["prop_short_exons_retained"] * ast["exons_len"]
                # Format as text
                ast["len_long_exons_retained"] = f"{ast['len_long_exons_retained']:.0f}"
                ast["len_short_exons_retained"] = f"{ast['len_short_exons_retained']:.0f}"
                ast["perc_exons_retained"] = f"{ast['prop_exons_retained'] * 100:.2f}"
                ast["perc_long_exons_retained"] = f"{ast['prop_long_exons_retained'] * 100:.2f}"
                ast["perc_short_exons_retained"] = f"{ast['prop_short_exons_retained'] * 100:.2f}"

            stats = [
                f"{curated_fasta_path}",  # Path to curated FASTA
                f"{curated_fasta_path.stem}",  # Locus name
                f"{ast['copies']}",  # Minimum number of copies in locus
                f"{ast['avg_copies']:.2f}",  # Average number of copies in locus
                f"{ast['sites']}",  # Alignment length
                f"{ast['gc']:.2f}",  # % GC content
                f"{ast['avg_pid']:.2f}",  # % Pairwise identity
                f"{ast['informative']}",  # Number of informative sites
                f"{ast['informativeness']:.2f}",  # % Informativeness
                f"{ast['missingness']:.2f}",  # % Missingness
                f"{ast['sequences']}",  # Number of sequences
                f"{ast['samples']}",  # Number of samples
                f"{ast['focal']}",  # Number of focal species
                f"{ast['outgroup']}",  # Number of outgroup species
                f"{ast['addons']}",  # Number of add-on samples
                f"{ast['species']}",  # Number of species
                f"{ast['genera']}",  # Number of genera
                f"{ast['cds_id']}",  # CDS ID
                f"{ast['exons_len']}",  # Total exon length of CDS
                f"{ast['len_long_exons_retained']}",  # Long exons retained in bp
                f"{ast['len_short_exons_retained']}",  # Short exons retained in bp
                f"{ast['perc_exons_retained']}",  # Proportion of CDS retained
                f"{ast['perc_long_exons_retained']}",  # Proportion of long exons retained
                f"{ast['perc_short_exons_retained']}",  # Proportion of short exons retained
            ]

            shared_aln_stats += ["\t".join(stats) + "\n"]
            dict_to_fasta(aln_trimmed, curated_fasta_path)

            message = f"'{curated_fasta_path.name}': curated [{elapsed_time(time.time() - start)}]"
        else:
            message = red(
                f"'{input_fasta_path.name}': FAILED (trimmed alignment had {len(aln_trimmed)} sequences)"
            )
    else:
        message = dim(f"'{curated_fasta_path.name}': SKIPPED (output FASTA already exists)")

    return message


def write_aln_stats(out_dir: Path, shared_aln_stats: list):
    stats_tsv_file = Path(out_dir, "captus-cluster_alignments.tsv")
    if not shared_aln_stats:
        if stats_tsv_file.exists() and not file_is_empty(stats_tsv_file):
            return stats_tsv_file
        else:
            return None
    else:
        with open(stats_tsv_file, "wt") as tsv_out:
            tsv_out.write(
                "\t".join(
                    [
                        "path",
                        "locus",
                        "copies",
                        "avg_copies",
                        "length",
                        "gc_content",
                        "avg_pid",
                        "informative_sites",
                        "informativeness",
                        "missingness",
                        "sequences",
                        "samples",
                        "focal_species",
                        "outgroup_species",
                        "addon_samples",
                        "species",
                        "genera",
                        "cds_id",
                        "cds_len",
                        "len_long_exons_retained",
                        "len_short_exons_retained",
                        "perc_exons_retained",
                        "perc_long_exons_retained",
                        "perc_short_exons_retained",
                    ]
                )
                + "\n"
            )
            tsv_out.writelines(sorted(shared_aln_stats))
        return stats_tsv_file
