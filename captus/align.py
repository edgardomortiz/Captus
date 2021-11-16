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

import json
import shutil
import subprocess
import time
from multiprocessing import Manager
from pathlib import Path

from tqdm import tqdm

from . import log
from . import settings_assembly as settings
from .bioformats import alignment_stats, dict_to_fasta, fasta_to_dict, pairwise_identity
from .misc import (bold, clipkit_path_version, dim, elapsed_time, format_dep_msg, is_dir_empty,
                   mafft_path_version, make_output_dir, python_library_check, quit_with_error, red,
                   set_ram, set_threads, successful_exit, tqdm_parallel_async_run, tqdm_serial_run)
from .version import __version__


def align(full_command, args):

    captus_start = time.time()
    out_dir, out_dir_msg = make_output_dir(args.out)
    log.logger = log.Log(Path(args.out, "captus-assembly_align.log"), stdout_verbosity_level=1)

    mar = 25  # Margin for aligning parameters and values

    ################################################################################################
    ############################################################################### STARTING SECTION
    log.log_section_header("Starting Captus-assembly: ALIGN", single_newline=False)
    log.log_explanation(
        "Welcome to the alignment step of Captus-assembly. In this step, Captus will collect all the"
        f" extracted markers across all samples in '{args.captus_extractions_dir}' and group them by"
        " marker. Then Captus will align each marker using MAFFT. If given, the reference loci are"
        " also included as alignment guides and removed after alignment", extra_empty_lines_after=0
    )
    log.log_explanation("For more information, please see https://github.com/edgardomortiz/Captus")

    log.log(f'{"Captus version":>{mar}}: {bold(f"v{__version__}")}')
    log.log(f'{"Command":>{mar}}: {bold(full_command)}')
    _, _, ram_GB, ram_GB_total = set_ram(args.ram)
    log.log(f'{"Max. RAM":>{mar}}: {bold(f"{ram_GB:.1f}GB")} {dim(f"(out of {ram_GB_total:.1f}GB)")}')
    threads_max, threads_total = set_threads(args.threads)
    log.log(f'{"Max. Threads":>{mar}}: {bold(threads_max)} {dim(f"(out of {threads_total})")}')
    log.log("")

    log.log(f'{"Dependencies":>{mar}}:')
    mafft_path, mafft_version, mafft_status = mafft_path_version(args.mafft_path)
    log.log(format_dep_msg(f'{"MAFFT":>{mar}}: ', mafft_version, mafft_status))
    clipkit_path, clipkit_version, clipkit_status = clipkit_path_version(args.clipkit_path)
    log.log(format_dep_msg(f'{"ClipKIT":>{mar}}: ', clipkit_version, clipkit_status))
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

    if args.redo_from:
        log.log(f'{"Redo from":>{mar}}: {bold(args.redo_from)}')
        prepare_redo(out_dir, args.redo_from)
        log.log("")

    skip_alignment = False
    if mafft_status == "not found":
        skip_alignment = True
        log.log(
            f"{bold('WARNING:')} MAFFT could not be found, the markers will be collected from"
            f" '{args.captus_extractions_dir}' but they will not be aligned. Please verify you have"
            " it installed or provide the full path to the program with '--mafft_path'"
        )
    skip_trimming = False
    if clipkit_status == "not found":
        skip_trimming = True
        log.log(
            f"{bold('WARNING:')} ClipKIT could not be found, the alignments will not be trimmed."
            "Please verify you have it installed or provide the full path to the program with"
            " '--clipkit_path'"
        )


    ################################################################################################
    ############################################################################# COLLECTING SECTION
    log.log_section_header("Collecting Extracted Markers")
    log.log_explanation(
        f"Now Captus will search within '{args.captus_extractions_dir}' for all the extracted"
        " markers across all samples, and group them in a single FASTA file per marker name, keeping"
        " them organized in subdirectories according to marker type and alignment format. The output"
        " FASTA files are not aligned yet. If provided, Captus will include the reference aminoacid/"
        "nucleotide sequences in the alignments. "
    )
    markers, markers_ignored = check_value_list(args.markers, settings.MARKER_DIRS)
    log.log(f'{"Markers to collect":>{mar}}: {bold(markers)} {dim(markers_ignored)}')
    formats, formats_ignored = check_value_list(args.formats, settings.FORMAT_DIRS)
    log.log(f'{"Alignment formats":>{mar}}: {bold(formats)} {dim(formats_ignored)}')
    log.log(f'{"Max. paralogs":>{mar}}: {bold(args.max_paralogs)}')
    log.log("")
    refs_paths = prepare_refs(args.captus_extractions_dir, mar)
    log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
    log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
    extracted_sample_dirs = find_extracted_sample_dirs(args.captus_extractions_dir)
    log.log(f'{"Samples to process":>{mar}}: {bold(len(extracted_sample_dirs))}')
    log.log("")
    log.log(make_output_dirtree(markers, formats, out_dir, settings.ALN_DIRS["UNAL"], mar))
    log.log("")
    collect_extracted_markers(markers, formats, args.max_paralogs, extracted_sample_dirs, out_dir,
                              settings.ALN_DIRS["UNAL"], refs_paths, args.overwrite, args.show_less)
    log.log("")


    ################################################################################################
    ############################################################################## ALIGNMENT SECTION
    log.log_section_header("Marker alignment with MAFFT")
    log.log_explanation(
        "Now Captus will align all collected markers using MAFFT. If you added the references to be"
        " used as alignment guides Captus will produce a directory with alignments including the"
        " references and a separate one with the references removed. "
    )
    if skip_alignment:
        quit_with_error(
            "MAFFT could not be found, markers were collected across samples but they will not be"
            " aligned. Verify you have MAFFT installed or provide the full path to the program with"
            " '--mafft_path'"
        )
    else:
        concurrent, threads_per_alignment = adjust_mafft_concurrency(args.concurrent, threads_max)
        log.log(f'{"Concurrent alignments":>{mar}}: {bold(concurrent)}')
        log.log(f'{"Threads per alignment":>{mar}}: {bold(threads_per_alignment)}')
        log.log("")
        log.log(f'{"Algorithm":>{mar}}: {bold(args.mafft_algorithm)}'
                f' {dim(settings.MAFFT_ALGORITHMS[args.mafft_algorithm]["aka"])}')
        log.log(f'{"Timeout":>{mar}}: {bold(args.mafft_timeout)}'
                f' {dim(f"[{elapsed_time(args.mafft_timeout)}]")}')
        log.log("")
        log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
        log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
        fastas_to_align = fastas_origs_dests(out_dir,
                                            settings.ALN_DIRS["UNAL"],
                                            Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["UNFI"]))
        log.log(f'{"FASTA files to align":>{mar}}: {bold(len(fastas_to_align))}')
        log.log("")
        log.log(make_output_dirtree(markers,
                                    formats,
                                    out_dir,
                                    Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["UNFI"]),
                                    mar))
        log.log("")

        mafft_params = []
        for fasta_orig in fastas_to_align:
            mafft_params.append((
                mafft_path,
                args.mafft_algorithm,
                threads_per_alignment,
                args.mafft_timeout,
                fasta_orig,
                fastas_to_align[fasta_orig],
                args.overwrite,
            ))

        if args.debug:
            tqdm_serial_run(mafft, mafft_params,
                            "Aligning with MAFFT", "MAFFT alignment completed",
                            "alignment", args.show_less)
        else:
            tqdm_parallel_async_run(mafft, mafft_params,
                                    "Aligning with MAFFT", "MAFFT alignment completed",
                                    "alignment", concurrent, args.show_less)
        log.log("")


    ################################################################################################
    ################################################################################ CLEANUP SECTION
    log.log_section_header("Paralog Filtering")
    log.log_explanation(
        "Now Captus will remove paralogs using the method(s) selected with '--filter_method'."
        " Afterwards, copies of the alignments without the reference sequences will also be created."
    )
    concurrent = threads_max
    filtering_refs = select_filtering_refs(refs_paths, markers, formats, args.filter_method)
    filter_method = args.filter_method.lower()
    if not filtering_refs:
        if args.filter_method == "careful":
            filter_method = None
        elif args.filter_method == "both":
            filter_method = "fast"

    log.log(f'{"Concurrent processes":>{mar}}: {bold(concurrent)}')
    log.log(f'{"Filtering method":>{mar}}: {bold(filter_method)}')
    log.log("")
    if filtering_refs:
        log.log(bold(f'{"References for filtering":>{mar}}:'))
        for marker in filtering_refs:
            log.log(f'{"  Marker":>{mar}}: {bold(marker)}')
            log.log(f'{"  Format":>{mar}}: {bold(filtering_refs[marker]["format_dir"][-2:])}')
            log.log(f'{"  Path":>{mar}}: {bold(filtering_refs[marker]["path"])}')
            log.log("")
    log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
    log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
    log.log("")

    if filter_method in ["fast", "both"]:
        fastas_to_filter = fastas_origs_dests(
            out_dir,
            Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["UNFI"]),
            Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["FAST"])
        )
        log.log("")
        log.log(bold(f'{"FAST paralog filtering":>{mar}}:'))
        log.log(f'{"FASTA files to process":>{mar}}: {bold(len(fastas_to_filter))}')
        log.log(make_output_dirtree(markers,
                                    formats,
                                    out_dir,
                                    Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["FAST"]),
                                    mar))
        log.log("")
        paralog_fast_filter(fastas_to_filter, args.overwrite, concurrent,
                            args.show_less, args.debug)
        log.log("")

    if filter_method in ["careful", "both"]:
        fastas_to_filter = fastas_origs_dests(
            out_dir,
            Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["UNFI"]),
            Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["CARE"])
        )
        log.log("")
        log.log(bold(f'{"CAREFUL paralog filtering":>{mar}}:'))
        log.log(f'{"FASTA files to process":>{mar}}: {bold(len(fastas_to_filter))}')
        log.log(make_output_dirtree(markers,
                                    formats,
                                    out_dir,
                                    Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["CARE"]),
                                    mar))
        log.log("")
        manager = Manager()
        shared_paralog_stats = manager.list()
        paralog_careful_filter(shared_paralog_stats, fastas_to_filter, filtering_refs,
                               concurrent, args.overwrite, args.show_less, args.debug)
        paralog_stats_tsv = write_paralog_stats(out_dir, shared_paralog_stats)
        log.log("")
        log.log(f'{"Paralog statistics":>{mar}}: {bold(paralog_stats_tsv)}')
        log.log("")


    ################################################################################################
    ###################################################################### REFERENCE REMOVAL SECTION
    log.log_section_header("Reference Sequences Removal")
    log.log_explanation(
        "Now Captus will create copies of the alignnments that will not include the reference"
        " sequences used as alignment guide and for paralog filtering"
    )
    try:
        remove_references = bool(any([path for marker in refs_paths
                                           for path in refs_paths[marker].values()]))
    except TypeError:
        remove_references = False

    if remove_references:
        if Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["UNFI"]).exists():
            fastas_to_rem_refs = fastas_origs_dests(
                out_dir,
                Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["UNFI"]),
                Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NREF"])
            )
            log.log("")
            log.log(f'{"FASTA files to process":>{mar}}: {bold(len(fastas_to_rem_refs))}')
            log.log(make_output_dirtree(markers,
                                        formats,
                                        out_dir,
                                        Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NREF"]),
                                        mar))
            log.log("")
            rem_refs(refs_paths, fastas_to_rem_refs, args.overwrite,
                     concurrent, args.show_less, args.debug)
            log.log("")

        if Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["FAST"]).exists():
            fastas_to_rem_refs = fastas_origs_dests(
                out_dir,
                Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["FAST"]),
                Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NRFA"])
            )
            log.log("")
            log.log(f'{"FASTA files to process":>{mar}}: {bold(len(fastas_to_rem_refs))}')
            log.log(make_output_dirtree(markers,
                                        formats,
                                        out_dir,
                                        Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NRFA"]),
                                        mar))
            log.log("")
            rem_refs(refs_paths, fastas_to_rem_refs, args.overwrite,
                     concurrent, args.show_less, args.debug)
            log.log("")

        if Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["CARE"]).exists():
            fastas_to_rem_refs = fastas_origs_dests(
                out_dir,
                Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["CARE"]),
                Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NRCA"])
            )
            log.log("")
            log.log(f'{"FASTA files to process":>{mar}}: {bold(len(fastas_to_rem_refs))}')
            log.log(make_output_dirtree(markers,
                                        formats,
                                        out_dir,
                                        Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NRCA"]),
                                        mar))
            log.log("")
            rem_refs(refs_paths, fastas_to_rem_refs, args.overwrite,
                     concurrent, args.show_less, args.debug)
            log.log("")


    ################################################################################################
    ############################################################################### TRIMMING SECTION
    log.log_section_header("Alignment Trimming with ClipKIT")
    log.log_explanation(
        "Now Captus will trim all the alignments using ClipKIT. The trimming strategy can be"
        " specified with the flag --clipkit_algorithm. Trimmed alignments are saved separately in a"
        " new directory called '03_aligned_trimmed'"
    )
    if skip_trimming:
        log.log(red(
            "ClipKIT could not be found, alignments cannot be trimmed. Verify you have ClipKIT"
            " installed or provide the full path to the program with '--clipkit_path'"
        ))
    else:
        log.log(f'{"Concurrent processes":>{mar}}: {bold(concurrent)}')
        log.log("")
        log.log(f'{"Algorithm":>{mar}}: {bold(args.clipkit_algorithm)}')
        log.log(f'{"Gaps threshold":>{mar}}: {bold(args.clipkit_gaps)}')
        log.log("")
        log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
        log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
        fastas_to_trim = fastas_origs_dests(out_dir,
                                            settings.ALN_DIRS["ALND"],
                                            settings.ALN_DIRS["TRIM"])
        log.log(f'{"FASTA files to trim":>{mar}}: {bold(len(fastas_to_trim))}')
        log.log("")
        log.log("Creating output directories:")
        if Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["UNFI"]).exists():
            log.log(make_output_dirtree(markers,
                                        formats,
                                        out_dir,
                                        Path(settings.ALN_DIRS["TRIM"], settings.ALN_DIRS["UNFI"]),
                                        mar))
            log.log("")
        if Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["FAST"]).exists():
            log.log(make_output_dirtree(markers,
                                        formats,
                                        out_dir,
                                        Path(settings.ALN_DIRS["TRIM"], settings.ALN_DIRS["FAST"]),
                                        mar))
            log.log("")
        if Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["CARE"]).exists():
            log.log(make_output_dirtree(markers,
                                        formats,
                                        out_dir,
                                        Path(settings.ALN_DIRS["TRIM"], settings.ALN_DIRS["CARE"]),
                                        mar))
            log.log("")
        if Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NREF"]).exists():
            log.log(make_output_dirtree(markers,
                                        formats,
                                        out_dir,
                                        Path(settings.ALN_DIRS["TRIM"], settings.ALN_DIRS["NREF"]),
                                        mar))
            log.log("")
        if Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NRFA"]).exists():
            log.log(make_output_dirtree(markers,
                                        formats,
                                        out_dir,
                                        Path(settings.ALN_DIRS["TRIM"], settings.ALN_DIRS["NRFA"]),
                                        mar))
            log.log("")
        if Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NRCA"]).exists():
            log.log(make_output_dirtree(markers,
                                        formats,
                                        out_dir,
                                        Path(settings.ALN_DIRS["TRIM"], settings.ALN_DIRS["NRCA"]),
                                        mar))
            log.log("")

        clipkit_params = []
        for fasta_orig in fastas_to_trim:
            clipkit_params.append((
                clipkit_path,
                args.clipkit_algorithm,
                args.clipkit_gaps,
                fasta_orig,
                fastas_to_trim[fasta_orig],
                args.overwrite,
            ))

        if args.debug:
            tqdm_serial_run(clipkit, clipkit_params,
                            "Trimming alignments with ClipKIT", "ClipKIT trimming completed",
                            "alignment", args.show_less)
        else:
            tqdm_parallel_async_run(clipkit, clipkit_params,
                                    "Trimming alignments with ClipKIT", "ClipKIT trimming completed",
                                    "alignment", concurrent, args.show_less)
    log.log("")


    ################################################################################################
    ################################################################################ CLEANUP SECTION
    log.log_section_header("Statistics Summarization and File Cleanup")
    log.log_explanation(
        "Empty directories will be removed as well as MAFFT and ClipKIT logs unless the flag"
        " '--keep_all' is enabled."
    )
    fastas_to_stats = list(Path(out_dir, settings.ALN_DIRS["ALND"]).rglob("*.f[an]a"))
    fastas_to_stats += list(Path(out_dir, settings.ALN_DIRS["TRIM"]).rglob("*.f[an]a"))
    manager = Manager()
    shared_aln_stats = manager.list()

    compute_aln_stats_params = []
    for fasta in fastas_to_stats:
        compute_aln_stats_params.append((shared_aln_stats, fasta))

    if args.debug:
        tqdm_serial_run(compute_aln_stats, compute_aln_stats_params,
                        "Computing alignment statistics",
                        "Alignment statistics completed",
                        "alignment", args.show_less)
    else:
        tqdm_parallel_async_run(compute_aln_stats, compute_aln_stats_params,
                                "Computing alignment statistics",
                                "Alignment statistics completed",
                                "alignment", concurrent, args.show_less)
    log.log("")

    aln_stats_tsv = write_aln_stats(out_dir, shared_aln_stats)
    if aln_stats_tsv:
        log.log(f'{"Alignment statistics":>{mar}}: {bold(aln_stats_tsv)}')
        log.log("")
        if all([numpy_found, pandas_found, plotly_found]):

            from .report import build_alignment_report

            log.log_explanation(
                "Generating Alignment Statistics report..."
            )
            aln_html_report, aln_html_msg = build_alignment_report(out_dir, aln_stats_tsv)
            log.log(f'{"Alignment report":>{mar}}: {bold(aln_html_report)}')
            log.log(f'{"":>{mar}}  {dim(aln_html_msg)}')
        else:
            log.log(
                f"{bold('WARNING:')} Captus uses 'numpy', 'pandas', and 'plotly' to generate  an HTML"
                " report based on the marker recovery statistics. At least one of these libraries could"
                " not be found, please verify these libraries are installed and available."
            )
    else:
        log.log(red("Skipping summarization step... (no alignment statistics files were produced)"))
    log.log("")

    if not args.keep_all:
        start = time.time()
        log.log("")
        log.log_explanation(
            "Deleting MAFFT's and ClipKIT's logs and other unnecessary files..."
        )
        reclaimed_bytes = 0
        files_to_delete = list(out_dir.resolve().rglob("*.mafft.log"))
        files_to_delete += list(out_dir.resolve().rglob("*.clipkit.log"))
        for del_file in files_to_delete:
            reclaimed_bytes += del_file.stat().st_size
            del_file.unlink()
        log.log("")
        log.log(
            f'    A total of {len(files_to_delete)} files'
            f' amounting to {reclaimed_bytes / 1024 ** 2:.2f}MB'
            f' were deleted in {elapsed_time(time.time() - start)}'
        )
    else:
        log.log(bold("No files were removed, '--keep_all' was enabled"))
    dirs_to_delete = [d for d in out_dir.resolve().rglob("*") if d.is_dir()]
    # Sort directories to show the deepest subdirectories first
    dirs_to_delete = sorted(dirs_to_delete, key=lambda x: len(x.parts), reverse=True)
    for del_dir in dirs_to_delete:
        try:
            del_dir.rmdir()
        except OSError:
            continue
    log.log("")


    ################################################################################################
    ################################################################################# ENDING SECTION
    successful_exit(
        "Captus-assembly: ALIGN -> successfully completed"
        f" [{elapsed_time(time.time() - captus_start)}]"
    )


def prepare_redo(out_dir, redo_from):
    alignment = [Path(out_dir, settings.ALN_DIRS["ALND"])]
    filtering = [Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["FAST"]),
                 Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["CARE"])]
    removal = [Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NREF"]),
               Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NRFA"]),
               Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NRCA"])]
    trimming = [Path(out_dir, settings.ALN_DIRS["TRIM"])]
    dirs_to_delete = []
    if redo_from == "alignment":
        dirs_to_delete = alignment + trimming
    elif redo_from == "filtering":
        dirs_to_delete = filtering + removal + trimming
    elif redo_from == "removal":
        dirs_to_delete = removal + trimming
    elif redo_from == "trimming":
        dirs_to_delete = trimming

    start = time.time()
    log.log("")
    log.log(bold(f"Deleting directories to redo from the '{redo_from}' stage:"))
    tqdm_cols = min(shutil.get_terminal_size().columns, 120)
    with tqdm(total=len(dirs_to_delete), ncols=tqdm_cols, unit="directory") as pbar:
        for del_dir in dirs_to_delete:
            shutil.rmtree(del_dir, ignore_errors=True)
            tqdm.write(f"'{del_dir}': deleted")
            pbar.update()
    log.log(bold(
        f" \u2514\u2500\u2192 Ready for redo from the '{redo_from}' stage,"
        f" {len(dirs_to_delete)} directories deleted [{elapsed_time(time.time() - start)}]"
    ))


def check_value_list(input_values: str, valid_values: dict):
    """
    Verify a list of comma-separated values against the keys of a dictionary. Return comma-separated
    list of accepted values and a commented comma-separated list of ignored values.
    """
    if "ALL" in input_values.upper():
        checked = ",".join([v for v in valid_values])
        ignored = ""
    else:
        checked = ",".join([v for v in input_values.upper().split(",") if v in valid_values])
        ignored = ",".join([v for v in input_values.split(",") if v.upper() not in valid_values])
    if ignored: ignored = f'(ignored values: {ignored})'
    return checked, ignored


def find_extracted_sample_dirs(captus_extractions_dir):
    """
    List sample extraction directories within 'captus_extractions_dir'. Exit with error if
    'captus_extractions_dir' is not valid or does not contain valid sample extraction directories
    """
    captus_extractions_dir = Path(captus_extractions_dir)
    if not captus_extractions_dir.is_dir():
        quit_with_error(
            f"'{captus_extractions_dir}' is not a valid directory, please provide a valid directory"
            " with '--captus_extractions_dir'"
        )
    extracted_sample_dirs = list(captus_extractions_dir.resolve().rglob("*__captus-ext"))
    if not extracted_sample_dirs:
        quit_with_error(
            f"Captus did not find valid sample directories within '{captus_extractions_dir}' please"
            " provide a valid directory with '--captus_extractions_dir"
        )
    return extracted_sample_dirs


def prepare_refs(captus_ext_dir, margin):
    json_path = Path(captus_ext_dir, settings.JSON_REFS)

    try:
        with open(json_path, "rt") as jin:
            refs_paths = json.load(jin)
            for marker in refs_paths:
                for info in refs_paths[marker]:
                    if refs_paths[marker][info]:
                        if info.endswith("_path"):
                            refs_paths[marker][info] = Path(refs_paths[marker][info])
    except FileNotFoundError:
        refs_paths = {
            "NUC" : {"AA_path": None, "AA_msg": "not used", "NT_path": None, "NT_msg": "not used"},
            "PTD" : {"AA_path": None, "AA_msg": "not used", "NT_path": None, "NT_msg": "not used"},
            "MIT" : {"AA_path": None, "AA_msg": "not used", "NT_path": None, "NT_msg": "not used"},
            "DNA" : {"AA_path": None, "AA_msg": None,       "NT_path": None, "NT_msg": "not used"},
        }

    log.log(bold(f'{"Reference datasets":>{margin}}:'))

    return refs_paths


def select_filtering_refs(refs_paths, markers, formats, method):
    filtering_refs = {}
    markers = markers.upper().split(",")
    formats = formats.upper().split(",")

    if method.lower() in ["careful", "both"]:
        if "NT" in formats:
            for marker in refs_paths:
                if marker == "DNA":
                    if refs_paths[marker]["NT_path"]:
                        filtering_refs[marker] = {"path": refs_paths[marker]["NT_path"],
                                                  "marker_dir": settings.MARKER_DIRS[marker],
                                                  "format_dir": settings.FORMAT_DIRS["MA"]}
                elif marker != "CLR":
                    if refs_paths[marker]["NT_path"]:
                        filtering_refs[marker] = {"path": refs_paths[marker]["NT_path"],
                                                  "marker_dir": settings.MARKER_DIRS[marker],
                                                  "format_dir": settings.FORMAT_DIRS["NT"]}
        elif "AA" in formats:
            for marker in refs_paths:
                if marker not in ["DNA", "CLR"]:
                    if refs_paths[marker]["AA_path"]:
                        filtering_refs[marker] = {"path": refs_paths[marker]["AA_path"],
                                                  "marker_dir": settings.MARKER_DIRS[marker],
                                                  "format_dir": settings.FORMAT_DIRS["AA"]}

    return filtering_refs


def make_output_dirtree(markers, formats, out_dir, base_dir, margin):
    """
    Create directory structure to receive extracted markers from all the samples in
    `captus_extractions_dir`, print directory tree created with explanations of the intended content
    of each directory
    """
    base_tree = Path(out_dir, base_dir).resolve()
    markers = sorted([(settings.MARKER_DIRS[m], m) for m in markers.upper().split(",")])
    formats = sorted([(settings.FORMAT_DIRS[f], f) for f in formats.upper().split(",")])
    valid_combos = sorted([(m[0], f[0]) for m in markers for f in formats
                           if (m[1], f[1]) in settings.VALID_MARKER_FORMAT_COMBO])
    corner = "\u2514\u2500\u2500 "
    branch = "\u251C\u2500\u2500 "
    space = "    "
    vline = "\u2502   "
    symbols = ["", corner]
    dirs = [Path(out_dir).resolve(), base_dir]
    margins = [f'{"Output directory tree":>{margin}}: ', f'{"":>{margin}}  ']
    for combo in valid_combos:
        make_output_dir(Path(base_tree, combo[0], combo[1]))
        if not combo[0] in dirs:
            if combo[0] == valid_combos[-1][0]:
                parent, fork = space, corner
            else:
                parent, fork = vline, branch
            dirs.append(combo[0])
            symbols.append(f"{space}{fork}")
            margins.append(f'{"":>{margin}}  ')
        dirs.append(combo[1])
        if symbols[-1] == f'{space}{parent}{corner}':
            symbols[-1] = f'{space}{parent}{branch}'
        symbols.append(f'{space}{parent}{corner}')
        margins.append(f'{"":>{margin}}  ')
    return "\n".join([f"{m}{bold(s)}{bold(d)}" for m, s, d in zip(margins, symbols, dirs)])


def collect_extracted_markers(
    markers, formats, max_paralogs, extracted_sample_dirs, out_dir, base_dir, refs_paths,
    overwrite, show_less
):
    # Collect markers from extraction folder
    source_files = [Path(settings.MARKER_DIRS[m], f"{m}{settings.FORMAT_SUFFIXES[f]}")
                   for m in markers.split(",") for f in formats.split(",")
                   if (m, f) in settings.VALID_MARKER_FORMAT_COMBO]
    out_dirs = [Path(out_dir, base_dir, settings.MARKER_DIRS[m], settings.FORMAT_DIRS[f])
                for m in markers.split(",") for f in formats.split(",")
                if (m, f) in settings.VALID_MARKER_FORMAT_COMBO]
    for o in out_dirs:
        if is_dir_empty(o) is True or overwrite is True:
            for fasta_file in o.glob("*"):
                fasta_file.unlink()

    write_fastas = []
    for sample_dir in extracted_sample_dirs:
        for source, destination in zip(source_files, out_dirs):
            source_file = Path(sample_dir, source)
            if source_file.exists():
                fasta_in = fasta_to_dict(source_file, ordered=True)
                for seq_name_full in fasta_in:
                    seq_name_parts = seq_name_full.split("|")
                    marker_name = seq_name_parts[1]
                    seq_name = seq_name_parts[0]
                    if len(seq_name_parts) == 3:
                        seq_name += f"|{seq_name_parts[2]}"
                    fasta_out = Path(destination, f"{marker_name}{source.suffix}")
                    if overwrite is True or not fasta_out.exists():
                        write_fastas.append(fasta_out)
    write_fastas = list(set(write_fastas))

    collect_sample_markers_params = []
    for sample_dir in extracted_sample_dirs:
        collect_sample_markers_params.append([
            write_fastas,
            max_paralogs,
            sample_dir,
            source_files,
            out_dirs,
            overwrite,
        ])
    tqdm_serial_run(collect_sample_markers, collect_sample_markers_params,
                    "Collecting extracted markers", "Collection of extracted markers finished",
                    "sample", show_less)

    # Erase the collected FASTAs with fewer than 'settings.MIN_SAMPLES_ALN' samples
    log.log("")
    start = time.time()
    collected_fastas = list(Path(out_dir, base_dir).rglob("*.f[an]a"))
    deleted_fastas = 0
    log.log(bold("Verifying collected FASTA files:"))
    tqdm_cols = min(shutil.get_terminal_size().columns, 120)
    with tqdm(total=len(collected_fastas), ncols=tqdm_cols, unit="file") as pbar:
        for coll_fasta in collected_fastas:
            coll_fasta_len = num_samples(fasta_to_dict(coll_fasta))
            if coll_fasta_len < settings.MIN_SAMPLES_ALN:
                coll_fasta.unlink()
                deleted_fastas += 1
            pbar.update()
    log.log(bold(
        f" \u2514\u2500\u2192 Collected FASTAs verified, {deleted_fastas} file(s) deleted for having"
        f" fewer than {settings.MIN_SAMPLES_ALN} samples [{elapsed_time(time.time() - start)}]"
    ))

    # Add references to all the possible collected FASTAs
    if refs_paths is not None:
        refs = [
            refs_paths["NUC"]["AA_path"], refs_paths["NUC"]["NT_path"],
            refs_paths["PTD"]["AA_path"], refs_paths["PTD"]["NT_path"],
            refs_paths["MIT"]["AA_path"], refs_paths["MIT"]["NT_path"],
            refs_paths["DNA"]["NT_path"], refs_paths["DNA"]["NT_path"],
        ]
        mrks = ["NUC", "NUC", "PTD", "PTD", "MIT", "MIT", "DNA", "DNA"]
        fmts = ["AA", "NT", "AA", "NT", "AA", "NT", "MA", "MF"]
        add_refs_params = []
        for r, m, f in zip(refs, mrks, fmts):
            if all([r, m in markers.upper().split(","), f in formats.upper().split(",")]):
                add_refs_params.append((
                    r, Path(out_dir, base_dir, settings.MARKER_DIRS[m], settings.FORMAT_DIRS[f])
                ))
        if bool(add_refs_params):
            log.log("")
            tqdm_serial_run(add_refs, add_refs_params,
                            "Adding reference markers", "Addition of reference markers finished",
                            "reference", show_less)


def collect_sample_markers(write_fastas, max_paralogs, sample_dir, source_files, out_dirs, overwrite):
    start = time.time()
    sample_name = sample_dir.parts[-1].replace("__captus-ext", "")
    markers_collected = []
    fastas_created = []
    messages = []

    for source, destination in zip(source_files, out_dirs):
        source_file = Path(sample_dir, source)
        if source_file.exists():
            fasta_in = fasta_to_dict(source_file, ordered=True)
            for seq_name_full in fasta_in:
                seq_name_parts = seq_name_full.split("|")
                marker_name = seq_name_parts[1]
                seq_name = seq_name_parts[0]
                if len(seq_name_parts) == 3:
                    seq_name += f"|{seq_name_parts[2]}"
                    paralog_rank = int(seq_name_parts[2])
                else:
                    paralog_rank = 0
                if paralog_rank < max_paralogs:
                    fasta_out = Path(destination, f"{marker_name}{source.suffix}")
                    if fasta_out in write_fastas:
                        markers_collected.append(marker_name)
                        fastas_created.append(str(fasta_out))
                        dict_to_fasta({seq_name: dict(fasta_in[seq_name_full])},
                                      fasta_out,
                                      append=True)
                    else:
                        messages.append(dim(
                            f"'{Path(*fasta_out.parts[-4:])}': skipped (output FASTA already exists)"
                        ))
    messages = sorted(list(set(messages)))
    messages.append(
        f"'{sample_name}': {len(set(fastas_created))} FASTA files created for"
        f" {len(set(markers_collected))} collected markers [{elapsed_time(time.time() - start)}]"
    )
    return "\n".join(messages)


def add_refs(ref_path, dest_dir):
    start = time.time()
    fastas_in_dest = list(Path(dest_dir).rglob("*.f[an]a"))
    ref_fasta = fasta_to_dict(ref_path)
    markers_in_ref = []
    fastas_found = []

    for fasta in fastas_in_dest:
        refs_needed = []
        fasta_in = fasta_to_dict(fasta)
        for seq_name in fasta_in:
            if "|query=" in fasta_in[seq_name]["description"] and not seq_name.endswith("|ref"):
                refs_needed.append(
                    fasta_in[seq_name]["description"].split("|query=")[1].split(":")[0]
                )
        for ref_name in set(refs_needed):
            name_parts = ref_name.split(settings.REFERENCE_CLUSTER_SEPARATOR)
            ref_out = f"{settings.REFERENCE_CLUSTER_SEPARATOR.join(name_parts[0:-1])}|ref"
            if ref_name in ref_fasta:
                markers_in_ref.append(name_parts[-1])
                fastas_found.append(ref_name)
            if ref_out not in fasta_in:
                fasta_in[ref_out] = ref_fasta[ref_name]
        dict_to_fasta(fasta_in, fasta)

    if fastas_found:
        message = (
            f"'{ref_path.name}': {len(set(fastas_found))} sequences from {len(set(markers_in_ref))}"
            f" different markers were appended to collected FASTA files"
            f" [{elapsed_time(time.time() - start)}]"
        )
    else:
        message = red(
            f"'{ref_path.name}': did not find collected FASTA files to which append references"
            f" [{elapsed_time(time.time() - start)}]"
        )
    return message


def num_samples(fasta_dict):
    """
    Determine number of different samples in alignment with sequence length greater than 0,
    and excluding sequences whose name ends in '|ref'
    """
    samples = [seq_name.split("|")[0] for seq_name in fasta_dict
               if not seq_name.endswith("|ref")
               and len(fasta_dict[seq_name]["sequence"].replace("-","")) > 0]
    return len(set(samples))


def adjust_mafft_concurrency(concurrent, threads_max):
    if concurrent == "auto":
        concurrent = max(threads_max // 2, 1)
    else:
        try:
            concurrent = int(concurrent)
        except ValueError:
            quit_with_error("Invalid value for '--concurrent', set it to 'auto' or use a number")
    if concurrent > threads_max: concurrent = min(threads_max, concurrent)
    threads_per_alignment = threads_max // concurrent
    return concurrent, threads_per_alignment


def fastas_origs_dests(dir_path: Path, orig_base_dir: str, dest_base_dir: str):
    """
    Search for FASTA files within `dir_path/orig_base_dir` with extensions `.faa` or `.fna` and
    return a dictionary with these paths (origins) as keys and the paths with `dest_base_dir`
    instead of `orig_base_dir` as values (destinations).
    """
    fastas_to_process = {}
    for path in list(Path(dir_path, orig_base_dir).rglob("*.f[an]a")):
        origin = path.resolve()
        # destination = Path(*[p if p != orig_base_dir else dest_base_dir for p in origin.parts])
        destination = Path(str(origin).replace(str(orig_base_dir), str(dest_base_dir)))
        fastas_to_process[origin] = destination
    return fastas_to_process


def mafft(
    mafft_path, mafft_algorithm, threads, mafft_timeout, fasta_in: Path, fasta_out: Path, overwrite
):
    start = time.time()
    fasta_out_short = Path(*fasta_out.parts[-3:])
    if num_samples(fasta_to_dict(fasta_in)) < settings.MIN_SAMPLES_ALN:
        message = dim(
            f"'{fasta_out_short}': skipped (input FASTA has fewer than"
            f" {settings.MIN_SAMPLES_ALN} samples)"
        )
        return message
    if fasta_out.exists():
        if len(fasta_to_dict(fasta_out)) == 0:
            fasta_out.unlink()
    if overwrite is True or not fasta_out.exists():
        mafft_cmd = [
            mafft_path,
            settings.MAFFT_ALGORITHMS[mafft_algorithm]["arg"],
            "--maxiterate", "1000",
            "--reorder",
            "--thread", f"{threads}",
            f"{fasta_in}",
        ]
        mafft_log_file = Path(fasta_out.parent, f"{fasta_out.stem}.mafft.log")
        with open(fasta_out, "w") as mafft_out:
            with open(mafft_log_file, "w") as mafft_log:
                mafft_log.write(f"Captus' MAFFT Command:\n  {' '.join(mafft_cmd)}\n\n\n")
            with open(mafft_log_file, "a") as mafft_log:
                try:
                    subprocess.run(mafft_cmd, stdout=mafft_out, stderr=mafft_log,
                                   timeout=mafft_timeout)
                    message = f"'{fasta_out_short}': aligned [{elapsed_time(time.time() - start)}]"
                except subprocess.TimeoutExpired:
                    message = (
                        f"'{red(fasta_out_short)}': '--max_time' exceeded"
                        f" [{elapsed_time(time.time() - start)}]"
                    )
                    fasta_out.unlink()
                    mafft_log.write(
                        "\n\nERROR: The alignment took too long to complete, increase"
                        " '--mafft_timeout' or switch to a faster '--mafft_algorithm' like 'retree1'"
                        " or 'retree2'\n"
                    )
    else:
        message = dim(f"'{fasta_out_short}': skipped (output FASTA already exists)")
    return message


def paralog_fast_filter(fastas_paths, overwrite, concurrent, show_less, debug):
    filter_paralogs_fast_params = []
    for fasta in fastas_paths:
        filter_paralogs_fast_params.append((
            fasta,
            fastas_paths[fasta],
            overwrite,
        ))
    if debug:
        tqdm_serial_run(filter_paralogs_fast, filter_paralogs_fast_params,
                        "Removing potential paralogs from alignments",
                        "Potential paralog removal completed",
                        "alignment", show_less)
    else:
        tqdm_parallel_async_run(filter_paralogs_fast, filter_paralogs_fast_params,
                                "Removing potential paralogs from alignments",
                                "Potential paralog removal completed",
                                "alignment", concurrent, show_less)


def filter_paralogs_fast(fasta_source, fasta_dest, overwrite):
    start = time.time()
    fasta_dest_short = Path(*fasta_dest.parts[-3:])
    fasta_with_paralogs, fasta_without_paralogs = fasta_to_dict(fasta_source), {}
    if overwrite is True or not fasta_dest.exists():
        for seq_name in fasta_with_paralogs:
            if ("hit=00" in fasta_with_paralogs[seq_name]["description"] or
                "|ref" in seq_name):
                seq_name_out = seq_name.replace("|00", "")
                fasta_without_paralogs[seq_name_out] = dict(fasta_with_paralogs[seq_name])
        if num_samples(fasta_without_paralogs) >= settings.MIN_SAMPLES_ALN:
            dict_to_fasta(fasta_without_paralogs, fasta_dest)
            message = f"'{fasta_dest_short}': paralogs removed [{elapsed_time(time.time() - start)}]"
        else:
            message = red(
                f"'{fasta_dest_short}': not saved (filtered FASTA had fewer than"
                f" {settings.MIN_SAMPLES_ALN} samples) [{elapsed_time(time.time() - start)}]"
            )
    else:
        message = dim(f"'{fasta_dest_short}': skipped (output FASTA already exists)")
    return message


def paralog_careful_filter(
    shared_paralog_stats, fastas_paths, filtering_refs, concurrent, overwrite, show_less, debug
):
    fastas = {}
    for marker in filtering_refs:
        for fasta in fastas_paths:
            if (fasta.parts[-2] == filtering_refs[marker]["format_dir"]
                and fasta.parts[-3] == filtering_refs[marker]["marker_dir"]):
                fastas[fasta.stem] = fasta

    filter_paralogs_careful_params = []
    for marker_name in fastas:
        fastas_marker = {}
        for fasta in fastas_paths:
            if fasta.stem == marker_name and fasta.parts[-3] == fastas[marker_name].parts[-3]:
                fastas_marker[fasta] = fastas_paths[fasta]
        filter_paralogs_careful_params.append((
            shared_paralog_stats,
            fastas[marker_name],
            fastas_marker,
            overwrite
        ))

    if debug:
        tqdm_serial_run(filter_paralogs_careful, filter_paralogs_careful_params,
                        "Removing potential paralogs from markers",
                        "Potential paralog removal completed",
                        "marker", show_less)
    else:
        tqdm_parallel_async_run(filter_paralogs_careful, filter_paralogs_careful_params,
                                "Removing potential paralogs from markers",
                                "Potential paralog removal completed",
                                "marker", concurrent, show_less)


def filter_paralogs_careful(shared_paralog_stats, fasta_model, fastas_paths, overwrite):

    start = time.time()

    aln = fasta_to_dict(fasta_model, ordered=True)
    fasta_model_marker = fasta_model.parts[-3][-3:]
    fasta_model_format = fasta_model.parts[-2][-2:].replace("es", "NT")

    refs = {}
    for seq in aln:
        if "query=" in aln[seq]["description"] and "hit=00" in aln[seq]["description"]:
            ref = aln[seq]["description"].split("query=")[1].split(":")[0]
            refs[ref] += 1 if ref in refs else 1
    best_ref_full_name = max(refs, key=refs.get)
    s = settings.REFERENCE_CLUSTER_SEPARATOR
    best_ref = s.join(best_ref_full_name.split(s)[:-1])
    best_ref_seq = aln[f"{best_ref}|ref"]["sequence"]

    tsv = []
    accepted = []
    for seq in aln:
        if "|" not in seq or seq.endswith("|ref"):
            accepted.append(seq)

    samples_with_paralogs = {}
    for seq in aln:
        if "|" in seq and not seq.endswith("|ref"):
            sample_name = "|".join(seq.split("|")[:-1])
            hit_num = seq.split("|")[-1]
            length_seq = len(aln[seq]["sequence"].replace("-", ""))
            lenght_ref = len(best_ref_seq.replace("-", ""))
            pid = pairwise_identity(best_ref_seq, aln[seq]["sequence"], fasta_model_format)
            paralog_score = (pid / 100) * (length_seq / lenght_ref)
            if sample_name in samples_with_paralogs:
                samples_with_paralogs[sample_name][seq] = paralog_score
            else:
                samples_with_paralogs[sample_name] = {seq: paralog_score}
            tsv.append([
                fasta_model_marker,     # [0]  marker type
                fasta_model_format,     # [1]  format used for filtering
                fasta_model.stem,       # [2]  locus name
                best_ref_full_name,     # [3]  reference name
                sample_name,            # [4]  sample name
                hit_num,                # [5]  hit ranking
                seq,                    # [6]  sequence name
                f"{length_seq}",        # [7]  ungapped sequence length
                f"{pid:.5f}",           # [8]  identity to reference
                f"{paralog_score:.5f}", # [9]  pid * (len(seq) / len(ref))
                f"{False}",             # [10] accepted as ortholog
            ])
    for sample in samples_with_paralogs:
        accepted.append(max(samples_with_paralogs[sample], key=samples_with_paralogs[sample].get))

    for row in tsv:
        if row[6] in accepted: row[10] = f"{True}"

    shared_paralog_stats += tsv
    messages = []
    fastas_saved = len(fastas_paths)
    for fasta in fastas_paths:
        fasta_with_paralogs, fasta_without_paralogs = fasta_to_dict(fasta), {}
        if overwrite is True or not fastas_paths[fasta].exists():
            for seq_name in fasta_with_paralogs:
                if seq_name in accepted:
                    if seq_name.endswith("|ref") or not "|" in seq_name:
                        fasta_without_paralogs[seq_name] = fasta_with_paralogs[seq_name]
                    else:
                        seq_name_out = "|".join(seq_name.split("|")[:-1])
                        fasta_without_paralogs[seq_name_out] = fasta_with_paralogs[seq_name]
            if num_samples(fasta_without_paralogs) >= settings.MIN_SAMPLES_ALN:
                dict_to_fasta(fasta_without_paralogs, fastas_paths[fasta])
            else:
                fastas_saved -= 1
        else:
            messages.append(dim(
                f"'{Path(*fastas_paths[fasta].parts[-3:])}': skipped (output FASTA already exists)"
            ))
    if fastas_saved > 0:
        messages.append((
            f"'{fasta_model.parts[-3]}-{fasta_model.stem}': paralogs removed from"
            f" {fastas_saved} files [{elapsed_time(time.time() - start)}]"
        ))
    else:
        messages.append(red(
            f"'{fasta_model.parts[-3]}-{fasta_model.stem}': not saved (filtered FASTAs had fewer"
            f" than {settings.MIN_SAMPLES_ALN} samples) [{elapsed_time(time.time() - start)}]"
        ))

    return "\n".join(messages)


def write_paralog_stats(out_dir, shared_paralog_stats):
    if not shared_paralog_stats:
        return red("No paralogs were found...")
    else:
        stats_tsv_file = Path(out_dir, "captus-assembly_align.paralog_stats.tsv")
        with open(stats_tsv_file, "wt") as tsv_out:
            tsv_out.write("\t".join(["marker_type",
                                     "format_filtered",
                                     "locus",
                                     "ref",
                                     "sample",
                                     "hit",
                                     "sequence",
                                     "length",
                                     "identity",
                                     "paralog_score",
                                     "accepted"]) + "\n")
            for record in sorted(shared_paralog_stats):
                tsv_out.write("\t".join(record)+"\n")
        return stats_tsv_file


def rem_refs(refs_paths, fastas_paths, overwrite, concurrent, show_less, debug):
    ref_names = []
    for marker in refs_paths:
        if refs_paths[marker]["NT_path"]:
            ref_path = refs_paths[marker]["NT_path"]
        elif refs_paths[marker]["AA_path"]:
            ref_path = refs_paths[marker]["AA_path"]
        else:
            continue
        for seq_name in fasta_to_dict(ref_path):
            name_parts = seq_name.split(settings.REFERENCE_CLUSTER_SEPARATOR)
            ref_name = f"{settings.REFERENCE_CLUSTER_SEPARATOR.join(name_parts[0:-1])}|ref"
            if not ref_name in ref_names:
                ref_names.append(ref_name)

    rem_refs_from_fasta_params = []
    for fasta in fastas_paths:
        rem_refs_from_fasta_params.append((
            fasta,
            fastas_paths[fasta],
            ref_names,
            overwrite,
        ))
    if debug:
        tqdm_serial_run(rem_refs_from_fasta, rem_refs_from_fasta_params,
                        "Removing references from alignments", "Removal of references completed",
                        "alignment", show_less)
    else:
        tqdm_parallel_async_run(rem_refs_from_fasta, rem_refs_from_fasta_params,
                                "Removing references from alignments", "References removal completed",
                                "alignment", concurrent, show_less)


def rem_refs_from_fasta(fasta_source: Path, fasta_dest: Path, ref_names: list, overwrite):
    start = time.time()
    fasta_dest_short = Path(*fasta_dest.parts[-4:])
    fasta_with_refs, fasta_without_refs = fasta_to_dict(fasta_source), {}
    if overwrite is True or not fasta_dest.exists():
        for seq_name in fasta_with_refs:
            if seq_name not in ref_names:
                fasta_without_refs[seq_name] = dict(fasta_with_refs[seq_name])
        if num_samples(fasta_without_refs) >= settings.MIN_SAMPLES_ALN:
            dict_to_fasta(fasta_without_refs, fasta_dest)
            message = (
                f"'{fasta_dest_short}': references removed [{elapsed_time(time.time() - start)}]"
            )
        else:
            message = red(
                f"'{fasta_dest_short}': not saved (filtered FASTA had fewer than"
                f" {settings.MIN_SAMPLES_ALN} samples) [{elapsed_time(time.time() - start)}]"
            )
    else:
        message = dim(f"'{fasta_dest_short}': skipped (output FASTA already exists)")
    return message


def clipkit(
    clipkit_path, clipkit_algorithm, clipkit_gaps, fasta_in: Path, fasta_out: Path, overwrite
):
    start = time.time()
    fasta_out_short = Path(*fasta_out.parts[-4:])
    if overwrite is True or not fasta_out.exists():
        clipkit_cmd = [
            clipkit_path,
            f"{fasta_in}",
            "--output", f"{fasta_out}",
            "--mode", f"{clipkit_algorithm}",
            "--gaps", f"{clipkit_gaps}",
            "--input_file_format", "fasta",
            "--output_file_format", "fasta",
            "--log",
        ]
        clipkit_log_file = Path(fasta_out.parent, f"{fasta_out.stem}.clipkit.log")
        with open(clipkit_log_file, "w") as clipkit_log:
            clipkit_log.write(f"Captus' ClipKIT Command:\n  {' '.join(clipkit_cmd)}\n\n\n")
        with open(clipkit_log_file, "a") as clipkit_log:
            subprocess.run(clipkit_cmd, stdout=clipkit_log, stderr=clipkit_log)
        clipkit_stats_file = Path(fasta_out.parent, f"{fasta_out}.log")
        with open(clipkit_log_file, "a") as clipkit_log:
            with open(clipkit_stats_file, "rt") as clipkit_stats:
                clipkit_log.write("#Position Status Type Gaps\n")
                for line in clipkit_stats:
                    clipkit_log.write(line)
        clipkit_stats_file.unlink()
        if num_samples(fasta_to_dict(fasta_out)) >= settings.MIN_SAMPLES_ALN:
            message = f"'{fasta_out_short}': trimmed [{elapsed_time(time.time() - start)}]"
        else:
            fasta_out.unlink()
            message = red(
                f"'{fasta_out_short}': not saved (trimmed FASTA had fewer than"
                f" {settings.MIN_SAMPLES_ALN} samples) [{elapsed_time(time.time() - start)}]"
            )
    else:
        message = dim(f"'{fasta_out_short}': skipped (output FASTA already exists)")
    return message


def compute_aln_stats(shared_aln_stats, fasta_path):

    start = time.time()

    fasta_path_parts = fasta_path.parts

    aln_marker, aln_format= "", ""
    for m in settings.MARKER_DIRS:
        if fasta_path_parts[-3] == settings.MARKER_DIRS[m]: aln_marker = m
    for f in settings.FORMAT_DIRS:
        if fasta_path_parts[-2] == settings.FORMAT_DIRS[f]: aln_format = f

    aln_stats = alignment_stats(fasta_path)

    aln_tsv = [[
        f"{fasta_path}",                                    # [0] alignment file location
        f'{bool(not "untrimmed" in fasta_path_parts[-5])}', # [1] alignment was trimmed
        f'{fasta_path_parts[-4].split("_")[1]}',            # [2] paralog filter applied
        f'{bool("no_refs" in fasta_path_parts[-4])}',       # [3] references removed
        aln_marker,                                         # [4] marker type
        aln_format,                                         # [5] alignment format
        fasta_path.stem,                                    # [6] locus name
        f'{aln_stats["sequences"]}',                        # [7] num sequences
        f'{aln_stats["sites"]}',                            # [8] num sites
        f'{aln_stats["informative"]}',                      # [9] num informative sites
        f'{aln_stats["uninformative"]}',                    # [10] num constant + singleton sites
        f'{aln_stats["constant"]}',                         # [11] num constant sites
        f'{aln_stats["singleton"]}',                        # [12] num singleton sites
        f'{aln_stats["patterns"]}',                         # [13] num unique columns
        f'{aln_stats["avg_pid"]}',                          # [14] average pairwise identity
        f'{aln_stats["missingness"]}',                      # [15] pct of gaps and Ns or Xs
    ]]
    shared_aln_stats += aln_tsv

    short_path = Path(*fasta_path.parts[-5:])
    message = f"'{short_path}': stats computed [{elapsed_time(time.time() - start)}]"

    return message


def write_aln_stats(out_dir, shared_aln_stats):
    if not shared_aln_stats:
        return red("No alignment filtering statistics found...")
    else:
        stats_tsv_file = Path(out_dir, "captus-assembly_align.stats.tsv")
        with open(stats_tsv_file, "wt") as tsv_out:
            tsv_out.write("\t".join(["path",
                                     "trimmed",
                                     "paralog_filter",
                                     "no_refs",
                                     "marker_type",
                                     "format",
                                     "locus",
                                     "seqs",
                                     "sites",
                                     "informative",
                                     "uninformative",
                                     "constant",
                                     "singleton",
                                     "patterns",
                                     "avg_pid",
                                     "missingness",]) + "\n")
            for record in sorted(shared_aln_stats):
                tsv_out.write("\t".join(record)+"\n")
        return stats_tsv_file
