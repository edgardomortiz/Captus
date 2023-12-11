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

import json
import shutil
import statistics
import subprocess
import time
from multiprocessing import Manager
from pathlib import Path

from tqdm import tqdm

from . import log, settings
from .bioformats import (alignment_stats, dict_to_fasta, fasta_to_dict, fasta_type,
                         pairwise_identity, sample_stats)
from .misc import (bold, clipkit_path_version, dim, dir_is_empty, elapsed_time, file_is_empty,
                   format_dep_msg, mafft_path_version, make_output_dir, muscle_path_version,
                   python_library_check, quit_with_error, red, set_ram, set_threads,
                   successful_exit, tqdm_parallel_async_run, tqdm_serial_run)
from .version import __version__


def align(full_command, args):

    captus_start = time.time()
    out_dir, out_dir_msg = make_output_dir(args.out)
    log.logger = log.Log(Path(args.out, "captus-assembly_align.log"), stdout_verbosity_level=1)

    mar = 26  # Margin for aligning parameters and values

    ################################################################################################
    ############################################################################### STARTING SECTION
    log.log_section_header("Starting Captus-assembly: ALIGN", single_newline=False)
    log.log_explanation(
        "Welcome to the alignment step of Captus-assembly. In this step, Captus will collect all the"
        f" extracted markers across all samples in '{args.captus_extractions_dir}' and group them by"
        " marker. Then Captus will align each marker using MAFFT or MUSCLE. If given, the reference"
        " loci are also included as alignment guides and removed after alignment",
        extra_empty_lines_after=0
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
    muscle_path, muscle_version, muscle_status = muscle_path_version(args.muscle_path)
    if args.align_method.startswith("mafft"):
        log.log(format_dep_msg(f'{"MAFFT":>{mar}}: ', mafft_version, mafft_status))
        log.log(format_dep_msg(f'{"MUSCLE":>{mar}}: ', "", "not used"))
    if args.align_method.startswith("muscle"):
        log.log(format_dep_msg(f'{"MAFFT":>{mar}}: ', "", "not used"))
        log.log(format_dep_msg(f'{"MUSCLE":>{mar}}: ', muscle_version, muscle_status))
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

    skip_collection = False
    skip_alignment  = False
    skip_filtering  = False
    skip_removal    = False
    if args.redo_from:
        log.log(f'{"Redo from":>{mar}}: {bold(args.redo_from)}')
        if args.redo_from.lower()   == "alignment":
            skip_collection = True
        elif args.redo_from.lower() == "filtering":
            skip_collection = True
            skip_alignment  = True
        elif args.redo_from.lower() == "removal":
            skip_collection = True
            skip_alignment  = True
            skip_filtering  = True
        elif args.redo_from.lower() == "trimming":
            skip_collection = True
            skip_alignment  = True
            skip_filtering  = True
            skip_removal    = True
        prepare_redo(out_dir, args.redo_from)
        log.log("")

    if mafft_status == "not found":
        log.log(
            f"{bold('WARNING:')} MAFFT could not be found, the markers will be collected from"
            f" '{args.captus_extractions_dir}' but they will not be aligned. Please verify you have"
            " it installed or provide the full path to the program with '--mafft_path'"
        )
    if muscle_status == "not found":
        log.log(
            f"{bold('WARNING:')} MUSCLE could not be found, the markers will be collected from"
            f" '{args.captus_extractions_dir}' but they will not be aligned. Please verify you have"
            " it installed or provide the full path to the program with '--muscle_path'"
        )
    if clipkit_status == "not found":
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
    concurrent = max(1, min(settings.MAX_HDD_READ_INSTANCES, threads_max))
    markers, markers_ignored = check_value_list(args.markers, settings.MARKER_DIRS)
    formats, formats_ignored = check_value_list(args.formats, settings.FORMAT_DIRS)
    show_less = not args.show_more
    refs_paths = prepare_refs(args.captus_extractions_dir, mar)
    if skip_collection:
        log.log(red(
            f"Skipping the marker collection step because you used '--redo_from {args.redo_from}'"
        ))
        log.log("")
    else:
        log.log(f'{"Markers to collect":>{mar}}: {bold(markers)} {dim(markers_ignored)}')
        log.log(f'{"Alignment formats":>{mar}}: {bold(formats)} {dim(formats_ignored)}')
        max_paralogs_msg = dim("(Collect all paralogs)") if args.max_paralogs == -1 else ""
        log.log(f'{"Max. paralogs to collect":>{mar}}: {bold(args.max_paralogs)} {max_paralogs_msg}')
        log.log(f'{"Min. samples to align":>{mar}}: {bold(args.min_samples)}')
        log.log("")
        log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
        log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
        extracted_sample_dirs, skipped_align = find_extracted_sample_dirs(args.captus_extractions_dir)
        log.log(f'{"Samples to process":>{mar}}: {bold(len(extracted_sample_dirs))}')
        log.log("")
        if skipped_align:
            log.log(f'{bold("WARNING:")} {len(skipped_align)} sample(s) will be skipped')
            for msg in skipped_align:
                log.log(msg)
        log.log("")

        log.log(make_output_dirtree(markers, formats, out_dir, settings.ALN_DIRS["UNAL"], mar))
        log.log("")

        collect_extracted_markers(markers, formats, args.max_paralogs, args.min_samples,
                                  extracted_sample_dirs, out_dir, settings.ALN_DIRS["UNAL"],
                                  refs_paths, args.overwrite, show_less)
        log.log("")



    ################################################################################################
    ############################################################################## ALIGNMENT SECTION
    log.log_section_header("Marker alignment with MAFFT or MUSCLE")
    log.log_explanation(
        "Now Captus will align all collected markers using MAFFT or MUSCLE. If you added the"
        " references to be used as alignment guides Captus will produce a directory with alignments"
        " including the references and a separate one with the references removed. "
    )
    if skip_alignment:
        log.log(red(
            f"Skipping the marker alignment step because you used '--redo_from {args.redo_from}'"
        ))
        log.log("")
    elif args.align_method.startswith("mafft") and mafft_status == "not found":
        quit_with_error(
            "MAFFT could not be found, markers were collected across samples but they will not be"
            " aligned. Verify you have MAFFT installed or provide the full path to the program with"
            " '--mafft_path'"
        )
    elif args.align_method.startswith("muscle") and muscle_status == "not found":
        quit_with_error(
            "MUSCLE could not be found, markers were collected across samples but they will not be"
            " aligned. Verify you have MUSCLE installed or provide the full path to the program with"
            " '--muscle_path'"
        )
    else:
        concurrent, threads_per_alignment = adjust_align_concurrency(args.concurrent, threads_max)
        log.log(f'{"Concurrent alignments":>{mar}}: {bold(concurrent)}')
        log.log(f'{"Threads per alignment":>{mar}}: {bold(threads_per_alignment)}')
        log.log("")
        log.log(f'{"Algorithm":>{mar}}: {bold(args.align_method)}'
                f' {dim(settings.ALIGN_ALGORITHMS[args.align_method]["aka"])}')
        log.log(f'{"Timeout":>{mar}}: {bold(args.timeout)}'
                f' {dim(f"[{elapsed_time(args.timeout)}]")}')
        log.log(f'{"Codon-align CDS":>{mar}}: {bold(not(args.disable_codon_align))}')
        log.log(f'{"Outgroup":>{mar}}: {bold(args.outgroup)}')
        log.log("")
        log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
        log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
        fastas_to_align = fastas_origs_dests(out_dir,
                                             settings.ALN_DIRS["UNAL"],
                                             Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["UNFR"]))
        log.log(f'{"FASTAs to align":>{mar}}: {bold(len(fastas_to_align))}')
        log.log("")
        log.log(make_output_dirtree(markers,
                                    formats,
                                    out_dir,
                                    Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["UNFR"]),
                                    mar))
        log.log("")

        aa_origs, aa_dests = [], []
        nt_origs, nt_dests = [], []
        if not args.disable_codon_align:
            for fasta_orig in sorted(fastas_to_align):
                if fasta_orig.parts[-2] == settings.FORMAT_DIRS["AA"]:
                    aa_origs.append(fasta_orig)
                    aa_dests.append(fastas_to_align[fasta_orig])
                    nt_equiv = Path(str(fasta_orig).replace(
                        settings.FORMAT_DIRS["AA"], settings.FORMAT_DIRS["NT"]).replace(".faa", ".fna"))
                    if nt_equiv in fastas_to_align:
                        nt_origs.append(nt_equiv)
                        nt_dests.append(fastas_to_align[nt_equiv])
                    else:
                        nt_origs.append("")
                        nt_dests.append("")
        other_origs, other_dests = [], []
        for fasta_orig in sorted(fastas_to_align):
            if fasta_orig not in aa_origs and fasta_orig not in nt_origs:
                other_origs.append(fasta_orig)
                other_dests.append(fastas_to_align[fasta_orig])

        align_params = []
        for aa_orig, aa_dest, nt_orig, nt_dest in zip(aa_origs, aa_dests, nt_origs, nt_dests):
            fasta_origs = list(filter(None, [aa_orig, nt_orig]))
            fasta_dests = list(filter(None, [aa_dest, nt_dest]))
            align_params.append((
                mafft_path,
                muscle_path,
                args.align_method,
                threads_per_alignment,
                args.timeout,
                args.outgroup,
                fasta_origs,
                fasta_dests,
                args.min_samples,
                args.overwrite,
            ))
        for fasta_orig, fasta_dest in zip(other_origs, other_dests):
            align_params.append((
                mafft_path,
                muscle_path,
                args.align_method,
                threads_per_alignment,
                args.timeout,
                args.outgroup,
                [fasta_orig],
                [fasta_dest],
                args.min_samples,
                args.overwrite,
            ))

        aligner = "MAFFT" if args.align_method.startswith("mafft") else "MUSCLE"

        if args.debug:
            tqdm_serial_run(msa, align_params,
                            f"Aligning with {aligner}", f"{aligner} alignment completed",
                            "alignment", show_less)
        else:
            tqdm_parallel_async_run(msa, align_params,
                                    f"Aligning with {aligner}", f"{aligner} alignment completed",
                                    "alignment", concurrent, show_less)
        log.log("")


    ################################################################################################
    ###################################################################### PARALOG FILTERING SECTION
    log.log_section_header("Paralog Filtering")
    log.log_explanation(
        "Now Captus will remove paralogs using the method(s) selected with '--filter_method'."
        " Afterwards, copies of the alignments without the reference sequences will also be created."
    )
    concurrent = threads_max
    if skip_filtering:
        log.log(red(
            f"Skipping the paralog filtering step because you used '--redo_from {args.redo_from}'"
        ))
        log.log("")
    else:
        filtering_refs = select_filtering_refs(refs_paths, markers, formats, args.filter_method)
        filter_method = args.filter_method.lower()
        if not filtering_refs:
            if args.filter_method == "informed":
                filter_method = None
            elif args.filter_method == "both":
                filter_method = "naive"
        if filter_method == "none":
            filter_method = None

        log.log(f'{"Concurrent processes":>{mar}}: {bold(concurrent)}')
        log.log(f'{"Filtering method":>{mar}}: {bold(filter_method)}')
        log.log("")
        if filtering_refs:
            log.log(bold(f'{"References for filtering":>{mar}}:'))
            for marker in filtering_refs:
                filtering_format = "AA" if filtering_refs[marker]["format_dir"] == "01_AA" else "NT"
                log.log(f'{"  Marker":>{mar}}: {bold(marker)}')
                log.log(f'{"  Format":>{mar}}: {bold(filtering_format)}')
                log.log(f'{"  Path":>{mar}}: {bold(filtering_refs[marker]["path"])}')
                log.log("")
        log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
        log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
        log.log("")

        if filter_method in ["naive", "both"]:
            fastas_to_filter = fastas_origs_dests(
                out_dir,
                Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["UNFR"]),
                Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NAIR"])
            )
            log.log("")
            log.log(bold(f'{"NAIVE paralog filtering":>{mar}}:'))
            log.log(f'{"FASTA files to process":>{mar}}: {bold(len(fastas_to_filter))}')
            log.log(make_output_dirtree(markers,
                                        formats,
                                        out_dir,
                                        Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NAIR"]),
                                        mar))
            log.log("")
            paralog_naive_filter(fastas_to_filter, args.min_samples, args.overwrite,
                                 concurrent, args.debug, show_less)
            log.log("")

        if filter_method in ["informed", "both"]:
            fastas_to_filter = fastas_origs_dests(
                out_dir,
                Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["UNFR"]),
                Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["INFR"])
            )
            log.log("")
            log.log(bold(f'{"INFORMED paralog filtering":>{mar}}:'))
            log.log(f'{"Filter tolerance":>{mar}}: {bold(args.tolerance)} {bold("STDEVS")}')
            log.log(f'{"FASTA files to process":>{mar}}: {bold(len(fastas_to_filter))}')
            log.log(make_output_dirtree(markers,
                                        formats,
                                        out_dir,
                                        Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["INFR"]),
                                        mar))
            log.log("")
            manager = Manager()
            shared_paralog_stats = manager.list()
            paralog_informed_filter(shared_paralog_stats, fastas_to_filter, filtering_refs,
                                    args.tolerance, args.min_samples, args.overwrite,
                                    concurrent, args.debug, show_less)
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
    if skip_removal:
        log.log(red(
            f"Skipping the reference removal step because you used '--redo_from {args.redo_from}'"
        ))
        log.log("")
    else:
        try:
            remove_references = bool(any([path for marker in refs_paths
                                          for path in refs_paths[marker].values()]))
        except TypeError:
            remove_references = False

        if remove_references:
            log.log("Creating output directories:")
            fastas_to_rem_refs = {}
            if Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["UNFR"]).exists():
                fastas_partial = fastas_origs_dests(
                    out_dir,
                    Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["UNFR"]),
                    Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["UNFI"])
                )
                fastas_to_rem_refs = {**fastas_to_rem_refs, **fastas_partial}
                log.log(make_output_dirtree(markers,
                                            formats,
                                            out_dir,
                                            Path(settings.ALN_DIRS["ALND"],
                                                 settings.ALN_DIRS["UNFI"]),
                                            mar))
                log.log("")

            if Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NAIR"]).exists():
                fastas_partial = fastas_origs_dests(
                    out_dir,
                    Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NAIR"]),
                    Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NAIV"])
                )
                fastas_to_rem_refs = {**fastas_to_rem_refs, **fastas_partial}
                log.log(make_output_dirtree(markers,
                                            formats,
                                            out_dir,
                                            Path(settings.ALN_DIRS["ALND"],
                                                 settings.ALN_DIRS["NAIV"]),
                                            mar))
                log.log("")

            if Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["INFR"]).exists():
                fastas_partial = fastas_origs_dests(
                    out_dir,
                    Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["INFR"]),
                    Path(settings.ALN_DIRS["ALND"], settings.ALN_DIRS["INFO"])
                )
                fastas_to_rem_refs = {**fastas_to_rem_refs, **fastas_partial}
                log.log(make_output_dirtree(markers,
                                            formats,
                                            out_dir,
                                            Path(settings.ALN_DIRS["ALND"],
                                                 settings.ALN_DIRS["INFO"]),
                                            mar))
                log.log("")
            log.log(f'{"FASTA files to process":>{mar}}: {bold(len(fastas_to_rem_refs))}')
            log.log("")
            rem_refs(refs_paths, fastas_to_rem_refs, args.min_samples,
                     args.overwrite, concurrent, args.debug, show_less)
            log.log("")


    ################################################################################################
    ############################################################################### TRIMMING SECTION
    log.log_section_header("Alignment Trimming with ClipKIT")
    log.log_explanation(
        "Now Captus will trim all the alignments using ClipKIT. The trimming strategy can be"
        " specified with the flag '--clipkit_method'. Once columns are trimmed by ClipKIT, Captus"
        " will remove sequences with less than '--min_coverage' than the mean length of sequences,"
        " ignoring gaps. Trimmed alignments are saved separately in a new directory called "
        "'03_aligned_trimmed'"
    )
    if clipkit_status == "not found":
        log.log(red(
            "ClipKIT could not be found, alignments cannot be trimmed. Verify you have ClipKIT"
            " installed or provide the full path to the program with '--clipkit_path'"
        ))
    else:
        log.log(f'{"Concurrent processes":>{mar}}: {bold(concurrent)}')
        log.log("")
        log.log(f'{"Algorithm":>{mar}}: {bold(args.clipkit_method)}')
        clipkit_gaps = "auto" if args.clipkit_method == "smart-gap" else args.clipkit_gaps
        clipkit_gaps = "auto" if args.min_data_per_column > 0 else args.clipkit_gaps
        log.log(f'{"Gaps threshold":>{mar}}: {bold(clipkit_gaps)}')
        if args.min_data_per_column > 0:
            log.log(f'{"Min. sites per column":>{mar}}: {bold(args.min_data_per_column)}')
        else:
            log.log(f'{"Min. sites per column":>{mar}}: {dim("NOT USED")}')
        log.log(f'{"Min. sequence coverage":>{mar}}: {bold(args.min_coverage)}')
        log.log(f'{"Min. samples to keep":>{mar}}: {bold(args.min_samples)}')
        log.log("")
        log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
        log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
        fastas_to_trim = fastas_origs_dests(out_dir,
                                            settings.ALN_DIRS["ALND"],
                                            settings.ALN_DIRS["TRIM"])
        log.log(f'{"FASTA files to trim":>{mar}}: {bold(len(fastas_to_trim))}')
        log.log("")
        log.log("Creating output directories:")
        if Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["UNFR"]).exists():
            log.log(make_output_dirtree(markers,
                                        formats,
                                        out_dir,
                                        Path(settings.ALN_DIRS["TRIM"], settings.ALN_DIRS["UNFR"]),
                                        mar))
            log.log("")
        if Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NAIR"]).exists():
            log.log(make_output_dirtree(markers,
                                        formats,
                                        out_dir,
                                        Path(settings.ALN_DIRS["TRIM"], settings.ALN_DIRS["NAIR"]),
                                        mar))
            log.log("")
        if Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["INFR"]).exists():
            log.log(make_output_dirtree(markers,
                                        formats,
                                        out_dir,
                                        Path(settings.ALN_DIRS["TRIM"], settings.ALN_DIRS["INFR"]),
                                        mar))
            log.log("")
        if Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["UNFI"]).exists():
            log.log(make_output_dirtree(markers,
                                        formats,
                                        out_dir,
                                        Path(settings.ALN_DIRS["TRIM"], settings.ALN_DIRS["UNFI"]),
                                        mar))
            log.log("")
        if Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NAIV"]).exists():
            log.log(make_output_dirtree(markers,
                                        formats,
                                        out_dir,
                                        Path(settings.ALN_DIRS["TRIM"], settings.ALN_DIRS["NAIV"]),
                                        mar))
            log.log("")
        if Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["INFO"]).exists():
            log.log(make_output_dirtree(markers,
                                        formats,
                                        out_dir,
                                        Path(settings.ALN_DIRS["TRIM"], settings.ALN_DIRS["INFO"]),
                                        mar))
            log.log("")

        clipkit_params = []
        for fasta_orig in fastas_to_trim:
            clipkit_params.append((
                clipkit_path,
                args.clipkit_method,
                args.clipkit_gaps,
                fasta_orig,
                fastas_to_trim[fasta_orig],
                args.min_data_per_column,
                args.min_coverage,
                args.min_samples,
                args.overwrite,
            ))

        if args.debug:
            tqdm_serial_run(clipkit, clipkit_params,
                            "Trimming alignments with ClipKIT", "ClipKIT trimming completed",
                            "alignment", show_less)
        else:
            tqdm_parallel_async_run(clipkit, clipkit_params,
                                    "Trimming alignments with ClipKIT", "ClipKIT trimming completed",
                                    "alignment", concurrent, show_less)
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
    fastas_to_stats = list(Path(out_dir, settings.ALN_DIRS["ALND"]).rglob("*.f[an]a"))
    fastas_to_stats += list(Path(out_dir, settings.ALN_DIRS["TRIM"]).rglob("*.f[an]a"))
    fastas_to_stats = sorted([file for file in fastas_to_stats
                              if not f"{file.name}".startswith(".")])
    manager = Manager()
    shared_sam_stats = manager.list()
    shared_aln_stats = manager.list()

    compute_stats_params = []
    for fasta in fastas_to_stats:
        compute_stats_params.append((
            shared_sam_stats,
            shared_aln_stats,
            fasta
        ))

    if args.debug:
        tqdm_serial_run(compute_stats, compute_stats_params,
                        "Computing alignment and sample statistics",
                        "Alignment and sample statistics completed",
                        "alignment", show_less)
    else:
        tqdm_parallel_async_run(compute_stats, compute_stats_params,
                                "Computing alignment and sample statistics",
                                "Alignment and sample statistics completed",
                                "alignment", concurrent, show_less)
    log.log("")

    aln_stats_tsv = write_aln_stats(out_dir, shared_aln_stats)
    sam_stats_tsv = write_sam_stats(out_dir, shared_sam_stats)
    if aln_stats_tsv:
        log.log(f'{"Alignment statistics":>{mar}}: {bold(aln_stats_tsv)}')
        log.log(f'{"Sample statistics":>{mar}}: {bold(sam_stats_tsv)}')
        log.log("")
        if all([numpy_found, pandas_found, plotly_found]):

            from .report import build_alignment_report

            log.log("")
            log.log_explanation(
                "Generating Alignment Statistics report..."
            )
            aln_html_report, aln_html_msg = build_alignment_report(out_dir,
                                                                   aln_stats_tsv,
                                                                   sam_stats_tsv)
            log.log(f'{"Alignment report":>{mar}}: {bold(aln_html_report)}')
            log.log(f'{"":>{mar}}  {dim(aln_html_msg)}')
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
        log.log_explanation(
            "Deleting MAFFT's and ClipKIT's logs and other unnecessary files..."
        )
        reclaimed_bytes = 0
        files_to_delete = list(out_dir.resolve().rglob("*.mafft.log"))
        files_to_delete += list(out_dir.resolve().rglob("*.muscle.log"))
        files_to_delete += list(out_dir.resolve().rglob("*.clipkit.log"))
        for del_file in files_to_delete:
            reclaimed_bytes += del_file.stat().st_size
            del_file.unlink()
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
    filtering = [Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NAIR"]),
                 Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["INFR"])]
    removal = [Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["UNFI"]),
               Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["NAIV"]),
               Path(out_dir, settings.ALN_DIRS["ALND"], settings.ALN_DIRS["INFO"])]
    trimming = [Path(out_dir, settings.ALN_DIRS["TRIM"])]
    dirs_to_delete = []
    if redo_from.lower() == "alignment":
        dirs_to_delete = alignment + trimming
    elif redo_from.lower() == "filtering":
        dirs_to_delete = filtering + removal + trimming
    elif redo_from.lower() == "removal":
        dirs_to_delete = removal + trimming
    elif redo_from.lower() == "trimming":
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
    if ignored:
        ignored = f'(ignored values: {ignored})'

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
    all_sample_dirs = sorted(list(captus_extractions_dir.resolve().rglob("*__captus-ext")))
    extracted_sample_dirs = []
    skipped = []
    for sample_dir in all_sample_dirs:
        if settings.SEQ_NAME_SEP in f"{sample_dir}".replace("__captus-ext", ""):
            sample_name = sample_dir.parts[-1].replace("__captus-ext", "")
            skipped.append(f"'{sample_dir.parts[-1]}': SKIPPED, pattern"
                           f" '{settings.SEQ_NAME_SEP}' not allowed in sample name"
                           f" '{sample_name}'")
        else:
            extracted_sample_dirs.append(sample_dir)
    if not extracted_sample_dirs:
        quit_with_error(
            f"Captus did not find valid sample directories within '{captus_extractions_dir}' please"
            " provide a valid directory with '--captus_extractions_dir"
        )

    return extracted_sample_dirs, skipped


def prepare_refs(captus_ext_dir, margin):
    json_path = Path(captus_ext_dir, settings.JSON_REFS)
    not_found = []
    found = []
    try:
        with open(json_path, "rt") as jin:
            refs_paths = json.load(jin)
            for marker in sorted(refs_paths):
                for info in refs_paths[marker]:
                    if refs_paths[marker][info]:
                        if info.endswith("_path"):
                            refs_paths[marker][info] = Path(refs_paths[marker][info])
                            msg = f'{marker:>{margin}}: {refs_paths[marker][info]} ({info})'
                            if not refs_paths[marker][info].exists():
                                not_found.append(msg)
                                refs_paths[marker][info] = None
                            else:
                                found.append(msg)
    except FileNotFoundError:
        refs_paths = {
            "NUC" : {"AA_path": None, "AA_msg": "not used", "NT_path": None, "NT_msg": "not used"},
            "PTD" : {"AA_path": None, "AA_msg": "not used", "NT_path": None, "NT_msg": "not used"},
            "MIT" : {"AA_path": None, "AA_msg": "not used", "NT_path": None, "NT_msg": "not used"},
            "DNA" : {"AA_path": None, "AA_msg": None,       "NT_path": None, "NT_msg": "not used"},
        }
    if not found and not not_found:
        log.log(f'{"WARNING":>{margin}}: No reference files found!')
    else:
        if found:
            log.log(bold(f'{"Reference datasets":>{margin}}:'))
            for ref_path in found:
                log.log(ref_path)
            log.log("")
        if not_found:
            log.log(red((f"WARNING: The following reference files were not found. If you moved the"
                          " reference files used for extraction from their original location you can"
                         f" change the path in the file '{json_path}':\n")))
            for ref_path in not_found:
                log.log(red(ref_path))
                log.log("")
            log.log("")

    return refs_paths


def select_filtering_refs(refs_paths, markers, formats, method):
    filtering_refs = {}
    markers = markers.upper().split(",")
    formats = formats.upper().split(",")
    if method.lower() in ["informed", "both"]:
        if "NT" in formats:
            for marker in markers:
                if (marker in ["NUC", "PTD", "MIT"]
                    and refs_paths[marker]["NT_path"]):
                    filtering_refs[marker] = {"path": refs_paths[marker]["NT_path"],
                                              "marker_dir": settings.MARKER_DIRS[marker],
                                              "format_dir": settings.FORMAT_DIRS["NT"]}
        # Prefer coding reference in nucleotides to aminoacids for paralog filtering
        if "AA" in formats:
            for marker in markers:
                if (marker in ["NUC", "PTD", "MIT"]
                    and marker not in filtering_refs
                    and refs_paths[marker]["AA_path"]):
                    filtering_refs[marker] = {"path": refs_paths[marker]["AA_path"],
                                              "marker_dir": settings.MARKER_DIRS[marker],
                                              "format_dir": settings.FORMAT_DIRS["AA"]}
        if "MA" in formats:
            for marker in markers:
                if (marker in ["DNA", "CLR"]
                    and marker not in filtering_refs
                    and refs_paths[marker]["NT_path"]):
                    filtering_refs[marker] = {"path": refs_paths[marker]["NT_path"],
                                              "marker_dir": settings.MARKER_DIRS[marker],
                                              "format_dir": settings.FORMAT_DIRS["MA"]}

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
        if combo[0] not in dirs:
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
    markers, formats, max_paralogs, min_samples, extracted_sample_dirs, out_dir, base_dir,
    refs_paths, overwrite, show_less
):
    source_files = [Path(settings.MARKER_DIRS[m], f"{m}{settings.FORMAT_SUFFIXES[f]}")
                   for m in markers.split(",") for f in formats.split(",")
                   if (m, f) in settings.VALID_MARKER_FORMAT_COMBO]

    out_dirs = [Path(out_dir, base_dir, settings.MARKER_DIRS[m], settings.FORMAT_DIRS[f])
                for m in markers.split(",") for f in formats.split(",")
                if (m, f) in settings.VALID_MARKER_FORMAT_COMBO]
    for od in out_dirs:
        if dir_is_empty(od) is True or overwrite is True:
            for fasta_file in od.glob("*"):
                fasta_file.unlink()

    # List per-sample input FASTA files available for marker collection
    sample_names = []
    fastas_per_sample = {}
    for sample_dir in extracted_sample_dirs:
        sample_name = sample_dir.parts[-1].replace("__captus-ext", "")
        sample_names.append(sample_name)
        for source, destination in zip(source_files, out_dirs):
            source_file = Path(sample_dir, source)
            if source_file.exists():
                fastas_per_sample[source_file] = {
                    "sample_name": sample_name,
                    "destination": destination,
                    "suffix": source.suffix,
                }
    if not fastas_per_sample:
        quit_with_error(
            f"No FASTA files found for marker type(s) '{markers}', verify your '--markers' and"
            " '--formats' parameters, perhaps these markers were not extracted in the previous step"
        )

    # Collect markers per sample's source FASTAs into a dictionary that can be updated in parallel
    fastas_per_marker = {}
    collect_marker_names_params = []
    for source_fasta_path in fastas_per_sample:
        collect_marker_names_params.append([
            fastas_per_marker,
            source_fasta_path,
            fastas_per_sample[source_fasta_path]["sample_name"],
            fastas_per_sample[source_fasta_path]["destination"],
            fastas_per_sample[source_fasta_path]["suffix"],
            max_paralogs,
        ])
    tqdm_serial_run(collect_sample_markers, collect_marker_names_params,
                    "Collecting extracted markers",
                    "Collection of extracted markers finished",
                    "source", show_less)

    # Write FASTAs compiled per marker when they have at least four samples
    log.log("")
    start = time.time()
    accepted = 0 # FASTA output was saved and has at least four samples
    rejected = 0 # FASTA output was not saved because has fewer than four samples
    skipped  = 0 # FASTA output skipped because already exists
    log.log(bold("Verifying and writing collected FASTA files:"))
    tqdm_cols = min(shutil.get_terminal_size().columns, 120)
    with tqdm(total=len(fastas_per_marker), ncols=tqdm_cols, unit="file") as pbar:
        for fasta in fastas_per_marker:
            if overwrite is True or not fasta.exists():
                if num_samples(fastas_per_marker[fasta]) >= min_samples:
                    dict_to_fasta(fastas_per_marker[fasta], fasta, sort=True)
                    accepted += 1
                else:
                    rejected += 1
            else:
                skipped += 1
            pbar.update()
    log.log(bold(
        f" \u2514\u2500\u2192 {len(fastas_per_marker)} collected FASTAs processed"
        f" [{elapsed_time(time.time() - start)}]"
    ))
    log.log(
        f"     {accepted} saved, {rejected} not saved for having fewer than {min_samples} samples,"
        f" {skipped} already existed and were skipped"
    )

    # Add references to all the possible collected FASTAs
    if refs_paths is not None:
        refs = [
            refs_paths["NUC"]["AA_path"], refs_paths["NUC"]["NT_path"],
            refs_paths["PTD"]["AA_path"], refs_paths["PTD"]["NT_path"],
            refs_paths["MIT"]["AA_path"], refs_paths["MIT"]["NT_path"],
            refs_paths["DNA"]["NT_path"], refs_paths["DNA"]["NT_path"],
            refs_paths["CLR"]["NT_path"], refs_paths["CLR"]["NT_path"],
        ]
        mrks = ["NUC", "NUC", "PTD", "PTD", "MIT", "MIT", "DNA", "DNA", "CLR", "CLR"]
        fmts = [ "AA",  "NT",  "AA",  "NT",  "AA",  "NT",  "MA",  "MF",  "MA",  "MF"]
        add_refs_params = []
        manager = Manager()
        shared_ref_names = manager.list()
        for r, m, f in zip(refs, mrks, fmts):
            if all([r, m in markers.upper().split(","), f in formats.upper().split(",")]):
                add_refs_params.append((
                    r,
                    Path(out_dir, base_dir, settings.MARKER_DIRS[m], settings.FORMAT_DIRS[f]),
                    shared_ref_names,
                ))
        if bool(add_refs_params):
            log.log("")
            tqdm_serial_run(add_refs, add_refs_params,
                            "Adding reference markers", "Addition of reference markers finished",
                            "reference", show_less)

    # Write ASTRAL-Pro sequence to sample equivalence tsv file
    astral_pro_tsv = write_astral_pro_seq_to_sam(out_dir,
                                                 max_paralogs,
                                                 shared_ref_names,
                                                 sample_names)
    if astral_pro_tsv:
        log.log("")
        log.log(f'{"ASTRAL-Pro seq-to-sample":>{26}}: {bold(astral_pro_tsv)}')


def collect_sample_markers(
    shared_fastas_per_marker, source_fasta_path, sample_name, destination, suffix, max_paralogs
):
    start = time.time()
    markers_collected = []

    fasta_in = fasta_to_dict(source_fasta_path)
    for seq_name_full in fasta_in:
        # Replace connecting 'n's only for alignment, otherwise MAFFT takes forever and opens
        # enormous gaps
        if Path(destination).parts[-1] not in [settings.FORMAT_DIRS["AA"], settings.FORMAT_DIRS["NT"]]:
            seq_no_ns = fasta_in[seq_name_full]["sequence"].replace("n", "-")
            fasta_in[seq_name_full]["sequence"] = seq_no_ns
        # Replace stop codons by X, MAFFT automatically removes the '*' symbol
        if Path(destination).parts[-1] == settings.FORMAT_DIRS["AA"]:
            seq_no_stops = fasta_in[seq_name_full]["sequence"].replace("*", "X")
            fasta_in[seq_name_full]["sequence"] = seq_no_stops
        seq_name_parts = seq_name_full.split(settings.SEQ_NAME_SEP)
        marker_name = seq_name_parts[1]
        fasta_out = Path(destination, f"{marker_name}{suffix}")
        seq_name = seq_name_parts[0]
        if fasta_out not in shared_fastas_per_marker:
            shared_fastas_per_marker[fasta_out] = {}
        if len(seq_name_parts) == 3:
            seq_name += f"{settings.SEQ_NAME_SEP}{seq_name_parts[2]}"
            paralog_rank = int(seq_name_parts[2])
        else:
            paralog_rank = 0
        if paralog_rank <= max_paralogs or max_paralogs == -1:
            markers_collected.append(marker_name)
            shared_fastas_per_marker[fasta_out][seq_name] = dict(fasta_in[seq_name_full])
    message = (
        f"'{sample_name}': {len(set(markers_collected)):,} markers collected from"
        f" source '{source_fasta_path.name}' [{elapsed_time(time.time() - start)}]"
    )

    return message


def add_refs(ref_path, dest_dir, shared_ref_names):
    start = time.time()
    fastas_in_dest = list(Path(dest_dir).rglob("*.f[an]a"))
    ref_fasta = fasta_to_dict(ref_path)
    markers_in_ref = []
    fastas_found = []

    for fasta in fastas_in_dest:
        refs_needed = []
        fasta_in = fasta_to_dict(fasta)
        for seq_name in fasta_in:
            if ("[query=" in fasta_in[seq_name]["description"]
                and not seq_name.endswith(f"{settings.SEQ_NAME_SEP}ref")):
                refs_needed.append(
                    fasta_in[seq_name]["description"].split("[query=")[1].split("]")[0]
                )
        for ref_name in set(refs_needed):
            name_parts = ref_name.split(settings.REF_CLUSTER_SEP)
            ref_out = (f"{settings.REF_CLUSTER_SEP.join(name_parts[0:-1])}"
                       f"{settings.SEQ_NAME_SEP}ref")
            if ref_out not in shared_ref_names:
                shared_ref_names.append(ref_out)
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
    samples = [seq_name.split(settings.SEQ_NAME_SEP)[0] for seq_name in fasta_dict
               if not seq_name.endswith(f"{settings.SEQ_NAME_SEP}ref")
               and len(fasta_dict[seq_name]["sequence"].replace("-","")) > 0]

    return len(set(samples))


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


def fastas_origs_dests(dir_path: Path, orig_base_dir: str, dest_base_dir: str):
    """
    Search for FASTA files within `dir_path/orig_base_dir` with extensions `.faa` or `.fna` and
    return a dictionary with these paths (origins) as keys and the paths with `dest_base_dir`
    instead of `orig_base_dir` as values (destinations).
    """
    fastas_to_process = {}
    for path in sorted(list(Path(dir_path, orig_base_dir).rglob("*.f[an]a"))):
        if not f"{path.name}".startswith("."):
            origin = path.resolve()
            # Troubleshoot later, needed if we want to take simlinks as input
            # destination = Path(dir_path,
            #                    dest_base_dir,
            #                    origin.parent.parts[-2],
            #                    origin.parent.parts[-1],
            #                    origin.name)
            destination = Path(str(origin).replace(str(orig_base_dir), str(dest_base_dir)))
            fastas_to_process[origin] = destination

    return fastas_to_process


def rehead_root_msa(fasta_in: Path, fasta_out: Path, outgroup: list):
    if outgroup:
        outgroup = outgroup.split(",")
    else:
        outgroup = []
    unaligned = fasta_to_dict(fasta_in)
    aligned = fasta_to_dict(fasta_out)
    reheaded = {}
    for sample_name in outgroup:
        for seq_name in sorted(aligned):
            if seq_name.split("__")[0] == sample_name:
                reheaded[seq_name] = {
                    "description": unaligned[seq_name]["description"],
                    "sequence": aligned[seq_name]["sequence"]
                }
    for seq_name in aligned:
        if seq_name.split("__")[0] not in outgroup:
            reheaded[seq_name] = {
                "description": unaligned[seq_name]["description"],
                "sequence": aligned[seq_name]["sequence"]
            }
    dict_to_fasta(reheaded, fasta_out)
    return


def msa(
    mafft_path, muscle_path, align_method, threads, timeout, outgroup, fastas_in: list,
    fastas_out: list, min_samples, overwrite
):

    start = time.time()

    # If two FASTAs are passed to 'fastas_in', they are in order: AA nd NT;
    # align the first FASTA, if it is AA it is used as model for NT equivalent
    # which is the same as codon-aligning the NT coding sequence
    fasta_in = fastas_in[0]
    fasta_out = fastas_out[0]

    fasta_out_short = Path(*fasta_out.parts[-3:])
    if num_samples(fasta_to_dict(fasta_in)) < min_samples:
        message = dim(
            f"'{fasta_out_short}': SKIPPED (input FASTA has fewer than"
            f" {min_samples} samples)"
        )
        return message
    if fasta_out.exists():
        if len(fasta_to_dict(fasta_out)) == 0:
            fasta_out.unlink()
    if overwrite is True or not fasta_out.exists():
        if align_method.startswith("mafft"):
            if f"{fasta_in}".lower().endswith(".faa"):
                mafft_cmd = [
                    mafft_path,
                    settings.ALIGN_ALGORITHMS[align_method]["arg"],
                    "--amino",
                    "--maxiterate", "1000",
                    "--reorder",
                    "--thread", f"{threads}",
                    f"{fasta_in}",
                ]
            else:
                mafft_cmd = [
                    mafft_path,
                    settings.ALIGN_ALGORITHMS[align_method]["arg"],
                    "--nuc",
                    "--maxiterate", "1000",
                    "--reorder",
                    "--thread", f"{threads}",
                    f"{fasta_in}",
                ]
            mafft_log_file = Path(fasta_out.parent, f"{fasta_out.stem}.mafft.log")
            with open(mafft_log_file, "w") as mafft_log:
                mafft_log.write(f"Captus' MAFFT Command:\n  {' '.join(mafft_cmd)}"
                                f" > {fasta_out}\n\n\n")
            with open(fasta_out, "w") as mafft_out:
                with open(mafft_log_file, "a") as mafft_log:
                    try:
                        subprocess.run(mafft_cmd, stdout=mafft_out, stderr=mafft_log,
                                       timeout=timeout)
                        if file_is_empty(fasta_out):
                            message = red(f"'{fasta_out_short}': FAILED alignment, empty output file")
                            fasta_out.unlink()
                        else:
                            rehead_root_msa(fasta_in, fasta_out, outgroup)
                            message = f"'{fasta_out_short}': aligned [{elapsed_time(time.time() - start)}]"
                    except subprocess.TimeoutExpired:
                        message = red(
                            f"'{fasta_out_short}': FAILED alignment, timeout"
                            f" exceeded [{elapsed_time(time.time() - start)}]"
                        )
                        fasta_out.unlink()
                        mafft_log.write(
                            "\n\nERROR: The alignment took too long to complete, increase"
                            " '--timeout' or switch to a faster '--align_method' like 'mafft_retree1'"
                            " or 'retree2'\n"
                        )
        elif align_method.startswith("muscle"):
            if f"{fasta_in}".lower().endswith(".faa"):
                muscle_cmd = [
                    muscle_path,
                    "-amino",
                    "-threads", f"{threads}",
                    settings.ALIGN_ALGORITHMS[align_method]["arg"], f"{fasta_in}",
                    "-output", f"{fasta_out}"
                ]
            else:
                muscle_cmd = [
                    muscle_path,
                    "-nt",
                    "-threads", f"{threads}",
                    settings.ALIGN_ALGORITHMS[align_method]["arg"], f"{fasta_in}",
                    "-output", f"{fasta_out}"
                ]
            muscle_log_file = Path(fasta_out.parent, f"{fasta_out.stem}.muscle.log")
            with open(muscle_log_file, "w") as muscle_log:
                muscle_log.write(f"Captus' MUSCLE Command:\n  {' '.join(muscle_cmd)}\n\n\n")
            with open(muscle_log_file, "a") as muscle_log:
                try:
                    subprocess.run(muscle_cmd, stdout=muscle_log, stderr=muscle_log,
                                   timeout=timeout)
                    if file_is_empty(fasta_out):
                        message = red(f"'{fasta_out_short}': FAILED alignment, empty output file")
                        fasta_out.unlink()
                    elif not fasta_out.exists():
                        message = red(f"'{fasta_out_short}': FAILED alignment, output not generated")
                    else:
                        rehead_root_msa(fasta_in, fasta_out, outgroup)
                        message = f"'{fasta_out_short}': aligned [{elapsed_time(time.time() - start)}]"
                except subprocess.TimeoutExpired:
                    message = (
                        f"'{red(fasta_out_short)}': FAILED alignment, timeout exceeded"
                        f" [{elapsed_time(time.time() - start)}]"
                    )
                    # No need to erase output is timeout is exceeded, MUSCLE doesn't create it
                    # fasta_out.unlink()
                    muscle_log.write(
                        "\n\nERROR: The alignment took too long to complete, increase"
                        " '--timeout' or switch to a faster '--align_method' like 'muscle_super5'\n"
                    )

        if len(fastas_in) == 2:
            codon_message = codon_align(mafft_path,
                                        muscle_path,
                                        align_method,
                                        threads,
                                        timeout,
                                        outgroup,
                                        fasta_out,
                                        fastas_in[1],
                                        fastas_out[1],
                                        min_samples,
                                        overwrite)
            message += f'\n{codon_message}'

    else:
        message = dim(f"'{fasta_out_short}': SKIPPED (output FASTA already exists)")

    return message


def codon_align(
    mafft_path, muscle_path, align_method, threads, timeout, outgroup, aa_aligned: Path,
    nt_orig: Path, nt_dest: Path, min_samples, overwrite
):
    if not aa_aligned.exists() or file_is_empty(aa_aligned):
        return msa(mafft_path, muscle_path, align_method, threads, timeout,
                   outgroup, [nt_orig], [nt_dest], min_samples, overwrite)
    else:
        start = time.time()
        fasta_out_short = Path(*nt_dest.parts[-3:])
        if num_samples(fasta_to_dict(nt_orig)) < min_samples:
            message = dim(
                f"'{fasta_out_short}': SKIPPED (input FASTA has fewer than"
                f" {min_samples} samples)"
            )
            return message
        if nt_dest.exists():
            if len(fasta_to_dict(nt_dest)) == 0:
                nt_dest.unlink()
        if overwrite is True or not nt_dest.exists():
            aa_aligned = fasta_to_dict(aa_aligned)
            nt_unaligned = fasta_to_dict(nt_orig)
            nt_aligned = {}
            for seq_name in aa_aligned:
                if seq_name in nt_unaligned:
                    seq_aa = aa_aligned[seq_name]["sequence"]
                    seq_nt = nt_unaligned[seq_name]["sequence"]
                    nt_start = 0
                    seq_nt_out = ""
                    for aa in seq_aa:
                        if aa == "-":
                            seq_nt_out += "-"*3
                        else:
                            seq_nt_out += seq_nt[nt_start:nt_start+3]
                            nt_start += 3
                    nt_aligned[seq_name] = {
                        "description": nt_unaligned[seq_name]["description"],
                        "sequence": seq_nt_out,
                    }
            dict_to_fasta(nt_aligned, nt_dest)
            message = f"'{fasta_out_short}': codon-aligned [{elapsed_time(time.time() - start)}]"
        else:
            message = dim(f"'{fasta_out_short}': SKIPPED (output FASTA already exists)")

        return message


def paralog_naive_filter(fastas_paths, min_samples, overwrite, concurrent, debug, show_less):
    filter_paralogs_naive_params = []
    for fasta in fastas_paths:
        filter_paralogs_naive_params.append((
            fasta,
            fastas_paths[fasta],
            min_samples,
            overwrite,
        ))
    if debug:
        tqdm_serial_run(filter_paralogs_naive, filter_paralogs_naive_params,
                        "Removing potential paralogs from alignments",
                        "Potential paralog removal completed",
                        "alignment", show_less)
    else:
        tqdm_parallel_async_run(filter_paralogs_naive, filter_paralogs_naive_params,
                                "Removing potential paralogs from alignments",
                                "Potential paralog removal completed",
                                "alignment", concurrent, show_less)


def filter_paralogs_naive(fasta_in: Path, fasta_out: Path, min_samples, overwrite):
    start = time.time()

    fasta_out_short = Path(*fasta_out.parts[-3:])
    if file_is_empty(fasta_in):
        return red(f"'{fasta_out_short}': FAILED paralog removal, input file was empty")

    if overwrite is True or not fasta_out.exists():
        fasta_with_paralogs, fasta_without_paralogs = fasta_to_dict(fasta_in), {}
        for seq_name in fasta_with_paralogs:
            if ("hit=00" in fasta_with_paralogs[seq_name]["description"] or
                f"{settings.SEQ_NAME_SEP}ref" in seq_name):
                seq_name_out = seq_name.replace(f"{settings.SEQ_NAME_SEP}00", "")
                fasta_without_paralogs[seq_name_out] = dict(fasta_with_paralogs[seq_name])
        if num_samples(fasta_without_paralogs) >= min_samples:
            dict_to_fasta(fasta_without_paralogs, fasta_out)
            message = f"'{fasta_out_short}': paralogs removed [{elapsed_time(time.time() - start)}]"
        else:
            message = red(
                f"'{fasta_out_short}': not saved (filtered FASTA had fewer than"
                f" {min_samples} samples) [{elapsed_time(time.time() - start)}]"
            )
    else:
        message = dim(f"'{fasta_out_short}': SKIPPED (output FASTA already exists)")

    return message


def paralog_informed_filter(
    shared_paralog_stats, fastas_paths, filtering_refs, tolerance, min_samples, overwrite,
    concurrent, debug, show_less
):
    fastas = {}
    for marker in filtering_refs:
        for fasta in fastas_paths:
            format_dir = filtering_refs[marker]["format_dir"]
            marker_dir = filtering_refs[marker]["marker_dir"]
            if fasta.parts[-2] == format_dir and fasta.parts[-3] == marker_dir:
                if marker_dir in fastas:
                    fastas[marker_dir][fasta.stem] = fasta
                else:
                    fastas[marker_dir] = {fasta.stem: fasta}

    filter_paralogs_informed_params = []
    for marker_dir in fastas:
        for marker_name in fastas[marker_dir]:
            fastas_marker = {}
            for fasta in fastas_paths:
                if (fasta.stem == marker_name
                    and fasta.parts[-3] == fastas[marker_dir][marker_name].parts[-3]):
                    fastas_marker[fasta] = fastas_paths[fasta]
            filter_paralogs_informed_params.append((
                shared_paralog_stats,
                fastas[marker_dir][marker_name],
                fastas_marker,
                tolerance,
                min_samples,
                overwrite
            ))

    if debug:
        tqdm_serial_run(filter_paralogs_informed, filter_paralogs_informed_params,
                        "Removing potential paralogs from markers",
                        "Potential paralog removal completed",
                        "marker", show_less)
    else:
        tqdm_parallel_async_run(filter_paralogs_informed, filter_paralogs_informed_params,
                                "Removing potential paralogs from markers",
                                "Potential paralog removal completed",
                                "marker", concurrent, show_less)


def filter_paralogs_informed(
    shared_paralog_stats, fasta_model, fastas_paths, tolerance, min_samples, overwrite
):

    start = time.time()

    output_missing = False
    for fasta in fastas_paths:
        if not fastas_paths[fasta].exists():
            output_missing = True
            break

    fasta_model_short = Path(*fastas_paths[fasta_model].parts[-3:])

    messages = []

    if overwrite is True or output_missing is True:
        if file_is_empty(fasta_model):
            return red(f"'{fasta_model_short}': FAILED paralog removal, input file was empty")

        aln = fasta_to_dict(fasta_model)
        fasta_model_marker = fasta_model.parts[-3][-3:]
        if fasta_model.parts[-2] == "01_AA":
            fasta_model_format = "AA"
        elif fasta_model.parts[-2] == "02_NT":
            fasta_model_format = "NT"
        elif fasta_model.parts[-2] == "01_matches":
            fasta_model_format = "NT"

        refs = {}
        for seq in aln:
            if not seq.endswith(f"{settings.SEQ_NAME_SEP}ref"):
                if "query=" in aln[seq]["description"] and "hit=00" in aln[seq]["description"]:
                    ref = aln[seq]["description"].split("query=")[1].split("]")[0]
                    if ref in refs:
                        refs[ref] += 1
                    else:
                        refs[ref] = 1
        best_ref_full_name = max(refs, key=refs.get)
        s = settings.REF_CLUSTER_SEP
        best_ref = s.join(best_ref_full_name.split(s)[:-1])
        best_ref_seq = aln[f"{best_ref}{settings.SEQ_NAME_SEP}ref"]["sequence"]

        tsv = []
        accepted = [seq for seq in aln if seq.endswith(f"{settings.SEQ_NAME_SEP}ref")]
        samples_with_paralogs = {}
        for seq in aln:
            if not seq.endswith(f"{settings.SEQ_NAME_SEP}ref"):
                if settings.SEQ_NAME_SEP in seq:
                    sample_name = settings.SEQ_NAME_SEP.join(seq.split(settings.SEQ_NAME_SEP)[:-1])
                    hit_num = seq.split(settings.SEQ_NAME_SEP)[-1]
                else:
                    sample_name = seq
                    hit_num = "00"
                length_seq = len(aln[seq]["sequence"].replace("-", ""))
                lenght_ref = len(best_ref_seq.replace("-", ""))
                pid = pairwise_identity(best_ref_seq, aln[seq]["sequence"],
                                        fasta_model_format,
                                        ignore_internal_gaps=True)
                paralog_score = (pid / 100) * (length_seq / lenght_ref)
                if sample_name in samples_with_paralogs:
                    samples_with_paralogs[sample_name][seq] = paralog_score
                else:
                    samples_with_paralogs[sample_name] = {seq: paralog_score}
                tsv.append([
                    fasta_model_marker,       # [0]  marker type
                    fasta_model_format,       # [1]  format used for filtering
                    fasta_model.stem,         # [2]  locus name
                    best_ref_full_name,       # [3]  reference name
                    sample_name,              # [4]  sample name
                    hit_num,                  # [5]  hit ranking
                    seq,                      # [6]  sequence name
                    f"{length_seq}",          # [7]  ungapped sequence length
                    f"{pid:.5f}",             # [8]  identity to reference
                    f"{paralog_score:.5f}",   # [9]  pid * (len(seq) / len(ref))
                    f"{False}",               # [10] accepted as ortholog
                ])
        for sample in samples_with_paralogs:
            accepted.append(max(samples_with_paralogs[sample], key=samples_with_paralogs[sample].get))

        pids = [float(row[8]) for row in tsv if row[6] in accepted]
        min_pid = statistics.mean(pids) - (tolerance * statistics.stdev(pids))
        for row in tsv:
            if row[6] in accepted:
                if float(row[8]) >= min_pid:
                    row[10] = f'{True}'
                else:
                    del accepted[accepted.index(row[6])]

        shared_paralog_stats += ["\t".join(row)+"\n" for row in tsv]
        fastas_saved = len(fastas_paths)
        for fasta in fastas_paths:
            fasta_with_paralogs, fasta_without_paralogs = fasta_to_dict(fasta), {}
            for seq_name in fasta_with_paralogs:
                if seq_name in accepted:
                    if (seq_name.endswith(f"{settings.SEQ_NAME_SEP}ref")
                        or settings.SEQ_NAME_SEP not in seq_name):
                        fasta_without_paralogs[seq_name] = fasta_with_paralogs[seq_name]
                    else:
                        seq_name_out = settings.SEQ_NAME_SEP.join(
                            seq_name.split(settings.SEQ_NAME_SEP)[:-1]
                        )
                        fasta_without_paralogs[seq_name_out] = fasta_with_paralogs[seq_name]
            if num_samples(fasta_without_paralogs) >= min_samples:
                dict_to_fasta(fasta_without_paralogs, fastas_paths[fasta])
            else:
                fastas_saved -= 1
        if fastas_saved > 0:
            messages.append((
                f"'{fasta_model.parts[-3]}-{fasta_model.stem}': paralogs removed from"
                f" {fastas_saved} files [{elapsed_time(time.time() - start)}]"
            ))
        else:
            messages.append(red(
                f"'{fasta_model.parts[-3]}-{fasta_model.stem}': not saved (filtered FASTAs had fewer"
                f" than {min_samples} samples) [{elapsed_time(time.time() - start)}]"
            ))
    else:
        messages.append(dim(
            f"'{fasta_model_short}': SKIPPED (output FASTA already exists)"
        ))

    return "\n".join(messages)


def write_paralog_stats(out_dir, shared_paralog_stats):
    if not shared_paralog_stats:
        return red("No paralogs were found...")
    else:
        stats_tsv_file = Path(out_dir, "captus-assembly_align.paralogs.tsv")
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
            tsv_out.writelines(sorted(shared_paralog_stats))
        return stats_tsv_file


def rem_refs(refs_paths, fastas_paths, min_samples, overwrite, concurrent, debug, show_less):
    ref_names = []
    for marker in refs_paths:
        if refs_paths[marker]["NT_path"]:
            ref_path = refs_paths[marker]["NT_path"]
        elif refs_paths[marker]["AA_path"]:
            ref_path = refs_paths[marker]["AA_path"]
        else:
            continue
        for seq_name in fasta_to_dict(ref_path):
            name_parts = seq_name.split(settings.REF_CLUSTER_SEP)
            ref_name = (f"{settings.REF_CLUSTER_SEP.join(name_parts[0:-1])}"
                        f"{settings.SEQ_NAME_SEP}ref")
            if ref_name not in ref_names:
                ref_names.append(ref_name)

    rem_refs_from_fasta_params = []
    for fasta in fastas_paths:
        rem_refs_from_fasta_params.append((
            fasta,
            fastas_paths[fasta],
            ref_names,
            min_samples,
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


def rem_refs_from_fasta(fasta_in: Path, fasta_out: Path, ref_names: list, min_samples, overwrite):
    start = time.time()

    fasta_out_short = Path(*fasta_out.parts[-4:])
    if file_is_empty(fasta_in):
        return red(f"'{fasta_out_short}': FAILED removal of references, input file was empty")

    if overwrite is True or not fasta_out.exists():
        fasta_with_refs, fasta_without_refs = fasta_to_dict(fasta_in), {}
        for seq_name in fasta_with_refs:
            if seq_name not in ref_names:
                fasta_without_refs[seq_name] = dict(fasta_with_refs[seq_name])
        if num_samples(fasta_without_refs) >= min_samples:
            dict_to_fasta(fasta_without_refs, fasta_out)
            message = (
                f"'{fasta_out_short}': references removed [{elapsed_time(time.time() - start)}]"
            )
        else:
            message = red(
                f"'{fasta_out_short}': not saved (filtered FASTA had fewer than"
                f" {min_samples} samples) [{elapsed_time(time.time() - start)}]"
            )
    else:
        message = dim(f"'{fasta_out_short}': SKIPPED (output FASTA already exists)")

    return message


def clipkit(
    clipkit_path, clipkit_method, clipkit_gaps, fasta_in: Path, fasta_out: Path, min_data_per_column,
    min_coverage, min_samples, overwrite
):
    start = time.time()

    fasta_out_short = Path(*fasta_out.parts[-4:])
    if file_is_empty(fasta_in):
        return red(f"'{fasta_out_short}': FAILED alignment trimming, input file was empty")

    if overwrite is True or not fasta_out.exists():
        if min_data_per_column > 0:
            num_seqs = len(fasta_to_dict(fasta_in))
            if num_seqs == min_data_per_column:
                clipkit_gaps = 1
            else:
                clipkit_gaps = 1 - (min(num_seqs, min_data_per_column) / num_seqs)
        clipkit_cmd = [
            clipkit_path,
            f"{fasta_in}",
            "--output", f"{fasta_out}",
            "--mode", f"{clipkit_method}",
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
        # Initially added because IQ-TREE fails when a sample has empty sequence:
        # Calculate mean ungapped length, remove individual sequences under 'min_coverage'
        # as proportion of the mean of ungapped lengths
        if fasta_out.exists():
            if file_is_empty(fasta_out):
                message = red(f"'{fasta_out_short}': FAILED trimming, check ClipKIT's log")
                fasta_out.unlink()
            else:
                fasta_trimmed = fasta_to_dict(fasta_out)
                ungapped_lengths = []
                aln_length = 0
                for seq_name in fasta_trimmed:
                    ungapped = len(fasta_trimmed[seq_name]["sequence"].replace("-", ""))
                    aln_length = len(fasta_trimmed[seq_name]["sequence"])
                    if ungapped > 0:
                        ungapped_lengths.append(ungapped)
                if ungapped_lengths:
                    mean_ungapped = statistics.mean(ungapped_lengths)
                else:
                    mean_ungapped = aln_length
                seqs_to_remove = []
                for seq_name in fasta_trimmed:
                    ungapped = len(fasta_trimmed[seq_name]["sequence"].replace("-", ""))
                    if ungapped / mean_ungapped < min_coverage:
                        seqs_to_remove.append(seq_name)
                if seqs_to_remove:
                    for seq_name in seqs_to_remove:
                        del fasta_trimmed[seq_name]
                    dict_to_fasta(fasta_trimmed, fasta_out)
                if num_samples(fasta_trimmed) >= min_samples:
                    message = f"'{fasta_out_short}': trimmed [{elapsed_time(time.time() - start)}]"
                else:
                    fasta_out.unlink()
                    message = red(
                        f"'{fasta_out_short}': not saved (trimmed FASTA had fewer than"
                        f" {min_samples} samples) [{elapsed_time(time.time() - start)}]"
                    )
        else:
            message = red(f"'{fasta_out_short}': FAILED trimming, check ClipKIT's log")
    else:
        message = dim(f"'{fasta_out_short}': SKIPPED (output FASTA already exists)")

    return message


def compute_stats(shared_sam_stats, shared_aln_stats, fasta_path):

    start = time.time()

    short_path = Path(*fasta_path.parts[-5:])
    if file_is_empty(fasta_path):
        return red(f"'{short_path}': FAILED statistics calculation, input file was empty")
    aln_type = fasta_type(fasta_path)
    if aln_type == "invalid":
        return red(f"'{short_path}': FAILED statistics calculation, invalid input file")

    fasta_path_parts = fasta_path.parts

    aln_marker, aln_format= "", ""
    for m in settings.MARKER_DIRS:
        if fasta_path_parts[-3] == settings.MARKER_DIRS[m]:
            aln_marker = m
    for f in settings.FORMAT_DIRS:
        if fasta_path_parts[-2] == settings.FORMAT_DIRS[f]:
            aln_format = f

    fasta_dict = fasta_to_dict(fasta_path)

    coding = bool(fasta_path_parts[-2] == "02_NT")
    aln_stats = alignment_stats(fasta_dict, aln_type, coding)
    aln_tsv = [[
        f"{fasta_path}",                                    # [0] alignment file location
        f'{bool("untrimmed" not in fasta_path_parts[-5])}', # [1] alignment was trimmed
        f'{fasta_path_parts[-4].split("_")[1]}',            # [2] paralog filter applied
        f'{bool("w_refs" in fasta_path_parts[-4])}',        # [3] still with references
        aln_marker,                                         # [4] marker type
        aln_format,                                         # [5] alignment format
        fasta_path.stem,                                    # [6] locus name
        f'{aln_stats["sequences"]}',                        # [7] num sequences
        f'{aln_stats["samples"]}',                          # [8] num samples
        f'{aln_stats["avg_copies"]}',                       # [9] avg num copies
        f'{aln_stats["sites"]}',                            # [10] num sites
        f'{aln_stats["informative"]}',                      # [11] num informative sites
        f'{aln_stats["informativeness"]}',                  # [12] pct of informative sites
        f'{aln_stats["uninformative"]}',                    # [13] num constant + singleton sites
        f'{aln_stats["constant"]}',                         # [14] num constant sites
        f'{aln_stats["singleton"]}',                        # [15] num singleton sites
        f'{aln_stats["patterns"]}',                         # [16] num unique columns
        f'{aln_stats["avg_pid"]}',                          # [17] average pairwise identity
        f'{aln_stats["missingness"]}',                      # [18] pct of gaps and Ns or Xs
        f'{aln_stats["gc"]}',                               # [19] GC content in %
        f'{aln_stats["gc_codon_p1"]}',                      # [20] GC content of pos1 in codon in %
        f'{aln_stats["gc_codon_p2"]}',                      # [21] GC content of pos2 in codon in %
        f'{aln_stats["gc_codon_p3"]}',                      # [22] GC content of pos3 in codon in %
    ]]
    shared_aln_stats += ["\t".join(row)+"\n" for row in aln_tsv]

    sam_stats = sample_stats(fasta_dict, aln_type, coding)
    sam_tsv = []
    for sample in sam_stats:
        sam_tsv.append([
            sample,                                         # [0] sample name
            f'{" / ".join(fasta_path.parts[-5:-1])}',       # [1] stage / marker / format
            fasta_path.stem,                                # [2] locus name
            f'{sam_stats[sample]["len_total"]}',            # [3] alignment length
            f'{sam_stats[sample]["len_gapped"]}',           # [4] length - terminal gaps
            f'{sam_stats[sample]["len_ungapped"]}',         # [5] length - all gaps
            f'{sam_stats[sample]["ambigs"]}',               # [6] num ambiguities
            f'{sam_stats[sample]["gc"]}',                   # [7] GC content in %
            f'{sam_stats[sample]["gc_codon_p1"]}',          # [8] GC content of pos1 in codon in %
            f'{sam_stats[sample]["gc_codon_p2"]}',          # [9] GC content of pos2 in codon in %
            f'{sam_stats[sample]["gc_codon_p3"]}',          # [10] GC content of pos3 in codon in %
            f'{sam_stats[sample]["num_copies"]}',           # [11] num copies
        ])
    shared_sam_stats += ["\t".join(row)+"\n" for row in sam_tsv]

    message = f"'{short_path}': stats computed [{elapsed_time(time.time() - start)}]"

    return message


def write_aln_stats(out_dir, shared_aln_stats):
    if not shared_aln_stats:
        return None
    else:
        stats_tsv_file = Path(out_dir, "captus-assembly_align.alignments.tsv")
        with open(stats_tsv_file, "wt") as tsv_out:
            tsv_out.write("\t".join(["path",
                                     "trimmed",
                                     "paralog_filter",
                                     "with_refs",
                                     "marker_type",
                                     "format",
                                     "locus",
                                     "seqs",
                                     "samples",
                                     "avg_copies",
                                     "sites",
                                     "informative",
                                     "informativeness",
                                     "uninformative",
                                     "constant",
                                     "singleton",
                                     "patterns",
                                     "avg_pid",
                                     "missingness",
                                     "gc",
                                     "gc_codon_p1",
                                     "gc_codon_p2",
                                     "gc_codon_p3",]) + "\n")
            tsv_out.writelines(sorted(shared_aln_stats))
        return stats_tsv_file


def write_sam_stats(out_dir, shared_sam_stats):
    if not shared_sam_stats:
        return None
    else:
        stats_tsv_file = Path(out_dir, "captus-assembly_align.samples.tsv")
        with open(stats_tsv_file, "wt") as tsv_out:
            tsv_out.write("\t".join(["sample",
                                     "stage_marker_format",
                                     "locus",
                                     "len_total",
                                     "len_gapped",
                                     "len_ungapped",
                                     "ambigs",
                                     "gc",
                                     "gc_codon_p1",
                                     "gc_codon_p2",
                                     "gc_codon_p3",
                                     "num_copies",]) + "\n")
            tsv_out.writelines(sorted(shared_sam_stats))
        return stats_tsv_file


def write_astral_pro_seq_to_sam(out_dir, max_paralogs, ref_names, sample_names):
    if not sample_names:
        return None
    else:
        seq_to_sam = []
        for name in sorted(ref_names):
            seq_to_sam.append(f'{name}\t{name}\n')
        for name in sorted(sample_names):
            seq_to_sam.append(f'{name}\t{name}\n')
            for p in range(max_paralogs+1):
                seq_to_sam.append(f'{name}{settings.SEQ_NAME_SEP}{p:02}\t{name}\n')
        stats_tsv_file = Path(out_dir, settings.ASTRAL_PRO_EQ)
        with open(stats_tsv_file, "wt") as tsv_out:
            tsv_out.writelines(seq_to_sam)
        return stats_tsv_file
