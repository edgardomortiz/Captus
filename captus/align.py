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

import shutil
import subprocess
import time
from pathlib import Path
from threading import Timer

from . import log
from . import settings_assembly as settings
from .bioformats import dict_to_fasta, fasta_to_dict, is_fasta_nt, translate_fasta_dict
from .misc import (bold, dim, elapsed_time, format_dep_msg, is_dir_empty, mafft_path_version,
                   make_output_dir, quit_with_error, red, set_ram, set_threads,
                   tqdm_parallel_async_run, tqdm_parallel_async_write, tqdm_serial_run)
from .version import __version__


def align(full_command, args):

    captus_start = time.time()
    out_dir, out_dir_msg = make_output_dir(args.out)
    log.logger = log.Log(Path(args.out, "captus-assembly_align.log"), stdout_verbosity_level=1)

    mar = 25  # Margin for aligning parameters and values

    ################################################################################################
    ############################################################################### STARTING SECTION
    log.log_section_header("Starting Captus-assembly: Align", single_newline=False)
    log.log_explanation(
        "Welcome to the alignment step of Captus-assembly. In this step, Captus will collect all the"
        f" extracted markers across all samples in '{args.captus_extractions_dir}' and group them by"
        " marker. Then Captus will align each marker using MAFFT. If given, the reference loci are"
        " also included as alignment guides and removed after alignment", extra_empty_lines_after=0
    )
    log.log_explanation("For more information, please see https://github.com/edgardomortiz/Captus")

    log.log(f'{"Captus version":>{mar}}: {bold(f"v{__version__}")}')
    log.log(f'{"Command":>{mar}}: {bold(full_command)}')
    ram_B, ram_MB, ram_GB, ram_GB_total = set_ram(args.ram)
    log.log(f'{"Max. RAM":>{mar}}: {bold(f"{ram_GB:.1f}GB")} {dim(f"(out of {ram_GB_total:.1f}GB)")}')
    threads_max, threads_total = set_threads(args.threads)
    log.log(f'{"Max. Threads":>{mar}}: {bold(threads_max)} {dim(f"(out of {threads_total})")}')
    log.log("")
    log.log(f'{"Dependencies":>{mar}}:')
    mafft_path, mafft_version, mafft_status = mafft_path_version(args.mafft_path)
    log.log(format_dep_msg(f'{"MAFFT":>{mar}}: ', mafft_version, mafft_status))
    log.log("")
    log.log(f'{"Output directory":>{mar}}: {bold(out_dir)}')
    log.log(f'{"":>{mar}}  {dim(out_dir_msg)}')
    log.log("")

    skip_alignment = False
    if mafft_status == "not found":
        skip_alignment = True
        log.log(
            f"{bold('WARNING:')} MAFFT could not be found, the markers will be collected from"
            f" '{args.captus_extractions_dir}' but ehy will not be aligned. Please verify you have"
            " it installed or provide the full path to the program with '--mafft_path'"
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
    log.log("")
    ref_paths = prepare_refs(args.nuc_refs,
                             args.ptd_refs,
                             args.mit_refs,
                             args.dna_refs,
                             args.clr_refs,
                             mar,
                             args.overwrite)
    log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
    log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
    extracted_sample_dirs = find_extracted_sample_dirs(args.captus_extractions_dir)
    log.log(f'{"Samples to process":>{mar}}: {bold(len(extracted_sample_dirs))}')
    log.log("")
    log.log(make_output_dirtree(markers, formats, out_dir, "01_unaligned", mar))
    log.log("")
    collect_extracted_markers(markers, formats, extracted_sample_dirs, out_dir,
                              "01_unaligned", ref_paths, args.overwrite, args.show_less)
    log.log("")


    ################################################################################################
    ############################################################################## ALIGNMENT SECTION
    log.log_section_header("Marker alignment with MAFFT")
    log.log_explanation(
        "Now Captus will align all collected markers using MAFFT. If you added the references to be"
        " used as alignment guides Captus will produce a directory with alignments including the"
        " references and a separate one with the references removed. "
    )
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
    fastas_to_align = fastas_origs_dests(out_dir, "01_unaligned", "02_aligned_unfiltered")
    log.log(f'{"FASTA files to align":>{mar}}: {bold(len(fastas_to_align))}')
    log.log("")
    log.log(make_output_dirtree(markers, formats, out_dir, "02_aligned_unfiltered", mar))
    log.log("")

    mafft_params = []
    for fasta_orig in fastas_to_align:
        mafft_params.append((
            mafft_path,
            args.mafft_algorithm,
            threads_per_alignment,
            fasta_orig,
            fastas_to_align[fasta_orig],
            args.mafft_timeout,
            args.overwrite,
        ))

    # Internal switch between parallel asynchronous and serial run (False for debugging)
    run_async = True
    if run_async:
        tqdm_parallel_async_run(mafft, mafft_params,
                                "Aligning with MAFFT", "MAFFT alignment completed",
                                "alignment", concurrent, args.show_less)
    else:
        tqdm_serial_run(mafft, mafft_params,
                        "Aligning with MAFFT", "MAFFT alignment completed",
                        "alignment", args.show_less)
    log.log("")


    ################################################################################################
    ################################################################################ CLEANUP SECTION
    log.log_section_header("Paralog Filtering")
    log.log_explanation(
        "Now Captus will remove paralogs using the method(s) selected with '--filter_method'."
        " Afterwards, copies of the alignments without the reference sequences will also be created."
    )
    concurrent = threads_max
    filtering_refs = {}
    if args.filter_method.lower() in ["careful", "both"]:
        for marker in ref_paths:
            if marker != "CLR":
                if ref_paths[marker]["NT_path"]:
                    filtering_refs[marker] = {"path": ref_paths[marker]["NT_path"],
                                              "marker_dir": settings.MARKER_DIRS[marker],
                                              "format_dir": settings.FORMAT_DIRS["NT"]}
                elif marker != "DNA" and ref_paths[marker]["AA_path"]:
                    filtering_refs[marker] = {"path": ref_paths[marker]["AA_path"],
                                              "marker_dir": settings.MARKER_DIRS[marker],
                                              "format_dir": settings.FORMAT_DIRS["AA"]}
    filter_method = args.filter_method.lower()
    if not filtering_refs:
        if args.filter_method == "careful":
            filter_method = None
        elif args.filiter_method == "both":
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
        fastas_to_filter = fastas_origs_dests(out_dir,
                                              "02_aligned_unfiltered",
                                              "03_aligned_fast_filter")
        log.log("")
        log.log(bold(f'{"FAST paralog filtering":>{mar}}:'))
        log.log(f'{"FASTA files to process":>{mar}}: {bold(len(fastas_to_filter))}')
        log.log(make_output_dirtree(markers, formats, out_dir, "03_aligned_fast_filter", mar))
        log.log("")
        paralog_fast_filter(fastas_to_filter, args.overwrite, concurrent,
                            args.show_less, run_async=True)
        log.log("")

    if filter_method in ["careful", "both"]:
        fastas_to_filter = fastas_origs_dests(out_dir,
                                              "02_aligned_unfiltered",
                                              "04_aligned_careful_filter")
        log.log("")
        log.log(bold(f'{"CAREFUL paralog filtering":>{mar}}:'))
        log.log(f'{"FASTA files to process":>{mar}}: {bold(len(fastas_to_filter))}')
        log.log(make_output_dirtree(markers, formats, out_dir, "04_aligned_careful_filter", mar))
        log.log("")
        tmp_dir, _ = make_output_dir(Path(out_dir, "tmp"))
        paralog_careful_filter(tmp_dir, fastas_to_filter, filtering_refs, args.overwrite,
                               concurrent, args.show_less, run_async=True)
        paralog_stats_tsv = collect_paralog_stats(out_dir, tmp_dir)
        log.log("")
        log.log(f'{"Paralog statistics":>{mar}}: {bold(paralog_stats_tsv)}')
        shutil.rmtree(tmp_dir, ignore_errors=True)
        log.log("")


    ################################################################################################
    ################################################################################ CLEANUP SECTION
    log.log_section_header("Reference Sequences Removal and File Cleanup")
    log.log_explanation(
        "Now Captus will create copies of the alignnments that will not include the reference"
        " sequences used as alignment guide and for paralog filtering. Empty directories will be"
        " removed as well as MAFFT logs unless the flag '--keep_all' is enabled."
    )
    try:
        remove_references = bool(any([path for marker in ref_paths
                                           for path in ref_paths[marker].values()]))
    except TypeError:
        remove_references = False

    if remove_references:
        if Path(out_dir, "02_aligned_unfiltered").exists():
            fastas_to_rem_refs = fastas_origs_dests(out_dir,
                                                    "02_aligned_unfiltered",
                                                    "05_aligned_unfiltered_no_refs")
            log.log("")
            log.log(f'{"FASTA files to process":>{mar}}: {bold(len(fastas_to_rem_refs))}')
            log.log(make_output_dirtree(markers, formats, out_dir,
                                        "05_aligned_unfiltered_no_refs", mar))
            log.log("")
            rem_refs(ref_paths, fastas_to_rem_refs, args.overwrite,
                     concurrent, args.show_less, run_async=True)
            log.log("")

        if Path(out_dir, "03_aligned_fast_filter").exists():
            fastas_to_rem_refs = fastas_origs_dests(out_dir,
                                                    "03_aligned_fast_filter",
                                                    "06_aligned_fast_filter_no_refs")
            log.log("")
            log.log(f'{"FASTA files to process":>{mar}}: {bold(len(fastas_to_rem_refs))}')
            log.log(make_output_dirtree(markers, formats, out_dir,
                                        "06_aligned_fast_filter_no_refs", mar))
            log.log("")
            rem_refs(ref_paths, fastas_to_rem_refs, args.overwrite,
                     concurrent, args.show_less, run_async=True)
            log.log("")

        if Path(out_dir, "04_aligned_careful_filter").exists():
            fastas_to_rem_refs = fastas_origs_dests(out_dir,
                                                    "04_aligned_careful_filter",
                                                    "07_aligned_careful_filter_no_refs")
            log.log("")
            log.log(f'{"FASTA files to process":>{mar}}: {bold(len(fastas_to_rem_refs))}')
            log.log(make_output_dirtree(markers, formats, out_dir,
                                        "07_aligned_careful_filter_no_refs", mar))
            log.log("")
            rem_refs(ref_paths, fastas_to_rem_refs, args.overwrite,
                     concurrent, args.show_less, run_async=True)
            log.log("")

    start = time.time()
    reclaimed_bytes = 0
    if not args.keep_all:
        files_to_delete = list(out_dir.resolve().rglob("*.mafft.log"))
        for del_file in files_to_delete:
            reclaimed_bytes += del_file.stat().st_size
            del_file.unlink()
        log.log("")
        log.log(
            f'A total of {len(files_to_delete)} files'
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
    log.log_section_header(
        "Captus-assembly: Align -> successfully completed"
        f" [{elapsed_time(time.time() - captus_start)}]"
    )
    log.log("")


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
    extracted_sample_dirs = list(captus_extractions_dir.resolve().rglob("*__captus-ext"))
    if not extracted_sample_dirs:
        quit_with_error(
            f"Captus did not find valid sample directories within '{captus_extractions_dir}' please"
            " provide a valid directory with '--captus_extractions_dir"
        )
    return extracted_sample_dirs


def prepare_refs(nuc_refs, ptd_refs, mit_refs, dna_refs, clr_refs, margin, overwrite):
    def get_protref_paths(refs, marker: str, margin, overwrite):
        aa_path = nt_path = ""
        aa_msg = nt_msg = dim("not used")
        if refs is None:
            return aa_path, aa_msg, nt_path, nt_msg
        refs = refs.split(",")
        if len(refs) == 1:
            refs = refs[0]
            if f"{refs}".lower() in settings.PROT_REFS[marker]:
                aa_path = settings.PROT_REFS[marker][f"{refs}".lower()]["AA"]
                aa_msg = f'{bold(refs)}\n{" " * (margin + 2)}{aa_path}'
                nt_path = settings.PROT_REFS[marker][f"{refs}".lower()]["NT"]
                nt_msg = f'{bold(refs)}\n{" " * (margin + 2)}{nt_path}'
            elif Path(refs).is_file() and is_fasta_nt(refs) is False:
                aa_path = Path(refs).resolve()
                aa_msg = bold(aa_path)
                nt_path = None
                nt_msg = "not provided"
            elif Path(refs).is_file() and is_fasta_nt(refs) is True:
                nt_path = Path(refs).resolve()
                nt_msg = bold(nt_path)
                aa_path, aa_msg = translate_ref(nt_path,
                                                settings.DEFAULT_GENETIC_CODES[marker]["id"],
                                                margin,
                                                overwrite)
            elif Path(refs).is_file() and is_fasta_nt(refs) == "not a FASTA":
                aa_msg = nt_msg = red("not a valid FASTA")
            else:
                aa_msg = nt_msg = red("file not found")
        elif len(refs) == 2:
            try:
                transtable = int(refs[1])
            except ValueError:
                transtable = False
            if bool(transtable):
                if Path(refs[0]).is_file() and is_fasta_nt(refs[0]) is True:
                    nt_path = Path(refs[0]).resolve()
                    nt_msg = bold(nt_path)
                    aa_path, aa_msg = translate_ref(nt_path, transtable, margin, overwrite)
                elif Path(refs[0]).is_file() and is_fasta_nt(refs[0]) is False:
                    aa_path = Path(refs[0]).resolve()
                    aa_msg = (
                        f'{bold(aa_path)}\n{" " * (margin + 2)}'
                        f'{dim("(This FASTA contained aminoacids, the Genetic Code was ignored)")}'
                    )
                elif Path(refs[0]).is_file() and is_fasta_nt(refs[0]) == "not a FASTA":
                    aa_msg = nt_msg = red("not a valid FASTA")
                else:
                    aa_msg = nt_msg = red("file not found")
            else:
                if all([Path(refs[0]).is_file(), is_fasta_nt(refs[0]) is False,
                        Path(refs[1]).is_file(), is_fasta_nt(refs[1]) is True]):
                    aa_path, nt_path = Path(refs[0]).resolve(), Path(refs[1]).resolve()
                    aa_msg, nt_msg = bold(aa_path), bold(nt_path)
                else:
                    aa_msg = red(f"'{refs[0]}' does not contain aminoacids, or it was not found.")
                    nt_msg = red(
                        f"'{refs[1]}' file does not contain nucleotides, or it was not found."
                    )
        return aa_path, aa_msg, nt_path, nt_msg

    def translate_ref(nt_path, transtable, margin, overwrite):
        start = time.time()
        suffix = settings.TRANSLATED_REF_SUFFIX
        if f"{nt_path}".endswith(".gz"):
            aa_path = Path(Path(nt_path).resolve().parent,
                           f'{Path(nt_path.replace(".gz", "")).stem}{suffix}')
        else:
            aa_path = Path(Path(nt_path).resolve().parent, f'{Path(nt_path).stem}{suffix}')
        if not aa_path.exists() or overwrite:
            nt_translated = translate_fasta_dict(fasta_to_dict(nt_path, ordered=True), transtable)
            dict_to_fasta(nt_translated, aa_path, wrap=80)
            aa_msg = (
                f'{bold(aa_path)}\n{" " * (margin + 2)}'
                f'{dim(f"(Translated using Genetic Code {transtable}")}'
                f'{dim(f" [{elapsed_time(time.time() - start)}])")}'
            )
        elif aa_path.exists():
            aa_msg = (
                f'{bold(aa_path)}\n{" " * (margin + 2)}'
                f'{dim(f"(Captus-translated FASTA found in the same location)")}'
            )
        return aa_path, aa_msg

    def get_dnaref_path(ref, margin):
        nt_path = ""
        nt_msg = dim("not used")
        if ref is None:
            return nt_path, nt_msg
        elif f"{ref}".lower() in settings.DNA_REFS:
            nt_path = settings.DNA_REFS[f"{ref}".lower()]
            nt_msg = f'{bold(ref)}\n{" " * (margin + 2)}{dim(nt_path)}'
        elif Path(ref).is_file() and is_fasta_nt(ref) is True:
            nt_path = Path(ref).resolve()
            nt_msg = bold(nt_path)
        else:
            nt_msg = red("not a valid FASTA")
        return nt_path, nt_msg

    if any([nuc_refs, ptd_refs, mit_refs, dna_refs, clr_refs]):
        log.log(bold(f'{"Reference datasets":>{margin}}:'))
        nuc_aa_path, nuc_aa_msg, nuc_nt_path, nuc_nt_msg = get_protref_paths(nuc_refs,
                                                                             "NUC",
                                                                             margin,
                                                                             overwrite)
        if bool(nuc_refs):
            log.log(f'{"Nuclear proteins AA":>{margin}}: {nuc_aa_msg}')
            log.log(f'{"Nuclear proteins NT":>{margin}}: {nuc_nt_msg}')
        ptd_aa_path, ptd_aa_msg, ptd_nt_path, ptd_nt_msg = get_protref_paths(ptd_refs,
                                                                             "PTD",
                                                                             margin,
                                                                             overwrite)
        if bool(ptd_refs):
            log.log(f'{"Plastidial proteins AA":>{margin}}: {ptd_aa_msg}')
            log.log(f'{"Plastidial proteins NT":>{margin}}: {ptd_nt_msg}')
        mit_aa_path, mit_aa_msg, mit_nt_path, mit_nt_msg = get_protref_paths(mit_refs,
                                                                             "MIT",
                                                                             margin,
                                                                             overwrite)
        if bool(mit_refs):
            log.log(f'{"Mitochondrial proteins AA":>{margin}}: {mit_aa_msg}')
            log.log(f'{"Mitochondrial proteins NT":>{margin}}: {mit_nt_msg}')
        dna_path, dna_msg = get_dnaref_path(dna_refs, margin)
        if bool(dna_refs):
            log.log(f'{"Miscellaneous DNA":>{margin}}: {dna_msg}')
        clr_path, clr_msg = get_dnaref_path(clr_refs, margin)
        if bool(clr_refs):
            log.log(f'{"Cluster-derived DNA":>{margin}}: {clr_msg}')
        log.log("")

        ref_paths = {
            "NUC": {"AA_path": nuc_aa_path, "NT_path": nuc_nt_path},
            "PTD": {"AA_path": ptd_aa_path, "NT_path": ptd_nt_path},
            "MIT": {"AA_path": mit_aa_path, "NT_path": mit_nt_path},
            "DNA": {"NT_path": dna_path},
            "CLR": {"NT_path": clr_path}
        }
        return ref_paths
    else:
        return None


def make_output_dirtree(markers, formats, out_dir, base_dir, margin):
    """
    Create directory structure to receive extracted markers from all the samples in
    `captus_extractions_dir`, print directory tree created with explanations of the intended content
    of each directory
    """
    base_tree = Path(out_dir, base_dir).resolve()
    markers = sorted([(settings.MARKER_DIRS[m], m) for m in markers.upper().split(",")])
    formats = sorted([(settings.FORMAT_DIRS[f], f) for f in formats.upper().split(",")])
    corner = "\u2514\u2500\u2500 "
    branch = "\u251C\u2500\u2500 "
    space = "    "
    vline = "\u2502   "
    symbols = ["", corner]
    dirs = [Path(out_dir).resolve(), base_dir]
    margins = [f'{"Output directory tree":>{margin}}: ', f'{"":>{margin}}  ']
    for m in markers:
        for f in formats:
            if (m[1], f[1]) in settings.VALID_MARKER_FORMAT_COMBO:
                make_output_dir(Path(base_tree, m[0], f[0]))
                if not m[0] in dirs:
                    if m[0] == markers[-1][0]:
                        parent, fork = space, corner
                    else:
                        parent, fork = vline, branch
                    dirs.append(m[0])
                    symbols.append(f"{space}{fork}")
                    margins.append(f'{"":>{margin}}  ')
                dirs.append(f[0])
                if symbols[-1] == f'{space}{parent}{corner}':
                    symbols[-1] = f'{space}{parent}{branch}'
                symbols.append(f'{space}{parent}{corner}')
                margins.append(f'{"":>{margin}}  ')
    return "\n".join([f"{m}{bold(s)}{bold(d)}" for m, s, d in zip(margins, symbols, dirs)])


def collect_extracted_markers(
    markers, formats, extracted_sample_dirs, out_dir, base_dir, ref_paths, overwrite, show_less
):
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
        elif is_dir_empty(o) is False and overwrite is False:
            quit_with_error(
                "Output directories are not empty, to replace previous results enable '--overwrite'"
                " or specify a different output directory with '--out'"
            )
    collect_sample_markers_params = []
    for sample_dir in extracted_sample_dirs:
        collect_sample_markers_params.append([
            sample_dir,
            source_files,
            out_dirs,
        ])
    tqdm_serial_run(collect_sample_markers, collect_sample_markers_params,
                    "Collecting extracted markers", "Collection of extracted markers finished",
                    "sample", show_less)
    if ref_paths is not None:
        refs = [
            ref_paths["NUC"]["AA_path"], ref_paths["NUC"]["NT_path"],
            ref_paths["PTD"]["AA_path"], ref_paths["PTD"]["NT_path"],
            ref_paths["MIT"]["AA_path"], ref_paths["MIT"]["NT_path"],
            ref_paths["DNA"]["NT_path"], ref_paths["DNA"]["NT_path"],
            ref_paths["CLR"]["NT_path"], ref_paths["CLR"]["NT_path"],
        ]
        mrks = ["NUC", "NUC", "PTD", "PTD", "MIT", "MIT", "DNA", "DNA", "CLR", "CLR"]
        fmts = ["AA", "NT", "AA", "NT", "AA", "NT", "MA", "MF", "MA", "MF"]
        exts = [".faa", ".fna", ".faa", ".fna", ".faa", ".fna", ".fna", ".fna", ".fna", ".fna"]
        add_refs_params = []
        for r, m, f, e in zip(refs, mrks, fmts, exts):
            if all([r, m in markers.upper().split(","), f in formats.upper().split(",")]):
                add_refs_params.append((
                    r,
                    Path(out_dir, base_dir, settings.MARKER_DIRS[m], settings.FORMAT_DIRS[f]),
                    e
                ))
        if bool(add_refs_params):
            log.log("")
            tqdm_serial_run(add_refs, add_refs_params,
                            "Adding reference markers", "Addition of reference markers finished",
                            "reference", show_less)


def collect_sample_markers(sample_dir, source_files, out_dirs):
    start = time.time()
    sample_name = sample_dir.parts[-1].replace("__captus-ext", "")
    markers_collected = []
    fastas_created = []
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
                markers_collected.append(marker_name)
                fasta_out = Path(destination, f"{marker_name}{source.suffix}")
                fastas_created.append(fasta_out)
                dict_to_fasta({seq_name: dict(fasta_in[seq_name_full])}, fasta_out, append=True)
    message = (
        f"'{sample_name}': {len(set(fastas_created))} FASTA files created for"
        f" {len(set(markers_collected))} collected markers [{elapsed_time(time.time() - start)}]"
    )
    return message


def add_refs(ref_path, dest_dir, extension):
    start = time.time()
    ref_fasta = fasta_to_dict(ref_path)
    markers_in_ref = []
    fastas_found = []
    for seq in ref_fasta:
        name_parts = seq.split(settings.REFERENCE_CLUSTER_SEPARATOR)
        seq_name = f"{settings.REFERENCE_CLUSTER_SEPARATOR.join(name_parts[0:-1])}|ref"
        marker_name = name_parts[-1]
        markers_in_ref.append(marker_name)
        fasta_out = Path(dest_dir, f"{marker_name}{extension}")
        if fasta_out.exists():
            with open(fasta_out, "rt") as fasta_to_check:
                for line in fasta_to_check:
                    if seq in line:
                        fastas_found.append(marker_name)
                        dict_to_fasta({seq_name: dict(ref_fasta[seq])}, fasta_out, append=True)
                        break
    if bool(fastas_found):
        message = (
            f"'{ref_path.name}': {len(set(fastas_found))} markers (out of"
            f" {len(set(markers_in_ref))}) were appended to collected FASTA files"
            f" [{elapsed_time(time.time() - start)}]"
        )
    else:
        message = red(
            f"'{ref_path.name}': did not find collected FASTA files to which append references"
            f" [{elapsed_time(time.time() - start)}]"
        )
    return message


def adjust_mafft_concurrency(concurrent, threads_max):
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
    for path in list(Path(dir_path, orig_base_dir).rglob("*.f[an]a")):
        origin = path.resolve()
        destination = Path(*[p if p != orig_base_dir else dest_base_dir for p in origin.parts])
        fastas_to_process[origin] = destination
    return fastas_to_process


def mafft(
    mafft_path, mafft_algorithm, threads, fasta_in: Path, fasta_out: Path, mafft_timeout, overwrite
):
    start = time.time()
    fasta_out_short = Path(*list(fasta_out.parts)[-3:])
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
            with open(mafft_log_file, "w") as mafft_log:
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


def rem_refs(ref_paths, fastas_paths, overwrite, concurrent, show_less, run_async):
    ref_names = []
    for marker in ref_paths:
        for fasta_path in ref_paths[marker]:
            if bool(ref_paths[marker][fasta_path]):
                for seq_name in fasta_to_dict(ref_paths[marker][fasta_path]):
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
    if run_async:
        tqdm_parallel_async_run(rem_refs_from_fasta, rem_refs_from_fasta_params,
                                "Removing references from alignments", "References removal completed",
                                "alignment", concurrent, show_less)
    else:
        tqdm_serial_run(rem_refs_from_fasta, rem_refs_from_fasta_params,
                        "Removing references from alignments", "Removal of references completed",
                        "alignment", show_less)


def rem_refs_from_fasta(fasta_source: Path, fasta_dest: Path, ref_names: list, overwrite):
    start = time.time()
    fasta_dest_short = Path(*list(fasta_dest.parts)[-3:])
    fasta_with_refs, fasta_without_refs = fasta_to_dict(fasta_source), {}
    if overwrite is True or not fasta_dest.exists():
        for seq_name in fasta_with_refs:
            if seq_name not in ref_names:
                fasta_without_refs[seq_name] = dict(fasta_with_refs[seq_name])
        dict_to_fasta(fasta_without_refs, fasta_dest)
        message = f"'{fasta_dest_short}': references removed [{elapsed_time(time.time() - start)}]"
    else:
        message = dim(f"'{fasta_dest_short}': skipped (output FASTA already exists)")
    return message


def paralog_fast_filter(fastas_paths, overwrite, concurrent, show_less, run_async):
    filter_paralogs_fast_params = []
    for fasta in fastas_paths:
        filter_paralogs_fast_params.append((
            fasta,
            fastas_paths[fasta],
            overwrite,
        ))
    if run_async:
        tqdm_parallel_async_run(filter_paralogs_fast, filter_paralogs_fast_params,
                                "Removing potential paralogs from alignments",
                                "Potential paralog removal completed",
                                "alignment", concurrent, show_less)
    else:
        tqdm_serial_run(filter_paralogs_fast, filter_paralogs_fast_params,
                        "Removing potential paralogs from alignments",
                        "Potential paralog removal completed",
                        "alignment", show_less)


def filter_paralogs_fast(fasta_source, fasta_dest, overwrite):
    start = time.time()
    fasta_dest_short = Path(*list(fasta_dest.parts)[-3:])
    fasta_with_paralogs, fasta_without_paralogs = fasta_to_dict(fasta_source), {}
    if overwrite is True or not fasta_dest.exists():
        for seq_name in fasta_with_paralogs:
            if ("hit=00" in fasta_with_paralogs[seq_name]["description"] or
                "|ref" in seq_name):
                seq_name_out = seq_name.replace("|00", "")
                fasta_without_paralogs[seq_name_out] = dict(fasta_with_paralogs[seq_name])
        dict_to_fasta(fasta_without_paralogs, fasta_dest)
        message = f"'{fasta_dest_short}': paralogs removed [{elapsed_time(time.time() - start)}]"
    else:
        message = dim(f"'{fasta_dest_short}': skipped (output FASTA already exists)")
    return message


def paralog_careful_filter(
    tmp_dir, fastas_paths, filtering_refs, overwrite, concurrent, show_less, run_async
):
    fastas = {}
    for marker in filtering_refs:
        for fasta in fastas_paths:
            if (fasta.parts[-2] == filtering_refs[marker]["format_dir"] and
                fasta.parts[-3] == filtering_refs[marker]["marker_dir"]):
                fastas[fasta.stem] = fasta

    filter_paralogs_careful_params = []
    for marker_name in fastas:
        fastas_marker = {}
        for fasta in fastas_paths:
            if fasta.stem == marker_name and fasta.parts[-3] == fastas[marker_name].parts[-3]:
                fastas_marker[fasta] = fastas_paths[fasta]
        filter_paralogs_careful_params.append((
            tmp_dir,
            fastas[marker_name],
            fastas_marker,
            overwrite
        ))

    if run_async:
        tqdm_parallel_async_run(filter_paralogs_careful, filter_paralogs_careful_params,
                                "Removing potential paralogs from markers",
                                "Potential paralog removal completed",
                                "marker", concurrent, show_less)
    else:
        tqdm_serial_run(filter_paralogs_careful, filter_paralogs_careful_params,
                        "Removing potential paralogs from markers",
                        "Potential paralog removal completed",
                        "marker", show_less)


def filter_paralogs_careful(tmp_dir, fasta_model, fastas_paths, overwrite):

    def calc_pid(s1, s2):
        s1_start, s2_start = len(s1) - len(s1.lstrip("-")), len(s2) - len(s2.lstrip("-"))
        s1_end, s2_end = len(s1.rstrip("-")), len(s2.rstrip("-"))
        overlap_length = min(s1_end, s2_end) - max(s1_start, s2_start)
        matches = 0
        if overlap_length > 0:
            for pos in range(max(s1_start, s2_start), min(s1_end, s2_end)):
                if s1[pos].upper() == s2[pos].upper():
                    matches += 1
            return matches * 100 / overlap_length
        else:
            return 0.00

    start = time.time()

    aln = fasta_to_dict(fasta_model, ordered=True)
    refs = {}
    for seq in aln:
        if "query=" in aln[seq]["description"] and "hit=00" in aln[seq]["description"]:
            ref = aln[seq]["description"].split("query=")[1].split(":")[0]
            if ref in refs:
                refs[ref] += 1
            else:
                refs[ref] = 1
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
            if sample_name in samples_with_paralogs:
                samples_with_paralogs[sample_name][seq] = calc_pid(best_ref_seq, aln[seq]["sequence"])
            else:
                samples_with_paralogs[sample_name] = {seq: calc_pid(best_ref_seq, aln[seq]["sequence"])}
            tsv.append([
                fasta_model.parts[-3][-3:],                   # marker type
                fasta_model.parts[-2][-2:],                   # format used for filtering
                fasta_model.stem,                             # locus name
                sample_name,                                  # sample name
                seq,                                          # sequence name
                hit_num,                                      # hit ranking
                best_ref_full_name,                           # reference name
                f"{samples_with_paralogs[sample_name][seq]}", # identity to reference
                f"{False}",                                   # accepted as ortholog
            ])
    for sample in samples_with_paralogs:
        accepted.append(max(samples_with_paralogs[sample], key=samples_with_paralogs[sample].get))

    for row in tsv:
        if row[4] in accepted:
            row[8] = f"{True}"

    with open(Path(tmp_dir, f"{best_ref}.tsv"), "wt") as tsv_out:
        for row in tsv:
            tsv_out.write("\t".join(row) + "\n")

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
            dict_to_fasta(fasta_without_paralogs, fastas_paths[fasta])
    message = (
        f"'{fasta_model.parts[-3]}-{fasta_model.stem}': paralogs removed from"
        f" {len(fastas_paths)} files [{elapsed_time(time.time() - start)}]"
    )
    return message


def collect_paralog_stats(out_dir, tmp_dir):
    tsv_files = sorted(list(Path(tmp_dir).resolve().glob("*.tsv")))
    if not tsv_files:
        return red("No paralog filtering statistics files found...")
    else:
        stats_tsv_file = Path(out_dir, "captus-assembly_paralog_filtering.tsv")
        with open(stats_tsv_file, "wt") as tsv_out:
            tsv_out.write("\t".join(["marker_type",
                                     "format_filtered",
                                     "locus",
                                     "sample",
                                     "sequence",
                                     "hit",
                                     "ref",
                                     "identity",
                                     "accepted"]) + "\n")
            for file in tsv_files:
                with open(file, "rt") as tsv_in:
                    for line in tsv_in:
                        tsv_out.write(line)
        return stats_tsv_file

