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

import platform
import shutil
import time
from pathlib import Path

from tqdm import tqdm

from . import log, settings
from .misc import (
    bold,
    dim,
    dir_is_empty,
    elapsed_time,
    file_is_empty,
    format_dep_msg,
    make_output_dir,
    python_library_check,
    quit_with_error,
    red,
    successful_exit,
)
from .version import __version__


def select(full_command, args):
    captus_start = time.time()
    out_dir, out_dir_msg = make_output_dir(args.out)
    log.logger = log.Log(Path(out_dir, "captus-select.log"), stdout_verbosity_level=1)

    mar = 33  # Margin for aligning parameters and values

    ################################################################################################
    ############################################################################### STARTING SECTION
    log.log_section_header("Starting Captus-design: SELECT", single_newline=False)
    log.log_explanation(
        "Welcome to the cluster selection step of Captus-design. In this step, Captus will make a"
        " copy of the clusters aligned and curated prudced by the previous step. The clusters can"
        " be selected according to several parameters, e.g., percentage of informative sites or total"
        " alignment length.",
        extra_empty_lines_after=0,
    )
    log.log_explanation("For more information, please see https://github.com/edgardomortiz/Captus")

    log.log(f"{'Captus version':>{mar}}: {bold(f'v{__version__}')}")
    log.log(f"{'Command':>{mar}}: {bold(full_command)}")
    log.log(f"{'OS':>{mar}}: {bold(platform.platform())}")
    log.log(f"{'Host':>{mar}}: {bold(platform.node())}")
    tsv_comment = f"#Captus v{__version__}\n#Command: {full_command}\n"
    log.log("")

    log.log(f"{'Python libraries':>{mar}}:")
    numpy_found, numpy_version, numpy_status = python_library_check("numpy")
    pandas_found, pandas_version, pandas_status = python_library_check("pandas")
    plotly_found, plotly_version, plotly_status = python_library_check("plotly")
    log.log(format_dep_msg(f"{'numpy':>{mar}}: ", numpy_version, numpy_status))
    log.log(format_dep_msg(f"{'pandas':>{mar}}: ", pandas_version, pandas_status))
    log.log(format_dep_msg(f"{'plotly':>{mar}}: ", plotly_version, plotly_status))
    log.log("")

    log.log(f"{'Overwrite files':>{mar}}: {bold(args.overwrite)}")
    log.log(f"{'Keep all files':>{mar}}: {bold(args.keep_all)}")
    log.log("")

    log.log(f"{'Output directory':>{mar}}: {bold(out_dir)}")
    log.log(f"{'':>{mar}}  {dim(out_dir_msg)}")
    log.log("")

    ################################################################################################
    ################################################################################ CLEANUP SECTION
    log.log_section_header("Loci selection")
    log.log_explanation(
        "Now Captus will copy the cluster alignments that fulfill the parameters chosen. It is"
        " recommended to enable '--dry_run' to test the parameters first."
    )

    aln_stats = load_aln_stats_tsv(args.captus_clusters_dir)
    aln_stats_filtered = filter_loci(
        aln_stats,
        args.avg_copies,
        args.length,
        args.pairwise_identity,
        args.gc_content,
        args.informative_sites,
        args.informativeness,
        args.missingness,
        args.num_sequences,
        args.num_samples,
        args.num_focal_species,
        args.num_outgroup_species,
        args.num_addon_samples,
        args.num_species,
        args.num_genera,
        args.cds_len,
        args.len_long_exons_retained,
        args.len_short_exons_retained,
        args.perc_total_cds_retained,
        args.perc_long_exons_retained,
        args.perc_short_exons_retained,
        mar,
    )
    log.log("")

    if args.dry_run:
        log.log("")
        log.log("")
        log.log(red("Disable '--dry_run' to copy the selected loci alignments..."))
    else:
        copy_loci(aln_stats_filtered, out_dir, args.overwrite, args.show_more)
        log.log("")
        start = time.time()
        log.log("")
        log.log_explanation("Writing selected alignments statistics...")
        aln_stats_tsv = write_aln_stats(out_dir, tsv_comment, aln_stats_filtered)
        if aln_stats_tsv:
            log.log(f"{'Alignment statistics':>{mar}}: {bold(aln_stats_tsv)}")
            log.log(f"{'':>{mar}}  {dim(f'File saved in {elapsed_time(time.time() - start)}')}")
            log.log("")
            if all([numpy_found, pandas_found, plotly_found]):
                from .report import build_design_report

                log.log("")
                log.log_explanation("Generating Alignment Statistics report...")
                aln_html_report, aln_html_msg = build_design_report(out_dir, aln_stats_tsv, "select")
                log.log(f"{'Alignment report':>{mar}}: {bold(aln_html_report)}")
                log.log(f"{'':>{mar}}  {dim(aln_html_msg)}")
            else:
                log.log(
                    f"{bold('WARNING:')} Captus uses 'numpy', 'pandas', and 'plotly' to generate an"
                    " HTML report based on the marker recovery statistics. At least one of these"
                    " libraries could not be found, please verify these libraries are installed and"
                    " available."
                )

    log.log("")

    ################################################################################################
    ################################################################################# ENDING SECTION
    successful_exit(
        f"Captus-design: SELECT -> successfully completed [{elapsed_time(time.time() - captus_start)}]"
    )


def load_aln_stats_tsv(clusters_dir: Path):
    start = time.time()
    aln_stats_tsv_path = Path(clusters_dir, "captus-cluster_alignments.tsv")
    if aln_stats_tsv_path.exists():
        aln_stats = {}
        with open(aln_stats_tsv_path, "rt") as stats:
            for line in stats:
                if not line.startswith("path") and not line.startswith("#"):
                    record = line.strip().split()
                    aln_stats[record[1]] = {
                        "path": Path(record[0]),
                        "copies": int(record[2]),
                        "avg_copies": float(record[3]),
                        "length": int(record[4]),
                        "gc_content": float(record[5]),
                        "avg_pid": float(record[6]),
                        "informative_sites": int(record[7]),
                        "informativeness": float(record[8]),
                        "missingness": float(record[9]),
                        "num_sequences": int(record[10]),
                        "num_samples": int(record[11]),
                        "num_focal_species": int(record[12]),
                        "num_outgroup_species": int(record[13]),
                        "num_addon_samples": int(record[14]),
                        "num_species": int(record[15]),
                        "num_genera": int(record[16]),
                        "cds_id": record[17],
                        "cds_len": record[18],
                        "len_long_exons_retained": record[19],
                        "len_short_exons_retained": record[20],
                        "perc_exons_retained": record[21],
                        "perc_long_exons_retained": record[22],
                        "perc_short_exons_retained": record[23],
                    }
        log.log(bold(f"Data from {aln_stats_tsv_path} loaded in {elapsed_time(time.time() - start)}"))
        return aln_stats
    else:
        quit_with_error(f"Alignment statistics file '{aln_stats_tsv_path}' not found.")


def filter_loci(
    aln_stats: dict,
    avg_copies: float,
    length: str,
    pairwise_identity: str,
    gc_content: str,
    informative_sites: str,
    informativeness: str,
    missingness: str,
    num_sequences: int,
    num_samples: int,
    num_focal_species: int,
    num_outgroup_species: int,
    num_addon_samples: int,
    num_species: int,
    num_genera: int,
    cds_len: str,
    len_long_exons_retained: str,
    len_short_exons_retained: str,
    perc_exons_retained: str,
    perc_long_exons_retained: str,
    perc_short_exons_retained: str,
    mar: int,
):
    log.log(f"{'':>{mar}}  {len(aln_stats)} loci found")
    log.log("")

    log.log(f"{'Avg. number of copies (--cop)':>{mar}}: {bold(avg_copies)}")
    try:
        min_par, max_par = (float(par) for par in avg_copies.split(","))
    except ValueError:
        quit_with_error(
            "Average number of copies must be given as two decimals separated by a comma without spaces"
        )
    aln_stats = {k: aln_stats[k] for k in aln_stats if min_par <= aln_stats[k]["avg_copies"] <= max_par}
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    log.log(f"{'Length (--len)':>{mar}}: {bold(length)} bp")
    try:
        min_par, max_par = (int(par) for par in length.split(","))
    except ValueError:
        quit_with_error("Length range must be given as two integers separated by a comma without spaces")
    aln_stats = {k: aln_stats[k] for k in aln_stats if min_par <= aln_stats[k]["length"] <= max_par}
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    log.log(f"{'Avg. pairwise identity (--pid)':>{mar}}: {bold(pairwise_identity)} %")
    try:
        min_par, max_par = (float(par) for par in pairwise_identity.split(","))
    except ValueError:
        quit_with_error(
            "Average pairwise identity range must be given as"
            " two decimals separated by a comma without spaces"
        )
    aln_stats = {k: aln_stats[k] for k in aln_stats if min_par <= aln_stats[k]["avg_pid"] <= max_par}
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    log.log(f"{'GC content (--gc)':>{mar}}: {bold(gc_content)} %")
    try:
        min_par, max_par = (float(par) for par in gc_content.split(","))
    except ValueError:
        quit_with_error(
            "GC content percentage range must be given as two"
            " decimals separated by a comma without spaces"
        )
    aln_stats = {k: aln_stats[k] for k in aln_stats if min_par <= aln_stats[k]["gc_content"] <= max_par}
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    log.log(f"{'Informative sites (--pis)':>{mar}}: {bold(informative_sites)}")
    try:
        min_par, max_par = (int(par) for par in informative_sites.split(","))
    except ValueError:
        quit_with_error(
            "Informative sites range must be given as two integers separated by a comma without spaces"
        )
    aln_stats = {
        k: aln_stats[k] for k in aln_stats if min_par <= aln_stats[k]["informative_sites"] <= max_par
    }
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    log.log(f"{'Informativeness (--inf)':>{mar}}: {bold(informativeness)} %")
    try:
        min_par, max_par = (float(par) for par in informativeness.split(","))
    except ValueError:
        quit_with_error(
            "Informativeness percentage range must be given as two"
            " decimals separated by a comma without spaces"
        )
    aln_stats = {
        k: aln_stats[k] for k in aln_stats if min_par <= aln_stats[k]["informativeness"] <= max_par
    }
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    log.log(f"{'Missingness (--mis)':>{mar}}: {bold(missingness)} %")
    try:
        min_par, max_par = (float(par) for par in missingness.split(","))
    except ValueError:
        quit_with_error(
            "Missingness percentage range must be given as two"
            " decimals separated by a comma without spaces"
        )
    aln_stats = {k: aln_stats[k] for k in aln_stats if min_par <= aln_stats[k]["missingness"] <= max_par}
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    log.log(f"{'Min. sequences (--seq)':>{mar}}: {bold(num_sequences)}")
    try:
        min_par = int(num_sequences)
    except ValueError:
        quit_with_error("Minimum number of sequences must be an integer")
    aln_stats = {k: aln_stats[k] for k in aln_stats if min_par <= aln_stats[k]["num_sequences"]}
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    log.log(f"{'Min. samples (--sam)':>{mar}}: {bold(num_samples)}")
    try:
        min_par = int(num_samples)
    except ValueError:
        quit_with_error("Minimum number of samples must be an integer")
    aln_stats = {k: aln_stats[k] for k in aln_stats if min_par <= aln_stats[k]["num_samples"]}
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    log.log(f"{'Min. focal species (--fos)':>{mar}}: {bold(num_focal_species)}")
    try:
        min_par = int(num_focal_species)
    except ValueError:
        quit_with_error("Minimum number of focal species must be an integer")
    aln_stats = {k: aln_stats[k] for k in aln_stats if min_par <= aln_stats[k]["num_focal_species"]}
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    log.log(f"{'Min. outgroup species (--ous)':>{mar}}: {bold(num_outgroup_species)}")
    try:
        min_par = int(num_outgroup_species)
    except ValueError:
        quit_with_error("Minimum number of outgroup species must be an integer")
    aln_stats = {k: aln_stats[k] for k in aln_stats if min_par <= aln_stats[k]["num_outgroup_species"]}
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    log.log(f"{'Min. add-on samples (--ads)':>{mar}}: {bold(num_addon_samples)}")
    try:
        min_par = int(num_addon_samples)
    except ValueError:
        quit_with_error("Minimum number of add-on samples must be an integer")
    aln_stats = {k: aln_stats[k] for k in aln_stats if min_par <= aln_stats[k]["num_addon_samples"]}
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    log.log(f"{'Min. species (--spp)':>{mar}}: {bold(num_species)}")
    try:
        min_par = int(num_species)
    except ValueError:
        quit_with_error("Minimum number of species must be an integer")
    aln_stats = {k: aln_stats[k] for k in aln_stats if min_par <= aln_stats[k]["num_species"]}
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    log.log(f"{'Min. genera (--gen)':>{mar}}: {bold(num_genera)}")
    try:
        min_par = int(num_genera)
    except ValueError:
        quit_with_error("Minimum number of genera must be an integer")
    aln_stats = {k: aln_stats[k] for k in aln_stats if min_par <= aln_stats[k]["num_genera"]}
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    log.log(f"{'Original CDS length (--cdl)':>{mar}}: {bold(cds_len)} bp")
    try:
        min_par, max_par = (int(par) for par in cds_len.split(","))
    except ValueError:
        quit_with_error(
            "Original CDS length range must be given as two integers separated by a comma without spaces"
        )
    aln_stats = {
        k: aln_stats[k]
        for k in aln_stats
        if (aln_stats[k]["cds_len"] == "NA" or min_par <= int(aln_stats[k]["cds_len"]) <= max_par)
    }
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    log.log(f"{'Len. long exons retained (--llr)':>{mar}}: {bold(len_long_exons_retained)} bp")
    try:
        min_par, max_par = (int(par) for par in len_long_exons_retained.split(","))
    except ValueError:
        quit_with_error(
            "Length of long exons retained range must be given as"
            " two integers separated by a comma without spaces"
        )
    aln_stats = {
        k: aln_stats[k]
        for k in aln_stats
        if (
            aln_stats[k]["len_long_exons_retained"] == "NA"
            or min_par <= int(aln_stats[k]["len_long_exons_retained"]) <= max_par
        )
    }
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    log.log(f"{'Len. short exons retained (--lsr)':>{mar}}: {bold(len_short_exons_retained)} bp")
    try:
        min_par, max_par = (int(par) for par in len_short_exons_retained.split(","))
    except ValueError:
        quit_with_error(
            "Length of short exons retained range must be given as"
            " two integers separated by a comma without spaces"
        )
    aln_stats = {
        k: aln_stats[k]
        for k in aln_stats
        if (
            aln_stats[k]["len_short_exons_retained"] == "NA"
            or min_par <= int(aln_stats[k]["len_short_exons_retained"]) <= max_par
        )
    }
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    log.log(f"{'Total CDS retained (--ptr)':>{mar}}: {bold(perc_exons_retained)} %")
    try:
        min_par, max_par = (float(par) for par in perc_exons_retained.split(","))
    except ValueError:
        quit_with_error(
            "Percentage of CDS retained range must be given as two"
            " decimals separated by a comma without spaces"
        )
    aln_stats = {
        k: aln_stats[k]
        for k in aln_stats
        if (
            aln_stats[k]["perc_exons_retained"] == "NA"
            or min_par <= float(aln_stats[k]["perc_exons_retained"]) <= max_par
        )
    }
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    log.log(f"{'Long exons retained (--plr)':>{mar}}: {bold(perc_long_exons_retained)} %")
    try:
        min_par, max_par = (float(par) for par in perc_long_exons_retained.split(","))
    except ValueError:
        quit_with_error(
            "Percentage of long exons retained range must be given as two"
            " decimals separated by a comma without spaces"
        )
    aln_stats = {
        k: aln_stats[k]
        for k in aln_stats
        if (
            aln_stats[k]["perc_long_exons_retained"] == "NA"
            or min_par <= float(aln_stats[k]["perc_long_exons_retained"]) <= max_par
        )
    }
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    log.log(f"{'Short exons retained (--psr)':>{mar}}: {bold(perc_short_exons_retained)} %")
    try:
        min_par, max_par = (float(par) for par in perc_short_exons_retained.split(","))
    except ValueError:
        quit_with_error(
            "Percentage of short exons retained range must be given as two"
            " decimals separated by a comma without spaces"
        )
    aln_stats = {
        k: aln_stats[k]
        for k in aln_stats
        if (
            aln_stats[k]["perc_short_exons_retained"] == "NA"
            or min_par <= float(aln_stats[k]["perc_short_exons_retained"]) <= max_par
        )
    }
    log.log(f"{'':>{mar}}  {len(aln_stats)} remaining loci")
    log.log("")

    footprint = 0
    for locus in aln_stats:
        if aln_stats[locus]["len_long_exons_retained"] != "NA":
            footprint += int(aln_stats[locus]["len_long_exons_retained"])
        else:
            footprint += aln_stats[locus]["length"]

    log.log("")
    log.log("")
    log.log(f"{'FINAL CAPTURE FOOTPRINT':>{mar}}: {bold(footprint)} bp in {bold(len(aln_stats))} loci")

    return aln_stats


def copy_loci(aln_stats: dict, out_dir: Path, overwrite: bool, show_more: bool):
    start = time.time()
    sel_dir, _ = make_output_dir(Path(out_dir, settings.DES_DIRS["AUT"]))
    _, _ = make_output_dir(Path(out_dir, settings.DES_DIRS["MAN"]))
    if overwrite or dir_is_empty(sel_dir):
        log.log("")
        log.log("")
        log.log(bold("Copying selected loci:"))
        tqdm_cols = min(shutil.get_terminal_size().columns, 120)
        with tqdm(total=len(aln_stats), ncols=tqdm_cols, unit="loci") as pbar:
            for locus in aln_stats:
                inner_start = time.time()
                src = aln_stats[locus]["path"]
                dst = Path(sel_dir, src.name)
                shutil.copy(f"{src}", f"{dst}")
                msg = f"'{src.name}': copied [{elapsed_time(time.time() - inner_start)}]"
                if show_more:
                    tqdm.write(msg)
                log.log(msg, print_to_screen=False)
                pbar.update()
        log.log(
            bold(
                f" \u2514\u2500\u2192 Successfully copied {len(aln_stats)} loci to"
                f" '{sel_dir}' [{elapsed_time(time.time() - start)}]"
            )
        )
    else:
        quit_with_error(f"'{sel_dir}': is not empty, enable '--overwrite' if needed")
    return


def write_aln_stats(out_dir: Path, tsv_comment: str, aln_stats_filtered: dict):
    stats_tsv_file = Path(out_dir, "captus-select_alignments.tsv")
    if not aln_stats_filtered:
        if stats_tsv_file.exists() and not file_is_empty(stats_tsv_file):
            return stats_tsv_file
        else:
            return None
    else:
        with open(stats_tsv_file, "wt") as tsv_out:
            tsv_out.write(tsv_comment)
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
            for locus in sorted(aln_stats_filtered):
                tsv_out.write(
                    "\t".join(
                        [
                            f"{aln_stats_filtered[locus]['path']}",
                            f"{locus}",
                            f"{aln_stats_filtered[locus]['copies']}",
                            f"{aln_stats_filtered[locus]['avg_copies']}",
                            f"{aln_stats_filtered[locus]['length']}",
                            f"{aln_stats_filtered[locus]['gc_content']}",
                            f"{aln_stats_filtered[locus]['avg_pid']}",
                            f"{aln_stats_filtered[locus]['informative_sites']}",
                            f"{aln_stats_filtered[locus]['informativeness']}",
                            f"{aln_stats_filtered[locus]['missingness']}",
                            f"{aln_stats_filtered[locus]['num_sequences']}",
                            f"{aln_stats_filtered[locus]['num_samples']}",
                            f"{aln_stats_filtered[locus]['num_focal_species']}",
                            f"{aln_stats_filtered[locus]['num_outgroup_species']}",
                            f"{aln_stats_filtered[locus]['num_addon_samples']}",
                            f"{aln_stats_filtered[locus]['num_species']}",
                            f"{aln_stats_filtered[locus]['num_genera']}",
                            f"{aln_stats_filtered[locus]['cds_id']}",
                            f"{aln_stats_filtered[locus]['cds_len']}",
                            f"{aln_stats_filtered[locus]['len_long_exons_retained']}",
                            f"{aln_stats_filtered[locus]['len_short_exons_retained']}",
                            f"{aln_stats_filtered[locus]['perc_exons_retained']}",
                            f"{aln_stats_filtered[locus]['perc_long_exons_retained']}",
                            f"{aln_stats_filtered[locus]['perc_short_exons_retained']}",
                        ]
                    )
                    + "\n"
                )
        return stats_tsv_file
