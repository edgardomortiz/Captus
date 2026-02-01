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

import argparse
import os
import sys
from pathlib import Path

# IQ-TREE parameters
SEED_NUM = 428571
BB_NUM = 1000

# FastTree parameters
SPR_NUM = 6
MLACC_NUM = 3

# Captus directories
CAPTUS_DIRS = {
    "unaligned": "01_unaligned",
    "untrimmed": "02_untrimmed",
    "trimmed": "03_trimmed",
    "unfiltered": "04_unfiltered",
    "naive": "05_naive",
    "informed": "06_informed",
    "unfiltered_w_refs": "01_unfiltered_w_refs",
    "naive_w_refs": "02_naive_w_refs",
    "informed_w_refs": "03_informed_w_refs",
    "NUC": "01_coding_NUC",
    "PTD": "02_coding_PTD",
    "MIT": "03_coding_MIT",
    "DNA": "04_misc_DNA",
    "CLR": "05_clusters",
    "AA": "01_AA",
    "NT": "02_NT",
    "GE": "03_genes",
    "GF": "04_genes_flanked",
    "MA": "01_matches",
    "MF": "02_matches_flanked",
}


def quit_with_error(message):
    """
    Displays the given message and ends the program's execution.
    """
    print(f"\nERROR: {message}\n")
    sys.exit(os.EX_SOFTWARE)


def create_iqtree_cmds(
    iqtree_path: str,
    fastas_paths: list,
    seq_type: str,
    phylogenies_dir: str,
    ufboot: int,
    alrt: int,
    force_bifurcating: bool,
    extra_options: str,
    threads: int,
):
    commands_list = []
    for fasta in fastas_paths:
        locus_name = fasta.stem
        cmd = [
            f"{iqtree_path}",
            "-s",
            f"{Path(fasta).resolve()}",
            "-st",
            f"{seq_type}",
            "-pre",
            f"{Path(Path(phylogenies_dir).resolve(), locus_name)}",
            "-nt",
            "AUTO",
            "-ntmax",
            f"{threads}",
            "-seed",
            f"{SEED_NUM}",
            "-m",
            "TEST",
        ]
        if ufboot > 0:
            cmd += ["-B", f"{ufboot}"]
        if alrt >= 0:
            cmd += ["-alrt", f"{alrt}"]
        if force_bifurcating is False:
            cmd += ["-czb"]
        if extra_options:
            cmd += extra_options.split()
        commands_list.append(" ".join(cmd))
    return commands_list


def create_fasttree_cmds(
    fasttree_path: str,
    fastas_paths: list,
    seq_type: str,
    phylogenies_dir: str,
    extra_options: str,
):
    commands_list = []
    for fasta in fastas_paths:
        locus_name = fasta.stem
        cmd = [
            f"{fasttree_path}",
            "-pseudo",
            "-spr",
            f"{SPR_NUM}",
            "-mlacc",
            f"{MLACC_NUM}",
            "-slownni",
        ]
        if seq_type == "DNA":
            cmd += [
                "-nt",
                "-gtr",
            ]
        if extra_options:
            cmd += extra_options.split()
        cmd += [
            f"{Path(fasta).resolve()}",
            f"1>{Path(Path(phylogenies_dir).resolve(), locus_name)}.treefile",
            f"2>{Path(Path(phylogenies_dir).resolve(), locus_name)}.fasttree.log",
        ]
        commands_list.append(" ".join(cmd))
    return commands_list


def main():
    parser = argparse.ArgumentParser(
        description="Create the commands to estimate phylogenies with IQ-TREE or FastTree"
        " corresponding to the alignments produced by Captus",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-a",
        "--captus_alignments",
        action="store",
        default="./04_alignments",
        dest="captus_alignments",
        help="Path to the directory that contains the output from the alignment step of Captus"
        " The path to a text file containing the list of paths to the alignments can also be"
        " provided, only alignments with extension .fna or .faa are accepted (if your alignments"
        " end in .faa use '-F AA')",
    )
    parser.add_argument(
        "-s",
        "--stage",
        action="store",
        default="trimmed",
        dest="stage",
        choices=["untrimmed", "trimmed"],
        help="Trimming stage",
    )
    parser.add_argument(
        "-f",
        "--filter",
        action="store",
        default="informed",
        dest="filter",
        choices=[
            "unfiltered",
            "naive",
            "informed",
            "unfiltered_w_refs",
            "naive_w_refs",
            "informed_w_refs",
        ],
        help="Paralog filter",
    )
    parser.add_argument(
        "-M",
        "--marker",
        action="store",
        default="NUC",
        dest="marker",
        choices=["NUC", "PTD", "MIT", "DNA", "CLR"],
        help="Marker type",
    )
    parser.add_argument(
        "-F",
        "--format",
        action="store",
        default="NT",
        dest="format",
        choices=["AA", "NT", "GE", "GF", "MA", "MF"],
        help="Alignment data format",
    )
    parser.add_argument(
        "-l",
        "--alignments_list",
        action="store",
        dest="alignments_list",
        help="Instead of searching within the Captus alignments directory you can provide a file"
        " list of paths to the alignments in FASTA format, one path per line (all the options above"
        " will be ignored except for '--format', i.e. if you want to process protein aligments use"
        " '--format AA')",
    )
    parser.add_argument(
        "-d",
        "--phylogenies_dir",
        action="store",
        default="./05_phylogenies",
        dest="phylogenies_dir",
        help="Destination directory for the estimated phylogenies",
    )
    parser.add_argument(
        "-o",
        "--out_dir",
        action="store",
        default="./phylo_commands",
        dest="out_dir",
        help="Output directory name",
    )
    parser.add_argument(
        "-c",
        "--out_commands",
        action="store",
        default="phylo_commands.txt",
        dest="out_commands",
        help="Prefix for output file, maybe use the original name of target file e.g. Allium353",
    )
    parser.add_argument(
        "-p",
        "--program",
        action="store",
        default="iqtree",
        dest="program",
        choices=["iqtree", "fasttree"],
        help="Phylogenetic software to use",
    )
    parser.add_argument(
        "-B",
        "--ufboot",
        action="store",
        default=0,
        type=int,
        dest="ufboot",
        help="Number of Ultrafast Bootstrap replicates for IQ-TREE, use values greater than 0 to"
        " enable",
    )
    parser.add_argument(
        "-A",
        "--alrt",
        action="store",
        default=-1,
        type=int,
        dest="alrt",
        help="Number of SH approximate likelihood replicates for IQ-TREE, 0 will perform the aLRT"
        " (Anisimova & Gascuel 2006) and values of at least 1 will perform the SH-alRT (Guindon et"
        " al. 2010)",
    )
    parser.add_argument(
        "--force_bifurcating",
        action="store_true",
        dest="force_bifurcating",
        help="Do not use the option -czb (Collapse Near Zero Branches) in IQ-TREE, to obtain fully"
        " bifurcating trees",
    )
    parser.add_argument(
        "-e",
        "--extra_options",
        action="store",
        default=None,
        dest="extra_options",
        help="Extra options to be passed to IQ-TREE or FastTree, enclose in quotations marks, e.g."
        ' "-wbtl -wsl"',
    )
    parser.add_argument(
        "--iqtree_path",
        action="store",
        default="iqtree",
        dest="iqtree_path",
        help="Path to IQ-TREE",
    )
    parser.add_argument(
        "--fasttree_path",
        action="store",
        default="fasttree",
        dest="fasttree_path",
        help="Path to FastTree",
    )
    parser.add_argument(
        "-t",
        "--threads",
        action="store",
        default=1,
        type=int,
        dest="threads",
        help="Threads to use for each IQTREE command. If FastTree is chosen, this value is ignored"
        " because FastTree is single-threaded",
    )
    args = parser.parse_args()

    fasta_ext = "fna"
    seq_type = "DNA"
    if args.format == "AA":
        fasta_ext = "faa"
        seq_type = "AA"
    fastas_paths = []
    if args.alignments_list is None:
        if not Path(args.captus_alignments).exists():
            quit_with_error(
                f"'{args.captus_alignments}' not found, verify this Captus alignment directory exists!"
            )
        elif Path(args.captus_alignments).is_dir():
            aln_dir = Path(
                args.captus_alignments,
                CAPTUS_DIRS[args.stage],
                CAPTUS_DIRS[args.filter],
                CAPTUS_DIRS[args.marker],
                CAPTUS_DIRS[args.format],
            )
            if not aln_dir.exists():
                quit_with_error(f"'{aln_dir}' not found, verify this Captus alignment directory exists!")
            fastas_paths = list(sorted(aln_dir.glob(f"*.{fasta_ext}")))
        elif Path(args.captus_alignments).is_file():
            with open(Path(args.captus_alignments), "rt") as paths_in:
                for line in paths_in:
                    fasta_path = Path(line.strip())
                    if fasta_path.exists() and fasta_path.suffix == f".{fasta_ext}":
                        fastas_paths.append(fasta_path)
    else:
        aln_paths_file = Path(args.alignments_list)
        if aln_paths_file.is_file():
            with open(aln_paths_file, "rt") as alns:
                for line in alns:
                    fasta_path = Path(line.strip())
                    if not fasta_path.is_file():
                        print(f"WARNING: file '{fasta_path}' not found, verify its location")
                    else:
                        fastas_paths.append(fasta_path)
        else:
            quit_with_error(f"'{aln_paths_file}' not found, verify its location")

    if not Path(args.out_dir).exists():
        try:
            Path(args.out_dir).mkdir(parents=True)
        except OSError:
            quit_with_error(f"Captus was unable to make the output directory {Path(args.out_dir)}")

    if not Path(args.phylogenies_dir).exists():
        try:
            Path(args.phylogenies_dir).mkdir(parents=True)
        except OSError:
            quit_with_error(
                f"Captus was unable to make the phylogenies directory {Path(args.phylogenies_dir)}"
            )

    if fastas_paths:
        if args.program.lower() == "iqtree":
            output_commands_list = create_iqtree_cmds(
                args.iqtree_path,
                fastas_paths,
                seq_type,
                args.phylogenies_dir,
                args.ufboot,
                args.alrt,
                args.force_bifurcating,
                args.extra_options,
                args.threads,
            )
        elif args.program.lower() == "fasttree":
            output_commands_list = create_fasttree_cmds(
                args.fasttree_path,
                fastas_paths,
                seq_type,
                args.phylogenies_dir,
                args.extra_options,
            )
        out_commands_path = Path(args.out_dir, args.out_commands)
        with open(Path(out_commands_path), "wt") as cmd:
            cmd.write("\n".join(output_commands_list) + "\n")
        print(f"A total of {len(output_commands_list)} commands saved to '{out_commands_path}'")
    else:
        print(
            f"No valid alignments or valid paths to alignments were found in '{args.captus_alignments}'"
        )


if __name__ == "__main__":
    main()
