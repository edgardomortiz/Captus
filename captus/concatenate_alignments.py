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
import gzip
import os
import random
import sys
from pathlib import Path


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
    "MF": "02_matches_flannked",
}

# Format extensions
EXT = {
    "fasta": "fas",
    "nexus": "nex",
    "phylip": "phy",
}


def fasta_to_dict(fasta_path):
    """
    Turns a FASTA file given with `fasta_path` into a dictionary. For example, for the sequence:
    ```text
    >k157_0 flag=1 multi=2.0000 len=274
    ATATTGATATTTCATAATAATAGTTTTTGAACTAAAAAGAAATTTTTCCTCCAATTATGTGGG
    ```
    Returns the dictionary:
    ```
    {'k157_0' : {
        'description': 'flag=1 multi=2.0000 len=274',
        'sequence': 'ATATTGATATTTCATAATAATAGTTTTTGAACTAAAAAGAAATTTTTCCTCCAAT...'
    }}
    ```
    """
    if f"{fasta_path}".endswith(".gz"):
        opener = gzip.open
    else:
        opener = open
    fasta_out = {}
    with opener(fasta_path, "rt") as fasta_in:
        seq = ""
        name = ""
        desc = ""
        for line in fasta_in:
            line = line.strip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if seq:
                    fasta_out[name] = {
                        "description": desc,
                        "sequence": seq,
                    }
                    seq = ""
                if len(line.split()) > 1:
                    name = line[1:].split()[0]
                    desc = " ".join(line.split()[1:])
                else:
                    name = line[1:].rstrip()
                    desc = ""
            else:
                seq += line
        if seq:
            fasta_out[name] = {
                "description": desc,
                "sequence": seq,
            }
    return fasta_out


def dict_to_fasta(
    in_fasta_dict, out_fasta_path, wrap=0, sort=False, shuffle=False, append=False, write_if_empty=False
):
    """
    Saves a `in_fasta_dict` from function `fasta_to_dict()` as a FASTA file to `out_fasta_path`
    """
    if f"{out_fasta_path}".endswith(".gz"):
        out_fasta_path = Path(f"{out_fasta_path}".strip(".gz"))
    action = "wt"
    if append is True:
        action = "at"
    if in_fasta_dict:
        if sort:
            in_fasta_dict = dict(sorted(in_fasta_dict.items(), key=lambda x: x[0]))
        if shuffle:
            in_fasta_dict_shuffled = list(in_fasta_dict.items())
            random.shuffle(in_fasta_dict_shuffled)
            in_fasta_dict = dict(in_fasta_dict_shuffled)
        with open(out_fasta_path, action) as fasta_out:
            for name in in_fasta_dict:
                header = f">{name} {in_fasta_dict[name]['description']}".strip()
                seq = in_fasta_dict[name]["sequence"]
                if wrap > 0:
                    seq_out = "\n".join([seq[i : i + wrap] for i in range(0, len(seq), wrap)])
                else:
                    seq_out = seq
                fasta_out.write(f"{header}\n{seq_out}\n")
    else:
        if write_if_empty:
            with open(out_fasta_path, action) as fasta_out:
                fasta_out.write("")

    return out_fasta_path


def dict_to_nexus(in_fasta_dict: dict, out_nexus_path: str, seq_type: str, wrap=0):
    ntax = len(in_fasta_dict)
    nchar = 0
    datatype = "DNA"
    if seq_type == "AA":
        datatype = "PROTEIN"
    max_seq_name_len = 0
    for seq_name in in_fasta_dict:
        if len(seq_name) > max_seq_name_len:
            max_seq_name_len = len(seq_name)
        nchar = len(in_fasta_dict[seq_name]["sequence"])
    interleave = ""
    if wrap > 0:
        interleave = "INTERLEAVE "

    with open(out_nexus_path, "wt") as nex:
        nex.write("#NEXUS\n")
        nex.write("BEGIN DATA;\n")
        nex.write(f"DIMENSIONS NTAX={ntax} NCHAR={nchar};\n")
        nex.write(f"FORMAT {interleave}DATATYPE={datatype} GAP=- MISSING=?;\n")
        nex.write("MATRIX\n")
        if wrap > 0:
            for i in range(0, nchar, wrap):
                for seq_name in in_fasta_dict:
                    padding = (max_seq_name_len - len(seq_name)) * " "
                    nex.write(
                        f"{seq_name}{padding} {in_fasta_dict[seq_name]['sequence'][i : i + wrap]}\n"
                    )
                nex.write("\n")
        else:
            for seq_name in in_fasta_dict:
                padding = (max_seq_name_len - len(seq_name)) * " "
                nex.write(f"{seq_name}{padding} {in_fasta_dict[seq_name]['sequence']}\n")
        nex.write(";\n")
        nex.write("END;\n")

    return out_nexus_path


def dict_to_phylip(in_fasta_dict: dict, out_phylip_path: str, wrap=0):
    ntax = len(in_fasta_dict)
    nchar = 0
    max_seq_name_len = 0
    for seq_name in in_fasta_dict:
        if len(seq_name) > max_seq_name_len:
            max_seq_name_len = len(seq_name)
        nchar = len(in_fasta_dict[seq_name]["sequence"])

    with open(out_phylip_path, "wt") as phy:
        phy.write(f"{ntax} {nchar}\n")
        if wrap > 0:
            for seq_name in in_fasta_dict:
                padding = (max_seq_name_len - len(seq_name)) * " "
                phy.write(f"{seq_name}{padding} {in_fasta_dict[seq_name]['sequence'][0:wrap]}\n")
            phy.write("\n")
            for i in range(wrap, nchar, wrap):
                for seq_name in in_fasta_dict:
                    phy.write(f"{in_fasta_dict[seq_name]['sequence'][i : i + wrap]}\n")
                phy.write("\n")
        else:
            for seq_name in in_fasta_dict:
                padding = (max_seq_name_len - len(seq_name)) * " "
                phy.write(f"{seq_name}{padding} {in_fasta_dict[seq_name]['sequence']}\n")

    return out_phylip_path


def quit_with_error(message):
    """
    Displays the given message and ends the program's execution.
    """
    print(f"\nERROR: {message}\n")
    sys.exit(os.EX_SOFTWARE)


def write_partitions(partitions: dict, seq_type: str, codon: str, out_path: str):
    """
    #raxml
    DNA, part1 = 1-100, 250-384
    DNA, part2 = 101-249\3, 102-249\3
    DNA, part3 = 103-249\3

    #nexus
    begin sets;
        charset part1 = 1-100 250-384;
        charset part2 = 101-249\3 102-249\3;
        charset part3 = 103-249\3;
    end;
    """
    out_path = Path(out_path)
    raxml_part_path = Path(out_path.parent, f"{out_path.stem}.raxml.part")
    nexus_part_path = Path(out_path.parent, f"{out_path.stem}.nexus.part")
    part_type = "DNA"
    if seq_type == "AA":
        part_type = "LG"
    raxml_parts = []
    nexus_parts = []
    acc_length = 0
    if codon != "123" and seq_type != "AA":
        codon = sorted([int(c) for c in list(codon)])
        for p in partitions:
            codon_coords = " ".join(
                [f"{acc_length + c}-{acc_length + partitions[p]['length']}\\3" for c in codon]
            )
            rax = f"{part_type}, {p} = {codon_coords}"
            raxml_parts.append(rax)
            nex = f"  CHARSET {p} = {codon_coords};"
            nexus_parts.append(nex)
            acc_length += partitions[p]["length"]
    else:
        for p in partitions:
            rax = f"{part_type}, {p} = {acc_length + 1}-{acc_length + partitions[p]['length']}"
            raxml_parts.append(rax)
            nex = f"  CHARSET {p} = {acc_length + 1}-{acc_length + partitions[p]['length']};"
            nexus_parts.append(nex)
            acc_length += partitions[p]["length"]

    with open(raxml_part_path, "wt") as rax:
        rax.write("\n".join(raxml_parts) + "\n")
    with open(nexus_part_path, "wt") as nex:
        nex.write("#NEXUS\n")
        nex.write("BEGIN SETS;\n")
        nex.write("\n".join(nexus_parts) + "\n")
        nex.write("END;\n")

    return raxml_part_path, nexus_part_path


def concat_alignments(
    fastas_paths: list,
):
    partitions = {}
    concatenated = {}
    for fasta in fastas_paths:
        fasta_in = fasta_to_dict(fasta)
        for seq_name in fasta_in:
            concatenated[seq_name] = {"sequence": [], "description": ""}
        for seq_name in fasta_in:
            part_len = len(fasta_in[seq_name]["sequence"])
            break
        partitions[fasta.stem] = {
            "length": part_len,
            "path": fasta,
        }

    for part in partitions:
        fasta_in = fasta_to_dict(partitions[part]["path"])
        for sample in concatenated:
            if sample in fasta_in:
                concatenated[sample]["sequence"].append(fasta_in[sample]["sequence"])
            else:
                concatenated[sample]["sequence"].append("-" * partitions[part]["length"])

    for sample in concatenated:
        concatenated[sample]["sequence"] = "".join(concatenated[sample]["sequence"])

    return partitions, concatenated


def main():
    parser = argparse.ArgumentParser(
        description="Concatenate the alignments contained in specific folders produced by Captus,"
        " save the supermatrix as FASTA, PHYLIP or NEXUS, including a partition file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-a",
        "--captus_alignments",
        action="store",
        default="./04_alignments",
        dest="captus_alignments",
        required=True,
        help="Path to the directory that contains the output from the alignment step of Captus"
        " The path to a text file containing the list of paths to the alignments can also be"
        " provided (if your alignments contain aminoacid sequences use '-F AA')",
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
        "-o",
        "--out_dir",
        action="store",
        default="./concatenated_alignments",
        dest="out_dir",
        help="Output directory name",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        action="store",
        default="supermatrix",
        dest="prefix",
        help="Prefix for output files",
    )
    parser.add_argument(
        "-e",
        "--out_format",
        action="store",
        default="fasta",
        dest="out_format",
        choices=["fasta", "nexus", "phylip"],
        help="Supermatrix format (will be used as file extension too)",
    )
    parser.add_argument(
        "-c",
        "--codon",
        action="store",
        default="123",
        dest="codon",
        choices=["123", "1", "2", "3", "12", "13", "23"],
        help="Codon or codon positions to include in the partition file, include all by default"
        " ignored for aminoacid alignments",
    )
    parser.add_argument(
        "-w",
        "--wrap",
        action="store",
        default=0,
        type=int,
        dest="wrap",
        help="Wrap sequence at this length in bp, 0 means no wrapping will be applied",
    )
    args = parser.parse_args()

    fasta_ext = "fna"
    seq_type = "DNA"
    if args.format == "AA":
        fasta_ext = "faa"
        seq_type = "AA"
    fastas_paths = []
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
                if fasta_path.exists():
                    fastas_paths.append(fasta_path)

    if not Path(args.out_dir).exists():
        try:
            Path(args.out_dir).mkdir(parents=True)
        except OSError:
            quit_with_error(f"Captus was unable to make the output directory {Path(args.out_dir)}")

    if fastas_paths:
        partitions, concatenated = concat_alignments(fastas_paths)

        out_supermatrix_path = Path(args.out_dir, f"{args.prefix}.{EXT[args.out_format]}")

        if args.out_format.lower() == "fasta":
            dict_to_fasta(concatenated, out_supermatrix_path, wrap=args.wrap)
        elif args.out_format.lower() == "nexus":
            dict_to_nexus(concatenated, out_supermatrix_path, seq_type, wrap=args.wrap)
        elif args.out_format.lower() == "phylip":
            dict_to_phylip(concatenated, out_supermatrix_path, wrap=args.wrap)

        raxml_part_path, nexus_part_path = write_partitions(
            partitions, seq_type, args.codon, out_supermatrix_path
        )

        print(f"Supermatrix file saved to '{out_supermatrix_path}'")
        print(f"RAxML-style partition file saved to '{raxml_part_path}'")
        print(f"NEXUS-style partition file saved to '{nexus_part_path}'")
    else:
        print(
            f"No valid alignments or valid paths to alignments were found in '{args.captus_alignments}'"
        )


if __name__ == "__main__":
    main()
