#!/usr/bin/env python3
"""
Copyright 2020-2024 Edgardo M. Ortiz (e.ortiz.v@gmail.com)
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
        seq  = ""
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
    in_fasta_dict, out_fasta_path, wrap=0, sort=False,
    shuffle=False, append=False, write_if_empty=False
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
                header = f'>{name} {in_fasta_dict[name]["description"]}'.strip()
                seq = in_fasta_dict[name]["sequence"]
                if wrap > 0:
                    seq_out = "\n".join([seq[i:i + wrap] for i in range(0, len(seq), wrap)])
                else:
                    seq_out = seq
                fasta_out.write(f'{header}\n{seq_out}\n')
    else:
        if write_if_empty:
            with open(out_fasta_path, action) as fasta_out:
                fasta_out.write("")

    return out_fasta_path


def quit_with_error(message):
    """
    Displays the given message and ends the program's execution.
    """
    print(f"\nERROR: {message}\n")
    sys.exit(os.EX_SOFTWARE)


def main():

    parser = argparse.ArgumentParser(description="Filter most common target per locus to create new reference target files")
    parser.add_argument(
        "-a", "--captus_alignments_dir",
        action="store",
        default="./04_alignments",
        dest="captus_alignments_dir",
        required=True,
        help="Path to the directory that contains the output from the alignment step of Captus",
    )
    parser.add_argument(
        "-o", "--out",
        action="store",
        default="./",
        dest="out",
        help="Output directory name",
    )
    parser.add_argument(
        "-s", "--suffix",
        action="store",
        default= "_single_target",
        dest="suffix",
        help="Suffix to add to the output reference target files",
    )
    args = parser.parse_args()

    aln_log_file = Path(args.captus_alignments_dir, "captus-assembly_align.log")
    paralog_file = Path(args.captus_alignments_dir, "captus-assembly_align.paralogs.tsv")

    targets_paths = {
        "NUC": {"AA_path": None, "AA_fasta": None, "AA_names": [], "NT_path": None, "NT_fasta": None, "NT_names": []},
        "PTD": {"AA_path": None, "AA_fasta": None, "AA_names": [], "NT_path": None, "NT_fasta": None, "NT_names": []},
        "MIT": {"AA_path": None, "AA_fasta": None, "AA_names": [], "NT_path": None, "NT_fasta": None, "NT_names": []},
        "DNA": {"NT_path": None, "NT_fasta": None, "NT_names": []},
        "CLR": {"NT_path": None, "NT_fasta": None, "NT_names": []},
    }

    if not Path(args.out).exists():
        try:
            Path(args.out).mkdir(parents=True)
        except OSError:
            quit_with_error(f"Captus was unable to make the output directory {Path(args.out)}")

    if not aln_log_file.is_file():
        quit_with_error(f"'{aln_log_file.name}' not found in '{args.out}', verify this is a valid Captus alignment directory")

    if not paralog_file.is_file():
        quit_with_error(f"'{paralog_file.name}' not found in '{args.out}', verify this is a valid Captus alignment directory")

    with open(aln_log_file, "rt") as log_in:
        for line in log_in:
            if ("                       NUC: " in line or
                "                       PTD: " in line or
                "                       MIT: " in line):
                marker = line.split()[0].strip(":")
                if "(AA_path)" in line:
                    target_path = Path(line.split()[1])
                    if target_path.is_file():
                        targets_paths[marker]["AA_path"] = target_path
                        targets_paths[marker]["AA_fasta"] = fasta_to_dict(target_path)
                elif "(NT_path)" in line:
                    target_path = Path(line.split()[1])
                    if target_path.is_file():
                        targets_paths[marker]["NT_path"] = target_path
                        targets_paths[marker]["NT_fasta"] = fasta_to_dict(target_path)
            elif ("                       DNA: " in line or
                  "                       CLR: " in line):
                marker = line.split()[0].strip(":")
                if "(NT_path)" in line:
                    target_path = Path(line.split()[1])
                    if target_path.is_file():
                        targets_paths[marker]["NT_path"] = target_path
                        targets_paths[marker]["NT_fasta"] = fasta_to_dict(target_path)

    with open(paralog_file, "rt") as par_in:
        for line in par_in:
            if not line.startswith("#") and not line.startswith("marker_type"):
                record = line.split()
                marker_type = record[0]
                format_filtered = record[1]
                target_name = record[3]
                if format_filtered == "NT":
                    targets_paths[marker_type]["NT_names"].append(target_name)
                elif format_filtered == "AA":
                    targets_paths[marker_type]["AA_names"].append(target_name)

    for marker in targets_paths:
        if marker in ["NUC", "PTD", "MIT"]:
            if len(targets_paths[marker]["AA_names"]) > 0:
                fasta_out = {}
                fasta_file_name = targets_paths[marker]["AA_path"]
                fasta_out_file_name = Path(args.out, f'{fasta_file_name.stem}{args.suffix}.faa')
                for target_name in targets_paths[marker]["AA_names"]:
                    try:
                        fasta_out[target_name] = targets_paths[marker]["AA_fasta"][target_name]
                    except KeyError:
                        print(f"'{target_name}' not found in '{fasta_file_name}'")
                dict_to_fasta(fasta_out, fasta_out_file_name)
                if fasta_file_name.is_file():
                    print(f"'{marker}' aminoacid reference target file saved to '{fasta_out_file_name}'")
        if marker in ["NUC", "PTD", "MIT", "DNA", "CLR"]:
            if len(targets_paths[marker]["NT_names"]) > 0:
                fasta_out = {}
                fasta_file_name = targets_paths[marker]["NT_path"]
                fasta_out_file_name = Path(args.out, f'{fasta_file_name.stem}{args.suffix}.fna')
                for target_name in targets_paths[marker]["NT_names"]:
                    try:
                        fasta_out[target_name] = targets_paths[marker]["NT_fasta"][target_name]
                    except KeyError:
                        print(f"'{target_name}' not found in '{fasta_file_name}'")
                dict_to_fasta(fasta_out, fasta_out_file_name)
                if fasta_file_name.is_file():
                    print(f"'{marker}' nucleotide reference target file saved to '{fasta_out_file_name}'")

    print("Done!")

if __name__ == "__main__":
    main()
