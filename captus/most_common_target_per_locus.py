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
import json
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


def quit_with_error(message):
    """
    Displays the given message and ends the program's execution.
    """
    print(f"\nERROR: {message}\n")
    sys.exit(os.EX_SOFTWARE)


def json_from_extractions(captus_extractions_dir: str):
    refs_json_file = Path(captus_extractions_dir, "captus-extract_refs.json")
    if refs_json_file.is_file():
        with open(refs_json_file, "rt") as jin:
            refs_json = json.load(jin)
    else:
        quit_with_error(
            f"The file '{refs_json_file}' was not found,"
            " please verify that your Captus extractions directory contains the file"
        )
    if refs_json:
        for marker_type in refs_json:
            refs_json[marker_type]["ref_names"] = {}
    else:
        quit_with_error(
            f"The file '{refs_json_file}' was empty,"
            " please verify that your Captus extractions directory contains a valid file"
        )

    stats_tsv_file = Path(captus_extractions_dir, "captus-extract_stats.tsv")
    if stats_tsv_file.is_file():
        with open(stats_tsv_file, "rt") as tsv:
            for line in tsv:
                if not line.startswith("#"):
                    record = line.strip().split("\t")
                    marker_type = record[1]
                    locus = record[2]
                    ref_name = record[3]
                    hit = record[7]
                    if hit == "00":
                        if locus in refs_json[marker_type]["ref_names"]:
                            if ref_name in refs_json[marker_type]["ref_names"][locus]:
                                refs_json[marker_type]["ref_names"][locus][ref_name] += 1
                            else:
                                refs_json[marker_type]["ref_names"][locus][ref_name] = 1
                        else:
                            refs_json[marker_type]["ref_names"][locus] = {ref_name: 1}
    else:
        quit_with_error(
            f"The file '{stats_tsv_file}' was not found,"
            " please verify that your Captus extractions directory contains the file"
        )

    return refs_json


def json_from_alignments(captus_alignments_dir: str):
    aln_log_file = Path(captus_alignments_dir, "captus-align.log")
    refs_json = {
        "NUC": {"AA_path": None, "NT_path": None, "ref_names": {}},
        "PTD": {"AA_path": None, "NT_path": None, "ref_names": {}},
        "MIT": {"AA_path": None, "NT_path": None, "ref_names": {}},
        "DNA": {"AA_path": None, "NT_path": None, "ref_names": {}},
        "CLR": {"AA_path": None, "NT_path": None, "ref_names": {}},
    }
    if aln_log_file.is_file():
        with open(aln_log_file, "rt") as log_in:
            for line in log_in:
                if (
                    "                       NUC: " in line
                    or "                       PTD: " in line
                    or "                       MIT: " in line
                    or "                       DNA: " in line
                    or "                       CLR: " in line
                ):
                    marker_type = line.split()[0].strip(":")
                    if "(AA_path)" in line:
                        ref_fasta_path = Path(line.split()[1])
                        if ref_fasta_path.is_file():
                            refs_json[marker_type]["AA_path"] = ref_fasta_path
                        else:
                            print(f"'{ref_fasta_path}' was not found, please verify the file location")
                    elif "(NT_path)" in line:
                        ref_fasta_path = Path(line.split()[1])
                        if ref_fasta_path.is_file():
                            refs_json[marker_type]["NT_path"] = ref_fasta_path
                        else:
                            print(f"'{ref_fasta_path}' was not found, please verify the file location")
    else:
        quit_with_error(
            f"The file '{aln_log_file}' was not found,"
            " please verify that your Captus alignments directory contains the file"
        )

    stats_tsv_file = Path(captus_alignments_dir, "captus-align_paralogs.tsv")
    if stats_tsv_file.is_file():
        with open(stats_tsv_file, "rt") as par_in:
            for line in par_in:
                if not line.startswith("#") and not line.startswith("marker_type"):
                    record = line.split()
                    marker_type = record[0]
                    locus = record[2]
                    ref_name = record[3]
                    hit = record[5]
                    if hit == "00":
                        if locus in refs_json[marker_type]["ref_names"]:
                            if ref_name in refs_json[marker_type]["ref_names"][locus]:
                                refs_json[marker_type]["ref_names"][locus][ref_name] += 1
                            else:
                                refs_json[marker_type]["ref_names"][locus][ref_name] = 1
                        else:
                            refs_json[marker_type]["ref_names"][locus] = {ref_name: 1}
    else:
        quit_with_error(
            f"The file '{stats_tsv_file}' was not found,"
            " please verify that your Captus alignments directory contains the file"
        )

    return refs_json


def save_mct(mct_fasta: dict, mct_fasta_path: Path, marker_type: str, seq_type: str, overwrite: bool):
    if mct_fasta_path.is_file():
        if overwrite:
            dict_to_fasta(mct_fasta, mct_fasta_path, write_if_empty=True)
            print(f"'{marker_type}' {seq_type} reference target file overwritten to '{mct_fasta_path}'")
        else:
            print(
                f"'{marker_type}' {seq_type} reference target file already exists"
                f" in '{mct_fasta_path}', use '--overwrite' to replace"
            )
    else:
        dict_to_fasta(mct_fasta, mct_fasta_path, write_if_empty=True)
        print(f"'{marker_type}' {seq_type} reference target file saved to '{mct_fasta_path}'")


def tally_json_and_save_mct(
    refs_json: dict, out_dir: str, suffix: str, min_samples: int, overwrite: bool
):
    for marker_type in refs_json:
        for locus in refs_json[marker_type]["ref_names"]:
            ref_names_count = refs_json[marker_type]["ref_names"][locus]
            num_samples = sum(refs_json[marker_type]["ref_names"][locus].values())
            if num_samples >= min_samples:
                refs_json[marker_type]["ref_names"][locus] = max(
                    ref_names_count, key=ref_names_count.get
                )
            else:
                refs_json[marker_type]["ref_names"][locus] = None
        if refs_json[marker_type]["AA_path"] is not None:
            ref_fasta_path = Path(refs_json[marker_type]["AA_path"])
            if ref_fasta_path.is_file():
                ref_fasta = fasta_to_dict(ref_fasta_path)
                ref_fasta_filtered = {}
                ref_fasta_filtered_path = Path(out_dir, f"{ref_fasta_path.stem}{suffix}.faa")
                for locus in sorted(refs_json[marker_type]["ref_names"]):
                    ref_name = refs_json[marker_type]["ref_names"][locus]
                    if ref_name is not None:
                        try:
                            ref_fasta_filtered[ref_name] = ref_fasta[ref_name]
                        except KeyError:
                            print(f"'{ref_name}' not found in '{ref_fasta_path}'")
                save_mct(
                    ref_fasta_filtered, ref_fasta_filtered_path, marker_type, "aminoacid", overwrite
                )
            else:
                print(f"'{ref_fasta_path}' was not found, please verify the file location")
        if refs_json[marker_type]["NT_path"] is not None:
            ref_fasta_path = Path(refs_json[marker_type]["NT_path"])
            if ref_fasta_path.is_file():
                ref_fasta = fasta_to_dict(ref_fasta_path)
                ref_fasta_filtered = {}
                ref_fasta_filtered_path = Path(out_dir, f"{ref_fasta_path.stem}{suffix}.fna")
                for locus in sorted(refs_json[marker_type]["ref_names"]):
                    ref_name = refs_json[marker_type]["ref_names"][locus]
                    if ref_name is not None:
                        try:
                            ref_fasta_filtered[ref_name] = ref_fasta[ref_name]
                        except KeyError:
                            print(f"'{ref_name}' not found in '{ref_fasta_path}'")
                save_mct(
                    ref_fasta_filtered, ref_fasta_filtered_path, marker_type, "nucleotide", overwrite
                )
            else:
                print(f"'{ref_fasta_path}' was not found, please verify the file location")


def main():
    parser = argparse.ArgumentParser(
        description="Filter most common target per locus to create new reference target files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-e",
        "--captus_extractions_dir",
        action="store",
        type=str,
        dest="captus_extractions_dir",
        help="Path to the directory that contains the output from the extraction step of Captus,"
        " tipically '03_extractions'",
    )
    parser.add_argument(
        "-a",
        "--captus_alignments_dir",
        action="store",
        type=str,
        dest="captus_alignments_dir",
        help="Path to the directory that contains the output from the alignment step of Captus,"
        " tipically '04_alignments'",
    )
    parser.add_argument(
        "-o",
        "--out_dir",
        action="store",
        default="./new_targets",
        type=str,
        dest="out_dir",
        help="Output directory name",
    )
    parser.add_argument(
        "-s",
        "--suffix",
        action="store",
        default="_MCT",
        type=str,
        dest="suffix",
        help="Suffix to add to the output reference target files (MCT = Most Common Target)",
    )
    parser.add_argument(
        "-m",
        "--min_samples",
        action="store",
        type=int,
        default=4,
        dest="min_samples",
        help="Only process loci found in at least this number of samples",
    )
    parser.add_argument(
        "--overwrite", action="store_true", dest="overwrite", help="Overwrite previous results"
    )
    args = parser.parse_args()

    if not Path(args.out_dir).exists():
        try:
            Path(args.out_dir).mkdir(parents=True)
        except OSError:
            quit_with_error(f"Captus was unable to make the output directory {Path(args.out_dir)}")

    if args.captus_extractions_dir is None and args.captus_alignments_dir is None:
        quit_with_error(
            "Please provide the path to a Captus extractions directory OR to a Captus alignments directory"
        )
    elif args.captus_extractions_dir is not None and args.captus_alignments_dir is not None:
        quit_with_error(
            "Please provide the path to a Captus extractions directory OR to a Captus alignments directory"
        )
    elif args.captus_extractions_dir:
        if Path(args.captus_extractions_dir).exists():
            refs_json = json_from_extractions(args.captus_extractions_dir)
        else:
            quit_with_error("Captus extractions directory not found, please verify the path is correct")
    elif args.captus_alignments_dir:
        if Path(args.captus_alignments_dir).exists():
            refs_json = json_from_alignments(args.captus_alignments_dir)
        else:
            quit_with_error("Captus alignments directory not found, please verify the path is correct")

    tally_json_and_save_mct(
        refs_json,
        args.out_dir,
        args.suffix,
        args.min_samples,
        args.overwrite,
    )


if __name__ == "__main__":
    main()
