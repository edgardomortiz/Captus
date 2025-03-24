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
import shutil
import subprocess
import sys
from pathlib import Path

# Locus name separator
REF_CLUSTER_SEP = "-"

# Retain cluster reps that are at least CLR_MIN_LENGTH_PROP * longest cluster representative seq
CLR_MIN_LENGTH_PROP = 0.66

# Retain cluster reps that are at least CLR_MIN_SIZE_PROP * deepest cluster
CLR_MIN_SIZE_PROP = 0.10


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


def split_mmseqs_clusters_file(all_seqs_file_path):
    """
    Takes the not-so-easy-to-parse output format from MMseqs2 where the start of each cluster is
    indicated by a repeated header and returns a list of clusters, each cluster in turn is formatted
    as a list where the even indices (starting with 0) contain the sequence headers and the odd
    indices contain the sequence. The list type was chosen because order matters within a cluster
    since the first sequence is the representative of the cluster.
    """
    if all_seqs_file_path.exists():
        mmseqs_clusters = []
        cluster = []
        last_line = ""
        with open(all_seqs_file_path, "rt") as clusters_in:
            for line in clusters_in:
                line = line.strip()

                # MMseqs2 modifies the sequence name, just check that two headers follow each other
                if all([last_line.startswith(">"), line.startswith(">"), len(cluster) > 2]):
                    del cluster[0]  # Erase repeated header at start of cluster
                    del cluster[-1]  # Erase present header that marks the start of next cluster
                    mmseqs_clusters.append(cluster)
                    cluster = [line]
                cluster.append(line)
                last_line = line
            del cluster[0]  # Erase repeated header at start of cluster
            mmseqs_clusters.append(cluster)
        return mmseqs_clusters
    else:
        return None


def select_refs_per_locus(
    fastas_paths: list,
    fasta_ext: str,
    min_identity: float,
    min_coverage: float,
    out_dir: Path,
    prefix: str,
    threads,
):
    # 1. Disalign sequences, replace "__" with "_", add locus to each seq, concatenate for clustering
    clust_input = {}
    for fasta_path in fastas_paths:
        fasta = fasta_to_dict(fasta_path)
        locus_name = fasta_path.name.replace(f".{fasta_ext}", "")
        for seq_name in fasta:
            new_seq_name = f"{seq_name.replace('__', '_')}-{locus_name}"
            new_seq = fasta[seq_name]["sequence"].replace("-", "").replace("n", "")
            clust_input[new_seq_name] = {
                "sequence": new_seq,
                "description": fasta[seq_name]["description"],
            }
    clust_input_path = Path(out_dir, "clust_input.fasta")
    dict_to_fasta(clust_input, clust_input_path)

    # 2. Cluster
    prefix = f"{prefix}_cl{min_identity:.2f}_cov{min_coverage:.2f}"
    tmp_path = Path(out_dir, "tmp")
    min_identity = float(min_identity / 100)
    min_coverage = float(min_coverage / 100)
    if threads == "auto":
        threads = os.cpu_count()
    else:
        threads = int(threads)
    mmseqs_cmd = [
        "mmseqs",
        "easy-cluster",
        f"{clust_input_path}",
        f"{Path(out_dir, prefix)}",
        f"{tmp_path}",
        "--spaced-kmer-mode",
        f"{0}",
        "-c",
        f"{min_coverage}",
        "--cov-mode",
        f"{1}",
        "--alignment-mode",
        f"{3}",
        "--min-seq-id",
        f"{min_identity}",
        "--seq-id-mode",
        f"{1}",
        "--cluster-mode",
        f"{2}",
        "--kmer-per-seq-scale",
        f"{0.3}",
        "--threads",
        f"{threads}",
        "--cluster-reassign",
        f"{1}",
        "-s",
        f"{7.5}",
        "--gap-extend",
        f"{1}",
    ]
    if fasta_ext == "fna":
        mmseqs_cmd += [
            "--gap-open",
            f"{3}",
        ]
    elif fasta_ext == "faa":
        mmseqs_cmd += [
            "--gap-open",
            f"{9}",
        ]
    mmseqs_log_file = Path(out_dir, f"{prefix}_mmseqs.log")
    with open(mmseqs_log_file, "w") as mmseqs_log:
        mmseqs_log.write(f"Captus' MMseqs2 Command:\n  {' '.join(mmseqs_cmd)}\n\n")
    with open(mmseqs_log_file, "a") as mmseqs_log:
        subprocess.run(mmseqs_cmd, stdout=mmseqs_log, stdin=mmseqs_log)

    all_seqs = Path(out_dir, f"{prefix}_all_seqs.fasta")
    cluster_tsv = Path(out_dir, f"{prefix}_cluster.tsv")
    rep_seq = Path(out_dir, f"{prefix}_rep_seq.fasta")

    # 3. Get max depth and length for each locus' clusters
    clusters = split_mmseqs_clusters_file(all_seqs)
    loci = {}
    for cluster in clusters:
        locus = cluster[0].split()[0].split(REF_CLUSTER_SEP)[-1]
        cluster_info = {
            "size": len(cluster) / 2,
            "length": len(cluster[1]),
        }
        if locus not in loci:
            loci[locus] = cluster_info
        else:
            if cluster_info["size"] > loci[locus]["size"]:
                loci[locus] = cluster_info
            elif cluster_info["size"] == loci[locus]["size"]:
                if cluster_info["length"] > loci[locus]["length"]:
                    loci[locus] = cluster_info

    # 4. Select and format cluster representatives
    loci_reps = {}
    for cluster in clusters:
        locus = cluster[0].split()[0].split(REF_CLUSTER_SEP)[-1]
        size = len(cluster) / 2
        length = len(cluster[1])
        if (
            size >= loci[locus]["size"] * CLR_MIN_SIZE_PROP
            and length >= loci[locus]["length"] * CLR_MIN_LENGTH_PROP
        ):
            seq_name = cluster[0].split()[0].replace(">", "")
            sequence = cluster[1]
            description = " ".join(cluster[0].split()[1:])
            description = f"[cluster_size={size:.0f}] {description}"
            if locus not in loci_reps:
                loci_reps[locus] = {
                    seq_name: {
                        "sequence": sequence,
                        "description": description,
                    }
                }
            else:
                loci_reps[locus][seq_name] = {
                    "sequence": sequence,
                    "description": description,
                }

    # 5. SAve new target file
    new_targets_path = Path(out_dir, f"{prefix}.fasta")
    for locus in sorted(loci_reps):
        dict_to_fasta(loci_reps[locus], new_targets_path, append=True, sort=True)

    # 6. Intermediate file cleanup
    clust_input_path.unlink()
    all_seqs.unlink()
    cluster_tsv.unlink()
    rep_seq.unlink()
    shutil.rmtree(tmp_path, ignore_errors=True)

    return new_targets_path


def main():
    parser = argparse.ArgumentParser(
        description="Filter most common target per locus to create new reference target files"
    )
    parser.add_argument(
        "-a",
        "--captus_alignments_dir",
        action="store",
        default="./04_alignments",
        dest="captus_alignments_dir",
        required=True,
        help="Path to the directory that contains the output from the alignment step of Captus",
    )
    parser.add_argument(
        "--stage_dir",
        action="store",
        default="03_trimmed",
        dest="stage_dir",
        choices=["01_unaligned", "02_untrimmed", "03_trimmed"],
        help="Alignment stage",
    )
    parser.add_argument(
        "--filter_dir",
        action="store",
        default="06_informed",
        dest="filter_dir",
        choices=["04_unfiltered", "05_naive", "06_informed"],
        help="Paralog filter",
    )
    parser.add_argument(
        "--marker_dir",
        action="store",
        default="01_coding_NUC",
        dest="marker_dir",
        choices=["01_coding_NUC", "02_coding_PTD", "03_coding_MIT", "04_misc_DNA", "05_clusters"],
        help="Marker type",
    )
    parser.add_argument(
        "--format_dir",
        action="store",
        default="02_NT",
        dest="format_dir",
        choices=["01_AA", "02_NT", "03_genes", "04_genes_flanked", "01_matches", "02_matches_flanked"],
        help="Alignment data format",
    )
    parser.add_argument(
        "-o",
        "--out",
        action="store",
        default="./new_targets",
        dest="out",
        help="Output directory name",
    )
    parser.add_argument(
        "--min_identity",
        action="store",
        default=75,
        type=float,
        dest="min_identity",
        help="Minimum identity percentage between sequences in a cluster",
    )
    parser.add_argument(
        "--min_coverage",
        action="store",
        default=33,
        type=float,
        dest="min_coverage",
        help="Any sequence in a cluster has to be at least this percent included in the length"
        " of the longest sequence in the cluster",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        action="store",
        default="new_targets",
        dest="prefix",
        help="Prefix for output files, maybe use the original name of target file e.g. Mega353",
    )
    parser.add_argument(
        "--threads",
        action="store",
        default="auto",
        dest="threads",
        help="Number of threads to use",
    )
    args = parser.parse_args()

    if not Path(args.out).exists():
        try:
            Path(args.out).mkdir(parents=True)
        except OSError:
            quit_with_error(f"Captus was unable to make the output directory {Path(args.out)}")

    aln_dir = Path(
        args.captus_alignments_dir,
        args.stage_dir,
        args.filter_dir,
        args.marker_dir,
        args.format_dir,
    )

    if not aln_dir.is_dir():
        quit_with_error(f"'{aln_dir}' not found, verify this is a valid Captus alignment directory")

    fasta_ext = "fna"
    if args.format_dir == "01_AA":
        fasta_ext = "faa"

    fastas_paths = list(aln_dir.glob(f"*.{fasta_ext}"))

    new_targets_path = select_refs_per_locus(
        fastas_paths,
        fasta_ext,
        args.min_identity,
        args.min_coverage,
        args.out,
        args.prefix,
        args.threads,
    )
    print(f"New target file saved to '{new_targets_path}'")

if __name__ == "__main__":
    main()
