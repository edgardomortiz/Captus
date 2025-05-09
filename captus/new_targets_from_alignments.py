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
import statistics
import subprocess
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

# Locus name separator in target files
REF_CLUSTER_SEP = "-"

# Separator between sample name and paralog number
SEQ_NAME_SEP = "__"


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
    out_dir: Path,
    prefix: str,
    include_references: bool,
    exclude_samples,
    min_seq_len: int,
    min_identity: float,
    min_coverage: float,
    wscore_proportion: float,
    coverage_proportion: float,
    length_proportion: float,
    size_proportion: float,
    best_only: bool,
    threads,
):
    prefix = f"{prefix}_i{min_identity:.2f}_c{min_coverage:.2f}_W{wscore_proportion:.2f}"
    prefix += f"_C{coverage_proportion:.2f}_L{length_proportion:.2f}_S{size_proportion:.2f}"
    clust_log_file = Path(out_dir, f"{prefix}.log")

    if exclude_samples is None:
        exclude_samples = []
    else:
        exclude_samples = [str(x) for x in exclude_samples.split(",")]

    # 1. Disalign sequences, replace "__" with "_", add locus to each seq, concatenate for clustering
    # keep track of highest wscore per locus in 'max_wscores'
    full_clust_input = {}
    all_targets_info = {}
    for fasta_path in fastas_paths:
        fasta = fasta_to_dict(fasta_path)
        locus = fasta_path.name.replace(f".{fasta_ext}", "")
        for seq_name in fasta:
            sample_name = seq_name.split(SEQ_NAME_SEP)[0]
            new_seq = fasta[seq_name]["sequence"].replace("-", "").replace("n", "")
            new_seq_len = len(new_seq)
            if sample_name not in exclude_samples and new_seq_len >= min_seq_len:
                if include_references is True:
                    new_seq_name = f"{seq_name.replace(SEQ_NAME_SEP, '_')}{REF_CLUSTER_SEP}{locus}"
                    new_seq = fasta[seq_name]["sequence"].replace("-", "").replace("n", "")
                    full_clust_input[new_seq_name] = {
                        "sequence": new_seq,
                        "description": fasta[seq_name]["description"],
                    }
                    if "wscore" in fasta[seq_name]["description"]:
                        wscore = float(fasta[seq_name]["description"].split("wscore=")[-1].split("]")[0])
                    else:
                        wscore = 1.0
                    if "cover" in fasta[seq_name]["description"]:
                        coverage = float(
                            fasta[seq_name]["description"].split("cover=")[-1].split("]")[0]
                        )
                    else:
                        coverage = 1.0
                    if "query" in fasta[seq_name]["description"]:
                        target = fasta[seq_name]["description"].split("query=")[-1].split("]")[0]
                    else:
                        target = new_seq_name.replace("__ref", "")
                    if target not in all_targets_info:
                        all_targets_info[target] = {
                            "wscore": wscore,
                            "coverage": coverage,
                            "count": 1,
                        }
                    else:
                        all_targets_info[target]["count"] += 1
                        if wscore > all_targets_info[target]["wscore"]:
                            all_targets_info[target]["wscore"] = wscore
                            all_targets_info[target]["coverage"] = coverage
                        elif wscore == all_targets_info[target]["wscore"]:
                            if coverage > all_targets_info[target]["coverage"]:
                                all_targets_info[target]["coverage"] = coverage
                else:
                    if seq_name.endswith("__ref"):
                        continue
                    else:
                        new_seq_name = f"{seq_name.replace('__', '_')}-{locus}"
                        new_seq = fasta[seq_name]["sequence"].replace("-", "").replace("n", "")
                        full_clust_input[new_seq_name] = {
                            "sequence": new_seq,
                            "description": fasta[seq_name]["description"],
                        }
                        if "wscore" in fasta[seq_name]["description"]:
                            wscore = float(
                                fasta[seq_name]["description"].split("wscore=")[-1].split("]")[0]
                            )
                        else:
                            wscore = 1.0
                        if "cover" in fasta[seq_name]["description"]:
                            coverage = float(
                                fasta[seq_name]["description"].split("cover=")[-1].split("]")[0]
                            )
                        else:
                            coverage = 1.0
                        if "query" in fasta[seq_name]["description"]:
                            target = fasta[seq_name]["description"].split("query=")[-1].split("]")[0]
                        else:
                            target = new_seq_name.replace("__ref", "")
                        if target not in all_targets_info:
                            all_targets_info[target] = {
                                "wscore": wscore,
                                "coverage": coverage,
                                "count": 1,
                            }
                        else:
                            all_targets_info[target]["count"] += 1
                            if wscore > all_targets_info[target]["wscore"]:
                                all_targets_info[target]["wscore"] = wscore
                                all_targets_info[target]["coverage"] = coverage
                            elif wscore == all_targets_info[target]["wscore"]:
                                if coverage > all_targets_info[target]["coverage"]:
                                    all_targets_info[target]["coverage"] = coverage

    # 2. Retain only most common target per locus and its data
    best_targets_info = {}
    for target in all_targets_info:
        locus = target.split(REF_CLUSTER_SEP)[-1]
        target_info = {
            "target": target,
            "wscore": all_targets_info[target]["wscore"],
            "coverage": all_targets_info[target]["coverage"],
            "count": all_targets_info[target]["count"],
        }
        if locus not in best_targets_info:
            best_targets_info[locus] = target_info
        else:
            if all_targets_info[target]["count"] > best_targets_info[locus]["count"]:
                best_targets_info[locus] = target_info
            elif all_targets_info[target]["count"] == best_targets_info[locus]["count"]:
                if all_targets_info[target]["wscore"] > best_targets_info[locus]["wscore"]:
                    best_targets_info[locus] = target_info
                elif all_targets_info[target]["wscore"] == best_targets_info[locus]["wscore"]:
                    if all_targets_info[target]["coverage"] > best_targets_info[locus]["coverage"]:
                        best_targets_info[locus] = target_info
    if not exclude_samples:
        exclude_samples = None
    else:
        exclude_samples = ",".join(exclude_samples)
    msg = (
        f"PREFIX: {prefix}\n"
        f"LOG: {clust_log_file}\n"
        f"SAMPLES TO EXCLUDE: {exclude_samples}\n"
        "\n"
        f"Starting with a total of {len(full_clust_input)} sequences in {len(best_targets_info)} loci\n"
        "\n"
    )
    print(msg)
    with open(clust_log_file, "w") as log:
        log.write(msg)

    # 3. Filter according to most common target, "WSCORE_PROPORTION", and "COVERAGE_PROPORTION"
    clust_input = {}
    for seq_name in full_clust_input:
        locus = seq_name.split(REF_CLUSTER_SEP)[-1]
        if "wscore" in full_clust_input[seq_name]["description"]:
            wscore = float(full_clust_input[seq_name]["description"].split("wscore=")[-1].split("]")[0])
        else:
            wscore = 1.0
        if "cover" in full_clust_input[seq_name]["description"]:
            coverage = float(full_clust_input[seq_name]["description"].split("cover=")[-1].split("]")[0])
        else:
            coverage = 1.0
        if "query" in full_clust_input[seq_name]["description"]:
            target = full_clust_input[seq_name]["description"].split("query=")[-1].split("]")[0]
        else:
            target = seq_name.replace("__ref", "")
        if target == best_targets_info[locus]["target"]:
            if (
                wscore >= best_targets_info[locus]["wscore"] * wscore_proportion
                and coverage >= best_targets_info[locus]["coverage"] * coverage_proportion
            ):
                clust_input[seq_name] = full_clust_input[seq_name]
        else:
            clust_input[seq_name] = full_clust_input[seq_name]
    clust_input_path = Path(out_dir, f"{prefix}_clust_input.fasta")
    dict_to_fasta(clust_input, clust_input_path)
    msg = (
        f"Retained {len(clust_input)} ({len(clust_input) / len(full_clust_input):.2%})"
        f" sequences after filtering by 'WSCORE_PROPORTION' of {wscore_proportion:.2f}"
        f" and 'COVERAGE_PROPORTION' of {coverage_proportion:.2f}"
        f" and 'MIN_SEQ_LEN' of {min_seq_len}\n"
        "\n"
        "CLUSTERING...\n"
        "\n"
    )
    print(msg)
    with open(clust_log_file, "a") as log:
        log.write(msg)

    # 4. Cluster
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
        "--mask",
        f"{1}",
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
    with open(clust_log_file, "a") as log:
        log.write(f"Captus' MMseqs2 Command:\n  {' '.join(mmseqs_cmd)}\n\n")
    with open(clust_log_file, "a") as log:
        subprocess.run(mmseqs_cmd, stdout=log, stdin=log)

    all_seqs = Path(out_dir, f"{prefix}_all_seqs.fasta")
    cluster_tsv = Path(out_dir, f"{prefix}_cluster.tsv")
    rep_seq = Path(out_dir, f"{prefix}_rep_seq.fasta")

    # 5. Get max depth and length for each locus' clusters
    clusters = split_mmseqs_clusters_file(all_seqs)
    cluster_sizes = []
    loci = {}
    for cluster in clusters:
        cluster_sizes.append(len(cluster) / 2)
        locus = cluster[0].split()[0].split(REF_CLUSTER_SEP)[-1]
        wscore = 1.0
        if "wscore" in cluster[0]:
            wscore = float(cluster[0].split("wscore=")[-1].split("]")[0])
        cluster_info = {
            "size": len(cluster) / 2,
            "length": len(cluster[1]),
            "wscore": wscore,
        }
        if locus not in loci:
            loci[locus] = cluster_info
        else:
            if cluster_info["size"] > loci[locus]["size"]:
                loci[locus] = cluster_info
            elif cluster_info["size"] == loci[locus]["size"]:
                if cluster_info["wscore"] > loci[locus]["wscore"]:
                    loci[locus] = cluster_info
                elif cluster_info["wscore"] == loci[locus]["wscore"]:
                    if cluster_info["length"] > loci[locus]["length"]:
                        loci[locus] = cluster_info
                    elif cluster_info["length"] == loci[locus]["length"]:
                        loci[locus] = cluster_info
    msg = (
        "\n"
        f"TOTAL NUMBER OF CLUSTERS: {len(cluster_sizes)}\n"
        f"MEAN CLUSTER SIZE: {statistics.mean(cluster_sizes):.2f}\n"
        f"MEDIAN CLUSTER SIZE: {statistics.median(cluster_sizes):.2f}\n"
        f"STD DEV CLUSTER SIZE: {statistics.stdev(cluster_sizes):.2f}\n"
        f"MAX CLUSTER SIZE: {max(cluster_sizes):.0f}\n"
        f"SINGLETON CLUSTERS: {cluster_sizes.count(1)}\n"
        "\n"
    )
    print(msg)
    with open(clust_log_file, "a") as log:
        log.write(msg)

    # 6. Select and format cluster representatives
    loci_reps = {}
    for cluster in clusters:
        locus = cluster[0].split()[0].split(REF_CLUSTER_SEP)[-1]
        size = len(cluster) / 2
        length = len(cluster[1])
        seq_name = cluster[0].split()[0].replace(">", "")
        sequence = cluster[1]
        description = " ".join(cluster[0].split()[1:])
        description = f"[cluster_size={size:.0f}] {description}"
        seq_desc = {
            "sequence": sequence,
            "description": description,
        }
        if best_only is False:
            if (
                size >= loci[locus]["size"] * size_proportion
                and length >= loci[locus]["length"] * length_proportion
            ):
                if locus not in loci_reps:
                    loci_reps[locus] = {seq_name: seq_desc}
                else:
                    loci_reps[locus][seq_name] = seq_desc
        else:
            wscore = 1.0
            if "wscore" in cluster[0]:
                wscore = float(cluster[0].split("wscore=")[-1].split("]")[0])
            locus_rep = {seq_name: seq_desc}
            if size > loci[locus]["size"]:
                loci_reps[locus] = locus_rep
            elif size == loci[locus]["size"]:
                if wscore > loci[locus]["wscore"]:
                    loci_reps[locus] = locus_rep
                elif wscore == loci[locus]["wscore"]:
                    if length > loci[locus]["length"]:
                        loci_reps[locus] = locus_rep
                    elif length == loci[locus]["length"]:
                        loci_reps[locus] = locus_rep
    reps_per_locus = []
    min4 = 0
    with open(clust_log_file, "a") as log:
        msg = "Loci with more than 3 representatives:"
        print(msg)
        log.write(f"{msg}\n")
        for locus in sorted(loci_reps):
            reps_per_locus.append(len(loci_reps[locus]))
            if len(loci_reps[locus]) > 3:
                min4 += 1
                print(f"{locus}:")
                log.write(f"{locus}:\n")
                for seq_name in loci_reps[locus]:
                    print(f"{seq_name} {loci_reps[locus][seq_name]['description']}")
                    log.write(f"{seq_name} {loci_reps[locus][seq_name]['description']}\n")
        msg = f"A total of {min4} loci had more than 3 representative sequences"
        print(msg)
        log.write(f"{msg}\n")
    msg = (
        "\n"
        f"TOTAL LOCI REMAINING: {len(reps_per_locus)}\n"
        f"TOTAL SEQS IN TARGET FILE: {sum(reps_per_locus)}\n"
        f"MEAN REPRESENTATIVES PER LOCUS: {statistics.mean(reps_per_locus):.2f}\n"
        f"MEDIAN REPRESENTATIVES PER LOCUS: {statistics.median(reps_per_locus):.2f}\n"
        f"STD DEV REPRESENTATIVES PER LOCUS: {statistics.stdev(reps_per_locus):.2f}\n"
        f"MAX REPRESENTATIVES PER LOCUS: {max(reps_per_locus):.0f}\n"
        f"LOCI WITH A SINGLE REPRESENTATIVE: {reps_per_locus.count(1)}\n"
        "\n"
    )
    print(msg)
    with open(clust_log_file, "a") as log:
        log.write(msg)

    # 7. Save new target file
    new_targets_path = Path(out_dir, f"{prefix}.fasta")
    for locus in sorted(loci_reps):
        dict_to_fasta(loci_reps[locus], new_targets_path, append=True, sort=True)
    msg = f"NEW TARGET FILE SAVED TO: '{new_targets_path}'\n"
    print(msg)
    with open(clust_log_file, "a") as log:
        log.write(msg)

    # 8. Intermediate file cleanup
    clust_input_path.unlink()
    all_seqs.unlink()
    cluster_tsv.unlink()
    rep_seq.unlink()
    shutil.rmtree(tmp_path, ignore_errors=True)

    return new_targets_path


def main():
    parser = argparse.ArgumentParser(
        description="Create a new target file from your own data by clustering the aligned sequences"
        " and choosing the best representatives per locus",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
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
        "-s",
        "--stage",
        action="store",
        default="untrimmed",
        dest="stage",
        choices=["unaligned", "untrimmed", "trimmed"],
        help="Alignment stage",
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
        help="Paralog filter, this is ignored when '--stage' is 'unaligned'",
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
        "--exclude_samples",
        action="store",
        dest="exclude_samples",
        help="Comma-separated list of samples to exclude (no spaces)",
    )
    parser.add_argument(
        "-o",
        "--out_dir",
        action="store",
        default="./new_targets",
        dest="out_dir",
        help="Output directory name",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        action="store",
        default="new_targets",
        dest="prefix",
        help="Prefix for output files, maybe use the original name of target file e.g. Allium353",
    )
    parser.add_argument(
        "--min_seq_len",
        action="store",
        default=33,
        type=int,
        dest="min_seq_len",
        help="Minimum sequence length to include in clustering",
    )
    parser.add_argument(
        "-i",
        "--min_identity",
        action="store",
        default=70,
        type=float,
        dest="min_identity",
        help="Minimum identity percentage between sequences in a cluster",
    )
    parser.add_argument(
        "-c",
        "--min_coverage",
        action="store",
        default=20,
        type=float,
        dest="min_coverage",
        help="Any sequence in a cluster has to be at least this percent included in the length"
        " of the longest sequence in the cluster",
    )
    parser.add_argument(
        "-W",
        "--wscore_proportion",
        action="store",
        default=0.55,
        type=float,
        dest="wscore_proportion",
        help="Cluster sequences with at least this proportion of the wscore of the highest wscore"
        " in the locus",
    )
    parser.add_argument(
        "-C",
        "--coverage_proportion",
        action="store",
        default=0.55,
        type=float,
        dest="coverage_proportion",
        help="Cluster sequences with at least this proportion of the coverage of the highest wscore"
        " sequence in the locus",
    )
    parser.add_argument(
        "-L",
        "--length_proportion",
        action="store",
        default=0.75,
        type=float,
        dest="length_proportion",
        help="Keep cluster representatives with at least this proportion of the length of the"
        " longest cluster representative in the locus",
    )
    parser.add_argument(
        "-S",
        "--size_proportion",
        action="store",
        default=0.55,
        type=float,
        dest="size_proportion",
        help="Keep cluster representatives with least this proportion of the size of the largest"
        " cluster in the locus",
    )
    parser.add_argument(
        "-B",
        "--best_only",
        action="store_true",
        dest="best_only",
        help="Enable to include a single representative per locus in the target file",
    )
    parser.add_argument(
        "-t",
        "--threads",
        action="store",
        default="auto",
        dest="threads",
        help="Number of threads to use",
    )
    args = parser.parse_args()

    if args.stage == "unaligned":
        aln_dir = Path(
            args.captus_alignments_dir,
            CAPTUS_DIRS[args.stage],
            CAPTUS_DIRS[args.marker],
            CAPTUS_DIRS[args.format],
        )
    else:
        aln_dir = Path(
            args.captus_alignments_dir,
            CAPTUS_DIRS[args.stage],
            CAPTUS_DIRS[args.filter],
            CAPTUS_DIRS[args.marker],
            CAPTUS_DIRS[args.format],
        )

    if not aln_dir.is_dir():
        quit_with_error(f"'{aln_dir}' not found, verify this Captus alignment directory exists!")

    if not Path(args.out_dir).exists():
        try:
            Path(args.out_dir).mkdir(parents=True)
        except OSError:
            quit_with_error(f"Captus was unable to make the output directory {Path(args.out_dir)}")

    fasta_ext = "fna"
    if args.format == "AA":
        fasta_ext = "faa"
    include_references = False
    if args.filter.endswith("_w_refs"):
        include_references = True

    fastas_paths = list(sorted(aln_dir.glob(f"*.{fasta_ext}")))

    select_refs_per_locus(
        fastas_paths,
        fasta_ext,
        args.out_dir,
        args.prefix,
        include_references,
        args.exclude_samples,
        args.min_seq_len,
        args.min_identity,
        args.min_coverage,
        args.wscore_proportion,
        args.coverage_proportion,
        args.length_proportion,
        args.size_proportion,
        args.best_only,
        args.threads,
    )


if __name__ == "__main__":
    main()
