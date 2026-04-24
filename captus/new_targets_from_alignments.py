#!/usr/bin/env python3
"""
Copyright 2020-2026 Edgardo M. Ortiz (e.ortiz.v@gmail.com)
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
import datetime
import gzip
import os
import random
import re
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
    "MF": "02_matches_flanked",
}

# Locus name separator in target files
REF_CLUSTER_SEP = "-"

# Separator between sample name and paralog number
SEQ_NAME_SEP = "__"

# Formatting
IND = " " * 2
END_FORMATTING = "\033[0m"
BOLD = "\033[1m"
DIM = "\033[2m"


def bold(text):
    return f"{BOLD}{text}{END_FORMATTING}"


def dim(text):
    return f"{DIM}{text}{END_FORMATTING}"


def remove_formatting(text):
    return re.sub(r"\033.*?m", r"", text)


def wlog(log_file: Path, message: str):
    with open(log_file, "at") as log:
        log.write(f"{remove_formatting(message)}\n")
    print(message)
    return


def now():
    now = dim(f"{datetime.datetime.now()}")
    return now


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


def msg_loci_stats(loci_fastas: dict, title: str):
    spacer = (26 - len(title)) * " "

    if not loci_fastas:
        msg = f"\n{IND}{spacer}{bold(title)}\n\n"
        msg += f"{IND}               Total loci: 0\n"
        return msg

    num_loci = len(loci_fastas)
    seqs_per_locus = [len(loci_fastas[locus]) for locus in loci_fastas]
    seqs_sd = 0
    if len(seqs_per_locus) > 1:
        seqs_sd = statistics.stdev(seqs_per_locus)

    msg = f"\n{IND}{spacer}{bold(title)}\n\n"
    msg += f"{IND}               Total loci: {num_loci}\n"
    msg += f"{IND}     Total singleton loci: {seqs_per_locus.count(1)}\n"
    msg += f"{IND}          Total sequences: {sum(seqs_per_locus)}\n\n"
    msg += f"{IND}    Median seqs per locus: {statistics.median(seqs_per_locus)}\n"
    msg += f"{IND}      Mean seqs per locus: {statistics.mean(seqs_per_locus):.3f}\n"
    msg += f"{IND}        SD seqs per locus: {seqs_sd:.3f}\n"
    msg += f"{IND}       Min seqs per locus: {min(seqs_per_locus)}\n"
    msg += f"{IND}       Max seqs per locus: {max(seqs_per_locus)}\n"

    return msg


def msg_samples_copies(loci_data: dict):
    if not loci_data:
        return ""

    num_samples = [loci_data[locus]["num_samples"] for locus in loci_data]
    samples_sd = 0
    if len(num_samples) > 1:
        samples_sd = statistics.stdev(num_samples)

    avg_copies = [loci_data[locus]["avg_copies"] for locus in loci_data]
    copies_sd = 0
    if len(avg_copies) > 1:
        copies_sd = statistics.stdev(avg_copies)

    msg = f"{IND} Median samples per locus: {statistics.median(num_samples)}\n"
    msg += f"{IND}   Mean samples per locus: {statistics.mean(num_samples):.3f}\n"
    msg += f"{IND}     SD samples per locus: {samples_sd:.3f}\n"
    msg += f"{IND}    Min samples per locus: {min(num_samples)}\n"
    msg += f"{IND}    Max samples per locus: {max(num_samples)}\n\n"
    msg += f"{IND}  Median copies per locus: {statistics.median(avg_copies):.3f}\n"
    msg += f"{IND}    Mean copies per locus: {statistics.mean(avg_copies):.3f}\n"
    msg += f"{IND}      SD copies per locus: {copies_sd:.3f}\n"
    msg += f"{IND}     Min copies per locus: {min(avg_copies):.3f}\n"
    msg += f"{IND}     Max copies per locus: {max(avg_copies):.3f}\n"

    return msg


def calc_per_locus_data(centroids: dict):
    centroids_data = {}
    for locus in centroids:
        num_seqs = 0
        samples = []
        inclusions = []
        for centroid in centroids[locus]:
            num_seqs += centroid["num_seqs"]
            samples += centroid["samples"]
            inclusions += centroid["includes"]
        num_samples = len(set(samples))
        centroids_data[locus] = {
            "num_seqs": num_seqs,
            "num_samples": num_samples,
            "avg_copies": num_seqs / num_samples,
            "includes": ",".join(sorted(list(set(inclusions)))),
        }

    return centroids_data


def get_wscore_coverage(description: str):
    wscore, coverage = None, None

    if "[wscore" in description and "[cover" in description:
        wscore = float(description.split("[wscore=")[1].split("]")[0])
        if "[cover=" in description:
            coverage = float(description.split("[cover=")[1].split("]")[0])
        elif "[coverage=" in description:
            coverage = float(description.split("[coverage=")[1].split("]")[0])

    return wscore, coverage


def prefilter_seqs(
    fastas_paths: list,
    out_dir: Path,
    prefix: str,
    exclude_samples,
    min_samples: int,
    min_seq_len: int,
    min_wscore_prop: float,
    min_coverage_prop: float,
    overwrite: bool,
    log: Path,
):

    def get_locus_name(fasta_path: str):
        fasta_name_parts = Path(fasta_path).name.split(".")

        if fasta_name_parts[-1].lower() == "gz":
            locus_name = ".".join(fasta_name_parts[:-2])
        else:
            locus_name = ".".join(fasta_name_parts[:-1])

        return locus_name

    if exclude_samples is None:
        exclude_samples = []
    else:
        exclude_samples = [str(x) for x in exclude_samples.split(",")]

    clust_input_file = Path(out_dir, f"{prefix}_clust_input.fasta")
    if clust_input_file.is_file():
        if overwrite is True:
            clust_input_file.unlink()
        else:
            quit_with_error(
                f"File '{clust_input_file}' already exists, use '--overwrite' to replace previous results"
            )

    loci_fastas = {}
    loci_data = {}
    loci_fastas_raw = {}

    for fasta_path in fastas_paths:
        samples = []
        locus_name = get_locus_name(fasta_path)
        locus_wscores = []
        locus_coverages = []
        locus_fasta_in = fasta_to_dict(fasta_path)
        loci_fastas_raw[locus_name] = locus_fasta_in
        locus_fasta_out = {}
        for seq_name in locus_fasta_in:
            sample_name = seq_name.split(SEQ_NAME_SEP)[0]
            samples.append(sample_name)
            seq = locus_fasta_in[seq_name]["sequence"].replace("-", "").replace("n", "")
            if len(seq) >= min_seq_len and sample_name not in exclude_samples:
                locus_fasta_out[f"{seq_name}{REF_CLUSTER_SEP}{locus_name}"] = {
                    "sequence": seq,
                    "description": locus_fasta_in[seq_name]["description"],
                }
                wscore, coverage = get_wscore_coverage(locus_fasta_in[seq_name]["description"])
                if wscore:
                    locus_wscores.append(wscore)
                if coverage:
                    locus_coverages.append(coverage)
        if locus_fasta_out:
            loci_fastas[locus_name] = locus_fasta_out
            if locus_wscores and locus_coverages:
                loci_data[locus_name] = {
                    "median_wscore": statistics.median(locus_wscores),
                    "median_coverage": statistics.median(locus_coverages),
                }
        if locus_name in loci_fastas and len(set(samples)) < min_samples:
            del loci_fastas[locus_name]
            del loci_data[locus_name]

    loci_fastas_raw_data = {}
    for locus in loci_fastas_raw:
        samples = []
        for seq_name in loci_fastas_raw[locus]:
            sample_name = seq_name.split(SEQ_NAME_SEP)[0]
            samples.append(sample_name)
        num_seqs = len(loci_fastas_raw[locus])
        num_samples = len(set(samples))
        loci_fastas_raw_data[locus] = {
            "num_seqs": num_seqs,
            "num_samples": num_samples,
            "avg_copies": num_seqs / num_samples,
        }
    msg = msg_loci_stats(loci_fastas_raw, "RAW STATS:")
    msg += "\n"
    msg += msg_samples_copies(loci_fastas_raw_data)
    wlog(log, msg)

    msg = (
        f"{now()} | Pre-filtering sequences "
        f"(exclude_samples={exclude_samples},"
        f" min_samples={min_samples},"
        f" min_seq_len={min_seq_len},"
        f" min_wscore_proportion={min_wscore_prop},"
        f" min_coverage_proportion={min_coverage_prop})"
    )
    wlog(log, msg)

    for locus_name in loci_fastas:
        locus_fasta_in = loci_fastas[locus_name]
        locus_fasta_out = {}
        for seq_name in locus_fasta_in:
            wscore, coverage = get_wscore_coverage(locus_fasta_in[seq_name]["description"])
            if wscore:
                if (
                    wscore >= loci_data[locus_name]["median_wscore"] * min_wscore_prop
                    and coverage >= loci_data[locus_name]["median_coverage"] * min_coverage_prop
                ):
                    locus_fasta_out[seq_name] = locus_fasta_in[seq_name]
            else:
                locus_fasta_out[seq_name] = locus_fasta_in[seq_name]
        if locus_fasta_out:
            loci_fastas[locus_name] = locus_fasta_out
        else:
            del loci_fastas[locus_name]

    if loci_fastas:
        for locus_name in sorted(loci_fastas):
            dict_to_fasta(loci_fastas[locus_name], clust_input_file, append=True)
        if clust_input_file.is_file():
            clust_input_data = {}
            for locus in loci_fastas:
                samples = []
                for seq_name in loci_fastas[locus]:
                    sample_name = seq_name.split(REF_CLUSTER_SEP)[0].split(SEQ_NAME_SEP)[0]
                    samples.append(sample_name)
                num_seqs = len(loci_fastas[locus])
                num_samples = len(set(samples))
                clust_input_data[locus] = {
                    "num_seqs": num_seqs,
                    "num_samples": num_samples,
                    "avg_copies": num_seqs / num_samples,
                }
            msg = msg_loci_stats(loci_fastas, "BEFORE CLUSTERING STATS:")
            msg += "\n"
            msg += msg_samples_copies(clust_input_data)
            wlog(log, msg)

            msg = f"{now()} | Pre-filtered sequences saved to '{clust_input_file}' for clustering"
            wlog(log, msg)

            return loci_fastas_raw_data, clust_input_data, clust_input_file
        else:
            quit_with_error(
                f"File '{clust_input_file}' was empty, relax your pre-filtering parameters and try again."
            )
    else:
        quit_with_error(
            "No sequences passed the initial filters, relax your pre-filtering parameters and try again."
        )


def cluster_seqs(
    clust_input_file: Path,
    clust_input_format: str,
    out_dir: Path,
    prefix: str,
    min_identity: float,
    min_coverage: float,
    threads,
    overwrite: bool,
    log: Path,
):
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
        f"{clust_input_file}",
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
    if clust_input_format == "NT":
        mmseqs_cmd += [
            "--gap-open",
            f"{3}",
        ]
    elif clust_input_format == "AA":
        mmseqs_cmd += [
            "--gap-open",
            f"{9}",
        ]

    clust_log_file = Path(out_dir, f"{prefix}_mmseqs.log")
    all_seqs = Path(out_dir, f"{prefix}_all_seqs.fasta")
    cluster_tsv = Path(out_dir, f"{prefix}_cluster.tsv")
    rep_seq = Path(out_dir, f"{prefix}_rep_seq.fasta")
    output_files = [clust_log_file, all_seqs, cluster_tsv, rep_seq]

    output_exists = any(of.is_file() for of in output_files)
    if output_exists:
        if overwrite is True:
            for of in output_files:
                if of.is_file():
                    of.unlink()
        else:
            quit_with_error(
                "Clustering output files exist, use '--overwrite' to replace previous results"
            )

    msg = f"{now()} | Clustering with MMseqs2 (min_identity={min_identity}, min_coverage={min_coverage})"
    wlog(log, msg)

    with open(clust_log_file, "a") as log:
        log.write(f"Captus' MMseqs2 Command:\n  {' '.join(mmseqs_cmd)}\n\n")
    with open(clust_log_file, "a") as log:
        subprocess.run(mmseqs_cmd, stdout=log, stdin=log)

    clust_output = split_mmseqs_clusters_file(all_seqs)

    if clust_output:
        clust_input_file.unlink()
        all_seqs.unlink()
        cluster_tsv.unlink()
        rep_seq.unlink()
        shutil.rmtree(tmp_path, ignore_errors=True)
        return clust_output
    else:
        quit_with_error("The clustering failed, delete the intermediate output files and try again.")


def select_targets(
    clust_output: dict,
    min_depth_prop: float,
    min_length_prop: float,
    max_target_copies: float,
    best_only: bool,
    split_paralogs: bool,
    min_avg_copies: float,
    min_samples_prop: float,
    log: Path,
):
    def keep_best_centroids(centroids: dict):
        best_centroids = {}
        for locus in centroids:
            best_centroids[locus] = [centroids[locus][0]]
        return best_centroids

    def remove_high_copy_centroids(centroids: dict, max_target_copies: float):
        low_copy_centroids = {}
        for locus in centroids:
            for centroid in centroids[locus]:
                if centroid["avg_copies"] <= max_target_copies:
                    if locus in low_copy_centroids:
                        low_copy_centroids[locus].append(centroid)
                    else:
                        low_copy_centroids[locus] = [centroid]
        return low_copy_centroids

    raw_centroids = {}

    for cluster in clust_output:
        centroid_locus = cluster[0].split()[0].split(REF_CLUSTER_SEP)[-1]
        centroid_sample = (
            cluster[0][1:].split()[0].replace(f"-{centroid_locus}", "").split(SEQ_NAME_SEP)[0]
        )
        includes = []
        samples = [centroid_sample]
        for i in range(2, len(cluster), 2):
            member_locus = cluster[i].split()[0].split(REF_CLUSTER_SEP)[-1]
            if member_locus != centroid_locus:
                includes.append(member_locus)
            sample_name = (
                cluster[i][1:].split()[0].replace(f"-{member_locus}", "").split(SEQ_NAME_SEP)[0]
            )
            samples.append(sample_name)

        header_pieces = cluster[0].split()
        seq_name = header_pieces[0].replace(">", "")
        num_seqs = int(len(cluster) / 2)
        samples = sorted(list(set(samples)))
        num_samples = len(samples)

        description = ""
        if len(header_pieces) > 1:
            description = " ".join(header_pieces[1:])
        description = f"[seqs={num_seqs:.0f}] [samples={num_samples:.0f}] {description}"

        wscore, coverage = get_wscore_coverage(description)

        centroid_info = {
            "seq_name": seq_name,
            "description": description,
            "sequence": cluster[1],
            "num_seqs": num_seqs,
            "num_samples": num_samples,
            "avg_copies": num_seqs / num_samples,
            "samples": samples,
            "length": len(cluster[1]),
            "wscore": wscore,
            "coverage": coverage,
            "includes": sorted(list(set(includes))),
        }

        if centroid_locus in raw_centroids:
            raw_centroids[centroid_locus].append(centroid_info)
        else:
            raw_centroids[centroid_locus] = [centroid_info]

    for locus in raw_centroids:
        raw_centroids[locus] = sorted(
            raw_centroids[locus], key=lambda x: (x["num_seqs"], x["length"], x["wscore"]), reverse=True
        )

    centroids = {}
    for locus in raw_centroids:
        centroids[locus] = [raw_centroids[locus][0]]
        if len(raw_centroids[locus]) > 1:
            top_centroid_depth = raw_centroids[locus][0]["num_seqs"]
            top_centroid_length = raw_centroids[locus][0]["length"]
            for i in range(1, len(raw_centroids[locus])):
                if (
                    raw_centroids[locus][i]["num_seqs"] >= top_centroid_depth * min_depth_prop
                    and raw_centroids[locus][i]["length"] >= top_centroid_length * min_length_prop
                ):
                    centroids[locus].append(raw_centroids[locus][i])

    clust_output_data = calc_per_locus_data(centroids)
    msg = msg_loci_stats(centroids, "AFTER CLUSTERING STATS:")
    msg += "\n"
    msg += msg_samples_copies(clust_output_data)
    wlog(log, msg)

    msg = (
        f"{now()} | Selecting new target per locus "
        f"(min_length_proportion={min_length_prop},"
        f" min_depth_proportion={min_depth_prop},"
        f" max_target_copies={max_target_copies},"
        f" best_only={best_only},"
        f" split_paralogs={split_paralogs})"
    )
    wlog(log, msg)

    if split_paralogs is True:
        msg = (
            f"{now()} | Paralog loci will be split if possible "
            f"(min_average_copies={min_avg_copies},"
            f" min_samples_proportion={min_samples_prop})"
        )
        wlog(log, msg)

    single_copy_centroids = {}
    multi_copy_centroids = {}
    if split_paralogs is True:
        for locus in clust_output_data:
            if clust_output_data[locus]["avg_copies"] <= min_avg_copies:
                single_copy_centroids[locus] = centroids[locus]
            else:
                if len(centroids[locus]) == 1:
                    multi_copy_centroids[locus] = centroids[locus]
                else:
                    paralog_suffix = 1
                    for centroid in centroids[locus]:
                        if (
                            centroid["num_samples"]
                            >= clust_output_data[locus]["num_samples"] * min_samples_prop
                        ):
                            new_locus = f"{locus}.PAR{paralog_suffix}"
                            centroid["seq_name"] = f"{centroid['seq_name']}.PAR{paralog_suffix}"
                            multi_copy_centroids[new_locus] = [centroid]
                            paralog_suffix += 1
                        else:
                            new_locus = f"{locus}.PAR0"
                            centroid["seq_name"] = f"{centroid['seq_name']}.PAR0"
                            if new_locus in multi_copy_centroids:
                                multi_copy_centroids[new_locus].append(centroid)
                            else:
                                multi_copy_centroids[new_locus] = [centroid]

        single_copy_centroids = remove_high_copy_centroids(single_copy_centroids, max_target_copies)
        multi_copy_centroids = remove_high_copy_centroids(multi_copy_centroids, max_target_copies)

        if best_only is True:
            single_copy_centroids = keep_best_centroids(single_copy_centroids)
            multi_copy_centroids = keep_best_centroids(multi_copy_centroids)

        single_copy_data = calc_per_locus_data(single_copy_centroids)
        multi_copy_data = calc_per_locus_data(multi_copy_centroids)

        if best_only is True:
            msg = bold("\n   AFTER SPLITTING PARALOGS")
            msg += bold("\n    AND RETAINING BEST ONLY:\n")
        else:
            msg = bold("\n   AFTER SPLITTING PARALOGS:\n")
        msg += msg_loci_stats(single_copy_centroids, "SINGLE-COPY LOCI:")
        msg += "\n"
        msg += msg_samples_copies(single_copy_data)
        msg += msg_loci_stats(multi_copy_centroids, "MULTI-COPY LOCI:")
        msg += "\n"
        msg += msg_samples_copies(multi_copy_data)
        wlog(log, msg)

        return single_copy_centroids, multi_copy_centroids
    else:
        centroids = remove_high_copy_centroids(centroids)
        if best_only is True:
            centroids = keep_best_centroids(centroids)

            centroids_data = calc_per_locus_data(centroids)

            msg = msg_loci_stats(centroids, "AFTER RETAINING BEST ONLY:")
            msg += "\n"
            msg += msg_samples_copies(centroids_data)
            wlog(log, msg)

        return centroids, None


def write_targets(
    raw_input_data: dict,
    clust_input_data: dict,
    single_copy_centroids: dict,
    multi_copy_centroids: dict,
    out_dir: Path,
    prefix: str,
    overwrite: bool,
    log: Path,
):
    target_file = Path(out_dir, f"{prefix}_targets.fasta")
    tsv_file = Path(out_dir, f"{prefix}_targets.tsv")
    output_files = [target_file, tsv_file]

    output_exists = any(of.is_file() for of in output_files)
    if output_exists:
        if overwrite is True:
            for of in output_files:
                if of.is_file():
                    of.unlink()
        else:
            quit_with_error("Target output files exist, use '--overwrite' to replace previous results")

    single_copy_data = calc_per_locus_data(single_copy_centroids)
    if multi_copy_centroids:
        multi_copy_data = calc_per_locus_data(multi_copy_centroids)

    target_data = {}
    loci_retained = []
    clu_fil_pcts = []
    clu_raw_pcts = []
    fil_raw_pcts = []

    for locus in single_copy_centroids:
        loci_retained.append(locus)
        clu_fil_samples = (
            single_copy_data[locus]["num_samples"] / clust_input_data[locus]["num_samples"]
        ) * 100
        clu_raw_samples = (
            single_copy_data[locus]["num_samples"] / raw_input_data[locus]["num_samples"]
        ) * 100
        fil_raw_samples = (
            clust_input_data[locus]["num_samples"] / raw_input_data[locus]["num_samples"]
        ) * 100
        clu_fil_pcts.append(clu_fil_samples)
        clu_raw_pcts.append(clu_raw_samples)
        fil_raw_pcts.append(fil_raw_samples)
        target_data[locus] = {
            "raw_locus": locus,
            "raw_seqs": raw_input_data[locus]["num_seqs"],
            "raw_samples": raw_input_data[locus]["num_samples"],
            "raw_avg_copies": f"{raw_input_data[locus]['avg_copies']:.3f}",
            "filtered_seqs": clust_input_data[locus]["num_seqs"],
            "filtered_samples": clust_input_data[locus]["num_samples"],
            "filtered_avg_copies": f"{clust_input_data[locus]['avg_copies']:.3f}",
            "locus": locus,
            "targets": len(single_copy_centroids[locus]),
            "sequences": single_copy_data[locus]["num_seqs"],
            "samples": single_copy_data[locus]["num_samples"],
            "avg_copies": f"{single_copy_data[locus]['avg_copies']:.3f}",
            "includes": single_copy_data[locus]["includes"],
            "clu_fil_samples": f"{clu_fil_samples:.3f}",
            "clu_raw_samples": f"{clu_raw_samples:.3f}",
            "fil_raw_samples": f"{fil_raw_samples:.3f}",
            "split": "FALSE",
        }

    if multi_copy_centroids:
        for locus in multi_copy_centroids:
            raw_locus = locus
            if len(locus.split(".")) > 1:
                suffix = locus.split(".")[-1]
                if suffix == "MINOR" or suffix.startswith("PAR"):
                    raw_locus = ".".join(locus.split(".")[:-1])
            loci_retained.append(raw_locus)
            clu_fil_samples = (
                multi_copy_data[locus]["num_samples"] / clust_input_data[raw_locus]["num_samples"]
            ) * 100
            clu_raw_samples = (
                multi_copy_data[locus]["num_samples"] / raw_input_data[raw_locus]["num_samples"]
            ) * 100
            fil_raw_samples = (
                clust_input_data[raw_locus]["num_samples"] / raw_input_data[raw_locus]["num_samples"]
            ) * 100
            clu_fil_pcts.append(clu_fil_samples)
            clu_raw_pcts.append(clu_raw_samples)
            fil_raw_pcts.append(fil_raw_samples)
            target_data[locus] = {
                "raw_locus": raw_locus,
                "raw_seqs": raw_input_data[raw_locus]["num_seqs"],
                "raw_samples": raw_input_data[raw_locus]["num_samples"],
                "raw_avg_copies": f"{raw_input_data[raw_locus]['avg_copies']:.3f}",
                "filtered_seqs": clust_input_data[raw_locus]["num_seqs"],
                "filtered_samples": clust_input_data[raw_locus]["num_samples"],
                "filtered_avg_copies": f"{clust_input_data[raw_locus]['avg_copies']:.3f}",
                "locus": locus,
                "targets": len(multi_copy_centroids[locus]),
                "sequences": multi_copy_data[locus]["num_seqs"],
                "samples": multi_copy_data[locus]["num_samples"],
                "avg_copies": f"{multi_copy_data[locus]['avg_copies']:.3f}",
                "includes": multi_copy_data[locus]["includes"],
                "clu_fil_samples": f"{clu_fil_samples:.3f}",
                "clu_raw_samples": f"{clu_raw_samples:.3f}",
                "fil_raw_samples": f"{fil_raw_samples:.3f}",
                "split": "FALSE",
            }

    loci_retained = list(set(loci_retained))

    for locus in clust_input_data:
        if locus not in loci_retained:
            fil_raw_samples = (
                clust_input_data[locus]["num_samples"] / raw_input_data[locus]["num_samples"]
            ) * 100
            # clu_fil_pcts.append(0)
            # clu_raw_pcts.append(0)
            fil_raw_pcts.append(fil_raw_samples)
            target_data[locus] = {
                "raw_locus": locus,
                "raw_seqs": raw_input_data[locus]["num_seqs"],
                "raw_samples": raw_input_data[locus]["num_samples"],
                "raw_avg_copies": f"{raw_input_data[locus]['avg_copies']:.3f}",
                "filtered_seqs": clust_input_data[locus]["num_seqs"],
                "filtered_samples": clust_input_data[locus]["num_samples"],
                "filtered_avg_copies": f"{clust_input_data[locus]['avg_copies']:.3f}",
                "locus": "NA",
                "targets": "NA",
                "sequences": "NA",
                "samples": "NA",
                "avg_copies": "NA",
                "includes": "NA",
                "clu_fil_samples": f"{0:.3f}",
                "clu_raw_samples": f"{0:.3f}",
                "fil_raw_samples": f"{fil_raw_samples:.3f}",
                "split": "FALSE",
            }

    for locus in raw_input_data:
        if locus not in clust_input_data:
            # clu_fil_pcts.append(0)
            # clu_raw_pcts.append(0)
            # fil_raw_pcts.append(0)
            target_data[locus] = {
                "raw_locus": locus,
                "raw_seqs": raw_input_data[locus]["num_seqs"],
                "raw_samples": raw_input_data[locus]["num_samples"],
                "raw_avg_copies": f"{raw_input_data[locus]['avg_copies']:.3f}",
                "filtered_seqs": "NA",
                "filtered_samples": "NA",
                "filtered_avg_copies": "NA",
                "locus": "NA",
                "targets": "NA",
                "sequences": "NA",
                "samples": "NA",
                "avg_copies": "NA",
                "includes": "NA",
                "clu_fil_samples": f"{0:.3f}",
                "clu_raw_samples": f"{0:.3f}",
                "fil_raw_samples": f"{0:.3f}",
                "split": "FALSE",
            }

    tsv_header = [
        "final_locus_name",
        "targets",
        "sequences",
        "samples",
        "avg_copies",
        "includes",
        "initial_locus_name",
        "initial_seqs",
        "initial_samples",
        "initial_avg_copies",
        "filter_seqs",
        "filter_samples",
        "filter_avg_copies",
        "filter_to_initial_samples",
        "final_to_initial_samples",
        "final_to_filter_samples",
        "split",
    ]

    target_fasta = {}
    num_split_loci = 0
    seq_numbers = []
    target_lengths = []

    with open(tsv_file, "wt") as tsv_out:
        tsv_out.write("\t".join(tsv_header) + "\n")
        for locus in sorted(target_data):
            record = [
                f"{target_data[locus]['locus']}",
                f"{target_data[locus]['targets']}",
                f"{target_data[locus]['sequences']}",
                f"{target_data[locus]['samples']}",
                f"{target_data[locus]['avg_copies']}",
                f"{target_data[locus]['includes']}",
                f"{target_data[locus]['raw_locus']}",
                f"{target_data[locus]['raw_seqs']}",
                f"{target_data[locus]['raw_samples']}",
                f"{target_data[locus]['raw_avg_copies']}",
                f"{target_data[locus]['filtered_seqs']}",
                f"{target_data[locus]['filtered_samples']}",
                f"{target_data[locus]['filtered_avg_copies']}",
                f"{target_data[locus]['fil_raw_samples']}",
                f"{target_data[locus]['clu_raw_samples']}",
                f"{target_data[locus]['clu_fil_samples']}",
                f"{target_data[locus]['split']}",
            ]

            if locus in single_copy_centroids:
                seq_numbers.append(len(single_copy_centroids[locus]))
                for centroid in single_copy_centroids[locus]:
                    new_seq_name = centroid["seq_name"].replace(SEQ_NAME_SEP, "_H")
                    target_fasta[new_seq_name] = {
                        "sequence": centroid["sequence"],
                        "description": centroid["description"],
                    }
                for centroid in single_copy_centroids[locus]:
                    target_lengths.append(len(centroid["sequence"]))
                    break
            if multi_copy_centroids and locus in multi_copy_centroids:
                seq_numbers.append(len(multi_copy_centroids[locus]))
                for centroid in multi_copy_centroids[locus]:
                    new_seq_name = centroid["seq_name"].replace(SEQ_NAME_SEP, "_H")
                    target_fasta[new_seq_name] = {
                        "sequence": centroid["sequence"],
                        "description": centroid["description"],
                    }
                for centroid in multi_copy_centroids[locus]:
                    target_lengths.append(len(centroid["sequence"]))
                    break
                if target_data[locus]["raw_locus"] != target_data[locus]["locus"]:
                    num_split_loci += 1
                    record[-1] = "TRUE"

            tsv_out.write("\t".join(record) + "\n")

    dict_to_fasta(target_fasta, target_file)

    msg = f"{now()} | New FASTA target file saved to '{target_file}'\n"
    msg += f"{now()} | Per-locus target information saved to '{tsv_file}'\n"
    wlog(log, msg)

    msg = bold("         % SAMPLES RETAINED:\n\n")
    msg += bold("             FILTERED / RAW:\n")
    msg += f"           FIL_RAW Median %: {statistics.median(fil_raw_pcts):.3f}\n"
    msg += f"             FIL_RAW Mean %: {statistics.mean(fil_raw_pcts):.3f}\n"
    msg += f"               FIL_RAW SD %: {statistics.stdev(fil_raw_pcts):.3f}\n"
    msg += f"              FIL_RAW Min %: {min(fil_raw_pcts):.3f}\n"
    msg += f"              FIL_RAW Max %: {max(fil_raw_pcts):.3f}\n"
    msg += "\n"
    msg += bold("                FINAL / RAW:\n")
    msg += f"           FIN_RAW Median %: {statistics.median(clu_raw_pcts):.3f}\n"
    msg += f"             FIN_RAW Mean %: {statistics.mean(clu_raw_pcts):.3f}\n"
    msg += f"               FIN_RAW SD %: {statistics.stdev(clu_raw_pcts):.3f}\n"
    msg += f"              FIN_RAW Min %: {min(clu_raw_pcts):.3f}\n"
    msg += f"              FIN_RAW Max %: {max(clu_raw_pcts):.3f}\n"
    msg += "\n"
    msg += bold("           FINAL / FILTERED:\n")
    msg += f"           FIN_FIL Median %: {statistics.median(clu_fil_pcts):.3f}\n"
    msg += f"             FIN_FIL Mean %: {statistics.mean(clu_fil_pcts):.3f}\n"
    msg += f"               FIN_FIL SD %: {statistics.stdev(clu_fil_pcts):.3f}\n"
    msg += f"              FIN_FIL Min %: {min(clu_fil_pcts):.3f}\n"
    msg += f"              FIN_FIL Max %: {max(clu_fil_pcts):.3f}\n"
    wlog(log, msg)

    msg = bold("\n          FINAL TARGET FILE:\n\n")
    msg += f"                 Total Loci: {len(seq_numbers)}\n"
    msg += f"       Total Singleton Loci: {seq_numbers.count(1)}\n"
    msg += f"           Total Split Loci: {num_split_loci}\n"
    msg += "\n"
    msg += f"            Total Sequences: {sum(seq_numbers)}\n"
    msg += f"      Median Seqs per Locus: {statistics.median(seq_numbers):.3f}\n"
    msg += f"        Mean Seqs per Locus: {statistics.mean(seq_numbers):.3f}\n"
    msg += f"          SD Seqs per Locus: {statistics.stdev(seq_numbers):.3f}\n"
    msg += f"         Min Seqs per Locus: {min(seq_numbers)}\n"
    msg += f"         Max Seqs per Locus: {max(seq_numbers)}\n"
    msg += "\n"
    msg += f"        Total Target Length: {sum(target_lengths)} bp\n"
    msg += f"       Median Target Length: {statistics.median(target_lengths):.3f} bp\n"
    msg += f"         Mean Target Length: {statistics.mean(target_lengths):.3f} bp\n"
    msg += f"           SD Target Length: {statistics.stdev(target_lengths):.3f} bp\n"
    msg += f"          Min Target Length: {min(target_lengths)} bp\n"
    msg += f"          Max Target Length: {max(target_lengths)} bp\n"
    wlog(log, msg)

    msg = f"{now()} | Finished successfully!\n"
    wlog(log, msg)

    return


def main():
    parser = argparse.ArgumentParser(
        description="Create a new target file from your own data by clustering the aligned sequences"
        " and choosing the best representatives per locus",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    captus_group = parser.add_argument_group("Captus alignments")
    captus_group.add_argument(
        "-a",
        "--captus_alignments_dir",
        action="store",
        dest="captus_alignments_dir",
        help="Path to the directory that contains the output from the alignment step of Captus",
    )
    captus_group.add_argument(
        "-s",
        "--stage",
        action="store",
        default="trimmed",
        dest="stage",
        choices=["unaligned", "untrimmed", "trimmed"],
        help="Alignment stage",
    )
    captus_group.add_argument(
        "-f",
        "--filter",
        action="store",
        default="unfiltered",
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
    captus_group.add_argument(
        "-M",
        "--marker",
        action="store",
        default="NUC",
        dest="marker",
        choices=["NUC", "PTD", "MIT", "DNA", "CLR"],
        help="Marker type",
    )
    captus_group.add_argument(
        "-F",
        "--format",
        action="store",
        default="NT",
        dest="format",
        choices=["AA", "NT", "GE", "GF", "MA", "MF"],
        help="Alignment data format",
    )

    alns_list_group = parser.add_argument_group("List of alignments")
    alns_list_group.add_argument(
        "-l",
        "--alignments_list",
        action="store",
        dest="alignments_list",
        help="Instead of searching within the Captus alignments directory you can provide a file"
        " with a list of paths to the alignments in FASTA format, one path per line (all the Captus-"
        "specific options like '--min_wscore_proportion' or '--min_coverage_proportion' will be"
        " ignored). These alignments are assumed to be in nucleotide, if you are providing peptide"
        " sequences please add '-F AA'",
    )

    prefiltering_group = parser.add_argument_group("Pre-filtering")
    prefiltering_group.add_argument(
        "--exclude_samples",
        action="store",
        dest="exclude_samples",
        help="Comma-separated list of samples to exclude (no spaces)",
    )
    prefiltering_group.add_argument(
        "--min_samples",
        action="store",
        default=4,
        type=int,
        dest="min_samples",
        help="Exclude locus if contains fewer samples than this number",
    )
    prefiltering_group.add_argument(
        "--min_seq_len",
        action="store",
        default=33,
        type=int,
        dest="min_seq_len",
        help="Minimum sequence length to include in clustering",
    )
    prefiltering_group.add_argument(
        "-W",
        "--min_wscore_proportion",
        action="store",
        default=0.66,
        type=float,
        dest="min_wscore_proportion",
        help="Only applies to alignments produced by Captus. The sequences contained within these"
        " alignments are annotated with a 'wscore' value. This filter removes sequences whose"
        " 'wscore' is not at least this proportion of the median 'wscore' in the alignment so they"
        " are not considered during clustering",
    )
    prefiltering_group.add_argument(
        "-C",
        "--min_coverage_proportion",
        action="store",
        default=0.66,
        type=float,
        dest="min_coverage_proportion",
        help="Only applies to alignments produced by Captus. The sequences contained within these"
        " alignments are annotated with a 'coverage' value. This filter removes sequences whose"
        " 'coverage' is not at least this proportion of the median 'coverage' in the alignment so"
        " they are not considered during clustering",
    )

    clustering_group = parser.add_argument_group("Clustering")
    clustering_group.add_argument(
        "-i",
        "--min_identity",
        action="store",
        default=70,
        type=float,
        dest="min_identity",
        help="Minimum identity percentage between sequences in a cluster",
    )
    clustering_group.add_argument(
        "-c",
        "--min_coverage",
        action="store",
        default=40,
        type=float,
        dest="min_coverage",
        help="Any sequence in a cluster has to be at least this percent included in the length"
        " of the longest sequence in the cluster",
    )

    targets_group = parser.add_argument_group("Target selection")
    targets_group.add_argument(
        "-L",
        "--min_length_proportion",
        action="store",
        default=0.75,
        type=float,
        dest="min_length_proportion",
        help="Keep cluster representatives with at least this proportion of the length of the"
        " longest cluster representative for the locus",
    )
    targets_group.add_argument(
        "-D",
        "--min_depth_proportion",
        action="store",
        default=0.55,
        type=float,
        dest="min_depth_proportion",
        help="Keep cluster representatives with least this proportion of the depth of the largest"
        " cluster for the locus, depth is defined as the number of sequences contained in a cluster",
    )
    targets_group.add_argument(
        "-B",
        "--best_only",
        action="store_true",
        dest="best_only",
        help="Enable to include a single representative per locus in the target file",
    )
    targets_group.add_argument(
        "-T",
        "--max_target_copies",
        action="store",
        default=1.33,
        type=float,
        dest="max_target_copies",
        help="Maximum average number of copies allowed per target in the final target file. Average"
        " copies are calculated dividing the total number of sequences by the total number of"
        " samples in the cluster that the target represents",
    )
    targets_group.add_argument(
        "-P", "--split_paralogs",
        action="store_true",
        dest="split_paralogs",
        help="Enable to split loci flagged as paralogous into separate loci. Only makes sense"
        " when using alignments that can have multiple copies per sample (multi-copy) like the"
        " 'unfiltered' alignments from Captus for example",
    )
    targets_group.add_argument(
        "-A",
        "--min_average_copies",
        action="store",
        default=1.33,
        type=float,
        dest="min_average_copies",
        help="Only applies to alignments with multiple sequences per sample (multi-copy). Average"
        " copies are calculated dividing the total number of sequences by the total number of"
        " samples. Alignments with more copies than this value will be flagged as paralogous",
    )
    targets_group.add_argument(
        "-S",
        "--min_samples_proportion",
        action="store",
        default=0.66,
        type=float,
        dest="min_samples_proportion",
        help="Only applies when '--split_paralogs' is enabled. For loci flagged as paralogous, only"
        " the cluster representatives that include at least this proportion of the total samples in"
        " the cluster will become the representative of a separate locus",
    )

    output_group = parser.add_argument_group("Output")
    output_group.add_argument(
        "-o",
        "--out_dir",
        action="store",
        default="./new_targets",
        dest="out_dir",
        help="Output directory name",
    )
    output_group.add_argument(
        "-p",
        "--prefix",
        action="store",
        default="new_targets",
        dest="prefix",
        help="Prefix for output files, e.g. Allium353",
    )

    other_group = parser.add_argument_group("Other")
    other_group.add_argument(
        "-t",
        "--threads",
        action="store",
        default="auto",
        dest="threads",
        help="Number of threads to use",
    )
    other_group.add_argument(
        "--overwrite",
        action="store_true",
        dest="overwrite",
        help="Enable to overwrite previous results",
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    if not Path(args.out_dir).exists():
        try:
            Path(args.out_dir).mkdir(parents=True)
        except OSError:
            quit_with_error(f"Captus was unable to make the output directory {Path(args.out_dir)}")

    log = Path(args.out_dir, f"{args.prefix}_ntfa.log")
    full_command = " ".join(sys.argv)
    msg = (
        f"{now()}\n"
        f"\n{bold('NEW TARGETS FROM ALIGNMENTS:')}\n"
        "Create a new target file selecting representative sequences from multiple loci alignments\n"
        f"\n{bold('COMMAND:')} {full_command}\n"
    )
    wlog(log, msg)

    fastas_paths = []

    if args.captus_alignments_dir:
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

        msg = f"{now()} | Verifying alignments in '{aln_dir}'..."
        wlog(log, msg)

        if not aln_dir.is_dir():
            quit_with_error(f"'{aln_dir}' not found, verify this Captus alignment directory exists!")

        fasta_ext = "fna"
        if args.format == "AA":
            fasta_ext = "faa"

        fastas_paths = list(sorted(aln_dir.glob(f"*.{fasta_ext}")))

    elif args.alignments_list:
        msg = f"{now()} | Verifying alignments in '{args.alignments_list}'..."
        wlog(log, msg)

        if not Path(args.alignments_list).exists():
            quit_with_error(f"'{args.alignments_list}' not found, verify this file exists!")

        with open(Path(args.alignments_list), "rt") as al:
            for line in al:
                fasta_path = Path(line.strip())
                if fasta_path.exists():
                    fastas_paths.append(fasta_path)
                else:
                    print(f"WARNING: the file '{fasta_path}' doesn't exist!")

    msg = f"{now()} | {len(fastas_paths)} alignments verified"
    wlog(log, msg)

    if not fastas_paths:
        quit_with_error(
            "No valid alignments were found, check your alignments directory or your list of alignments"
        )

    raw_input_data, clust_input_data, clust_input_file = prefilter_seqs(
        fastas_paths,
        Path(args.out_dir),
        args.prefix,
        args.exclude_samples,
        args.min_samples,
        args.min_seq_len,
        args.min_wscore_proportion,
        args.min_coverage_proportion,
        args.overwrite,
        log,
    )

    clust_output = cluster_seqs(
        clust_input_file,
        args.format,
        Path(args.out_dir),
        args.prefix,
        args.min_identity,
        args.min_coverage,
        args.threads,
        args.overwrite,
        log,
    )

    single_copy_centroids, multi_copy_centroids = select_targets(
        clust_output,
        args.min_depth_proportion,
        args.min_length_proportion,
        args.max_target_copies,
        args.best_only,
        args.split_paralogs,
        args.min_average_copies,
        args.min_samples_proportion,
        log,
    )

    write_targets(
        raw_input_data,
        clust_input_data,
        single_copy_centroids,
        multi_copy_centroids,
        Path(args.out_dir),
        args.prefix,
        args.overwrite,
        log,
    )


if __name__ == "__main__":
    main()
