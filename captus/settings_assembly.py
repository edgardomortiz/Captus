#!/usr/bin/env python3
"""
Copyright 2020 Edgardo M. Ortiz (e.ortiz.v@gmail.com)
https://github.com/edgardomortiz/Captus

This module contains hard-coded settings for Captus-assembly

This file is part of Captus. Captus is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Captus is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Captus. If
not, see <http://www.gnu.org/licenses/>.
"""

import platform
from pathlib import Path

# FASTA valid filename extensions:
FASTA_VALID_EXTENSIONS = [".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz"]

# FASTQ valid filename extensions:
FASTQ_VALID_EXTENSIONS = [".fq", ".fastq", ".fq.gz", ".fastq.gz"]

# Fraction of total RAM available to Captus when using 'auto' in --ram
RAM_FRACTION = 0.98

# Path to this file:
SETTINGS_ASSEMBLY_PATH = Path(__file__).resolve()

# Data directory
DATA_DIR = Path(SETTINGS_ASSEMBLY_PATH.parent.parent, "data")

# FASTA file of Illumina adaptors
ILLUMINA_ADAPTORS = Path(DATA_DIR, "adaptors_illumina.fasta")

# FASTA file of BGISEQ, DNBSEQ, MGISEQ adaptors
BGISEQ_ADAPTORS = Path(DATA_DIR, "adaptors_bgi-dnb-mgi.fasta")

# FASTA file of BGISEQ, DNBSEQ, MGISEQ adaptors
POLYA_ADAPTORS = Path(DATA_DIR, "polyA.fasta")

# FASTA reference genome of the phiX virus
PHIX_REF_GENOME = Path(DATA_DIR, "phix174_ill.ref.fa.gz")

# FASTA file of sequencing artifacts
SEQUENCING_ARTIFACTS = Path(DATA_DIR, "sequencing_artifacts.fasta")

# R Markdown to generate Sequence Quality Reports
CAPTUS_QC_RMD_REPORT = Path(DATA_DIR, "captus-assembly_clean.Rmd")

# Keep reads with a minimum length of this value after trimming, for the adaptor removal stage and
# the quality filtering/trimming stage
BBDUK_MIN_LENGTH = 21

# K-mer size for first round of adaptor cleaning
BBDUK_ADAPTOR_ROUND1_K = 21

# Minimum K-mer size for first round of adaptor cleaning
BBDUK_ADAPTOR_ROUND1_MINK = 11

# Hamming distance for first round of adaptor cleaning
BBDUK_ADAPTOR_ROUND1_HDIST = 2

# K-mer size for second round of adaptor cleaning
BBDUK_ADAPTOR_ROUND2_K = 19

# Hamming distance for second round of adaptor cleaning
BBDUK_ADAPTOR_ROUND2_HDIST = 1

# Minimum K-mer size for second round of adaptor cleaning
BBDUK_ADAPTOR_ROUND2_MINK = 9

# K-mer size for contaminant filtering during qulity filtering and trimming
BBDUK_QUALITY_K = 31

# Hamming distance for contaminant filtering during qulity filtering and trimming
BBDUK_QUALITY_HDIST = 1

# Trim reads' low-quality nucleotides from these directions l=left, r=right, lr=both
BBDUK_QUALITY_QTRIM = "lr"

# Maximum number of Ns allowed in a read
BBDUK_QUALITY_MAXNS = 5

# Compression level for final FASTQ files, lower is faster
BBDUK_QUALITY_ZIPLEVEL = 5

# Maximum number of concurrent FastQC instances, there may not be increase on speed after around 16
# simultaneous runs of FastQC because of hard drive limitations (may improve with SSDs)
FASTQC_MAX_INSTANCES = 16

# Compression level for subsampled FASTQ files
REFORMAT_ZIPLEVEL = 5

# Caculate the average read length of a FASTQ with this many reads
NUM_READS_TO_CALCULATE_MEAN_READ_LENGTH = 100000

# Include largest kmer sizes in MEGAHIT's 'k_list' if they are at most this larger than the mean
# read length of the input FASTQs
DELTA_MEAN_READ_LENGTH_TO_MAX_KMER_SIZE = 40

# MEGAHIT's '--bubble-level' parameter, apparently reducing the merging level of the bubbles can
# lead to longer contigs, MEGAHIT's default is 2, we try experimentally 1
MEGAHIT_BUBBLE_LEVEL = 1

# MEGAHIT's '--prune-level' parameter, since the data is not metagenomic we should use 3 to recover
# longer contigs
MEGAHIT_PRUNE_LEVEL = 3

# MEGAHIT's '--prune-depth' parameter, the default is sqrt(median_of_kmer_depth) when '--prune-level'
# is set to 3, this results in low coverage contigs being removed too aggresively, to recover as
# much assembled sequence as possible set to 1
MEGAHIT_PRUNE_DEPTH = 1

# Minimum RAM in bytes for a MEGAHIT assembly
MEGAHIT_MIN_RAM_B = 4 * 1024 ** 3

# Minimum number of threads for a MEGAHIT assembly
MEGAHIT_MIN_THREADS = 4

# Minimum kmer size to attempt FASTA to FASTG conversion, usually kmer sizes in the 20s produce too
# complex graphs that are hard to visualize
MIN_KMER_SIZE_FOR_FASTG = 31

# Minimum RAM reserved a BLAT instance in bytes (2GB)
BLAT_MIN_RAM_B = 2 * 1024 ** 3

# Separator in between contigs from the same Scipio's hit
SCIPIO_CONTIG_SEPARATOR = "n" * 50

# Separator between sample and protein cluster for Scipio's reference proteins. For example,
# Angiosperms353 uses a '-', meaning that anything before the '-' is the sample name and after the
# '-' is the protein cluster identifier
# For example, this is how the Angiopsmerms353 formats the FASTA headers in their reference proteins
# dataset:
# >WIGA-7577
# >WXVX-7577
# >XBKS-7577
# >XFJG-7577
# >AEPI-7583
# >AXNH-7583
# >Ambtr-7583
# >Arath-7583
# By the presence of the '-' we can tell that the first four (4) represent protein sequences from
# the same protein cluster '7577' that come from different species ('WIGA', 'WXVX', 'XBKS', 'XFJG').
# The same goes for the last four (4) which belong to protein cluster '7583'
REFERENCE_CLUSTER_SEPARATOR = "-"

# Bundled sets of protein references
PROT_REFS = {
    "NUC": {
        "angiosperms353": {
            "AA": Path(DATA_DIR, "Angiosperms353.FAA"),
            "NT": "",
        },
    },
    "PTD": {
        "angiospermsptd": {
            "AA": Path(DATA_DIR, "AngiospermsPTD.FAA"),
            "NT": "",
        },
    },
    "MIT": {
        "angiospermsmit": {
            "AA": Path(DATA_DIR, "AngiospermsMIT.FAA"),
            "NT": "",
        },
    }
}

# Bundled miscellaneous DNA references
DNA_REFS = {
    "": Path("")
}

# Bundled dependencies directory
DEPS_DIR = Path(Path(__file__).resolve().parent.parent, "dependencies")

# Bundled Scipio v1.4.1 path
BUNDLED_SCIPIO = Path(DEPS_DIR, "scipio-1.4", "scipio.1.4.1.pl")

# Bundled Scipio v1.4.1 path
BUNDLED_SCIPIO_YAML2GFF = Path(DEPS_DIR, "scipio-1.4", "yaml2gff.1.4.pl")

# Bundled BLAT >= v36x7 path
os_type = platform.system()
if os_type == "Darwin":
    BUNDLED_BLAT = Path(DEPS_DIR, "blat", "blat_36x9_mac")
elif os_type == "Linux":
    BUNDLED_BLAT = Path(DEPS_DIR, "blat", "blat_36x9_linux")
else:
    BUNDLED_BLAT = "unknown system"

# Extraction directories per marker type
MARKER_DIRS = {
    "NUC": "01_coding_NUC",         # Nuclear protein markers
    "PTD": "02_coding_PTD",         # Plastidial protein markers
    "MIT": "03_coding_MIT",         # Mitochondrial protein markers
    "DNA": "04_misc_DNA",           # Miscellaneous DNA markers
    "CLR": "05_clusters",           # Cluster-derived DNA markers
}

# Filename suffixes for the different extraction FASTA formats
FORMAT_SUFFIXES = {
    "AA": "_coding_AA.faa",         # coding sequences in aminoacids
    "NT": "_coding_NT.fna",         # coding sequences in nucleotides
    "GE": "_genes.fna",             # complete gene sequences in nucleotides
    "GF": "_genes_flanked.fna",     # complete gene sequences with additional flanking bp
    "MA": "_matches.fna",           # matching DNA sequences from miscellaneous DNA extraction
    "MF": "_matches_flanked.fna",   # matching DNA sequences with additional flanking bp
}

# Directory names for the different extraction FASTA formats
FORMAT_DIRS = {
    "AA": "01_AA",                  # coding sequences in aminoacids
    "NT": "02_NT",                  # coding sequences in nucleotides
    "GE": "03_genes",               # complete gene sequences in nucleotides
    "GF": "04_genes_flanked",       # complete gene sequences with additional flanking bp
    "MA": "01_matches",             # matching DNA sequences from miscellaneous DNA extraction
    "MF": "02_matches_flanked",     # matching DNA sequences with additional flanking bp
}

# Translated protein reference suffix
TRANSLATED_REF_SUFFIX = ".captus.faa"

# Valid combinations of marker directories and format directories
VALID_MARKER_FORMAT_COMBO = [(m, f) for m in ["NUC","PTD","MIT"] for f in ["AA","NT","GE","GF"]]
VALID_MARKER_FORMAT_COMBO += [(m, f) for m in ["DNA","CLR"] for f in ["MA","MF"]]

# Scipio's initial round will run with 'min_score' multiplied by this factor
SCIPIO_SCORE_INITIAL_FACTOR = 0.9

# Scipio's accepted intron penalty controls what types of intron borders are accepted:
# 1.0 = GT--AG, GC--AG
# 1.1 = GT--AG, GC--AG, AT--AC
# 1.2 = GT--AG, GC--AG, AT--AC, GG--AG, GA--AG
SCIPIO_ACCEPTED_INTRON_PENALTY = 1.1

# Minimum BLAT score for Scipio's runs when number of reference proteins is <= 'max_loci_scipio2x'
SCIPIO_2X_BLAT_MIN_SCORE = 15

# Minimum BLAT score for Scipio's only run when number of reference proteins is > 'max_loci_scipio2x'
SCIPIO_1X_BLAT_MIN_SCORE = 30

# BLAT identity as proportion of Scipio's 'min_identity' parameter
SCIPIO_BLAT_IDENTITY_FACTOR = 0.9

# Default Genetic Codes to set Scipio's --transtable
DEFAULT_GENETIC_CODES = {
    "NUC": {"id": 1, "name": "Standard"},
    "PTD": {"id": 11, "name": "Bacterial, Archaeal and Plant Plastid"},
    "MIT": {"id": 1, "name": "Standard"}
}

# Divergent genetic codes that will need a lower 'SCIPIO_BLAT_IDENTITY_FACTOR' for BLAT to find hits
# these genetic codes have aminoacid substitutions that may make it harder for BLAT to find a hit
DIVERGENT_GENETIC_CODES = [2, 3, 5, 9, 12, 13, 14, 21, 24, 26, 33]

# Scipio's maximum identity when using a divergent genetic code
SCIPIO_MAX_IDENTITY_DIV_CODE = 66

# Extra settings for the final round of Scipio according to the genome of the genes
SCIPIO_GENOME_SETTINGS = {
    # Change here the settings for nuclear genes:
    "NUC": [
        "--max_assemble_size=75000",
        "--region_size=1000",
        # "--blat_params=-oneOff=1",
        # "--blat_tilesize=6",
        "--exhaust_align_size=500",
        "--exhaust_gap_size=21",
        "--min_dna_coverage=0",
        "--max_move_exon=6",
        "--gap_to_close=21", # keep <=21 or it breaks reconstruction across several contigs
    ],
    # Change here the settings for plastidial genes:
    "PTD": [
        "--max_assemble_size=9000",
        "--region_size=0",
        "--blat_params=-oneOff=1",
        "--blat_tilesize=6",
        "--exhaust_align_size=9000",
        "--exhaust_gap_size=21",
        "--min_dna_coverage=0.2",
        "--max_move_exon=6",
        "--gap_to_close=21",
    ],
    # Change here the settings for mitochondrial genes:
    "MIT": [
        "--max_assemble_size=9000",
        "--region_size=0",
        "--blat_params=-oneOff=1",
        "--blat_tilesize=6",
        "--exhaust_align_size=9000",
        "--exhaust_gap_size=21",
        "--min_dna_coverage=0.2",
        "--max_move_exon=6",
        "--gap_to_close=21",
    ],
}

# Tolerance proportion for determining if two hits are compatible in their 'match_id' or percentage
# of identity to the reference, or determining if they have an acceptable margin of overlap. These
# two conditions are used to determine compatible pairs of hits during the assembly of partial hits
# during the extraction of non-coding markers
DNA_TOLERANCE_PROP = 0.05

# Absolute maximum overlap in bp between two adjacent partial hits, this value will truncate
# 'DNA_TOLERANCE_PROP'*'q_size' to not allow large overlaps when reference sequence is very long
# If the contig was assembled wthin Captus, the kmer size was recorded in the contig name and will
# be used as the absolute maximum overlap between contigs
DNA_MAX_OVERLAP_BP = 99

# Minimum coverage for non-coding hits before attempting greedy assembly
DNA_MIN_COVERAGE_BEFORE_ASSEMBLY = 10

# During non-coding sequence extraction, a version of the sequence matched with extra upstream and
# downstream nucleotides is extracted in case they help in alignment or provide useful extra data
# This value determines how many extra nucleotides to extract up- and/or downstream
DNA_UP_DOWN_STREAM_BP = 1000

# Maximum number of FASTA files to rehead and write simultaneously. Too many may hurt performance
# if HDD, it may improve with SSDs
MAX_WRITING_INSTANCES = 16

# Clustering identity percentage, if '--cl_min_identity' is left as 'auto' it becomes 90% of the
# '--dna_min_identity' value
MMSEQS2_BLAT_DNA_IDENTITY_FACTOR = 0.98

# Minimum clustering identity when set to 'auto'
MMSEQS_MIN_AUTO_MIN_IDENTITY = 75

# MMseqs sensitivity set the same as default 4 (experiment increasing a little?)
MMSEQS2_SENSITIVITY = 4.0

# Coverage mode for MMseqs2, parameter '--cov-mode':
# 0: coverage of query and target
# 1: coverage of target
# 2: coverage of query
# 3: target seq. length has to be at least x% of query length
# 4: query seq. length has to be at least x% of target length
# 5: short seq. needs to be at least x% of the other seq. length
#                  --cov-mode
#                  0    1    2
# Q: MAVGTACRPA  60%  IGN  60%
# T: -AVGTAC---  60% 100%  IGN
#        -c 0.7    -    +    -
#        -c 0.6    +    +    +
MMSEQS2_COV_MODE = 1

# From MMseqs2 help:
# Scale k-mer per sequence based on sequence length as kmer-per-seq val + scale x seqlen
# kmer-per-seq is 21 by default, so if we set this value to 0.3 for a 500 bp sequence, MMseqs2 will
# sample (500*0.3)+21 = 171 kmers from that sequence
MMSEQS2_KMER_PER_SEQ_SCALE = 0.3

# GFF feature colors, Okabe & Ito palette
GFF_COLORS = {
    "NUC": "#e69f00",  # orange
    "PTD": "#009e73",  # green
    "MIT": "#56B4e9",  # skyblue
    "DNA": "#d55e00",  # vermillion
    "CLR": "#cc79a7",  # purple
}

# MAFFT algorithms choices to real MAFFT syntax translation
MAFFT_ALGORITHMS = {
    "genafpair": {
        "arg": "--genafpair",
        "aka": "(a.k.a. E-INS-i)"
    },
    "localpair": {
        "arg": "--localpair",
        "aka": "(a.k.a. L-INS-i)"
    },
    "globalpair": {
        "arg": "--globalpair",
        "aka": "(a.k.a. G-INS-i)"
    },
    "retree1": {
        "arg": "--retree 1",
        "aka": "(a.k.a. FFT-NS-1)"
    },
    "retree2": {
        "arg": "--retree 2",
        "aka": "(a.k.a. FFT-NS-2)"
    },
}