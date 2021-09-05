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
RAM_FRACTION = 0.99

# Path to this file:
SETTINGS_ASSEMBLY_PATH = Path(__file__).resolve()

# Data directory
DATA_DIR = Path(SETTINGS_ASSEMBLY_PATH.parent.parent, "data")

# Code directory
CODE_DIR = Path(SETTINGS_ASSEMBLY_PATH.parent.parent, "captus")

# FASTA file of Illumina adaptors
ILLUMINA_ADAPTORS = Path(DATA_DIR, "adaptors_illumina.fasta")

# FASTA file of BGISEQ, DNBSEQ, MGISEQ adaptors
BGISEQ_ADAPTORS = Path(DATA_DIR, "adaptors_bgi-dnb-mgi.fasta")

# Combined FASTA file of Illumina, BGISEQ, DNBSEQ, MGISEQ adaptors
COMBINED_ADAPTORS = Path(DATA_DIR, "adaptors_combined.fasta")

# FASTA file of BGISEQ, DNBSEQ, MGISEQ adaptors
POLYA_ADAPTORS = Path(DATA_DIR, "polyA.fasta")

# FASTA reference genome of the phiX virus
PHIX_REF_GENOME = Path(DATA_DIR, "phix174_ill.ref.fa.gz")

# FASTA file of sequencing artifacts
SEQUENCING_ARTIFACTS = Path(DATA_DIR, "sequencing_artifacts.fasta")

# Adaptor list for FastQC or Falco
QC_ADAPTORS_LIST = Path(DATA_DIR, "qc_adaptors_list.txt")

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

# K-mer size for contaminant filtering during quality filtering and trimming
BBDUK_QUALITY_K = 31

# Hamming distance for contaminant filtering during quality filtering and trimming
BBDUK_QUALITY_HDIST = 1

# Trim reads' low-quality nucleotides from these directions l=left, r=right, lr=both
BBDUK_QUALITY_QTRIM = "lr"

# Maximum number of Ns allowed in a read
BBDUK_QUALITY_MAXNS = 5

# Compression level for final FASTQ files, lower is faster
BBDUK_QUALITY_ZIPLEVEL = 5

# Maximum number of concurrent FastQC or Falco instances, there may not be increase on speed after
# 16 simultaneous runs of FastQC or Falco because of hard drive limitations (may improve with SSDs)
QC_STATS_MAX_INSTANCES = 16

# Compression level for subsampled FASTQ files
REFORMAT_ZIPLEVEL = 5

# File names for QC extra tables with statistics
QC_FILES = {
    "REBA": "reads_bases.tsv",
    "ADR1": "adaptors_round1.tsv",
    "ADR2": "adaptors_round2.tsv",
    "CONT": "contaminants.tsv",
    "PBSQ": "per_base_seq_qual.tsv",
    "PSQS": "per_seq_qual_scores.tsv",
    "PBSC": "per_base_seq_content.tsv",
    "PSGC": "per_seq_gc_content.tsv",
    "SLEN": "seq_len_dist.tsv",
    "SDUP": "seq_dup_levels.tsv",
    "ADCO": "adaptor_content.tsv",
}

# Caculate the average read length of a FASTQ with this many reads
NUM_READS_TO_CALCULATE_MEAN_READ_LENGTH = 100000

# Include largest kmer sizes in MEGAHIT's 'k_list' if they are at most this larger than the mean
# read length of the input FASTQs
DELTA_MEAN_READ_LENGTH_TO_MAX_KMER_SIZE = 31

# Minimum RAM in bytes for a MEGAHIT assembly
MEGAHIT_MIN_RAM_B = 4 * 1024 ** 3

# Defaults for MEGAHIT, optimized for hybridization capture or genome skimming
MEGAHIT_K_LIST = "31,39,47,63,79,95,111,127,143,159,175"
MEGAHIT_MIN_COUNT = 2
MEGAHIT_PRUNE_LEVEL = 2

# Presets for MEGAHIT assemblies of RNAseq and WGS
MEGAHIT_PRESETS = {
    "RNA": {
        "k_list": "27,47,67,87,107,127,147,167",
        "min_count": 2,
        "prune_level": 2,
        "min_ram_B": 8 * 1024 ** 3, # 8GB
    },
    "WGS": {
        "k_list": "31,39,51,71,91,111,131,151,171",
        "min_count": 3,
        "prune_level": 2,
        "min_ram_B": 8 * 1024 ** 3, # 8GB
    },
}

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
            "NT": Path(DATA_DIR, "Angiosperms353.FNA"),
        },
        "mega353": {
            "AA": Path(DATA_DIR, "Mega353_centroids80.FAA"),
            "NT": Path(DATA_DIR, "Mega353_centroids80.FNA"),
        },
    },
    "PTD": {
        "seedplantsptd": {
            "AA": Path(DATA_DIR, "SeedPlantsPTD.FAA"),
            "NT": Path(DATA_DIR, "SeedPlantsPTD.FNA"),
        },
    },
    "MIT": {
        "seedplantsmit": {
            "AA": Path(DATA_DIR, "SeedPlantsMIT.FAA"),
            "NT": Path(DATA_DIR, "SeedPlantsMIT.FNA"),
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

# JSON with paths to references filename
JSON_REFS = "captus-assembly_extract.refs.json"

# Valid combinations of marker directories and format directories
VALID_MARKER_FORMAT_COMBO =  [(m, f) for m in ["NUC","PTD","MIT"] for f in ["AA","NT","GE","GF"]]
VALID_MARKER_FORMAT_COMBO += [(m, f) for m in ["DNA","CLR"] for f in ["MA","MF"]]

# Scipio's initial round will run with 'min_score' multiplied by this factor
# Empirically, we observed that a score of 0.14 in the first round gave good results, so it is no
# longer necessary to use a lower socre for the first round of Scipio
SCIPIO_SCORE_INITIAL_FACTOR = 1.0

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
    "NUC": {"id":  1, "name": "Standard"},
    "PTD": {"id": 11, "name": "Bacterial, Archaeal and Plant Plastid"},
    "MIT": {"id":  1, "name": "Standard"}
}

# Divergent genetic codes that will need a lower 'SCIPIO_BLAT_IDENTITY_FACTOR' for BLAT to find hits
# these genetic codes have aminoacid substitutions that may make it harder for BLAT to find a hit
DIVERGENT_GENETIC_CODES = [2, 3, 5, 9, 12, 13, 14, 21, 24, 26, 33]

# Scipio's maximum identity when using a divergent genetic code
SCIPIO_MAX_IDENTITY_DIV_CODE = 66

# Basic Scipio setting that are genome-specific, more importantly make BLAT parameters more
# restrictive in total assembly size and maximum intron sizes
SCIPIO_GENOME_BASIC_SETTINGS = {
    # Change here the settings for nuclear genes:
    "NUC": [
    ],
    # Change here the settings for plastidial genes:
    "PTD": [
        "--region_size=0",
        "--blat_params=-maxIntron=2000",
        "--max_assemble_size=9000",
    ],
    # Change here the settings for mitochondrial genes:
    "MIT": [
        "--region_size=0",
        "--blat_params=-maxIntron=9000",
        "--max_assemble_size=50000",
    ],
}

# Extra settings for the final round of Scipio according to the genome of the genes
# Always keep 'gap_to_close' <= 21 or reconstruction of genes across several contigs breaks down
SCIPIO_GENOME_EXTRA_SETTINGS = {
    # Change here the settings for nuclear genes:
    "NUC": [
        "--blat_params=-oneOff=1",
        "--blat_tilesize=6",
        "--exhaust_align_size=5000",
        "--exhaust_gap_size=500",
        "--max_move_exon=10",
    ],
    # Change here the settings for plastidial genes:
    "PTD": [
        "--region_size=0",
        "--blat_params=-oneOff=1 -maxIntron=2000",
        "--blat_tilesize=6",
        "--exhaust_align_size=2000",
        "--exhaust_gap_size=900",
        "--max_assemble_size=9000",
        "--min_intron_len=500",
        "--max_move_exon=10",
        "--gap_to_close=120",
    ],
    # Change here the settings for mitochondrial genes:
    "MIT": [
        "--region_size=0",
        "--blat_params=-oneOff=1 -maxIntron=9000",
        "--blat_tilesize=6",
        "--exhaust_align_size=9000",
        "--exhaust_gap_size=900",
        "--max_assemble_size=50000",
        "--max_move_exon=10",
        "--gap_to_close=120",
    ],
}

# Maximum insertion allowed expressed as proportion of the reference sequence length
DNA_MAX_INSERT_PROP = 0.5

# Absolute maximum insertion allowed in bp
DNA_MAX_INSERT_SIZE = 1000

# Tolerance proportion for determining if two hits are compatible in their percentage of 'identity'
# to the reference, or determining if they have an acceptable margin of overlap. These two
# conditions are used to determine compatible pairs of hits during the assembly of partial hits
# during the extraction of non-coding markers
DNA_TOLERANCE_PROP = 0.05

# Minimum coverage for non-coding hits before attempting greedy assembly
DNA_MIN_COVERAGE_BEFORE_ASSEMBLY = 10

# When reads are assembled within Captus using MEGAHIT, the maximum overlap between adjacent hits is
# the size of the largest kmer used for the assembly. If the assembly comes from a different source
# then the maximum overlap is fixed to this value
DNA_MAX_OVERLAP_BP = 99

# During non-coding sequence extraction, a version of the sequence matched with extra upstream and
# downstream nucleotides is extracted in case they help in alignment or provide useful extra data
# This value determines how many extra nucleotides to extract up- and/or downstream
DNA_UP_DOWN_STREAM_BP = 1000

# Maximum number of FASTA files to rehead and write simultaneously. Too many may hurt performance
# if HDD, it may improve with SSDs
MAX_WRITING_INSTANCES = 16

# Clustering identity percentage, if '--cl_min_identity' is left as 'auto' it becomes 98% of the
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

# A cluster must have at least this proportion of the total number of samples to be used as
# reference
CLR_MIN_SAMPLE_PROP = 0.3

# GFF feature colors, Okabe & Ito palette
GFF_COLORS = {
    "NUC": "#e69f00",  # orange
    "PTD": "#009e73",  # green
    "MIT": "#56B4e9",  # skyblue
    "DNA": "#d55e00",  # vermillion
    "CLR": "#cc79a7",  # purple
}

# Alignment output folders
ALN_DIRS = {
    "UNAL": "01_unaligned",
    "ALND": "02_aligned_untrimmed",
    "TRIM": "03_aligned_trimmed",
    "UNFI": "01_unfiltered",
    "FAST": "02_fast",
    "CARE": "03_careful",
    "NREF": "04_unfiltered_no_refs",
    "NRFA": "05_fast_no_refs",
    "NRCA": "06_careful_no_refs",
}

# MAFFT algorithms choices to real MAFFT syntax translation
MAFFT_ALGORITHMS = {
    "auto": {
        "arg": "--auto",
        "aka": "(a.k.a. Automatic)"
    },
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

# Minimum of sequences allowed in an alignment
MIN_SAMPLES_ALN = 3
