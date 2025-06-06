#!/usr/bin/env python3
"""
Copyright 2020-2025 Edgardo M. Ortiz (e.ortiz.v@gmail.com)
https://github.com/edgardomortiz/Captus

This module contains hard-coded settings for Captus

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

# Seed for randomization
RANDOM_SEED = 142857

# The default recursion limit in Python is 1000, but the complexity of the Needleman-Wunsch
# and Smith-Waterman alignment algorithms is len(seq1) * len(seq2), considering that we used these
# ony to align segments of unmatched proteins, a protein fragment length of 2000 is much more than
# enough
RECURSION_LIMIT = 2000 * 2000

# FASTA valid filename extensions:
FASTA_VALID_EXTENSIONS = [".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz"]

# FASTQ valid filename extensions:
FASTQ_VALID_EXTENSIONS = [".fq", ".fastq", ".fq.gz", ".fastq.gz"]

# GFF valid filename extensions:
GFF_VALID_EXTENSIONS = [".gff", ".gff3", ".gtf", ".gff.gz", ".gff3.gz", ".gtf.gz"]

# Fraction of total RAM available to Captus when using 'auto' in --ram
RAM_FRACTION = 0.99

# Fraction of total RAM available to Captus when running Java applications
RAM_FRACTION_JAVA = 0.70

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

# Maximum number of concurrent reading processes from the HDD. For example, if too many FastQC or
# Falco instances are launched performance suffers due to hard drive limitations (might improve
# in the future with SSDs)
MAX_HDD_READ_INSTANCES = 4

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

# Caculate the min, max, average read length of a FASTQ with this many reads
NUM_READS_TO_CALCULATE_STATS = 100000

# Maximum difference between percetange of A and T or C and G in the last base of a read, if the
# difference is larger thatn this value the last base of the reads is trimmed (Illumina sequencing
# artifact)
MAX_DELTA_AT_GC = 10

# Include largest kmer sizes in MEGAHIT's 'k_list' if they are at most this larger than the mean
# read length of the input FASTQs
DELTA_MEAN_READ_LENGTH_TO_MAX_KMER_SIZE = 31

# Minimum RAM in bytes for a MEGAHIT assembly
MEGAHIT_MIN_RAM_B = 4 * 1024**3

# Presets for MEGAHIT assemblies of RNAseq and WGS
MEGAHIT_PRESETS = {
    "CAPSKIM": {
        "k_list": "31,39,47,63,79,95,111,127,143,159,175",
        "min_count": 2,
        "prune_level": 2,
        "min_ram_B": 6 * 1024**3,  # 4GB
        "min_threads": 8,
        "extra_options": None,
    },
    "RNA": {
        "k_list": "27,47,67,87,107,127,147,167",
        "min_count": 2,
        "prune_level": 2,
        "min_ram_B": 8 * 1024**3,  # 8GB
        "min_threads": 8,
        "extra_options": None,
    },
    "WGS": {
        "k_list": "31,39,49,69,89,109,129,149,169",
        "min_count": 3,
        "prune_level": 2,
        "min_ram_B": 12 * 1024**3,  # 12GB
        "min_threads": 8,
        "extra_options": "--no-mercy",
    },
}

# Minimum number of threads for a MEGAHIT assembly
MEGAHIT_MIN_THREADS = 4

# Minimum kmer size to attempt FASTA to FASTG conversion, usually kmer sizes in the 20s produce too
# complex graphs that are hard to visualize
MIN_KMER_SIZE_FOR_FASTG = 31

# Name of depht of coverage file for MEGAHIT and Salmon
CONTIGS_DEPTH = "contigs_depth.tsv"

# Salmon output directories names:
SALMON_INDEX_DIR = "00_salmon_index"
SALMON_QUANT_DIR = "01_salmon_quant"

# BLAT RAM multiplier for protein searches
BLAT_PROT_RAM_FACTOR = 4

# BLAT RAM multiplier for nuecleotide searches
BLAT_DNA_RAM_FACTOR = 2

# BLAT/Scipio are now parallelized, is good to split the load across threads
EXTRACT_MIN_THREADS = 4

# Hidden option: Switch to True to fill protein gaps with X
FILL_GAP_WITH_X = False

# Minimum stretch of missing aminoacids in recovered protein to be filled by "X"
SCIPIO_MIN_GAP_LEN_TO_X = 5

# Maximum gap size as proportion of reference length. For example, 0.5 means that only gaps that are
# at most as long as half the length of the reference will be translated and aligned
SCIPIO_MAX_GAP_AS_REF_PROP = 0.3

# Use a more relaxed identity between the translated gap and the unmatched segment in the reference
# by subtracting this number from the 'min_identity' chosen for extraction.
# minimum gap identity = min identity - SCIPIO_MAX_GAP_DELTA_IDENTITY
SCIPIO_MAX_GAP_DELTA_IDENTITY = 0.3

# When aligning the translations from the three reading frames to a protein segment, penalize the
# match rate by this number as many times as stop codons are found in the translation corresponding
# to each reading frame
SCIPIO_STOP_PENALTY = 0.5

# Multiply 'wscore' by this much, for each frameshift in the protein
# e.g. for 3 frameshifts = wscore * 0.975 * 0.975 * 0.975
SCIPIO_FRAMESHIFT_PENALTY = 0.975

# Multiply 'wscore' by this much, for each additional contig used in the assembly
# e.g. for 3 contigs = wscore * 0.995 * 0.995
EXTRA_CONTIG_PENALTY = 0.999

# Maximum number of mismatches between recovered protein as given by Scipio and the new translation
# performed by Captus after checking and fixing the gene model
SCIPIO_MAX_MISMATCHES = 3

# Separator in between contigs from the same Scipio's hit
SCIPIO_CONTIG_SEP = "n" * 50

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
REF_CLUSTER_SEP = "-"

# Bundled sets of protein references
PROT_REF_TARGETS = {
    "NUC": {
        "angiosperms353": {
            "AA": Path(DATA_DIR, "Angiosperms353.FAA"),
            "NT": Path(DATA_DIR, "Angiosperms353.FNA"),
            "transtable": 1,
            "notes": "Clustered at 0.98 identity and 0.1 coverage with mmseqs2, isolated Xs removed",
        },
        "mega353": {
            "AA": Path(DATA_DIR, "Mega353.FAA"),
            "NT": Path(DATA_DIR, "Mega353.FNA"),
            "transtable": 1,
            "notes": "Clustered at 0.78 identity and 0.1 coverage with mmseqs2, isolated Xs removed",
        },
    },
    "PTD": {
        "seedplantsptd": {
            "AA": Path(DATA_DIR, "SeedPlantsPTD.FAA"),
            "NT": Path(DATA_DIR, "SeedPlantsPTD.FNA"),
            "transtable": 11,
        },
    },
    "MIT": {
        "seedplantsmit": {
            "AA": Path(DATA_DIR, "SeedPlantsMIT.FAA"),
            "NT": Path(DATA_DIR, "SeedPlantsMIT.FNA"),
            "transtable": 1,
        },
    },
}

# Bundled miscellaneous DNA references
DNA_REF_TARGETS = {
    "angiosperms353": Path(DATA_DIR, "Angiosperms353.FNA"),
    "mega353": Path(DATA_DIR, "Mega353.FNA"),
    "seedplantsptd": Path(DATA_DIR, "SeedPlantsPTD.FNA"),
    "seedplantsmit": Path(DATA_DIR, "SeedPlantsMIT.FNA"),
}

# Bundled dependencies directory
DEPS_DIR = Path(Path(__file__).resolve().parent.parent, "dependencies")

# Bundled Scipio v1.4.1 path
BUNDLED_SCIPIO = Path(DEPS_DIR, "scipio-1.4", "scipio.1.4.1.pl")

# Bundled BLAT >= v36x7 path
os_type = platform.system()
if os_type == "Darwin":
    BUNDLED_BLAT = Path(DEPS_DIR, "blat", "blat_37x1_mac")
elif os_type == "Linux":
    BUNDLED_BLAT = Path(DEPS_DIR, "blat", "blat_37x1_linux")
else:
    BUNDLED_BLAT = f"No support for {os_type} systems at this time."

# Ddirectory for storing imported reference target files
REF_TARGETS_DIR = "00_ref_targets"

# Temporary directory for splitting reference target files
REF_TARGETS_SPLIT_DIR = "01_ref_targets_split"

# Create groups of roughly this amount of sequences when splitting the reference file
REF_SPLIT_CHUNK_SIZE = 1000

# Extraction directories per marker type
MARKER_DIRS = {
    "NUC": "01_coding_NUC",  # Nuclear protein markers
    "PTD": "02_coding_PTD",  # Plastidial protein markers
    "MIT": "03_coding_MIT",  # Mitochondrial protein markers
    "DNA": "04_misc_DNA",  # Miscellaneous DNA markers
    "CLR": "05_clusters",  # Cluster-derived DNA markers
}

# Extraction directories per marker type
DIR_MARKERS = {
    "01_coding_NUC": "NUC",  # Nuclear protein markers
    "02_coding_PTD": "PTD",  # Plastidial protein markers
    "03_coding_MIT": "MIT",  # Mitochondrial protein markers
    "04_misc_DNA": "DNA",  # Miscellaneous DNA markers
    "05_clusters": "CLR",  # Cluster-derived DNA markers
}

# Filename suffixes for the different extraction FASTA formats
FORMAT_SUFFIXES = {
    "AA": "_coding_AA.faa",  # coding sequences in aminoacids
    "NT": "_coding_NT.fna",  # coding sequences in nucleotides
    "GE": "_genes.fna",  # complete gene sequences in nucleotides
    "GF": "_genes_flanked.fna",  # complete gene sequences with additional flanking bp
    "MA": "_matches.fna",  # matching DNA sequences from miscellaneous DNA extraction
    "MF": "_matches_flanked.fna",  # matching DNA sequences with additional flanking bp
}

# Directory names for the different extraction FASTA formats
FORMAT_DIRS = {
    "AA": "01_AA",  # coding sequences in aminoacids
    "NT": "02_NT",  # coding sequences in nucleotides
    "GE": "03_genes",  # complete gene sequences in nucleotides
    "GF": "04_genes_flanked",  # complete gene sequences with additional flanking bp
    "MA": "01_matches",  # matching DNA sequences from miscellaneous DNA extraction
    "MF": "02_matches_flanked",  # matching DNA sequences with additional flanking bp
}

# Translated protein reference suffix
TRANSLATED_REF_SUFFIX = ".captus.faa"

# JSON with paths to references filename
JSON_REFS = "captus-extract_refs.json"

# Valid combinations of marker directories and format directories
VALID_MARKER_FORMAT_COMBO = [(m, f) for m in ["NUC", "PTD", "MIT"] for f in ["AA", "NT", "GE", "GF"]]
VALID_MARKER_FORMAT_COMBO += [(m, f) for m in ["DNA", "CLR"] for f in ["MA", "MF"]]

# Scipio's initial round will run with 'min_score' multiplied by this factor
# Empirically, we observed that a score of 0.14 in the first round gave good results, so it is no
# longer necessary to use a lower socre for the first round of Scipio
SCIPIO_SCORE_INITIAL_FACTOR = 1.0

# Scipio's accepted intron penalty controls what types of intron borders are accepted:
# 1.0 = GT--AG, GC--AG
# 1.1 = GT--AG, GC--AG, AT--AC
# 1.2 = GT--AG, GC--AG, AT--AC, GG--AG, GA--AG
SCIPIO_ACCEPTED_INTRON_PENALTY = 1.1

# Default Genetic Codes to set Scipio's --transtable
DEFAULT_GENETIC_CODES = {
    "NUC": {"id": 1, "name": "Standard"},
    "PTD": {"id": 11, "name": "Bacterial, Archaeal and Plant Plastid"},
    "MIT": {"id": 1, "name": "Standard"},
}

# Divergent genetic codes that will need a lower 'SCIPIO_BLAT_IDENTITY_FACTOR' for BLAT to find hits
# these genetic codes have aminoacid substitutions that may make it harder for BLAT to find a hit
DIVERGENT_GENETIC_CODES = [2, 3, 5, 9, 12, 13, 14, 21, 24, 26, 33]

# Scipio's maximum identity when using a divergent genetic code
SCIPIO_MAX_IDENT_DIV_CODE = 66

# Basic Scipio setting that are genome-specific, more importantly make BLAT parameters more
# restrictive in total assembly size and maximum intron sizes
SCIPIO_ROUND1_SETTINGS = {
    # Change here the initial settings for nuclear genes:
    "NUC": {
        "scipio": [
            "--blat_params=-maxIntron=50000",
            "--min_dna_coverage=0.2",
        ],
        # They keys and values of BLAT settings must match the string in Scipio's '--blat_params'
        "blat": {
            "tile_size": None,
            "one_off": None,
            "min_score": 15,
            "max_intron": 50000,
        },
    },
    # Change here the initial settings for plastidial genes:
    "PTD": {
        "scipio": [
            "--region_size=0",
            "--blat_params=-maxIntron=1300",
            "--max_assemble_size=9000",
        ],
        # They keys and values of BLAT settings must match the string in Scipio's '--blat_params'
        "blat": {
            "tile_size": None,
            "one_off": None,
            "min_score": 15,
            "max_intron": 1300,
        },
    },
    # Change here the initial settings for mitochondrial genes:
    "MIT": {
        "scipio": [
            "--region_size=0",
            "--blat_params=-maxIntron=9000",
            "--max_assemble_size=50000",
        ],
        # They keys and values of BLAT settings must match the string in Scipio's '--blat_params'
        "blat": {
            "tile_size": None,
            "one_off": None,
            "min_score": 15,
            "max_intron": 9000,
        },
    },
}

# Extra settings for the final round of Scipio according to the genome of the genes
# Always keep 'gap_to_close' <= 21 or reconstruction of genes across several contigs breaks down
SCIPIO_ROUND2_SETTINGS = {
    # Change here the final settings for nuclear genes:
    "NUC": {
        "scipio": [
            "--blat_params=-maxIntron=50000",
            "--blat_tilesize=6",
            "--exhaust_align_size=30000",
            "--exhaust_gap_size=21",
            "--min_dna_coverage=0.2",
        ],
        # They keys and values of BLAT settings must match the string in Scipio's '--blat_params'
        "blat": {
            "tile_size": 6,
            "one_off": None,
            "min_score": 15,
            "max_intron": 50000,
        },
    },
    # Change here the final settings for plastidial genes:
    "PTD": {
        "scipio": [
            "--region_size=0",
            "--blat_params=-oneOff=1 -maxIntron=1300",
            "--exhaust_align_size=1300",
            "--max_assemble_size=9000",
            "--min_intron_len=500",
            "--gap_to_close=84",  # Don't set >90, genes found in 1st round and lost in 2nd
        ],
        # They keys and values of BLAT settings must match the string in Scipio's '--blat_params'
        "blat": {
            "tile_size": None,
            "one_off": 1,
            "min_score": 15,
            "max_intron": 1300,
        },
    },
    # Change here the final settings for mitochondrial genes:
    "MIT": {
        "scipio": [
            "--region_size=0",
            "--blat_params=-oneOff=1 -maxIntron=9000",
            "--exhaust_align_size=9000",
            "--max_assemble_size=50000",
            "--gap_to_close=84",  # Don't set >90, genes found in 1st round and lost in 2nd
        ],
        # They keys and values of BLAT settings must match the string in Scipio's '--blat_params'
        "blat": {
            "tile_size": None,
            "one_off": 1,
            "min_score": 15,
            "max_intron": 9000,
        },
    },
}

# BLAT identity as proportion of Scipio's 'min_identity' parameter
PROT_BLAT_IDENT_PROP = 0.90

# BLAT identity proportion for DNA to DNA searches
SEARCH_IDENT_PROP = 0.99

# BLAT identity proportion for DNA to DNA searches
DNA_BLAT_MIN_SCORE = 20

# Maximum gap between adjacent partial BLAT hits (-maxIntron), to limit the hit assembly size
# when searching in chromosome level assemblies (should be increased for synteny analysis?)
DNA_BLAT_MAX_INSERTION = 50000

# Maximum organellar genome size, avoid looking for organellar genes in larger contigs
MAX_ORGANELLE_SIZE = {
    "PTD": 550_000,
    "MIT": 1_000_000,
}

# Minimum wscore for a PSL records
MIN_ABS_WSCORE = 0.001

# Minimum wscore as proportio of the best wscore for the target sequence in PSL records
MIN_PROP_WSCORE = 0.02

# BLAT size multipliers per marker type
SIZE_MUL = {
    "NUC": 3,
    "PTD": 3,
    "MIT": 3,
    "DNA": 1,
    "CLR": 1,
}

# Genes that partially overlap
VALID_OVERLAPS = [
    ["atpB", "atpE"],
    ["ndhC", "ndhK"],
    ["ndhF", "ycf1"],
    ["psbC", "psbD"],
]

# Maximum percetange of overlap between two hits, nuclear genes should not overlap, but we allow
# a small percentage of overlap
HIT_MAX_PCT_OVERLAP = 5.0

# During prefiltering of BLAT hits, keep this top N best targets per locus
BEST_N_TARGETS_GLOBAL = 100

# For extractions with two rounds, keep this many targets after the first round (before was 1)
BEST_N_TARGETS_INITIAL = 7

# Minimum contig size to be considered a full chromosome in bp
MIN_CHROM_SIZE = 10_000_000

# Maximum number of cross-loci comparisons to make when evaluating overlaps
MAX_CROSS_LOCI_COMP = 400_000_000

# Minimum length of terminal exon for organellar proteins
SCIPIO_MIN_LEN_FINAL_EXON = 21

# Maximum insertion allowed expressed as proportion of the reference sequence length
DNA_MAX_INSERT_PROP = 0.2

# Absolute maximum insertion allowed in bp
DNA_MAX_INSERT_SIZE = 200

# Tolerance proportion for determining if two hits are compatible in their margin of overlap. Used
# during the assembly of partial hits during the extraction of non-coding markers
DNA_TOLERANCE_LEN = 0.0500

# Tolerance proportion for determining if two hits are compatible in their percentage of 'identity'
# to the reference. Used during the assembly of partial hits during the extraction of non-coding
# markers (3.33% seems to perform OK across several conditions)
DNA_TOLERANCE_PID = 0.0333

# Minimum coverage for non-coding hits before attempting greedy assembly
DNA_MIN_COVERAGE_BEFORE_ASSEMBLY = 10

# When reads are assembled within Captus using MEGAHIT, the maximum overlap between adjacent hits is
# the size of the largest kmer used for the assembly. If the assembly comes from a different source
# then the maximum overlap is fixed to this value
DNA_MAX_OVERLAP_BP = 99

# Smallest reliable overlap for two contigs to be concatenated during assembly of flanked hits
DNA_MIN_OVERLAP_BP = 21

# Separator for gaps of undetermined lengths during assemble of flanked hits
DNA_CONTIG_SEPARATOR = "n" * 10

# During non-coding sequence extraction, a version of the sequence matched with extra upstream and
# downstream nucleotides is extracted in case they help in alignment or provide useful extra data
# This value determines how many extra nucleotides to extract up- and/or downstream
DNA_UP_DOWN_STREAM_BP = 1000

# Maximum number of FASTA files to rehead and write simultaneously. Too many may hurt performance
# if HDD, it may improve with SSDs
MAX_HDD_WRITE_INSTANCES = 24

# Clustering identity percentage, if '--cl_min_identity' is left as 'auto' it becomes 99% of the
# '--dna_min_identity' value
MMSEQS_BLAT_DNA_IDENTITY_FACTOR = 0.99

# Minimum clustering identity when set to 'auto'
MMSEQS_MIN_AUTO_MIN_IDENTITY = 75

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
MMSEQS_COV_MODE = 1

# Penalty for opening a gap when aligning sequences during clustering. The lower the value the
# slower clustering becomes. Not recommended to go lower than 3, minimum possible value is 1
MMSEQS_GAP_OPEN = 3

# Penalty for extending a gap when aligning sequences during clustering. The lower the value the
# slower clustering becomes. Minimum possible value is 1"
MMSEQS_GAP_EXTEND = 1

# From MMseqs2 help:
# Scale k-mer per sequence based on sequence length as kmer-per-seq val + scale x seqlen
# kmer-per-seq is 21 by default, so if we set this value to 0.3 for a 500 bp sequence, MMseqs2 will
# sample (500*0.3)+21 = 171 kmers from that sequence
MMSEQS_KMER_PER_SEQ_SCALE = 0.3

# MMseqs sensitivity set at the max 7.5 (If too slow for a dataset, this can be decreased via options)
MMSEQS_SENSITIVITY = 7.5

# Memory proportion to allocate to MMseqs
MMSEQS_RAM_FRACTION = 0.99

# MMask low complexity regions prior to clustering 1 = True, 0 = False
MMSEQS_MASK_LOW_COMPLEXITY = 1

# A cluster must have at least this proportion of the total number of samples to be used as
# reference
CLR_MIN_SAMPLE_PROP = 0.66

# Minimum proportion of longest cluster representative to be retained for secondary cluster
# representatives
CLR_REP_MIN_LEN_PROP = 0.8

# Minimum proportion of total samples allowed in a secondary cluster representative
CLR_REP_MIN_SAMPLE_PROP = 0.01

# Extraction statistics table header
EXT_STATS_HEADER = [
    "sample_name",
    "marker_type",
    "locus",
    "ref_name",
    "ref_coords",
    "ref_type",
    "ref_len_matched",
    "hit",
    "pct_recovered",
    "pct_identity",
    "score",
    "wscore",
    "hit_len",
    "cds_len",
    "intron_len",
    "flanks_len",
    "frameshifts",
    "hit_contigs",
    "hit_l50",
    "hit_l90",
    "hit_lg50",
    "hit_lg90",
    "ctg_names",
    "ctg_strands",
    "ctg_coords",
]

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
    "ALND": "02_untrimmed",
    "TRIM": "03_trimmed",
    "UNFR": "01_unfiltered_w_refs",
    "NAIR": "02_naive_w_refs",
    "INFR": "03_informed_w_refs",
    "UNFI": "04_unfiltered",
    "NAIV": "05_naive",
    "INFO": "06_informed",
}

# MAFFT/MUSCLE algorithms choices to real MAFFT/MUSCLE syntax translation
ALIGN_ALGORITHMS = {
    "mafft_auto": {"arg": "--auto", "aka": "(a.k.a. MAFFT's Automatic)"},
    "mafft_genafpair": {"arg": "--genafpair", "aka": "(a.k.a. MAFFT's E-INS-i)"},
    "mafft_localpair": {"arg": "--localpair", "aka": "(a.k.a. MAFFT's L-INS-i)"},
    "mafft_globalpair": {"arg": "--globalpair", "aka": "(a.k.a. MAFFT's G-INS-i)"},
    "mafft_retree1": {"arg": "--retree 1", "aka": "(a.k.a. MAFFT's FFT-NS-1)"},
    "mafft_retree2": {"arg": "--retree 2", "aka": "(a.k.a. MAFFT's FFT-NS-2)"},
    "muscle_align": {"arg": "-align", "aka": "(a.k.a. MUSCLE 5's PPP)"},
    "muscle_super5": {"arg": "-super5", "aka": "(a.k.a. MUSCLE 5's Super5)"},
}

# TAPER default parameters
TAPER_PARAMS = [
    {"k": 5, "p": 0.25, "q": 0.1, "L": 30},
    {"k": 9, "p": 0.25, "q": 0.25, "L": 54},
    {"k": 17, "p": 0.1, "q": 0.5, "L": float("inf")},
]

# TAPER default mask character
TAPER_MASK = "-"

# Default MAFFT algorithm for Captus Design, use keys from ALIGN_ALGORITHMS (above)
DESIGN_ALIGN_ALGORITHM = "mafft_genafpair"

# Separator used in output sequences names, for example used to distinguish sample name from gene
# name or to indicate copy number
SEQ_NAME_SEP = "__"

# File name for sequence-to-sample equivalence table used by ASTRAL-Pro to analyze trees that
# include paralogs
ASTRAL_PRO_EQ = "captus-align_astral-pro.tsv"

# Import data for clustering file names
DES_SUFFIXES = {
    "CDS": "_markers.fasta",  # Full CDS, concatenated exons
    "LONG": "_exons_long.fasta",  # Individual exons >= bait_length
    "SHORT": "_exons_short.fasta",  # Individual exons < bait_length
    "DATA": "_exons_data.tsv",  # Exon and intron data, tab-separated values table
    "MARKERS": "_markers.fasta",  # Any sequences, when no GFF is found
    "FILTER": "_markers_filtered.fasta",  # Sequences filtered according to '--min_seq_len'
    "DEDUPED": "_markers_deduped.fasta",  # Deduplicated sequences, ready for clustering
}

# Design output folders
DES_DIRS = {
    "CAT": "01_concatenated",  # inside 01_clustered_markers
    "CLR": "02_clustered",  # inside 01_clustered_markers
    "ALN": "03_aligned",  # inside 01_clustered_markers
    "CUR": "04_curated",  # inside 01_clustered_markers
    "AUT": "01_automatic",  # inside 02_selected_markers
    "MAN": "02_manual",  # inside 02_selected_markers
    "RAW": "01_raw",  # inside 03_baits
    "FIL": "02_filtered",  # inside 03_baits
    "BAI": "03_clustered",  # inside 03_baits
    "TAR": "04_baits_and_targets",  # inside 03_baits
}

# Design output files
DES_FILES = {
    "BFEX": "baits_full_exons.fasta",
    "BFNE": "baits_full_no_exons.fasta",
    "BDEX": "baits_exons.fasta",
    "BDNE": "baits_no_exons.fasta",
    "LONG": "long_exons.fasta",
    "BEXM": "baits_exons_mapped.fasta",
    "BCAT": "baits_concat.fasta",
    "BACC": "baits_accepted.fasta",
    "BREJ": "baits_rejected.fasta",
}

# Automatic maximum sequence length for deduplication and clustering according to clustering program
MAX_SEQ_LEN = {
    "mmseqs": 65535,
    "vsearch": 5000,
}

# Split bait file if chunks of this size
BAITS_SPLIT_SIZE = 500000

# Binary chunk size
CHUNK_SIZE = 8192

# Minimum proportion of the longest sequence to be considered as long in Captus design
MIN_SEQ_LEN_PROP = 0.5

# Max indel allowed when mapping to long exons as proportion of bait length
MAX_INDEL_BAIT_PROP = 0.1
