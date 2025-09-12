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

import copy
import gzip
import math
import random
import re
import shutil
import statistics
import subprocess
import sys
import tarfile
import time
import urllib
from pathlib import Path

from . import settings
from .misc import ElapsedTimeThread, bold, elapsed_time, gzip_compress, pigz_compress

# Regular expression to match MEGAHIT headers assembled within Captus
CAPTUS_MEGAHIT_HEADER_REGEX = r"^NODE_\d+_length_\d+_cov_\d+\.\d+_k_\d{2,3}_flag_\d"

# Regular expression to match MEGAHIT headers
MEGAHIT_HEADER_REGEX = r"^k\d{2,3}_\d+\sflag=\d\smulti=\d+\.\d+\slen=\d+$"

# Regular expression to match SKESA headers
SKESA_HEADER_REGEX = r"^Contig_\d+_\d+\.\d+$|^Contig_\d+_\d+\.\d+_Circ$"

# Regular expression to match Spades headers
SPADES_HEADER_REGEX = r"^NODE_\d+_length_\d+_cov_\d+\.\d+|^EDGE_\d+_length_\d+_cov_\d+\.\d+"

# Set of valid nucleotides, including IUPAC ambiguities
NT_IUPAC = {
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "T": ["T"],
    "M": ["A", "C"],
    "R": ["G", "A"],
    "W": ["A", "T"],
    "S": ["G", "C"],
    "Y": ["C", "T"],
    "K": ["G", "T"],
    "V": ["A", "G", "C"],
    "H": ["A", "C", "T"],
    "D": ["G", "A", "T"],
    "B": ["G", "C", "T"],
    "N": ["A", "C", "G", "T"],
    "-": ["-"],
}

# Calculate pairwise identities for nucleotides
NT_PIDS = {}
i = 0
for a in sorted(NT_IUPAC):
    for b in sorted(NT_IUPAC)[i:]:
        expand = NT_IUPAC[a] + NT_IUPAC[b]
        combos = len(NT_IUPAC[a]) * len(NT_IUPAC[b])
        matches = 0
        for nuc in set(expand):
            matches += (expand.count(nuc) * (expand.count(nuc) - 1)) / 2
        NT_PIDS[f"{a}{b}"] = matches / combos
    i += 1
NT_PIDS["--"] = 0.0

# fmt: off
# # Set of valid aminoacids, including IUPAC ambiguities
AA_IUPAC = {
    "A": ["A"],
    "B": ["D", "N"],
    "C": ["C"],
    "D": ["D"],
    "E": ["E"],
    "F": ["F"],
    "G": ["G"],
    "H": ["H"],
    "I": ["I"],
    "J": ["I", "L"],
    "K": ["K"],
    "L": ["L"],
    "M": ["M"],
    "N": ["N"],
    "O": ["O"],
    "P": ["P"],
    "Q": ["Q"],
    "R": ["R"],
    "S": ["S"],
    "T": ["T"],
    "U": ["U"],
    "V": ["V"],
    "W": ["W"],
    "Y": ["Y"],
    "Z": ["E", "Q"],
    "X": ["A","C","D","E","F","G","H","I","K","L","M","N","O","P","Q","R","S","T","U","V","W","Y"],
    "*": ["*"],
    "-": ["-"],
}
# fmt: on

# Calculate pairwise identities for aminoacids
AA_PIDS = {}
i = 0
for a in sorted(AA_IUPAC):
    for b in sorted(AA_IUPAC)[i:]:
        expand = AA_IUPAC[a] + AA_IUPAC[b]
        combos = len(AA_IUPAC[a]) * len(AA_IUPAC[b])
        matches = 0
        for nuc in set(expand):
            matches += (expand.count(nuc) * (expand.count(nuc) - 1)) / 2
        AA_PIDS[a + b] = matches / combos
    i += 1
AA_PIDS["--"] = 0.0

# Ambiguous aminoacids used in fuzzy translation
AA_AMBIGS = {
    ("D", "N"): "N",  # "B"
    ("E", "Q"): "Q",  # "Z"
    ("I", "L"): "L",  # "J"
}

# Since the list of nucleotides is a subset of the list of aminoacids we need a list of only the
# aminoacid letters that are not also nucleotide letters
AA_NOT_IN_NT = set(AA_IUPAC) - set(NT_IUPAC)

# Reverse complement dictionary
REV_COMP_DICT = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "g",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "r": "y",
    "y": "r",
    "s": "s",
    "w": "w",
    "K": "M",
    "M": "K",
    "B": "V",
    "V": "B",
    "k": "m",
    "m": "k",
    "b": "v",
    "v": "b",
    "D": "H",
    "H": "D",
    "N": "N",
    "d": "h",
    "h": "d",
    "n": "n",
    ".": ".",
    "-": "-",
    "?": "?",
}

# Codons corresponding to aminoacid strings in GENETIC_CODES
CODONS = {
    "base1": "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
    "base2": "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
    "base3": "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG",
}

# fmt: off
# NCBI Genetic Codes (modified from ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt)
GENETIC_CODES = {
    1: {
        "name": "Standard",
        "aa": "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "ss": "---M------**--*----M---------------M----------------------------",
    },
    2: {
        "name": "Vertebrate Mitochondrial",
        "aa": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        "ss": "----------**--------------------MMMM----------**---M------------",
    },
    3: {
        "name": "Yeast Mitochondrial",
        "aa": "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "ss": "----------**----------------------MM---------------M------------",
    },
    4: {
        "name": "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma",
        "aa": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "ss": "--MM------**-------M------------MMMM---------------M------------",
    },
    5: {
        "name": "Invertebrate Mitochondrial",
        "aa": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        "ss": "---M------**--------------------MMMM---------------M------------",
    },
    6: {
        "name": "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear",
        "aa": "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "ss": "--------------*--------------------M----------------------------",
    },
    9: {
        "name": "Echinoderm Mitochondrial; Flatworm Mitochondrial",
        "aa": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "ss": "----------**-----------------------M---------------M------------",
    },
    10: {
        "name": "Euplotid Nuclear",
        "aa": "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "ss": "----------**-----------------------M----------------------------",
    },
    11: {
        "name": "Bacterial, Archaeal and Plant Plastid",
        "aa": "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "ss": "---M------**--*----M------------MMMM---------------M------------",
    },
    12: {
        "name": "Alternative Yeast Nuclear",
        "aa": "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "ss": "----------**--*----M---------------M----------------------------",
    },
    13: {
        "name": "Ascidian Mitochondrial",
        "aa": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
        "ss": "---M------**----------------------MM---------------M------------",
    },
    14: {
        "name": "Alternative Flatworm Mitochondrial",
        "aa": "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "ss": "-----------*-----------------------M----------------------------",
    },
    15: {
        "name": "Blepharisma Macronuclear",
        "aa": "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "ss": "----------*---*--------------------M----------------------------",
    },
    16: {
        "name": "Chlorophycean Mitochondrial",
        "aa": "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "ss": "----------*---*--------------------M----------------------------",
    },
    21: {
        "name": "Trematode Mitochondrial",
        "aa": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "ss": "----------**-----------------------M---------------M------------",
    },
    22: {
        "name": "Scenedesmus obliquus Mitochondrial",
        "aa": "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "ss": "------*---*---*--------------------M----------------------------",
    },
    23: {
        "name": "Thraustochytrium Mitochondrial",
        "aa": "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "ss": "--*-------**--*-----------------M--M---------------M------------",
    },
    24: {
        "name": "Rhabdopleuridae Mitochondrial",
        "aa": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        "ss": "---M------**-------M---------------M---------------M------------",
    },
    25: {
        "name": "Candidate Division SR1 and Gracilibacteria",
        "aa": "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "ss": "---M------**-----------------------M---------------M------------",
    },
    26: {
        "name": "Pachysolen tannophilus Nuclear",
        "aa": "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "ss": "----------**--*----M---------------M----------------------------",
    },
    27: {
        "name": "Karyorelict Nuclear",
        "aa": "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "ss": "--------------*--------------------M----------------------------",
    },
    28: {
        "name": "Condylostoma Nuclear",
        "aa": "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "ss": "----------**--*--------------------M----------------------------",
    },
    29: {
        "name": "Mesodinium Nuclear",
        "aa": "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "ss": "--------------*--------------------M----------------------------",
    },
    30: {
        "name": "Peritrich Nuclear",
        "aa": "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "ss": "--------------*--------------------M----------------------------",
    },
    31: {
        "name": "Blastocrithidia Nuclear",
        "aa": "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "ss": "----------**-----------------------M----------------------------",
    },
    32: {
        "name": "Balanophoraceae Plastid",
        "aa": "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "ss": "---M------*---*----M------------MMMM---------------M------------",
    },
    33: {
        "name": "Cephalodiscidae Mitochondrial",
        "aa": "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        "ss": "---M-------*-------M---------------M---------------M------------",
    },
}
# fmt: on


def score_matrix_to_dict(matrix):
    """
    Convert a matrix given in tabular format into a dictionary for alignment functions

    Parameters
    ----------
    matrix : list
        Two-dimensional matrix defining alignment scores
    """
    matrix_as_dict = {}
    col_labels = matrix[0]
    for i in range(1, len(matrix)):
        row_label = matrix[i][0]
        for j in range(1, len(matrix[i])):
            key = f"{row_label}{col_labels[j]}"
            matrix_as_dict[key] = matrix[i][j]
    matrix_as_dict["gap"] = min(matrix_as_dict.values())

    return matrix_as_dict


# fmt: off
# PAM250 matrix for scoring protein alignments
PAM250 = score_matrix_to_dict(
    [
        [".","A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","B","J","Z","X","*"],
        ["A", +2, -2, +0, +0, -2, +0, +0, +1, -1, -1, -2, -1, -1, -3, +1, +1, +1, -6, -3, +0, +0, -1, +0, -1, -8],
        ["R", -2, +6, +0, -1, -4, +1, -1, -3, +2, -2, -3, +3, +0, -4, +0, +0, -1, +2, -4, -2, -1, -3, +0, -1, -8],
        ["N", +0, +0, +2, +2, -4, +1, +1, +0, +2, -2, -3, +1, -2, -3, +0, +1, +0, -4, -2, -2, +2, -3, +1, -1, -8],
        ["D", +0, -1, +2, +4, -5, +2, +3, +1, +1, -2, -4, +0, -3, -6, -1, +0, +0, -7, -4, -2, +3, -3, +3, -1, -8],
        ["C", -2, -4, -4, -5,+12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3, +0, -2, -8, +0, -2, -4, -5, -5, -1, -8],
        ["Q", +0, +1, +1, +2, -5, +4, +2, -1, +3, -2, -2, +1, -1, -5, +0, -1, -1, -5, -4, -2, +1, -2, +3, -1, -8],
        ["E", +0, -1, +1, +3, -5, +2, +4, +0, +1, -2, -3, +0, -2, -5, -1, +0, +0, -7, -4, -2, +3, -3, +3, -1, -8],
        ["G", +1, -3, +0, +1, -3, -1, +0, +5, -2, -3, -4, -2, -3, -5, +0, +1, +0, -7, -5, -1, +0, -4, +0, -1, -8],
        ["H", -1, +2, +2, +1, -3, +3, +1, -2, +6, -2, -2, +0, -2, -2, +0, -1, -1, -3, +0, -2, +1, -2, +2, -1, -8],
        ["I", -1, -2, -2, -2, -2, -2, -2, -3, -2, +5, +2, -2, +2, +1, -2, -1, +0, -5, -1, +4, -2, +3, -2, -1, -8],
        ["L", -2, -3, -3, -4, -6, -2, -3, -4, -2, +2, +6, -3, +4, +2, -3, -3, -2, -2, -1, +2, -3, +5, -3, -1, -8],
        ["K", -1, +3, +1, +0, -5, +1, +0, -2, +0, -2, -3, +5, +0, -5, -1, +0, +0, -3, -4, -2, +1, -3, +0, -1, -8],
        ["M", -1, +0, -2, -3, -5, -1, -2, -3, -2, +2, +4, +0, +6, +0, -2, -2, -1, -4, -2, +2, -2, +3, -2, -1, -8],
        ["F", -3, -4, -3, -6, -4, -5, -5, -5, -2, +1, +2, -5, +0, +9, -5, -3, -3, +0, +7, -1, -4, +2, -5, -1, -8],
        ["P", +1, +0, +0, -1, -3, +0, -1, +0, +0, -2, -3, -1, -2, -5, +6, +1, +0, -6, -5, -1, -1, -2, +0, -1, -8],
        ["S", +1, +0, +1, +0, +0, -1, +0, +1, -1, -1, -3, +0, -2, -3, +1, +2, +1, -2, -3, -1, +0, -2, +0, -1, -8],
        ["T", +1, -1, +0, +0, -2, -1, +0, +0, -1, +0, -2, +0, -1, -3, +0, +1, +3, -5, -3, +0, +0, -1, -1, -1, -8],
        ["W", -6, +2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4, +0, -6, -2, -5,+17, +0, -6, -5, -3, -6, -1, -8],
        ["Y", -3, -4, -2, -4, +0, -4, -4, -5, +0, -1, -1, -4, -2, +7, -5, -3, -3, +0,+10, -2, -3, -1, -4, -1, -8],
        ["V", +0, -2, -2, -2, -2, -2, -2, -1, -2, +4, +2, -2, +2, -1, -1, -1, +0, -6, -2, +4, -2, +2, -2, -1, -8],
        ["B", +0, -1, +2, +3, -4, +1, +3, +0, +1, -2, -3, +1, -2, -4, -1, +0, +0, -5, -3, -2, +3, -3, +2, -1, -8],
        ["J", -1, -3, -3, -3, -5, -2, -3, -4, -2, +3, +5, -3, +3, +2, -2, -2, -1, -3, -1, +2, -3, +5, -2, -1, -8],
        ["Z", +0, +0, +1, +3, -5, +3, +3, +0, +2, -2, -3, +0, -2, -5, +0, +0, -1, -6, -4, -2, +2, -2, +3, -1, -8],
        ["X", -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -8],
        ["*", -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, +1],
    ]
)
# fmt: on


def get_read_stats(fastq_path, num_reads):
    """
    Determine mean read length to adjust MEGAHIT's '--k-list' and '--min-contig-len' accordingly
    """
    line_count = 0
    read_lengths = []
    if f"{fastq_path}".endswith(".gz"):
        opener = gzip.open
    else:
        opener = open
    try:
        with opener(fastq_path, "rt") as fastq:
            for line in fastq:
                line_count += 1
                if line_count % 4 == 0:
                    read_lengths.append(len(line.strip("\n")))
                if line_count == num_reads * 4:
                    break
    except gzip.BadGzipFile:
        return False

    read_stats = {
        "min_read_length": min(read_lengths),
        "max_read_length": max(read_lengths),
        "mean_read_length": math.ceil(statistics.mean(read_lengths)),
    }

    return read_stats


def genetic_code(genetic_code_id):
    """
    Given the `genetic_code_id` returns the code-specific list of codon to aminoacid equivalences,
    list of alternative start codons, and list of alternative stop codons

    Parameters
    ----------
    genetic_code_id : int
        Genetic Code's numeric ID, see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

    Returns
    -------
    dict
        A dictionary with three lists: aminoacids, alternative starts, and alternative stops
    """
    aminos = {}
    starts = {}
    stops = {}
    for b1, b2, b3, a, s in zip(
        CODONS["base1"],
        CODONS["base2"],
        CODONS["base3"],
        GENETIC_CODES[genetic_code_id]["aa"],
        GENETIC_CODES[genetic_code_id]["ss"],
    ):
        codon = f"{b1}{b2}{b3}"
        aminos[codon] = a
        if s == "M":
            starts[codon] = s
        if s == "*":
            stops[codon] = s

    return {"aminos": aminos, "starts": starts, "stops": stops}


def translate(seq, genetic_code: dict, frame=1, start_as_M=False):
    """
    Translate nucleotide sequence `seq` given a specific reading `frame` and the dictionary
    produced by the function `genetic_code()`

    Parameters
    ----------
    seq : str
        Nucleotide sequence
    genetic_code : dict
        A dictionary with three lists: aminoacids, alternative starts, and alternative stops
    frame : int, optional
        Reading frame, accepted values are: 1, 2, 3, -1, -2, -3, by default 1

    Returns
    -------
    str
        Aminoacid sequence
    """

    def resolve_codon(codon):
        """
        If `codon` contains IUPAC nucleotide ambiguity codes, solve them and return a list of all
        possible codons to perform fuzzy translation
        """
        if "-" in codon or codon.count("N") > 1:
            return []
        else:
            return [
                f"{p1}{p2}{p3}"
                for p1 in NT_IUPAC[codon[0]]
                for p2 in NT_IUPAC[codon[1]]
                for p3 in NT_IUPAC[codon[2]]
            ]

    def unresolve_aminoacid(aminos: list):
        """
        Because of fuzzy translation, a codon may produce several aminoacids. This function reduces
        the possibilities to a single residue. It considers the ambiguous aminoacids in `AA_AMBIGS`.
        """
        if len(aminos) == 1:
            return aminos[0]
        aminos = tuple(sorted(set(aminos)))
        if len(aminos) == 1:
            return aminos[0]
        elif len(aminos) == 2 and aminos in AA_AMBIGS:
            return AA_AMBIGS[aminos]
        else:
            return "X"

    seq = seq.replace(" ", "").replace("-", "").upper()
    if frame < 0:
        seq = reverse_complement(seq)

    codons = [seq[p : p + 3] for p in range(abs(frame) - 1, len(seq), 3) if len(seq[p : p + 3]) == 3]

    if len(codons) < 1:
        return ""

    seq_AA = ""

    starts_fuzzy = []
    for triplet in resolve_codon(codons[0]):
        if triplet in genetic_code["starts"] and start_as_M:
            starts_fuzzy.append(genetic_code["starts"][triplet])
        elif triplet in genetic_code["aminos"]:
            starts_fuzzy.append(genetic_code["aminos"][triplet])
        else:
            starts_fuzzy.append("X")
    seq_AA += unresolve_aminoacid(starts_fuzzy)

    if len(codons) == 1:
        return seq_AA

    for codon in codons[1:-1]:
        aminos_fuzzy = []
        for triplet in resolve_codon(codon):
            if triplet in genetic_code["aminos"]:
                aminos_fuzzy.append(genetic_code["aminos"][triplet])
            else:
                aminos_fuzzy.append("X")
        seq_AA += unresolve_aminoacid(aminos_fuzzy)

    stops_fuzzy = []
    for triplet in resolve_codon(codons[-1]):
        if triplet in genetic_code["stops"]:
            stops_fuzzy.append(genetic_code["stops"][triplet])
        elif triplet in genetic_code["aminos"]:
            stops_fuzzy.append(genetic_code["aminos"][triplet])
        else:
            stops_fuzzy.append("X")
    seq_AA += unresolve_aminoacid(stops_fuzzy)

    return seq_AA


def translate_fasta_dict(in_fasta_dict: dict, genetic_code_id: int, frame="guess", start_as_M=False):
    """
    Translates `in_fasta_dict` with `genetic_code_id`. If `frame` is not specified, every sequence
    is translated in the 6 reading frames and the one with the fewest stop codons is returned, when
    two or more are tied, the one with the fewest Xs is chosen.

    Parameters
    ----------
    in_fasta_dict : dict
        Dictionary produced by the function `fasta_to_dict` after processing a FASTA file
    genetic_code_id : int
        Genetic Code's numeric ID, see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    frame : str, int optional
        Reading frame, accepted values are: 1, 2, 3, -1, -2, -3, by default "guess"

    Returns
    -------
    dict
        A `fasta_to_dict` structure with the translated sequences
    """

    def guess_frame(seq, genetic_code, start_as_M):
        translations, counts_Xs, counts_stops = [], [], []
        for f in [1, 2, 3, -1, -2, -3]:
            translation = translate(seq, genetic_code, f, start_as_M)
            Xs = translation.count("X")
            if translation[-1] == "X":
                Xs -= 1
            counts_Xs.append(Xs)
            stops = translation.count("*")
            if translation[-1] == "*":
                stops -= 1
            counts_stops.append(stops)
            translations.append(translation)
        min_stops_idxs = [i for i in range(len(counts_stops)) if counts_stops[i] == min(counts_stops)]
        if len(min_stops_idxs) == 1:
            return translations[min_stops_idxs[0]]
        else:
            min_Xs = min([counts_Xs[i] for i in min_stops_idxs])
            for i in min_stops_idxs:
                if counts_Xs[i] == min_Xs:
                    return translations[i]

    gc = genetic_code(genetic_code_id)
    translated_fasta = {}
    if frame == "guess":
        for seq_name in in_fasta_dict:
            translated_fasta[seq_name] = {
                "description": in_fasta_dict[seq_name]["description"],
                "sequence": guess_frame(in_fasta_dict[seq_name]["sequence"], gc, start_as_M),
            }
    else:
        for seq_name in in_fasta_dict:
            translated_fasta[seq_name] = {
                "description": in_fasta_dict[seq_name]["description"],
                "sequence": translate(in_fasta_dict[seq_name]["sequence"], gc, frame, start_as_M),
            }

    return translated_fasta


def fix_premature_stops(in_fasta_dict: dict):
    """
    Returns 'None" if no premature stops were found or a fasta_dict with the premature stops
    replaced by 'X'

    Parameters
    ----------
    in_fasta_dict : dict
        Input fasta_dict has to be in aminoacids
    """
    has_premature_stops = False
    out_fasta_dict = {}
    for seq_name in in_fasta_dict:
        seq_out = (
            f"{in_fasta_dict[seq_name]['sequence'][:-1].replace('*', 'X')}"
            f"{in_fasta_dict[seq_name]['sequence'][-1]}"
        )
        seq_out = seq_out.replace(" ", "").replace("-", "")
        if in_fasta_dict[seq_name]["sequence"] != seq_out:
            has_premature_stops = True
        out_fasta_dict[seq_name] = {
            "sequence": seq_out,
            "description": in_fasta_dict[seq_name]["description"],
        }
    if has_premature_stops:
        return out_fasta_dict
    else:
        return None


def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence
    """
    return "".join([REV_COMP_DICT[n] for n in seq[::-1]])


sys.setrecursionlimit(settings.RECURSION_LIMIT)


def align_prots(s1, s2, method, scoring_matrix=PAM250):
    def gapless(s1, s2, aln):
        len_s1 = len(s1)
        len_s2 = len(s2)
        max_len = max(len(s1), len(s2))
        min_len = min(len(s1), len(s2))
        aln["s1_aln"] = s1
        aln["s2_aln"] = s2
        if len(s1) < len(s2):
            aln["s1_aln"] += "-" * (len(s2) - len(s1))
        elif len(s1) > len(s2):
            aln["s2_aln"] += "-" * (len(s1) - len(s2))
        for i in range(max_len - min_len + 1):
            matches, mismatches = 0, []
            s1_aln, s1_start, s1_end = s1, 0, len_s1
            s2_aln, s2_start, s2_end = s2, 0, len_s2
            if len_s2 < len_s1:
                s2_aln, s2_start, s2_end = list("-" * max_len), i, i + len_s2
            else:
                s1_aln, s1_start, s1_end = list("-" * max_len), i, i + len_s1
            for j in range(min_len):
                if len_s2 < len_s1:
                    pid = AA_PIDS["".join(sorted(f"{s1[i + j]}{s2[j]}"))]
                    matches += pid
                    if pid < 0.5:
                        mismatches.append(i + j)
                    s2_aln[i + j] = s2[j]
                else:
                    pid = AA_PIDS["".join(sorted(f"{s1[j]}{s2[i + j]}"))]
                    matches += pid
                    if pid < 0.5:
                        mismatches.append(i + j)
                    s1_aln[i + j] = s1[j]
            try:
                match_rate = matches / min_len
            except ZeroDivisionError:
                match_rate = 0.0
            if match_rate > aln["match_rate"]:
                aln["matches"] = matches
                aln["mismatches"] = mismatches
                aln["match_rate"] = match_rate
                aln["s1_aln"] = "".join(s1_aln)
                aln["s1_start"] = s1_start
                aln["s1_end"] = s1_end
                aln["s2_aln"] = "".join(s2_aln)
                aln["s2_start"] = s2_start
                aln["s2_end"] = s2_end

        return aln

    def best_score(val_matrix, row, col, nt_left, nt_top, method, scoring_matrix):
        left_score = val_matrix[row][col - 1] + scoring_matrix["gap"]
        diag_score = val_matrix[row - 1][col - 1] + scoring_matrix[f"{nt_left}{nt_top}"]
        up_score = val_matrix[row - 1][col] + scoring_matrix["gap"]
        if diag_score >= left_score and diag_score >= up_score:
            score = diag_score
            direction = 3
        elif left_score > diag_score and left_score >= up_score:
            score = left_score
            direction = 2
        else:
            score = up_score
            direction = 1
        if method == "nw":
            return score, direction
        if method == "sw":
            return max(score, 0), direction

    def direction_matrix(s1, s2, method, scoring_matrix):
        seq_left = f" {s1}"
        seq_top = f" {s2}"
        val_matrix = [[0] * len(seq_top) for i in range(len(seq_left))]
        dir_matrix = [[0] * len(seq_top) for i in range(len(seq_left))]
        for i in range(len(val_matrix)):
            dir_matrix[i][0] = 1
        for i in range(len(val_matrix[0])):
            dir_matrix[0][i] = 2
        dir_matrix[0][0] = 0
        if method == "nw":
            for i in range(len(val_matrix)):
                val_matrix[i][0] = i * scoring_matrix["gap"]
            for i in range(len(val_matrix[0])):
                val_matrix[0][i] = i * scoring_matrix["gap"]
            for L in range(1, len(seq_left)):
                for t in range(1, len(seq_top)):
                    val_matrix[L][t], dir_matrix[L][t] = best_score(
                        val_matrix, L, t, seq_left[L], seq_top[t], "nw", scoring_matrix
                    )
            return dir_matrix
        if method == "sw":
            top_score = 0
            opt_l = 0
            opt_t = 0
            for L in range(1, len(seq_left)):
                for t in range(1, len(seq_top)):
                    val_matrix[L][t], dir_matrix[L][t] = best_score(
                        val_matrix, L, t, seq_left[L], seq_top[t], "sw", scoring_matrix
                    )
                    if val_matrix[L][t] > top_score:
                        top_score = val_matrix[L][t]
                        opt_l = L
                        opt_t = t
            return dir_matrix, opt_l, opt_t

    def traceback_nw(dir_matrix, seq_left, seq_top):
        seq_left = f" {seq_left}"
        seq_top = f" {seq_top}"
        traceback_nw.align_left = ""
        traceback_nw.align_top = ""

        def stepback(dir_matrix, seq_left, seq_top, L, t):
            if L < 0 or t < 0:
                return
            if dir_matrix[L][t] == 3:
                stepback(dir_matrix, seq_left, seq_top, L - 1, t - 1)
                traceback_nw.align_top += seq_top[t]
                traceback_nw.align_left += seq_left[L]
            elif dir_matrix[L][t] == 2:
                stepback(dir_matrix, seq_left, seq_top, L, t - 1)
                traceback_nw.align_top += seq_top[t]
                traceback_nw.align_left += "-"
            elif dir_matrix[L][t] == 1:
                stepback(dir_matrix, seq_left, seq_top, L - 1, t)
                traceback_nw.align_top += "-"
                traceback_nw.align_left += seq_left[L]
            else:
                return

        stepback(dir_matrix, seq_left, seq_top, len(seq_left) - 1, len(seq_top) - 1)
        return traceback_nw.align_left, traceback_nw.align_top

    def traceback_sw(dir_matrix, opt_l, opt_t, seq_left, seq_top):
        seq_left = f" {seq_left}"
        seq_top = f" {seq_top}"
        traceback_sw.align_left = ""
        traceback_sw.align_top = ""

        def stepback(dir_matrix, seq_left, seq_top, L, t):
            if L < 0 or t < 0:
                return
            if dir_matrix[L][t] == 3:
                stepback(dir_matrix, seq_left, seq_top, L - 1, t - 1)
                traceback_sw.align_top += seq_top[t]
                traceback_sw.align_left += seq_left[L]
            elif dir_matrix[L][t] == 2:
                stepback(dir_matrix, seq_left, seq_top, L, t - 1)
                traceback_sw.align_top += seq_top[t]
                traceback_sw.align_left += "-"
            elif dir_matrix[L][t] == 1:
                stepback(dir_matrix, seq_left, seq_top, L - 1, t)
                traceback_sw.align_top += "-"
                traceback_sw.align_left += seq_left[L]
            else:
                return

        stepback(dir_matrix, seq_left, seq_top, opt_l, opt_t)
        add_top = ""
        add_left = ""
        if len(seq_left) > len(seq_top):
            add_left = seq_left[opt_l + 1 :]
            add_top = seq_top[opt_t + 1 :]
            add_top += (len(add_left) - len(add_top)) * "-"
        else:
            add_top = seq_top[opt_t + 1 :]
            add_left = seq_left[opt_l + 1 :]
            add_left += (len(add_top) - len(add_left)) * "-"
        return f"{traceback_sw.align_left}{add_left}", f"{traceback_sw.align_top}{add_top}"

    def needleman_wunsch(s1, s2, aln, scoring_matrix):
        s1_aln, s2_aln = traceback_nw(direction_matrix(s1, s2, "nw", scoring_matrix), s1, s2)
        aln["s1_start"] = len(s1_aln) - len(s1_aln.lstrip("-"))
        aln["s1_end"] = len(s1_aln.rstrip("-"))
        aln["s1_aln"] = s1_aln
        aln["s2_start"] = len(s2_aln) - len(s2_aln.lstrip("-"))
        aln["s2_end"] = len(s2_aln.rstrip("-"))
        aln["s2_aln"] = s2_aln
        aln_start = max(aln["s1_start"], aln["s2_start"])
        aln_end = min(aln["s1_end"], aln["s2_end"])
        for i in range(aln_start, aln_end):
            if "-" not in [s1_aln[i], s2_aln[i]]:
                pid = AA_PIDS["".join(sorted(f"{s1_aln[i]}{s2_aln[i]}"))]
                aln["matches"] += pid
                if pid < 0.5:
                    aln["mismatches"].append(i)
        try:
            aln["match_rate"] = aln["matches"] / (aln_end - aln_start)
        except ZeroDivisionError:
            aln["match_rate"] = 0.0

        return aln

    def smith_waterman(s1, s2, aln, scoring_matrix):
        dir_mat, opt_l, opt_t = direction_matrix(s1, s2, "sw", scoring_matrix)
        s1_aln, s2_aln = traceback_sw(dir_mat, opt_l, opt_t, s1, s2)
        aln["s1_start"] = len(s1_aln) - len(s1_aln.lstrip("-"))
        aln["s1_end"] = len(s1_aln.rstrip("-"))
        aln["s1_aln"] = s1_aln
        aln["s2_start"] = len(s2_aln) - len(s2_aln.lstrip("-"))
        aln["s2_end"] = len(s2_aln.rstrip("-"))
        aln["s2_aln"] = s2_aln
        aln_start = max(aln["s1_start"], aln["s2_start"])
        aln_end = min(aln["s1_end"], aln["s2_end"])
        for i in range(aln_start, aln_end):
            if "-" not in [s1_aln[i], s2_aln[i]]:
                pid = AA_PIDS["".join(sorted(f"{s1_aln[i]}{s2_aln[i]}"))]
                aln["matches"] += pid
                if pid < 0.5:
                    aln["mismatches"].append(i)
        try:
            aln["match_rate"] = aln["matches"] / (aln_end - aln_start)
        except ZeroDivisionError:
            aln["match_rate"] = 0.0

        return aln

    valid_methods = ["gapless", "nw", "sw"]
    if method not in valid_methods:
        return False

    # Alignment data template:
    aln = {
        "matches": 0,
        "mismatches": [],
        "match_rate": 0.0,
        "s1_start": 0,
        "s1_end": 0,
        "s1_aln": "",
        "s2_start": 0,
        "s2_end": 0,
        "s2_aln": "",
    }

    # Capitalize sequences and replace unusual
    s1 = s1.upper().replace("U", "C").replace("O", "X")
    s2 = s2.upper().replace("U", "C").replace("O", "X")

    if method == "gapless":
        return gapless(s1, s2, aln)
    elif method == "nw":
        return needleman_wunsch(s1, s2, aln, scoring_matrix)
    elif method == "sw":
        return smith_waterman(s1, s2, aln, scoring_matrix)
    else:
        return False


def pairwise_identity(seq1: str, seq2: str, seq_type: str, ignore_internal_gaps=False):
    """
    Given a pair of aligned sequences 'seq1' and 'seq2', calculate the identity of the overlapping
    region

    Parameters
    ----------
    seq1 : str
        Sequence 1
    seq2 : str
        Sequence 2
    seq_type : str
        'AA' if aminoacid, 'NT' if nucleotide
    ignore_internal_gaps : bool, optional
        Wheter to include the length of internal gaps in calculation or not, by default False

    Returns
    -------
    float
        Pairwise sequence identity as percentage
    """

    PIDS = {}
    if seq_type == "NT":
        PIDS = NT_PIDS
    elif seq_type == "AA":
        PIDS = AA_PIDS

    seq1 = seq1.upper()
    seq2 = seq2.upper()
    seq1_start, seq2_start = len(seq1) - len(seq1.lstrip("-")), len(seq2) - len(seq2.lstrip("-"))
    seq1_end, seq2_end = len(seq1.rstrip("-")), len(seq2.rstrip("-"))
    overlap_length = min(seq1_end, seq2_end) - max(seq1_start, seq2_start)
    matches = 0
    aligned_length = 0
    if overlap_length > 0:
        for pos in range(max(seq1_start, seq2_start), min(seq1_end, seq2_end)):
            pair = "".join(sorted(f"{seq1[pos]}{seq2[pos]}"))
            if pair != "--":
                if ignore_internal_gaps:
                    if "-" not in pair:
                        try:
                            matches += PIDS[pair]
                            aligned_length += 1
                        except KeyError:
                            continue
                else:
                    try:
                        matches += PIDS[pair]
                        aligned_length += 1
                    except KeyError:
                        continue
        try:
            return (matches / aligned_length) * 100
        except ZeroDivisionError:
            return 0.00
    else:
        return 0.00


def site_pairwise_identity(site: str, seq_type: str):
    """
    Calculates the mean pairwise identity of an alignment site (column), gap vs non-gap comparisons
    are counted as mismatches, but gap vs gap comparisons are excluded. Fractional matches are also
    possible when comparing ambiguities, for example R vs A is 0.5 matches because R is A or G.
    If terminal gaps are represented by '#' they are removed from the site before the calculation.

    Parameters
    ----------
    site : str
        Transposed alignment site (column)
    seq_type : str
        'AA' if aminoacid, 'NT' if nucleotide

    Returns
    -------
    float
        Mean pairwise identity of the site as percentage
    """

    PIDS = {}
    if seq_type == "NT":
        PIDS = NT_PIDS
    elif seq_type == "AA":
        PIDS = AA_PIDS

    site = site.replace("#", "").upper()
    len_site = len(site)
    if len_site < 2:
        return None
    else:
        ss_site = sorted(set(site))
        counts = {r: site.count(r) for r in ss_site}
        combos = (len_site * (len_site - 1)) / 2
        if "-" in counts:
            combos -= (counts["-"] * (counts["-"] - 1)) / 2
        if combos == 0:
            return [0, 0]
        matches = 0
        for ra in counts:
            try:
                matches += ((counts[ra] * (counts[ra] - 1)) / 2) * PIDS[ra + ra]
            except KeyError:
                continue
        i = 0
        for ra in ss_site[i:]:
            for rb in ss_site[i + 1 :]:
                try:
                    matches += counts[ra] * counts[rb] * PIDS[ra + rb]
                except KeyError:
                    continue
            i += 1
        return [(matches / combos) * 100, combos]


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
    in_fasta_dict,
    out_fasta_path,
    wrap=0,
    sort=False,
    shuffle=False,
    append=False,
    write_if_empty=False,
):
    """
    Saves a `in_fasta_dict` from function `fasta_to_dict()` as a FASTA file to `out_fasta_path`
    """
    compress = False
    if f"{out_fasta_path}".endswith(".gz"):
        compress = True
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
    if out_fasta_path.exists() and compress is True:
        if shutil.which("pigz"):
            pigz_compress(out_fasta_path, 2)
        elif shutil.which("gzip"):
            gzip_compress(out_fasta_path)
        out_fasta_path_gz = Path(out_fasta_path.parent, f"{out_fasta_path.name}.gz")
        if out_fasta_path_gz.exists():
            out_fasta_path = out_fasta_path_gz
    return out_fasta_path


def fasta_type(fasta_path):
    """
    Verify FASTA format and that sequence is only nucleotides, returns 'NT' when the file contains
    only nucleotides or nucleotide IUPAC ambiguity codes and 'AA' when it contains aminoacids.
    When other characters are found in the sequence or the files doesn't have headers it returns
    'invalid'
    """
    has_headers = False
    has_aminos = False
    if f"{fasta_path}".endswith(".gz"):
        opener = gzip.open
    else:
        opener = open
    with opener(fasta_path, "rt") as fasta_to_check:
        for line in fasta_to_check:
            line = line.strip("\n")
            if not line:
                continue
            if line.startswith(">"):
                has_headers = True
            else:
                if has_headers is True:
                    line_set = set(line.replace(" ", "").upper().replace("U", "T"))
                    if line_set - set(NT_IUPAC):
                        if not line_set - set(AA_IUPAC):
                            has_aminos = True
                            break
                else:
                    return "invalid"
    if has_headers and has_aminos:
        return "AA"
    elif has_headers and not has_aminos:
        return "NT"
    else:
        return "invalid"


def alignment_stats(fasta_dict, aln_type, coding: bool):
    """
    Tally number of sequences, sites (total, singletons, informative, constant), patterns (unique
    sites), average pairwise identity and percentage of missingness (including gaps and ambiguities)
    """

    def aligned_length(fasta_dict):
        if len(fasta_dict) < 2:
            return False
        length = -1
        for seq in fasta_dict:
            if length == -1:
                length = len(fasta_dict[seq]["sequence"])
            else:
                if len(fasta_dict[seq]["sequence"]) != length:
                    print(seq, fasta_dict[seq])
                    return False
        if length > 0:
            return length
        else:
            return False

    def num_samples(fasta_dict):
        sample_names = []
        for seq in fasta_dict:
            sample_names.append(seq.split(settings.SEQ_NAME_SEP)[0])

        return len(set(sample_names))

    def mark_terminal_gaps(fasta_dict):
        fasta_list = []
        for seq_name in fasta_dict:
            seq = fasta_dict[seq_name]["sequence"]
            left_gaps = "#" * (len(seq) - len(seq.lstrip("-")))
            right_gaps = "#" * (len(seq) - len(seq.rstrip("-")))
            fasta_list.append(f"{left_gaps}{seq.strip('-')}{right_gaps}")

        return fasta_list

    def transpose_aln(fasta_list: list, num_sites):
        sites = [""] * num_sites
        for seq in fasta_list:
            for pos in range(num_sites):
                sites[pos] += seq[pos]

        return sites

    def clean_patterns(sites: list, missing_chars: list):
        """
        Replace missing data with '-' for determination of pattern type. Lists of missing data
        symbols take from:
        http://www.iqtree.org/doc/Frequently-Asked-Questions#how-does-iq-tree-treat-gapmissingambiguous-characters
        We added '#" to represent terminal gaps.

        Parameters
        ----------
        sites : list
            Transposed alignment sites (columns)

        Returns
        -------
        list
            Transposed alignment sites with missing data replaced by '-'
        """
        clean = []
        for site in sites:
            site_out = ""
            for r in site:
                if r in missing_chars:
                    site_out += "-"
                else:
                    site_out += r
            clean.append(site_out)

        return clean

    def pattern_type(pattern, seq_type):
        pattern = pattern.replace("-", "")
        if not pattern:
            return "constant"
        r_sets = []
        if seq_type == "NT":
            r_sets = [set(NT_IUPAC[r]) for r in pattern]
        elif seq_type == "AA":
            r_sets = [set(AA_IUPAC[r]) for r in pattern]
        set_intersect = r_sets[0].intersection(*r_sets[1:])
        if set_intersect:
            return "constant"
        else:
            pat_expand = "".join(["".join(s) for s in r_sets])
            shared_chars = [c for c in set(pat_expand) if pat_expand.count(c) >= 2]
        if len(shared_chars) >= 2:
            return "informative"
        else:
            return "singleton"

    def gc_content(fasta_dict, seq_type, coding):
        if seq_type != "NT":
            return "NA", "NA", "NA", "NA"
        else:
            gc_bp = []
            for seq_name in fasta_dict:
                for r in fasta_dict[seq_name]["sequence"].replace("-", ""):
                    if r in NT_IUPAC:
                        g_prop = NT_IUPAC[r].count("G") / len(NT_IUPAC[r])
                        c_prop = NT_IUPAC[r].count("C") / len(NT_IUPAC[r])
                        gc_bp.append(g_prop + c_prop)
            gc_total = round(sum(gc_bp) / len(gc_bp) * 100, 2)
            if coding:
                gc_bp_p1 = [gc_bp[i] for i in range(len(gc_bp)) if i % 3 == 0]
                gc_bp_p2 = [gc_bp[i] for i in range(len(gc_bp)) if i % 3 == 1]
                gc_bp_p3 = [gc_bp[i] for i in range(len(gc_bp)) if i % 3 == 2]
                gc_codon_p1 = round(sum(gc_bp_p1) / len(gc_bp_p1) * 100, 2)
                gc_codon_p2 = round(sum(gc_bp_p2) / len(gc_bp_p2) * 100, 2)
                gc_codon_p3 = round(sum(gc_bp_p3) / len(gc_bp_p3) * 100, 2)
                return gc_total, gc_codon_p1, gc_codon_p2, gc_codon_p3
            else:
                return gc_total, "NA", "NA", "NA"

    stats = {
        "sequences": "NA",
        "samples": "NA",
        "avg_copies": 0.0,
        "sites": "NA",
        "informative": 0,
        "informativeness": 0.0,
        "uninformative": 0,
        "constant": 0,
        "singleton": 0,
        "patterns": "NA",
        "avg_pid": "NA",
        "missingness": "NA",
        "gc": "NA",
        "gc_codon_p1": "NA",
        "gc_codon_p2": "NA",
        "gc_codon_p3": "NA",
    }

    if aln_type == "NT":
        missing_chars = ["?", ".", "~", "O", "N", "X", "#"]
    elif aln_type == "AA":
        missing_chars = ["?", ".", "~", "*", "X", "#"]

    for seq_name in fasta_dict:
        seq_upper = fasta_dict[seq_name]["sequence"].upper()
        fasta_dict[seq_name]["sequence"] = seq_upper
    stats["sequences"] = len(fasta_dict)
    stats["samples"] = num_samples(fasta_dict)
    stats["avg_copies"] = round(stats["sequences"] / stats["samples"], 2)

    num_sites = aligned_length(fasta_dict)
    if not num_sites:
        return stats
    stats["sites"] = num_sites

    sites = transpose_aln(mark_terminal_gaps(fasta_dict), num_sites)
    ids_combos = [site_pairwise_identity(site, aln_type) for site in sites]
    ids_combos = [idc for idc in ids_combos if idc is not None]
    identities = 0
    combos = 0
    for idc in ids_combos:
        identities += idc[0] * idc[1]
        combos += idc[1]
    stats["avg_pid"] = round(identities / combos, 2)

    sites = clean_patterns(sites, missing_chars)
    for site in sites:
        stats[pattern_type(site, aln_type)] += 1
    stats["informativeness"] = round(stats["informative"] / stats["sites"] * 100, 2)
    stats["uninformative"] = stats["constant"] + stats["singleton"]

    stats["patterns"] = len(set(sites))

    sites_concat = "".join(sites)
    if len(sites_concat):
        missing_data = sites_concat.count("-")
        for c in missing_chars:
            missing_data += sites_concat.count(c)
        stats["missingness"] = round(missing_data / len(sites_concat) * 100, 2)

    gc_total, gc_codon_p1, gc_codon_p2, gc_codon_p3 = gc_content(fasta_dict, aln_type, coding)
    stats["gc"] = gc_total
    stats["gc_codon_p1"] = gc_codon_p1
    stats["gc_codon_p2"] = gc_codon_p2
    stats["gc_codon_p3"] = gc_codon_p3

    return stats


def sample_stats(fasta_dict, aln_type, coding):
    """
    Compile sample-wise statistic across alignments, last stage of `align` module

    Parameters
    ----------
    fasta_dict : dict
        FASTA as dictionary, from fasta_to_dict()
    aln_type : str
        NT for nucleotides or AA for aminoacids
    coding : bool
        Nucleotide alignment is coding or not

    Returns
    -------
    stats
        A dictionary with statistics of the sample within an alignment
    """
    stats = {}
    for seq_name in fasta_dict:
        if seq_name.endswith(f"{settings.SEQ_NAME_SEP}ref"):
            sam_name = seq_name
        else:
            sam_name = seq_name.split(settings.SEQ_NAME_SEP)[0]
        len_seq = len(fasta_dict[seq_name]["sequence"])
        len_gapped = len(fasta_dict[seq_name]["sequence"].strip("-"))
        ungapped_seq = fasta_dict[seq_name]["sequence"].replace("-", "")
        len_ungapped = len(ungapped_seq)
        ambigs = 0
        gc_bp = []
        gc_total = "NA"
        gc_codon_p1 = "NA"
        gc_codon_p2 = "NA"
        gc_codon_p3 = "NA"
        if aln_type == "AA":
            for r in ungapped_seq.upper():
                if r in list("BZJX?"):
                    ambigs += 1
        elif aln_type == "NT" and len_ungapped > 0:
            for r in ungapped_seq.upper():
                if r in list("MRWSYKVHDBN?"):
                    ambigs += 1
                if r in NT_IUPAC:
                    g_prop = NT_IUPAC[r].count("G") / len(NT_IUPAC[r])
                    c_prop = NT_IUPAC[r].count("C") / len(NT_IUPAC[r])
                    gc_bp.append(g_prop + c_prop)
            gc_total = round(sum(gc_bp) / len(gc_bp) * 100, 2)
            if coding:
                gc_bp_p1 = [gc_bp[i] for i in range(len(gc_bp)) if i % 3 == 0]
                gc_bp_p2 = [gc_bp[i] for i in range(len(gc_bp)) if i % 3 == 1]
                gc_bp_p3 = [gc_bp[i] for i in range(len(gc_bp)) if i % 3 == 2]
                gc_codon_p1 = round(sum(gc_bp_p1) / len(gc_bp_p1) * 100, 2)
                gc_codon_p2 = round(sum(gc_bp_p2) / len(gc_bp_p2) * 100, 2)
                gc_codon_p3 = round(sum(gc_bp_p3) / len(gc_bp_p3) * 100, 2)
        if sam_name not in stats:
            stats[sam_name] = {
                "len_total": len_seq,
                "len_gapped": [len_gapped],
                "len_ungapped": [len_ungapped],
                "ambigs": [ambigs],
                "gc": [gc_total],
                "gc_codon_p1": [gc_codon_p1],
                "gc_codon_p2": [gc_codon_p2],
                "gc_codon_p3": [gc_codon_p3],
                "num_copies": 1,
            }
        else:
            stats[sam_name]["len_gapped"].append(len_gapped)
            stats[sam_name]["len_ungapped"].append(len_ungapped)
            stats[sam_name]["ambigs"].append(ambigs)
            stats[sam_name]["gc"].append(gc_total)
            stats[sam_name]["gc_codon_p1"].append(gc_codon_p1)
            stats[sam_name]["gc_codon_p2"].append(gc_codon_p2)
            stats[sam_name]["gc_codon_p3"].append(gc_codon_p3)
            stats[sam_name]["num_copies"] += 1
    for sam_name in stats:
        if len(stats[sam_name]["len_gapped"]) > 1:
            for k in stats[sam_name]:
                if isinstance(stats[sam_name][k], list):
                    try:
                        stats[sam_name][k] = f"{statistics.mean(stats[sam_name][k]):.2f}"
                    except TypeError:
                        stats[sam_name][k] = "NA"
        else:
            for k in stats[sam_name]:
                if isinstance(stats[sam_name][k], list):
                    stats[sam_name][k] = stats[sam_name][k][0]

    return stats


def fasta_headers_to_spades(fasta_dict):
    """
    Given a `fasta_dict` from the function `fasta_to_dict()`, transform MEGAHIT or SKESA headers
    into Spades-like (Velvet has identical header format to Spades, afaik) headers like this :
    ```text
        '>NODE_0_length_274_cov_2.0[_k_157_flag_1]' [*optional]
    ```
    This makes it easier to parse and show in Bandage, also Scipio behaves better when FASTA headers
    don't contain spaces.
    The conversion must be done AFTER `megahit_toolkit contig2fastg` has been run on the original
    MEGAHIT FASTA file. If not, MEGAHIT can't process it and obtain an assembly graph.
    """

    assembler = "unknown"
    for name in fasta_dict:
        spades_header_regex = re.compile(SPADES_HEADER_REGEX)
        megahit_header_regex = re.compile(MEGAHIT_HEADER_REGEX)
        skesa_header_regex = re.compile(SKESA_HEADER_REGEX)
        header = f"{name} {fasta_dict[name]['description']}".strip()
        if spades_header_regex.match(header):
            assembler = "spades-like"
            return fasta_dict, assembler
        elif megahit_header_regex.match(header):
            assembler = "megahit"
        elif skesa_header_regex.match(header):
            assembler = "skesa"
        break

    fasta_dict_spades_headers = {}

    if assembler == "megahit":
        for name in fasta_dict:
            header = f"{name} {fasta_dict[name]['description']}".strip().split()

            # MEGAHIT header example: 'k157_0 flag=1 multi=2.0000 len=274'
            new_header = "_".join(
                [
                    "NODE",
                    header[0].split("_")[1],
                    "length",
                    header[3].split("=")[1],
                    "cov",
                    header[2].split("=")[1],
                    "k",
                    header[0].split("_")[0][1:],
                    "flag",
                    header[1].split("=")[1],
                ]
            )
            fasta_dict_spades_headers[new_header] = {
                "description": "",
                "sequence": fasta_dict[name]["sequence"],
            }
        return fasta_dict_spades_headers, assembler

    elif assembler == "skesa":
        for name in fasta_dict:
            header = name.strip().split("_")

            # SKESA header example: 'Contig_5_7.08923', if circular: 'Contig_5_7.08923_Circ'
            new_header = [
                "NODE",
                header[1],
                "length",
                f"{len(fasta_dict[name]['sequence'])}",
                "cov",
                header[2],
            ]
            if "_Circ" in header:
                new_header += ["circular"]
            new_header = "_".join(new_header)
            fasta_dict_spades_headers[new_header] = {
                "description": "",
                "sequence": fasta_dict[name]["sequence"],
            }
        return fasta_dict_spades_headers, assembler

    else:
        return fasta_dict, assembler


def scipio_yaml_to_dict(
    yaml_path,
    min_score,
    min_identity,
    min_coverage,
    marker_type,
    transtable,
    max_paralogs,
    paralog_tolerance,
    predict,
):
    """
    Process Scipio's YAML output, verify translation of each model, add extra aminoacid in gaps if
    they can be translated

    Parameters
    ----------
    yaml_path : Path
        Path to Scipio's YAML path
    min_identity : int
        Filter model with at leas this percentage of identity to the reference protein
    min_coverage : int
        Filter model with at leas this percentage of coverage to the reference protein
    marker_type : str
        Protein origin: NUC, PTD, MIT
    transtable : int
        Translation table number
    predict : bool
        Attempt translation of introns flagged as dubious by Scipio. If translation does not
        introduce preamture stop codons a new coding segment is added to the recovered protein

    Returns
    -------
    [type]
        [description]
    """

    def load_scipio_yaml(yaml_path):
        """
        Load YAML output from Scipio, returns a dictionary

        Parameters
        ----------
        yaml_path : Path or str
            Path to the .yaml file from Scipio

        Returns
        -------
        dict
            A dictionary containing the YAML data
        """

        nones = ["''", "~", "[]"]
        seqshift = {
            "nucl_start": None,
            "nucl_end": None,
            "dna_start": None,
            "dna_end": None,
            "prot_start": None,
            "prot_end": None,
        }
        skeys = {
            "          - nucl_start": "nucl_start",
            "            nucl_end": "nucl_end",
            "            dna_start": "dna_start",
            "            dna_end": "dna_end",
            "            prot_start": "prot_start",
            "            prot_end": "prot_end",
        }
        matching = {
            "type": None,
            "nucl_start": None,
            "nucl_end": None,
            "dna_start": None,
            "dna_end": None,
            "seq": "",
            "prot_start": None,
            "prot_end": None,
            "seqshifts": [],
            "mismatchlist": [],
            "undeterminedlist": [],
            "inframe_stopcodons": [],
            "translation": "",
            "overlap": None,
        }
        mkeys = {
            "      - type": "type",
            "        nucl_start": "nucl_start",
            "        nucl_end": "nucl_end",
            "        dna_start": "dna_start",
            "        dna_end": "dna_end",
            "        seq": "seq",
            "        prot_start": "prot_start",
            "        prot_end": "prot_end",
            "        seqshifts": "seqshifts",
            "        mismatchlist": "mismatchlist",
            "        undeterminedlist": "undeterminedlist",
            "        inframe_stopcodons": "inframe_stopcodons",
            "        translation": "translation",
            "        overlap": "overlap",
        }
        hit = {
            "ID": None,
            "status": None,
            "reason": "",
            "prot_len": None,
            "prot_start": None,
            "prot_end": None,
            "prot_seq": None,
            "target": None,
            "target_len": None,
            "strand": None,
            "dna_start": None,
            "dna_end": None,
            "matches": None,
            "mismatches": None,
            "undetermined": None,
            "unmatched": None,
            "additional": None,
            "score": None,
            "upstream": None,
            "upstream_gap": None,
            "matchings": [],
            "stopcodon": None,
            "downstream": None,
            "downstream_gap": None,
        }
        hkeys = {
            "  - ID": "ID",
            "    status": "status",
            "    reason": "reason",
            "    prot_len": "prot_len",
            "    prot_start": "prot_start",
            "    prot_end": "prot_end",
            "    prot_seq": "prot_seq",
            "    target": "target",
            "    target_len": "target_len",
            "    strand": "strand",
            "    dna_start": "dna_start",
            "    dna_end": "dna_end",
            "    matches": "matches",
            "    mismatches": "mismatches",
            "    undetermined": "undetermined",
            "    unmatched": "unmatched",
            "    additional": "additional",
            "    score": "score",
            "    upstream": "upstream",
            "    upstream_gap": "upstream_gap",
            "    matchings": "matchings",
            "    stopcodon": "stopcodon",
            "    downstream": "downstream",
            "    downstream_gap": "downstream_gap",
        }

        ints = [
            "    prot_len",
            "    prot_start",
            "    prot_end",
            "    target_len",
            "    dna_start",
            "    dna_end",
            "    matches",
            "    mismatches",
            "    undetermined",
            "    unmatched",
            "    additional",
            "    score",
            "        nucl_start",
            "        nucl_end",
            "        dna_start",
            "        dna_end",
            "        prot_start",
            "        prot_end",
            "        overlap",
            "          - nucl_start",
            "            nucl_end",
            "            dna_start",
            "            dna_end",
            "            prot_start",
            "            prot_end",
        ]

        yaml = {}
        with open(yaml_path, "rt") as yaml_in:
            protein = ""
            hit_num = ""
            for line in yaml_in:
                line = line.strip("\n")
                if all([not line.startswith("#"), line != "---", line]):
                    if not line.startswith(" "):
                        protein = line.split("_(")[0] if "_(" in line else line.split(":")[0]
                        hit_num = int(line.split("_(")[1].replace("):", "")) if "_(" in line else 0
                    elif ": " in line:
                        k, v = line.split(": ")[0], ": ".join(line.split(": ")[1:])
                        if v in nones:
                            continue
                        if k in ints:
                            try:
                                v = int(v)
                            except ValueError:
                                pass
                        else:
                            if k == "    target" and " " in v:
                                v = v.split()[0].strip("'").strip('"')
                            if v.startswith("'"):
                                v = v.strip("'").replace(" ", "")
                            if v.startswith("["):
                                try:
                                    v = [int(n) for n in v[1:-1].split(", ")]
                                except ValueError:
                                    pass
                        if k == "  - ID":
                            if protein not in yaml:
                                yaml[protein] = {hit_num: [copy.deepcopy(hit)]}
                            else:
                                if hit_num not in yaml[protein]:
                                    yaml[protein][hit_num] = [copy.deepcopy(hit)]
                                else:
                                    yaml[protein][hit_num].append(copy.deepcopy(hit))
                            yaml[protein][hit_num][-1][hkeys[k]] = v
                        elif k in hkeys:
                            yaml[protein][hit_num][-1][hkeys[k]] = v
                        elif k == "      - type":
                            yaml[protein][hit_num][-1]["matchings"].append(copy.deepcopy(matching))
                            yaml[protein][hit_num][-1]["matchings"][-1][mkeys[k]] = v
                        elif k in mkeys:
                            yaml[protein][hit_num][-1]["matchings"][-1][mkeys[k]] = v
                        elif k == "          - nucl_start":
                            yaml[protein][hit_num][-1]["matchings"][-1]["seqshifts"].append(
                                copy.deepcopy(seqshift)
                            )
                            yaml[protein][hit_num][-1]["matchings"][-1]["seqshifts"][-1][skeys[k]] = v
                        elif k in skeys:
                            yaml[protein][hit_num][-1]["matchings"][-1]["seqshifts"][-1][skeys[k]] = v
                    elif line.startswith("          - "):
                        mm = int(line.replace("          - ", ""))
                        yaml[protein][hit_num][-1]["matchings"][-1]["mismatchlist"].append(mm)

        return yaml

    def insert_shifts(mat: dict):
        """
        Scipio provides a list of sequence shifts, these include insertions of new aminoacids by
        adding 3n nucleotides to specific positions, as well as insertion of single-Ns or double-Ns
        in order to repair frameshifts. This functions inserts the proper number of Ns to
        fix errors in translation due to shifts in reading frame.

        Parameters
        ----------
        mat : dict
            Contains all the needed data for an exon (protein coordinates, contig coordinates,
            translation, etc.), it follows Scipio's structure.

        Returns
        -------
        str, str
            First element of the tuple is the sequence with inserted shifts in nucleotides,
            second element is its corresponding translation as it comes from Scipio.
        """
        translation_out = mat["translation"].strip("'")
        seq_out = ""
        inserts = {}
        if mat["seqshifts"]:
            for shift in mat["seqshifts"]:
                shift_len = shift["dna_end"] - shift["dna_start"]
                outside_codons = shift_len % 3
                if outside_codons > 0:
                    pos = shift["dna_start"] - mat["dna_start"] + shift_len - 1
                    inserts[pos] = "N" * (3 - outside_codons)
            for pos in range(len(mat["seq"])):
                seq_out += mat["seq"][pos]
                if pos in inserts:
                    seq_out += inserts[pos]
        else:
            seq_out = mat["seq"]
        if seq_out == "":
            translation_out = ""

        return seq_out, translation_out

    def concat_cds(mod: dict):
        """
        Concatenate coding sequence in the model

        Parameters
        ----------
        mod : dict
            Gene model formatted by function parse_model()

        Returns
        -------
        str, str
            First element in tuple is concatenated sequence in nucleotide, second tuple item is
            concatenated sequence in aminoacid as it comes from Scipio
        """
        cds_nt, cds_aa = "", ""
        for i in range(len(mod["mat_types"])):
            if mod["mat_types"][i] == "exon":
                cds_nt += mod["mat_nt"][i]
                cds_aa += mod["mat_aa"][i]

        return cds_nt, cds_aa

    def align_exon_nt_to_scipio_aa(seq_nt: str, gencode: dict, seq_aa: str):
        """
        Given the exon sequence in nucleotides, translate in the three positive reading frames and
        align each translation to Scipio's provided translation to determine the exon's correct
        reading frame as well as the number of leading and trailing bases outside of the reading
        frame

        Parameters
        ----------
        seq_nt : str
            Exon sequence in nucleotides
        gencode : dict
            Genetic code, output from function genetic_code()
        seq_aa : str
            Scipio's provided translation for the exon

        Returns
        -------
        dict
            Dictionary with the alignment information, including reading frame and number of
            leading and trailing untranslated nucleotides
        """
        rfs = {i: translate(seq_nt, gencode, frame=i, start_as_M=False) for i in [1, 2, 3]}
        alns = {i: align_prots(rfs[i], seq_aa, "gapless") for i in rfs}
        rf, aln = max(alns.items(), key=(lambda x: x[1]["match_rate"]))  # By GS
        aln["rf"] = rf
        aln["lead"] = rf - 1
        aln["trail"] = (len(seq_nt) - aln["lead"]) % 3

        return aln

    def remove_leading(mod: dict, i: int, lead: int, aln: dict):
        """
        Remove leading untranslated nucleotides from exonic sequence, add them to the end of the
        previous non-coding segment, and adjust match starting coordinate if needed

        Parameters
        ----------
        mod : dict
            Gene model formatted by function parse_model()
        i : int
            Current exon index in the model
        lead : int
            Number of leading nucleotides to remove
        aln : dict
            Alignment output from function align_nt_to_aa()

        Returns
        -------
        dict
            Modified gene model
        """
        non_coding = ["intron", "intron?", "gap", "upstream"]
        if lead > 0:
            # Remove leading bases from exonic nucleotide sequence
            mod["hit_starts"][i] += lead
            mod["mat_nt"][i] = mod["mat_nt"][i][lead:]
            # Append removed bases to previous non-coding segment if in same contig
            if (
                i > 0
                and mod["mat_types"][i - 1] in non_coding
                and mod["hit_ids"][i] == mod["hit_ids"][i - 1]
            ):
                mod["mat_nt"][i - 1] += mod["mat_nt"][i][:lead]
        # Remove initial extra aminoacid from translation if needed
        if aln["s1_start"] > 0:
            mod["ref_starts"][i] += aln["s1_start"]
            mod["mat_aa"][i] = mod["mat_aa"][i][aln["s1_start"] :]

        return mod

    def remove_trailing(mod: dict, i: int, trail: int, aln: dict, mod_len: int):
        """
        Remove trailing untranslated nucleotides from exonic sequence, add them to the beggining of
        the next non-coding segment, and adjust match ending coordinate if needed

        Parameters
        ----------
        mod : dict
            Gene model formatted by function parse_model()
        i : int
            Current exon index in the model
        trail : int
            Number of trailing nucleotides to remove
        aln : dict
            Alignment output from function align_nt_to_aa()
        mod_len : int
            Length of mod["hit_ids"]

        Returns
        -------
        dict
            Modified gene model
        """
        non_coding = ["intron", "intron?", "gap", "downstream"]
        if trail > 0:
            # Remove trailing bases from exonic nucleotide sequence
            mod["hit_ends"][i] -= trail
            mod["mat_nt"][i] = mod["mat_nt"][i][:-trail]
            # Prepend removed bases to next non-coding segment if in same contig
            if (
                i < mod_len - 1
                and mod["mat_types"][i + 1] in non_coding
                and mod["hit_ids"][i] == mod["hit_ids"][i + 1]
            ):
                mod["mat_nt"][i + 1] = f"{mod['mat_nt'][i][-trail:]}{mod['mat_nt'][i + 1]}"
        # Remove final extra aminoacid from translation if needed
        if aln["s1_end"] < aln["s2_end"]:
            mod["ref_ends"][i] -= aln["s2_end"] - aln["s1_end"]
            mod["mat_aa"][i] = mod["mat_aa"][i][: aln["s1_end"]]

        return mod

    def remove_short_terminal_exon(mod: dict):
        mod_len = len(mod["mat_types"])
        try:
            exon_indexes = [i for i in range(mod_len) if mod["mat_types"][i] == "exon"]
        except ValueError:
            return mod
        if len(exon_indexes) > 1:
            last_exon_idx = exon_indexes[-1]
            if (
                mod["mat_types"][last_exon_idx - 1] == "intron"
                and len(mod["mat_nt"][last_exon_idx]) < settings.SCIPIO_MIN_LEN_FINAL_EXON
            ):
                last_intron_idx = last_exon_idx - 1
                mod["ref_starts"] = mod["ref_starts"][:last_intron_idx]
                mod["ref_ends"] = mod["ref_ends"][:last_intron_idx]
                mod["hit_ids"] = mod["hit_ids"][:last_intron_idx]
                mod["hit_starts"] = mod["hit_starts"][:last_intron_idx]
                mod["hit_ends"] = mod["hit_ends"][:last_intron_idx]
                mod["mat_types"] = mod["mat_types"][:last_intron_idx]
                mod["mat_notes"] = mod["mat_notes"][:last_intron_idx]
                mod["mat_nt"] = mod["mat_nt"][:last_intron_idx]
                mod["mat_aa"] = mod["mat_aa"][:last_intron_idx]
                ref_ends = list(filter(None, mod["ref_ends"]))
                mod["mismatches"] = list(filter(lambda x: x <= max(ref_ends), mod["mismatches"]))
                if len(set(mod["hit_ids"])) < len(mod["hit_contigs"]):
                    mod["hit_contigs"] = mod["hit_contigs"][:-1]
                return mod
            else:
                return mod
        else:
            return mod

    def fix_model(mod: dict, gencode: dict, marker_type: str):
        """
        When translation of concatenated exons does not match the translation given by Scipio, we
        need to check every exon for inconsistencies like extra dangling nucleotides that break the
        reading frame

        Parameters
        ----------
        mod : dict
            Gene model formatted by function parse_model()
        gencode : dict
            Genetic code, output from function genetic_code()
        marker_type : str
            "NUC", "PTD", or "MIT"

        Returns
        -------
        dict
            Modified gene model
        """

        # Scipio rarely finds a distant, ver short terminal exon for some proteins like psbB, we
        # need to remove that exon and fix the rest of the model (remove intron before the short
        # terminal exon, remove stop codon sequence if found after the short exon)
        if marker_type in ["PTD", "MIT"]:
            mod = remove_short_terminal_exon(mod)

        # Simply exit if model has no exons
        if "exon" not in mod["mat_types"]:
            return None

        mod_len = len(mod["mat_types"])

        # 1. Best case scenario: gene model translates exactly as Scipio and doesn't need fixing
        cds_nt, cds_aa = concat_cds(mod)
        if translate(cds_nt, gencode, frame=1, start_as_M=False) == cds_aa and len(cds_nt) % 3 == 0:
            return mod

        # 2. Fix cases where Scipio's translation is longer than possible given the exon seq
        mat_alns = [None] * mod_len
        for i in range(mod_len):
            if mod["mat_types"][i] == "exon":
                aln = align_exon_nt_to_scipio_aa(mod["mat_nt"][i], gencode, mod["mat_aa"][i])
                if aln["s1_start"] > 0:
                    start_aa = aln["s1_start"] - 1 if aln["lead"] > 1 else aln["s1_start"]
                    mod["mat_aa"][i] = mod["mat_aa"][i][start_aa:]
                if aln["s2_end"] - aln["s1_end"] > 0:
                    end_aa = aln["s1_end"] + 1 if aln["trail"] > 1 else aln["s1_end"]
                    mod["mat_aa"][i] = mod["mat_aa"][i][:end_aa]
                mat_alns[i] = align_exon_nt_to_scipio_aa(mod["mat_nt"][i], gencode, mod["mat_aa"][i])
        cds_nt, cds_aa = concat_cds(mod)
        if translate(cds_nt, gencode, frame=1, start_as_M=False) == cds_aa and len(cds_nt) % 3 == 0:
            return mod

        # 3. Fix dangling nucleotides in exons when gaps are found, fixes coordinates as well
        for i in range(mod_len):
            if mod["mat_types"][i] == "exon":
                if "gap before" in mod["mat_notes"][i]:
                    mod = remove_leading(mod, i, mat_alns[i]["lead"], mat_alns[i])
                    mat_alns[i] = align_exon_nt_to_scipio_aa(mod["mat_nt"][i], gencode, mod["mat_aa"][i])
                if "gap after" in mod["mat_notes"][i]:
                    mod = remove_trailing(mod, i, mat_alns[i]["trail"], mat_alns[i], mod_len)
                    mat_alns[i] = align_exon_nt_to_scipio_aa(mod["mat_nt"][i], gencode, mod["mat_aa"][i])
        cds_nt, cds_aa = concat_cds(mod)
        if translate(cds_nt, gencode, frame=1, start_as_M=False) == cds_aa and len(cds_nt) % 3 == 0:
            return mod

        # 4. Remove extra dangling bases that break reading frame
        prev_i = None
        last_exon_i = [i for i in range(mod_len) if mod["mat_types"][i] == "exon"][-1]
        for i in range(mod_len):
            if mod["mat_types"][i] == "exon":
                # Remove extra leading nucleotides/aminoacids from first exon
                if prev_i is None:
                    if mat_alns[i]["lead"] > 0:
                        mod = remove_leading(mod, i, mat_alns[i]["lead"], mat_alns[i])
                # Continue with remaining exons
                else:
                    if mat_alns[prev_i]["trail"] + mat_alns[i]["lead"] in [0, 3]:
                        pass
                    elif mat_alns[prev_i]["trail"] == 0 and mat_alns[i]["lead"] in [1, 2]:
                        mod = remove_leading(mod, i, mat_alns[i]["lead"], mat_alns[i])
                    elif mat_alns[prev_i]["trail"] == 1 and mat_alns[i]["lead"] in [0, 1]:
                        mod = remove_trailing(mod, prev_i, 1, mat_alns[prev_i], mod_len)
                        mod = remove_leading(mod, i, mat_alns[i]["lead"], mat_alns[i])
                    elif mat_alns[prev_i]["trail"] == 2:
                        if mat_alns[i]["lead"] == 0:
                            mod = remove_trailing(mod, prev_i, 2, mat_alns[prev_i], mod_len)
                        elif mat_alns[i]["lead"] == 2:
                            # If last exon has an extra terminal aminoacid
                            if mat_alns[prev_i]["s1_end"] < mat_alns[prev_i]["s2_end"]:
                                # Trim a single base from the start of current exon
                                mod = remove_leading(mod, i, 1, mat_alns[i])
                            # If last exon has two extra bp but not extra amino
                            else:
                                # Trim a single base from the end of last exon
                                mod = remove_trailing(mod, prev_i, 1, mat_alns[prev_i], mod_len)
                # Remove extra trailing nucleotides/aminoacids from last exon
                if i == last_exon_i and mat_alns[i]["trail"] > 0:
                    mod = remove_trailing(mod, i, mat_alns[i]["trail"], mat_alns[i], mod_len)
                mat_alns[i] = align_exon_nt_to_scipio_aa(mod["mat_nt"][i], gencode, mod["mat_aa"][i])
                prev_i = i

        cds_nt, cds_aa = concat_cds(mod)
        cds_nt_translated = translate(cds_nt, gencode, frame=1, start_as_M=False)
        if cds_nt_translated == cds_aa:
            return mod
        else:
            aln = align_prots(cds_nt_translated, cds_aa, "gapless")
            mismatches = len(aln["mismatches"]) - sum(aln["s2_aln"][i] == "X" for i in aln["mismatches"])
            if mismatches <= settings.SCIPIO_MAX_MISMATCHES:
                return mod
            else:
                return None

    def translate_between_exons(mod: dict, i: int, gencode: dict, predict: bool, min_identity: float):
        if (
            len(mod["mat_nt"][i]) > 2
            and mod["mat_types"][i - 1] == mod["mat_types"][i + 1] == "exon"
            and mod["hit_ids"][i - 1] == mod["hit_ids"][i] == mod["hit_ids"][i + 1]
        ):
            # Align previous exon to its translation to find out its trailing bases
            prev_aln = align_exon_nt_to_scipio_aa(mod["mat_nt"][i - 1], gencode, mod["mat_aa"][i - 1])
            # Align next exon to to its translation to find out its leading bases
            next_aln = align_exon_nt_to_scipio_aa(mod["mat_nt"][i + 1], gencode, mod["mat_aa"][i + 1])
            # Take trailing bases from previous and leading bases from next and add to current chunk
            current_chunk = ""
            if prev_aln["trail"] > 0:
                current_chunk += mod["mat_nt"][i - 1][-prev_aln["trail"] :]
            current_chunk += mod["mat_nt"][i]
            current_chunk += mod["mat_nt"][i + 1][: next_aln["lead"]]

            # 'seq_chunks', 'lead', and 'trail' only get modified when translation is acceptable
            seq_chunks = None
            lead = None
            trail = None

            if mod["mat_types"][i] == "gap":
                # Attempt translation of current chunk and then alignment to unmatched segment of
                # reference protein ('prot_gap'):
                # ref_seq[0] = protein sequence, ref_seq[1] = starting coordinate
                ref_seq = mod["ref_seqs"][mod["hit_ids"][i]]
                prot_gap = ref_seq[0][
                    max(mod["ref_ends"][i - 1] - ref_seq[1], 0) : min(
                        mod["ref_starts"][i + 1] - ref_seq[1], len(ref_seq[0])
                    )
                ]
                len_gap = len(prot_gap)
                max_gap_size = mod["ref_size"] * settings.SCIPIO_MAX_GAP_AS_REF_PROP
                min_gap_identity = (min_identity / 100) - settings.SCIPIO_MAX_GAP_DELTA_IDENTITY

                if len_gap > 0:
                    len_chunk = len(current_chunk) // 3
                    overlap = min(len_gap, len_chunk) / max(len_gap, len_chunk)
                    # Don't align if 'current_chunk' is too long compared to the unmatched protein
                    if len_chunk <= max_gap_size and len_gap * len_chunk <= settings.RECURSION_LIMIT:
                        rfs = {
                            i: translate(current_chunk, gencode, frame=i, start_as_M=False)
                            for i in [1, 2, 3]
                        }
                        # Global alignment if they overlap at least 80%
                        if overlap >= 0.8:
                            alns = {i: align_prots(rfs[i], prot_gap, "nw") for i in rfs}
                        # Local alignment if the overlap is smaller
                        else:
                            alns = {i: align_prots(rfs[i], prot_gap, "sw") for i in rfs}
                        # Prioritize filling gaps without skipping leading or trailing bases
                        if (
                            len(current_chunk) % 3 == 0
                            and "*" not in rfs[1]
                            and alns[1]["match_rate"] >= min_gap_identity
                        ):
                            lead, trail = 0, 0
                            seq_chunks = ["", mod["mat_nt"][i].upper()]
                            mod["mismatches"] += [
                                pos + mod["ref_ends"][i - 1] for pos in alns[1]["mismatches"]
                            ]
                        # Rank reading frames by their match rate penalized by the number of stops
                        else:
                            for rf in alns:
                                alns[rf]["match_rate"] *= settings.SCIPIO_STOP_PENALTY ** alns[rf][
                                    "s1_aln"
                                ].count("*")
                            rf, aln = max(alns.items(), key=(lambda x: x[1]["match_rate"]))
                            if aln["match_rate"] >= min_gap_identity and "*" not in aln["s1_aln"]:
                                lead = rf - 1
                                trail = (len(current_chunk) - lead) % 3
                                if trail > 0:
                                    seq_chunks = [
                                        f"{mod['mat_nt'][i][:lead].lower()}",
                                        f"{mod['mat_nt'][i][lead:-trail].upper()}",
                                        f"{mod['mat_nt'][i][-trail:].lower()}",
                                    ]
                                else:
                                    seq_chunks = [
                                        f"{mod['mat_nt'][i][:lead].lower()}",
                                        f"{mod['mat_nt'][i][lead:].upper()}",
                                    ]
                                mod["mismatches"] += [
                                    pos + mod["ref_ends"][i - 1] for pos in aln["mismatches"]
                                ]

            if predict:
                if mod["mat_types"][i] == "intron?" or mod["mat_types"][i] == "gap":
                    rf1 = translate(current_chunk, gencode, frame=1, start_as_M=False)
                    if len(current_chunk) % 3 == 0 and "*" not in rf1:
                        lead, trail = 0, 0
                        seq_chunks = ["", mod["mat_nt"][i].upper()]

                # # When predicting we think we should only add the intervening segment if it can
                # # be completely translated (i.e. divisible by 3 and without stop codons), the
                # # following code would add translation even if it includes frameshifts
                # rfs = {i: translate(current_chunk, gencode, frame=i, start_as_M=False)
                #        for i in [1,2,3]}
                # no_stops = {i: rfs[i] for i in rfs if not "*" in rfs[i]}
                # if no_stops:
                #     if 1 in no_stops:
                #         if len(current_chunk) % 3 == 0:
                #             lead, trail = 0, 0
                #             seq_chunks = ["", mod["mat_nt"][i].upper()]
                #     else:
                #         rf, aa = max(no_stops.items(), key=(lambda x: len(x[1])))
                #         lead = rf - 1
                #         trail = (len(current_chunk) - lead) % 3
                #         if trail > 0:
                #             seq_chunks = [f'{mod["mat_nt"][i][:lead].lower()}',
                #                           f'{mod["mat_nt"][i][lead:-trail].upper()}',
                #                           f'{mod["mat_nt"][i][-trail:].lower()}']
                #         else:
                #             seq_chunks = [f'{mod["mat_nt"][i][:lead].lower()}',
                #                           f'{mod["mat_nt"][i][lead:].upper()}']

            if seq_chunks:
                mod["seq_flanked"] += "".join(seq_chunks)
                mod["seq_gene"] += "".join(seq_chunks)
                mod["seq_nt"] += seq_chunks[1]
                mod["ref_starts"][i] = mod["ref_ends"][i - 1]
                mod["ref_ends"][i] = mod["ref_starts"][i + 1]
                mod["hit_starts"][i] = mod["hit_ends"][i - 1] + lead
                mod["hit_ends"][i] = mod["hit_starts"][i + 1] - trail
                mod["mat_types"][i] += "/exon"
                if "gap" in mod["mat_types"][i]:
                    mod["mat_notes"][i] += "/filled"
                elif "intron?" in mod["mat_types"][i]:
                    mod["mat_notes"][i] += "/predicted"
            else:
                mod["seq_flanked"] += mod["mat_nt"][i].lower()
                mod["seq_gene"] += mod["mat_nt"][i].lower()

        else:
            mod["seq_flanked"] += mod["mat_nt"][i].lower()
            mod["seq_gene"] += mod["mat_nt"][i].lower()

        return mod

    def merge_adjoining_exons(mod: dict):
        last_exon = None
        for i in range(len(mod["mat_types"])):
            if mod["mat_types"][i].endswith("exon"):
                if last_exon is None:
                    last_exon = i
                else:
                    if (
                        mod["hit_ids"][last_exon] == mod["hit_ids"][i]
                        and mod["hit_starts"][i] <= mod["hit_ends"][last_exon]
                    ):
                        mod["ref_ends"][last_exon] = mod["ref_ends"][i]
                        mod["hit_ends"][last_exon] = mod["hit_ends"][i]
                        mod["ref_starts"][i] = None
                        mod["ref_ends"][i] = None
                        mod["hit_starts"][i] = None
                        mod["hit_ends"][i] = None
                    else:
                        last_exon = i

        return mod

    def check_gaps_and_concat_seqs(mod: dict, gencode: dict, predict: bool, min_identity: float):
        # If the model doesn't end with 'downstream' or 'exon' trim it until last piece is an 'exon'
        if mod["mat_types"][-1] not in ["exon", "downstream"]:
            for i in reversed(range(len(mod["mat_types"]))):
                if mod["mat_types"][i] == "exon":
                    for data in [
                        "ref_starts",
                        "ref_ends",
                        "hit_ids",
                        "hit_starts",
                        "hit_ends",
                        "mat_types",
                        "mat_notes",
                        "mat_nt",
                        "mat_aa",
                    ]:
                        mod[data] = mod[data][: i + 1]
                    break

        # Proceed with the checkup and concatenation
        last_exon = None
        for i in range(len(mod["mat_types"])):
            if mod["mat_types"][i] == "upstream":
                mod["seq_flanked"] += mod["mat_nt"][i].lower()
            if mod["mat_types"][i] == "exon":
                if settings.FILL_GAP_WITH_X:
                    if last_exon is not None:  # if this is not the first exon
                        gap_len = mod["ref_starts"][i] - mod["ref_ends"][last_exon]
                        if (
                            "gap before" in mod["mat_notes"][i]
                            and "gap after" in mod["mat_notes"][last_exon]
                        ) or gap_len >= settings.SCIPIO_MIN_GAP_LEN_TO_X:
                            if "/exon" not in mod["mat_types"][i - 1]:
                                mod["seq_nt"] += "n" * 3 * gap_len
                mod["seq_flanked"] += mod["mat_nt"][i].upper()
                mod["seq_gene"] += mod["mat_nt"][i].upper()
                mod["seq_nt"] += mod["mat_nt"][i].upper()
                last_exon = i
            if mod["mat_types"][i] in ["gap", "intron?"]:
                mod = translate_between_exons(mod, i, gencode, predict, min_identity)
            if mod["mat_types"][i] in ["intron", "separator"]:
                mod["seq_flanked"] += mod["mat_nt"][i].lower()
                mod["seq_gene"] += mod["mat_nt"][i].lower()
            if mod["mat_types"][i] == "stopcodon":
                mod["seq_flanked"] += mod["mat_nt"][i].upper()
                mod["seq_gene"] += mod["mat_nt"][i].upper()
            if mod["mat_types"][i] == "downstream":
                if (
                    mod["mat_types"][i - 1] == "stopcodon"
                    and mod["mat_nt"][i - 1] == mod["mat_nt"][i][:3]
                ):
                    mod["seq_flanked"] += mod["mat_nt"][i][3:].lower()
                else:
                    mod["seq_flanked"] += mod["mat_nt"][i].lower()
        mod["seq_aa"] = translate(mod["seq_nt"], gencode, frame=1, start_as_M=False)

        return mod

    def concat_coords(ids: list, starts: list, ends: list):
        coords = ""
        last_exon = None
        for i in range(len(ids)):
            if starts[i] is not None:
                coord = f"{abs(starts[i])}-{abs(ends[i])}"
                if last_exon is None:
                    coords += coord
                else:
                    if ids[i] == ids[last_exon]:
                        coords += f",{coord}"
                    else:
                        coords += f"\n{coord}"
                last_exon = i

        return coords

    def calc_prot_len_matched(mod: dict):
        starts = [
            mod["ref_starts"][i]
            for i in range(len(mod["mat_types"]))
            if ("predicted" not in mod["mat_notes"][i] and "exon" in mod["mat_types"][i])
        ]
        ends = [
            mod["ref_ends"][i]
            for i in range(len(mod["mat_types"]))
            if ("predicted" not in mod["mat_notes"][i] and "exon" in mod["mat_types"][i])
        ]
        prot_len_matched = sum([j - i for i, j in zip(starts, ends)])
        # In very rare cases Scipio returns overlaps, especially when hit is spread over multiple
        # contigs, discount the overlap so it doesn't influence 'wscore'
        for i in range(len(starts) - 1):
            overlap = starts[i + 1] - ends[i]
            if overlap < 0:
                prot_len_matched -= abs(overlap) + 1

        return prot_len_matched

    def parse_model(
        yaml_mod: dict,
        protein_name: str,
        marker_type: str,
        gencode: dict,
        predict: bool,
        min_identity: float,
    ):
        """
        A gene model is composed by several sequence matchings (e.g. exons, introns, etc.) which can
        be found across several contigs. This function organizes all relevant matching across
        multiple contigs in a dictionary of lists for easy manipulation.

        Parameters
        ----------
        mod : dict
            Gene model directly imported from Scipio's YAML output
        protein_name : str
            Name of protein sequence used as reference
        marker_type : str
            NUC, PTD, or PTD

        Returns
        -------
        [type]
            Model data parsed as a dictionary
        """

        # This is the dictionary structure for containing the gene model:
        mod = {
            "ref_name": "",  # full name of protein sequence used as reference
            "ref_seqs": {},  # segments of reference protein matched, with start and end coords
            "ref_size": 0,  # length of protein sequence used as reference
            "ref_coords": "",  # reference protein matching intervals
            "ref_starts": [],  # 'prot_start' per accepted exon matching
            "ref_ends": [],  # 'prot_end' per accepted exon matching
            "hit_ids": [],  # NUC, PTD, or MIT followed by Scipio's hit 'ID'(s)
            "hit_contigs": [],  # accepted 'target' names (= contig names in assembly)
            "hit_coords": "",  # matching intervals per 'target'
            "hit_starts": [],  # 'dna_start' for each accepted exon
            "hit_ends": [],  # 'dna_end' for each accepted exon
            "mat_types": [],  # 'type' of each 'matching' across 'target'(s)
            "mat_notes": [],  # info about the 'matching', e.g.: exon followed by gap
            "mat_nt": [],  # 'seq' in nucleotides for each accepted 'matching'
            "mat_aa": [],  # 'translation' in aminoacids for each accepted exon
            "strand": [],  # 'strand' matched per 'target'
            "matches": [],  # number of aminoacid matches per 'target'
            "mismatches": [],  # number of aminoacid mismatches per 'target'
            "coverage": 0.0,  # (matches + mismatches) / ref_size * 100
            "identity": 0.0,  # matches / (matches + mismatches) * 100
            "score": 0.0,  # (matches - mismatches) / ref_size
            "wscore": 0.0,  # score * (length AA / max length AA across refs)
            "gapped": False,  # recovered protein has gaps with respect to the reference
            "seq_flanked": "",  # concatenation of 'mat_nt'
            "seq_gene": "",  # 'seq_flanked' without 'upstream' or 'downstream' nucleotides
            "seq_nt": "",  # concatenation of 'mat_nt' only for exons
            "seq_aa": "",  # concatenation of 'mat_aa'
            "match_len": 0,  # number of reference aminoacids matched
        }
        add_separator = False

        for ctg in yaml_mod:  # ctg = contig (gene models can span several contigs)
            if ctg["prot_start"] < ctg["prot_end"]:  # only accept protein stretches longer than 0 bp
                ref_starts = []
                ref_ends = []
                hit_ids = []
                hit_starts = []
                hit_ends = []
                mat_types = []
                mat_notes = []
                mat_nt = []
                mat_aa = []
                mismatches = []
                has_overlap = False
                hit_id = f"{marker_type}{ctg['ID']}"
                # Sometimes target (contig) names can have spaces and be in between '', fix those
                if " " in ctg["target"]:
                    ctg["target"] = ctg["target"].strip("'").split()[0]
                if add_separator:
                    ref_starts.append(None)
                    ref_ends.append(None)
                    hit_ids.append(None)
                    hit_starts.append(None)
                    hit_ends.append(None)
                    mat_types.append("separator")
                    mat_notes.append("separator")
                    mat_nt.append(settings.SCIPIO_CONTIG_SEP)
                    mat_aa.append(None)
                if ctg["upstream"]:
                    ref_starts.append(None)
                    ref_ends.append(None)
                    hit_ids.append(hit_id)
                    hit_starts.append(None)
                    hit_ends.append(None)
                    mat_types.append("upstream")
                    mat_notes.append("upstream")
                    mat_nt.append(ctg["upstream"])
                    mat_aa.append(None)
                exons_in_ctg = sum(mat["type"] == "exon" for mat in ctg["matchings"])
                current_exon = 0
                for mat in ctg["matchings"]:
                    if mat["type"] == "exon":
                        mat["seq"], mat["translation"] = insert_shifts(mat)
                        current_exon += 1
                        note = ""
                        if current_exon == 1 and "gap to querystart" in ctg["reason"]:
                            note += "gap before/"
                        if current_exon == 1 and "gap to previous hit" in ctg["reason"]:
                            note += "gap before/"
                        if current_exon == exons_in_ctg and "gap to next hit" in ctg["reason"]:
                            note += "gap after/"
                        if current_exon == exons_in_ctg and "gap to queryend" in ctg["reason"]:
                            note += "gap after/"
                        if mat["seq"]:
                            ref_starts.append(mat["prot_start"])
                            ref_ends.append(mat["prot_end"])
                            hit_ids.append(hit_id)
                            hit_starts.append(mat["dna_start"])
                            hit_ends.append(mat["dna_end"])
                            mat_types.append("exon")
                            mat_notes.append(note)
                            mat_nt.append(mat["seq"])
                            mat_aa.append(mat["translation"])
                            mismatches += mat["mismatchlist"]
                            if mat["overlap"]:
                                has_overlap = True
                    else:
                        ref_starts.append(None)
                        ref_ends.append(None)
                        hit_ids.append(hit_id)
                        hit_starts.append(None)
                        hit_ends.append(None)
                        mat_types.append(mat["type"])
                        mat_notes.append(mat["type"])
                        mat_nt.append(mat["seq"])
                        mat_aa.append(None)
                if ctg["stopcodon"]:
                    ref_starts.append(None)
                    ref_ends.append(None)
                    hit_ids.append(hit_id)
                    hit_starts.append(None)
                    hit_ends.append(None)
                    mat_types.append("stopcodon")
                    mat_notes.append("stopcodon")
                    mat_nt.append(ctg["stopcodon"])
                    mat_aa.append(None)
                if ctg["downstream"]:
                    ref_starts.append(None)
                    ref_ends.append(None)
                    hit_ids.append(hit_id)
                    hit_starts.append(None)
                    hit_ends.append(None)
                    mat_types.append("downstream")
                    mat_notes.append("downstream")
                    mat_nt.append(ctg["downstream"])
                    mat_aa.append(None)
                if "exon" in mat_types:
                    mod["ref_name"] = protein_name
                    mod["ref_seqs"][hit_id] = (ctg["prot_seq"], ctg["prot_start"], ctg["prot_end"])
                    mod["ref_size"] = ctg["prot_len"]
                    mod["ref_starts"] += ref_starts
                    mod["ref_ends"] += ref_ends
                    mod["hit_ids"] += hit_ids
                    mod["hit_contigs"] += [ctg["target"]]
                    mod["hit_starts"] += hit_starts
                    mod["hit_ends"] += hit_ends
                    mod["mat_types"] += mat_types
                    mod["mat_notes"] += mat_notes
                    mod["mat_nt"] += mat_nt
                    mod["mat_aa"] += mat_aa
                    mod["strand"] += [ctg["strand"]]
                    mod["matches"] += [ctg["matches"]]
                    mod["mismatches"] += mismatches
                    add_separator = False if has_overlap else True

        mod = fix_model(mod, gencode, marker_type)
        if not mod:
            return None

        mod = check_gaps_and_concat_seqs(mod, gencode, predict, min_identity)
        prot_len_matched = calc_prot_len_matched(mod)
        mod = merge_adjoining_exons(mod)
        mod["ref_coords"] = concat_coords(mod["hit_ids"], mod["ref_starts"], mod["ref_ends"])
        mod["hit_coords"] = concat_coords(mod["hit_ids"], mod["hit_starts"], mod["hit_ends"])
        mod["hit_ids"] = "\n".join(list(dict.fromkeys(list(filter(None, mod["hit_ids"])))))
        mod["hit_contigs"] = "\n".join(mod["hit_contigs"])
        mod["strand"] = "\n".join(mod["strand"])
        mismatches = len(set(mod["mismatches"]))
        matches = prot_len_matched - mismatches
        mod["coverage"] = prot_len_matched / mod["ref_size"] * 100
        mod["identity"] = matches / prot_len_matched * 100
        mod["score"] = (matches - mismatches) / mod["ref_size"]
        mod["gapped"] = bool("gap" in "".join(mod["mat_notes"]))
        mod["match_len"] = min(prot_len_matched, mod["ref_size"])

        return mod

    gencode = genetic_code(transtable)
    yaml = load_scipio_yaml(yaml_path)
    unfilter_models = {}
    for prot in yaml:  # prot = protein (reference protein)
        for yaml_model in yaml[prot]:  # yaml_mod = gene model (hit, or paralog)
            model = parse_model(
                yaml[prot][yaml_model], prot, marker_type, gencode, predict, min_identity
            )
            if model:
                if prot in unfilter_models:
                    unfilter_models[prot][yaml_model] = model
                else:
                    unfilter_models[prot] = {yaml_model: model}

    # Separate reference protein names formatted like in the Angiosperms353.FAA file to get
    # the name of the protein cluster only when EVERY reference protein has the
    # 'settings.REFERENCE_CLUSTER_SEPARATOR'
    refs_have_separators = True
    for prot in unfilter_models:
        if settings.REF_CLUSTER_SEP not in prot:
            refs_have_separators = False
            break

    # Keep track of longest recovered AA length per locus
    max_len_aa_recov = {}
    for prot in unfilter_models:
        for model in unfilter_models[prot]:
            if refs_have_separators:
                ref_cluster = prot.split(settings.REF_CLUSTER_SEP)[-1]
            else:
                ref_cluster = prot
            if ref_cluster in max_len_aa_recov:
                if unfilter_models[prot][model]["match_len"] > max_len_aa_recov[ref_cluster]:
                    max_len_aa_recov[ref_cluster] = unfilter_models[prot][model]["match_len"]
            else:
                max_len_aa_recov[ref_cluster] = unfilter_models[prot][model]["match_len"]

    # Filter models by 'min_score', 'min_identity', and 'min_coverage', calculate 'wscore'
    # penalize 'wscore' by number of frameshifts, remove paralogs according to 'paralog_tolerance'
    filter_models = {}
    for prot in unfilter_models:
        accepted_models = []
        for model in unfilter_models[prot]:
            if refs_have_separators:
                ref_cluster = prot.split(settings.REF_CLUSTER_SEP)[-1]
            else:
                ref_cluster = prot
            frameshifts = len(
                set(
                    [
                        math.ceil((p + 1) / 3)
                        for p in range(len(unfilter_models[prot][model]["seq_nt"]))
                        if unfilter_models[prot][model]["seq_nt"][p] == "N"
                    ]
                )
            )
            contigs = len(set(unfilter_models[prot][model]["hit_contigs"].split("\n")))
            unfilter_models[prot][model]["wscore"] = (
                (
                    unfilter_models[prot][model]["score"]
                    * (unfilter_models[prot][model]["match_len"] / max_len_aa_recov[ref_cluster])
                )
                * (settings.SCIPIO_FRAMESHIFT_PENALTY**frameshifts)
                * (settings.EXTRA_CONTIG_PENALTY ** (contigs - 1))
            )
            if (
                unfilter_models[prot][model]["score"] >= min_score
                and unfilter_models[prot][model]["identity"] >= min_identity
                and unfilter_models[prot][model]["coverage"] >= min_coverage
            ):
                accepted_models.append(unfilter_models[prot][model])
        if accepted_models:
            accepted_models = sorted(accepted_models, key=lambda i: i["wscore"], reverse=True)
            max_wscore = accepted_models[0]["wscore"]
            accepted_models = [
                model for model in accepted_models if model["wscore"] >= max_wscore / paralog_tolerance
            ]
            if max_paralogs > -1:
                accepted_models = accepted_models[: max_paralogs + 1]
            filter_models[prot] = accepted_models
    unfilter_models = None

    # Keep only the best hit 'models[protein][0]' with highest score and its paralogs for each
    # protein cluster. Check 'settings_assembly.py' for a more detailed description of how to
    # format your protein reference files under 'REFERENCE_CLUSTER_SEPARATOR'
    if filter_models:
        best_models = {}
        for prot in filter_models:
            if refs_have_separators:
                ref_cluster = prot.split(settings.REF_CLUSTER_SEP)[-1]
            else:
                ref_cluster = prot
            if ref_cluster not in best_models:
                best_models[ref_cluster] = filter_models[prot]
            else:
                if filter_models[prot][0]["wscore"] > best_models[ref_cluster][0]["wscore"]:
                    best_models[ref_cluster] = filter_models[prot]
                elif filter_models[prot][0]["wscore"] == best_models[ref_cluster][0]["wscore"]:
                    top_so_far = statistics.mean(
                        [
                            best_models[ref_cluster][i]["wscore"]
                            for i in range(len(best_models[ref_cluster]))
                        ]
                    )
                    tied = statistics.mean(
                        [filter_models[prot][i]["wscore"] for i in range(len(filter_models[prot]))]
                    )
                    if tied > top_so_far:
                        best_models[ref_cluster] = filter_models[prot]
        filter_models = None
        return best_models
    else:
        return None


def parse_psl_record(psl_line):
    record = psl_line.strip().split()
    try:
        psl = {
            "matches": int(record[0]),
            "mismatches": int(record[1]),
            "rep_matches": int(record[2]),
            "n_count": int(record[3]),
            "q_num_insert": int(record[4]),
            "q_base_insert": int(record[5]),
            "t_num_insert": int(record[6]),
            "t_base_insert": int(record[7]),
            "strand": record[8],
            "q_name": record[9],
            "q_size": int(record[10]),
            "q_start": int(record[11]),
            "q_end": int(record[12]),
            "t_name": record[13],
            "t_size": int(record[14]),
            "t_start": int(record[15]),
            "t_end": int(record[16]),
            "block_count": int(record[17]),
            "block_sizes": [int(x) for x in record[18].strip(",").split(",")],
            "q_starts": [int(x) for x in record[19].strip(",").split(",")],
            "t_starts": [int(x) for x in record[20].strip(",").split(",")],
        }
        if len(record) == 22:
            psl["wscore"] = float(record[21])
        return psl
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        print(f"Malformed PSL record: {psl_line.strip()}")
        return None


def calculate_psl_identity(
    matches: int,
    rep_matches: int,
    mismatches: int,
    q_num_insert: int,
    t_num_insert: int,
    strand: str,
    q_start: int,
    q_end: int,
    t_size: int,
    t_start: int,
    t_end: int,
    block_sizes: list,
    t_starts: list,
    q_type="unknown",
    q_is_mrna=False,
):
    """
    Adapted from:
    https://genome-source.gi.ucsc.edu/gitlist/kent.git/raw/master/src/utils/pslScore/pslScore.pl
    Calculates pct_identity from blat record in .psl
    query and target refer to BLAT nomenclature
    BLAT query = Captus target loci
    BLAT target = Captus assembled contigs

    Args:
        matches (int): number of matches that are not repeats
        rep_matches (int): number of matches that are part of repeats
        mismatches (int): number of mismatches
        q_num_insert (int): number of inserts in query
        t_num_insert (int): number of inserts in target
        strand (str): query strand, a second + or - refers to the target strand when query is protein
        q_start (int): alignment start position in query
        q_end (int): alignment end position in query
        t_size (int): target sequence size
        t_start (int): alignment start position in target
        t_end (int): alignment end position in target
        block_sizes (list): list of integers corresponding to each block size
        t_starts (list): list of integers corresponding to start positions of each block in target
        q_type (str, optional): dna or protein, will determine if unknown. Defaults to "unknown".
        q_is_mrna (bool, optional): if the query is mRNA. Defaults to False.
    """

    def psl_is_protein(strand, t_size, t_start, t_end, block_sizes, t_starts):
        """Returns a size multiplier of 1 if the PSL is DNA or 3 if the PSL is PROTEIN"""
        size_mul = 1
        if len(strand) > 1:
            direction = strand[-1]
            t_end_last_block = t_starts[-1] + (3 * block_sizes[-1])
            if direction == "+":
                if t_end == t_end_last_block:
                    size_mul = 3
            elif direction == "-":
                if t_start == (t_size - t_end_last_block):
                    size_mul = 3
        return size_mul

    if q_type == "protein":
        size_mul = 3
    elif q_type == "dna":
        size_mul = 1
    else:
        size_mul = psl_is_protein(strand, t_size, t_start, t_end, block_sizes, t_starts)

    milli_bad = 0
    q_ali_size = size_mul * (q_end - q_start)
    t_ali_size = t_end - t_start
    ali_size = q_ali_size
    if t_ali_size < q_ali_size:
        ali_size = t_ali_size
    if ali_size <= 0:
        return milli_bad
    size_dif = q_ali_size - t_ali_size
    if size_dif < 0:
        if q_is_mrna is True:
            size_dif = 0
        else:
            size_dif = -size_dif
    insert_factor = q_num_insert
    # Let's not penalize insertions in the contigs by disabling the next two lines
    # if is_mrna is False:
    #     insert_factor += t_num_insert
    total = size_mul * (matches + rep_matches + mismatches)
    if total != 0:
        round_away_from_zero = 3 * math.log(1 + size_dif)
        if round_away_from_zero < 0:
            round_away_from_zero = int(round_away_from_zero - 0.5)
        else:
            round_away_from_zero = int(round_away_from_zero + 0.5)
        milli_bad = (1000 * (mismatches * size_mul + insert_factor + round_away_from_zero)) / total

    return 100.0 - milli_bad * 0.10


def blat_misc_dna_psl_to_dict(
    psl_path,
    target_dict,
    min_identity,
    min_coverage,
    marker_type,
    disable_stitching,
    max_paralogs,
    paralog_tolerance,
):
    """
    Parse .psl from BLAT, assemble greedily the partial hits, and return the best set of hits if
    the reference contains more than a single sequence of the same type (analogous to the reference
    file Angiosperms353.FAA)
    """

    def merge_blocks(q_size: int, strand: str, block_sizes: list, q_starts: list, t_starts: list):
        """
        BLAT split alignments in block intercalated by insertions, we merge adjacent blocks
        when an intervening insertion is not bigger than DNA_MAX_INSERT_SIZE or
        qSize * DNA_MAX_INSERT_PROP
        query and target refer to BLAT nomenclature
        BLAT query = Captus target loci
        BLAT target = Captus assembled contigs

        Args:
            q_size (int): size of query sequence
            strand (str): + or -, applies to query
            block_sizes (list): length in bp of each subaligment block
            q_starts (list): start coordinates of query
            t_starts (list): start coordinates of target

        Returns:
            q_starts_mb, q_ends_mb, t_starts_mb, t_ends_mb: list of starts and ends for query and
                                                            target, reordered for the query when
                                                            strand is negative, adjacent blocks
                                                            merged when the insertions sizes are
                                                            accepted (see above)
        """
        q_ends = [q_starts[i] + block_sizes[i] for i in range(len(block_sizes))]
        t_ends = [t_starts[i] + block_sizes[i] for i in range(len(block_sizes))]
        q_starts_mb, t_starts_mb = [q_starts[0]], [t_starts[0]]
        q_ends_mb, t_ends_mb = [], []
        for i in range(len(block_sizes) - 1):
            if t_starts[i + 1] - t_ends[i] >= min(
                q_size * settings.DNA_MAX_INSERT_PROP, settings.DNA_MAX_INSERT_SIZE
            ):
                q_starts_mb.append(q_starts[i + 1])
                t_starts_mb.append(t_starts[i + 1])
                q_ends_mb.append(q_ends[i])
                t_ends_mb.append(t_ends[i])
        q_ends_mb.append(q_ends[-1])
        t_ends_mb.append(t_ends[-1])
        if strand == "-":
            q_starts_mb_neg = [q_size - x for x in q_ends_mb[::-1]]
            q_ends_mb_neg = [q_size - x for x in q_starts_mb[::-1]]
            return q_starts_mb_neg, q_ends_mb_neg, t_starts_mb, t_ends_mb
        else:
            return q_starts_mb, q_ends_mb, t_starts_mb, t_ends_mb

    def determine_matching_region(q_size, q_start, q_end, t_size, t_start, t_end, strand):
        """
        Determine if a contig matches entirely the query, or if it is a partial hit. In cases of
        partial hits determine if it is a split hit (proximal, middle, or distal) or if the hit
        is partial and subsumed within a larger stretch of sequence unrelated to the query
        """
        region = ""
        if q_end - q_start >= q_size * (1 - settings.DNA_TOLERANCE_LEN):
            region = "full"
        elif (
            q_start <= q_size * settings.DNA_TOLERANCE_LEN
            and (q_size - q_end) <= q_size * settings.DNA_TOLERANCE_LEN
        ):
            region = "full"
        elif (
            q_start >= q_size * settings.DNA_TOLERANCE_LEN
            and (q_size - q_end) >= q_size * settings.DNA_TOLERANCE_LEN
        ):
            region = "middle"
        elif q_start <= q_size * settings.DNA_TOLERANCE_LEN:
            region = "proximal"
        elif (q_size - q_end) <= q_size * settings.DNA_TOLERANCE_LEN:
            region = "distal"
        else:  # hit is only partial and surrounded by a large proportion of unmatched sequence
            region = "wedged"
        # Now check if the flanks in the contig are too long to be part of a multi-part hit
        left_flank_too_long = t_start / (t_end - t_start) < settings.DNA_TOLERANCE_LEN
        right_flank_too_long = (t_size - t_end) / (t_end - t_start) < settings.DNA_TOLERANCE_LEN
        if region == "middle":
            if left_flank_too_long or right_flank_too_long:
                region = "wedged"
        elif (region == "proximal" and strand == "+") or (region == "distal" and strand == "-"):
            if right_flank_too_long:
                region = "wedged"
        elif (region == "proximal" and strand == "-") or (region == "distal" and strand == "+"):
            if left_flank_too_long:
                region = "wedged"
        return region

    def determine_max_overlap(contig_name):
        template = re.compile(CAPTUS_MEGAHIT_HEADER_REGEX)
        if template.match(contig_name):
            max_overlap_bp = int(contig_name.split("_")[7])
        else:
            max_overlap_bp = settings.DNA_MAX_OVERLAP_BP
        return max_overlap_bp

    def greedy_assembly_partial_hits(hits_list, max_overlap_bp, max_paralogs: int):
        """
        Partial hits are assembled greedily starting by the hits with the highest 'wscore',
        producing as many paralogs as necessary. Paralogs are merged from partial hits when all
        the hits within the paralog are compatible (they have overlap and identity within tolerated
        limits)
        """
        full_hits, part_hits = [], []

        for hit in hits_list:
            if hit["region"] == "full":
                full_hits.append(hit)
            else:
                part_hits.append(hit)

        paths = []
        if part_hits:
            # Sort partial hits by their 'wscore'
            part_hits = sorted(part_hits, key=lambda i: i["wscore"], reverse=True)
            # Make an empty list as long as the number of hits, in the worst-case scenario no hit is
            # compatible with any other
            paths = [[] for h in part_hits]
            # Fill the first path with the first hit in 'part_hits'
            paths[0] = [part_hits[0]]

        # Check the remaining hits, only accept hits in a path when they are compatible with all the
        # other hits in that path
        if len(part_hits) > 1:
            for h2 in part_hits[1:]:
                for path in paths:
                    compatible = True
                    if path:
                        for h1 in path:
                            if not pair_is_compatible(h1, h2, max_overlap_bp):
                                compatible = False
                                break
                    if compatible:
                        path.append(h2)
                        break

        # Sort each path by starting query coordinate
        sorted_paths = []
        for path in paths:
            if path:
                if len(path) > 1:
                    sorted_paths.append(sorted(path, key=lambda i: i["q_start"][0]))
                else:
                    sorted_paths.append(path)
        del paths

        assembly_paths = [[full_hits[f]] for f in range(len(full_hits))] + sorted_paths

        # Search FASTA assembly input ('target_dict') and extract/stitch needed sequence fragments
        assembly = extract_and_stitch_edges(assembly_paths, max_paralogs)
        assembly_paths = None
        return assembly

    def pair_is_compatible(h1, h2, max_overlap_bp):
        if abs(h1["identity"] - h2["identity"]) > (settings.DNA_TOLERANCE_PID * 100.0):
            return False
        overlap = min(h1["q_end"][-1], h2["q_end"][-1]) - max(h1["q_start"][0], h2["q_start"][0])
        if overlap < 1:
            return True
        if (
            h1["match_len"] * (1 - settings.DNA_TOLERANCE_LEN) <= overlap
            or h2["match_len"] * (1 - settings.DNA_TOLERANCE_LEN) <= overlap
        ):
            return False
        return bool(overlap <= max_overlap_bp)

    def extract_and_stitch_edges(assembly_paths, max_paralogs: int):
        """
        Returns a list of stitched sequences with metadata, sorted by relevance: 'full' hits first
        sorted by 'wscore', followed by assembled hits sorted by 'wscore', and finally
        partial unassembled hits sorted by 'wscore'.
        Assembly: For overlaps follow the coordinate system of the query and remove the overlap from
        the partial hit with the lower 'identity', for non-overlapped partial hits, concatenate hits
        intercalating them with as many Ns as indicated by gap in query
        """
        raw_assembly = []
        for path in assembly_paths:
            len_path = len(path)
            if len_path > 1:
                for i in range(1, len_path):
                    if path[i]["region"] == "proximal":
                        path[i]["region"] = "middle"
                for i in range(len_path - 1):
                    if path[i]["region"] == "distal":
                        path[i]["region"] = "middle"
            asm_hit = {
                "ref_name": path[0]["ref_name"],  # full name of non-coding reference sequence
                "ref_size": path[0]["ref_size"],  # length of non-coding reference sequence
                "ref_coords": join_coords(path[0]["q_start"], path[0]["q_end"]),  # coords in ref
                "hit_ids": path[0]["hit_id"],  # 'DNA##' or 'CLR##'
                "hit_contigs": path[0]["hit_contig"],  # contig name(s) used in assembly
                "hit_coords": join_coords(path[0]["t_start"], path[0]["t_end"]),  # coords in match
                "strand": path[0]["strand"],  # contigs strand(s)
                "matches": path[0]["matches"],  # accumulated matches across targets
                "mismatches": path[0]["mismatches"],  # accumulated mismatches across targets
                "coverage": path[0]["coverage"],  # ((matches + mismatches) / ref_size) * 100
                "identity": path[0]["identity"],  # (matches / (matches + mismatches)) * 100
                "score": path[0]["score"],  # Scipio-like score as (matches - mismatches) / ref_size
                "wscore": path[0]["wscore"],  # Scipio-like * (len matched / locus max len matched)
                "gapped": path[0]["gapped"],  # set to True when is assembly of partial hits
                "region": path[0]["region"],  # full, proximal, middle, distal with respect to ref
                # assembled sequence match plus upstream and downstream buffer
                "seq_flanked": extract_psl_sequence(
                    target_dict, path[0], settings.DNA_UP_DOWN_STREAM_BP, flanked=True
                ),
                # assembled sequence match
                "seq_gene": extract_psl_sequence(target_dict, path[0], 0),
                "seq_nt": "",  # not used
                "seq_aa": "",  # not used
                "match_len": path[0]["q_end"][-1] - path[0]["q_start"][0],  # Reference length matched
            }

            if len(path) > 1:
                ref_starts = [list(path[0]["q_start"])]
                ref_ends = [list(path[0]["q_end"])]
                hit_ids = [path[0]["identity"]]
                sum_matches = path[0]["matches"]
                sum_mismatches = path[0]["mismatches"]
                for h in range(len(path) - 1):
                    asm_hit["hit_ids"] += f"\n{path[h + 1]['hit_id']}"
                    asm_hit["hit_contigs"] += f"\n{path[h + 1]['hit_contig']}"
                    asm_hit["hit_coords"] += (
                        f"\n{join_coords(path[h + 1]['t_start'], path[h + 1]['t_end'])}"
                    )
                    asm_hit["strand"] += f"\n{path[h + 1]['strand']}"
                    hit_ids.append(path[h + 1]["identity"])
                    asm_hit["region"] += f",{path[h + 1]['region']}"
                    next_seq_flanked = extract_psl_sequence(
                        target_dict, path[h + 1], settings.DNA_UP_DOWN_STREAM_BP, flanked=True
                    )
                    next_seq_gene = extract_psl_sequence(target_dict, path[h + 1], 0)
                    sum_matches += path[h + 1]["matches"]
                    sum_mismatches += path[h + 1]["mismatches"]

                    overlap = path[h]["q_end"][-1] - path[h + 1]["q_start"][0]
                    # Negative 'overlap' is a gap that has to be filled with 'n's
                    if overlap <= 0:
                        gap_len = abs(overlap) if overlap < 0 else 0
                        asm_hit["seq_flanked"] = stitch_contigs(
                            asm_hit["seq_flanked"],
                            asm_hit["hit_contigs"].split("\n")[-1],
                            next_seq_flanked,
                            path[h + 1]["hit_contig"],
                            gap_len,
                        )
                        asm_hit["seq_gene"] += f"{'n' * gap_len}{next_seq_gene}"
                        ref_starts.append(list(path[h + 1]["q_start"]))
                        ref_ends.append(list(path[h + 1]["q_end"]))
                    else:
                        # Ignore overlapped portion from the hit with smaller 'identity' to the ref
                        if path[h]["identity"] >= path[h + 1]["identity"]:
                            asm_hit["seq_flanked"] += next_seq_flanked[overlap:]
                            asm_hit["seq_gene"] += next_seq_gene[overlap:]
                            ref_starts.append(list(path[h + 1]["q_start"]))
                            ref_starts[-1][0] += overlap
                            ref_ends.append(list(path[h + 1]["q_end"]))
                        else:
                            asm_hit["seq_flanked"] = (
                                f"{asm_hit['seq_flanked'][:-overlap]}{next_seq_flanked}"
                            )
                            asm_hit["seq_gene"] = f"{asm_hit['seq_gene'][:-overlap]}{next_seq_gene}"
                            ref_ends[-1][-1] -= overlap
                            ref_starts.append(list(path[h + 1]["q_start"]))
                            ref_ends.append(list(path[h + 1]["q_end"]))

                # Reformat the matched query coordinates for decorating the FASTA description line
                asm_hit["ref_coords"] = "\n".join(
                    [join_coords(s, e) for s, e in zip(ref_starts, ref_ends)]
                )
                # To avoid inflating the 'score' and 'coverage' of hits with insertions we need to
                # find out first the 'matched_len' discounting gaps and insertions to ref
                match_len = ref_ends[-1][-1] - ref_starts[0][0] - asm_hit["seq_gene"].count("n")
                asm_hit["coverage"] = match_len / asm_hit["ref_size"] * 100.0
                # Calculate the mean 'identity' of all the partial hits used in the assembled path
                asm_hit["identity"] = statistics.mean(hit_ids)
                # Recalculate the 'score' and 'wscore' using sum of matches/mismatches from all
                # partial hits used in the assemble path
                matches = (sum_matches * match_len) / (sum_matches + sum_mismatches)
                mismatches = (sum_mismatches * match_len) / (sum_matches + sum_mismatches)
                asm_hit["score"] = (matches - mismatches) / asm_hit["ref_size"]
                full_len = len(asm_hit["seq_gene"].replace("n", ""))
                asm_hit["wscore"] = asm_hit["score"] * (full_len / asm_hit["ref_size"])
                asm_hit["gapped"] = bool("n" in asm_hit["seq_gene"])
                asm_hit["match_len"] = match_len

            # Append hits to the global assembly 'raw_assembly'
            raw_assembly.append(dict(asm_hit))

        # Sorting hits by 'wscore'
        raw_assembly = sorted(raw_assembly, key=lambda i: i["wscore"], reverse=True)
        # Filter hit by 'min_coverage' and 'min_identity'
        assembly = []
        for hit in raw_assembly:
            if hit["identity"] >= min_identity and hit["coverage"] >= min_coverage:
                assembly.append(hit)
        # Filter by 'max_paralogs'
        if max_paralogs > -1:
            assembly = assembly[: max_paralogs + 1]
        raw_assembly = None
        return assembly

    def join_coords(starts, ends):
        """
        Given a list of starts [1,30,50] and a list of ends [12,38,100] produce a string for
        decorating the sequence description as '1-12,30-38,50-100', don't make them 1-based yet,
        that is done by the functions 'write_gff3' and 'write_fastas_and_report'
        """
        return ",".join([f"{s}-{e}" for s, e in zip(starts, ends)])

    def extract_psl_sequence(fasta_dict, hit_dict, up_down_stream_bp, flanked=False):
        """
        Extract sequence from a 'fasta_dict' object using BLAT's PSL coordinate style
        """
        contig = hit_dict["hit_contig"]
        strand = hit_dict["strand"]
        starts = list(hit_dict["t_start"])
        ends = list(hit_dict["t_end"])
        ## region = hit_dict["region"]
        sequence = ""
        if flanked:
            for i in range(len(starts)):
                # When we extract the flanked matches and there are intervening blocks of unmatched
                # sequence we must extract those as well, for example matching CDS against a genome will
                # match the exons, in the flanked version we recover the introns as well
                sequence += fasta_dict[contig]["sequence"][starts[i] : ends[i]].upper()
                try:
                    sequence += fasta_dict[contig]["sequence"][ends[i] : starts[i + 1]].lower()
                except IndexError:
                    continue
            # Now we extract the flanks, being aware of strand
            seq_len = len(fasta_dict[contig]["sequence"])
            start_flank = max((starts[0] - up_down_stream_bp), 0)
            end_flank = min((ends[-1] + up_down_stream_bp), seq_len)
            left_flank = fasta_dict[contig]["sequence"][start_flank : starts[0]].lower()
            right_flank = fasta_dict[contig]["sequence"][ends[-1] : end_flank].lower()
            sequence = f"{left_flank}{sequence}{right_flank}"
            # if region == "full":
            #     sequence = f"{left_flank}{sequence}{right_flank}"
            # elif ((region  == "proximal" and strand == "+")
            #       or (region == "distal" and strand == "-")):
            #     sequence = f"{left_flank}{sequence}"
            # elif ((region  == "proximal" and strand == "-")
            #       or (region == "distal" and strand == "+")):
            #     sequence = f"{sequence}{right_flank}"
        else:
            for s, e in zip(starts, ends):
                sequence += fasta_dict[contig]["sequence"][s:e].upper()
        if strand == "+":
            return sequence
        elif strand == "-":
            return reverse_complement(sequence)

    def stitch_contigs(current_ctg_seq, current_ctg_name, next_ctg_seq, next_ctg_name, gap_len):
        current_ctg_cov = 0
        if "_cov_" in current_ctg_name:
            current_ctg_cov = float(current_ctg_name.split("_cov_")[-1].split("_")[0])
        current_trail = ""
        for i in reversed(range(len(current_ctg_seq))):
            if current_ctg_seq[i].isupper():
                current_trail = current_ctg_seq[i + 1 :]
                break

        next_ctg_cov = 0
        if "_cov_" in next_ctg_name:
            next_ctg_cov = float(next_ctg_name.split("_cov_")[-1].split("_")[0])
        next_lead = ""
        for i in range(len(next_ctg_seq)):
            if next_ctg_seq[i].isupper():
                next_lead = next_ctg_seq[:i]

        check_len = min(len(current_trail), len(next_lead))
        scaffold_overlap = 0
        # Only check for potential overlaps > DNA_MIN_OVERLAP_BP (21 bp)
        if check_len > settings.DNA_MIN_OVERLAP_BP:
            while check_len >= settings.DNA_MIN_OVERLAP_BP:
                if overlap_matches(current_trail[-check_len:], next_lead[0:check_len]):
                    scaffold_overlap = check_len
                    break
                check_len -= 1
        scaffold = ""
        sep_len = gap_len - (len(current_trail) + len(next_lead))
        separator = "n" * sep_len if sep_len > 0 else settings.DNA_CONTIG_SEPARATOR

        if scaffold_overlap == 0:
            scaffold = f"{current_ctg_seq}{separator}{next_ctg_seq}"
        else:
            if current_ctg_cov >= next_ctg_cov:
                scaffold = f"{current_ctg_seq}{next_ctg_seq[scaffold_overlap:]}"
            else:
                scaffold = f"{current_ctg_seq[:-scaffold_overlap]}{next_ctg_seq}"

        return scaffold

    def overlap_matches(seq1, seq2):
        matches = 0
        seq_len = len(seq1)
        for n1, n2 in zip(seq1.upper(), seq2.upper()):
            matches += NT_PIDS["".join(sorted([n1, n2]))]
        if matches == seq_len:
            return True
        if matches / seq_len >= (1 - settings.DNA_TOLERANCE_PID):
            return True
        if seq_len - matches <= 1.0 and seq_len <= 10:
            return True
        return False

    raw_dna_hits = {}
    hit = {}
    hit_num = 1
    with open(psl_path) as psl_in:
        for line in psl_in:
            p = parse_psl_record(line)
            if p is None:
                continue
            q_starts_mb, q_ends_mb, t_starts_mb, t_ends_mb = merge_blocks(
                p["q_size"],
                p["strand"],
                p["block_sizes"],
                p["q_starts"],
                p["t_starts"],
            )
            coverage = (p["matches"] + p["rep_matches"] + p["mismatches"]) / p["q_size"]
            score = (p["matches"] + p["rep_matches"] - p["mismatches"]) / p["q_size"]
            wscore = score * coverage
            pct_coverage = coverage * 100
            pct_identity = calculate_psl_identity(
                p["matches"],
                p["rep_matches"],
                p["mismatches"],
                p["q_num_insert"],
                p["t_num_insert"],
                p["strand"],
                q_starts_mb[0],
                q_ends_mb[-1],
                p["t_size"],
                t_starts_mb[0],
                t_ends_mb[-1],
                p["block_sizes"],
                p["t_starts"],
                q_type="dna",
            )

            if disable_stitching:  # prevent locus assembly across multiple contigs
                region = "full"
            else:
                region = determine_matching_region(
                    p["q_size"],
                    q_starts_mb[0],
                    q_ends_mb[-1],
                    p["t_size"],
                    t_starts_mb[0],
                    t_ends_mb[-1],
                    p["strand"],
                )
            gapped = not bool(region == "full")

            # Ignore hits with not enough pct_coverage of the reference, or partial hits immersed in
            # larger segments of unrelated sequence
            if not (pct_coverage < settings.DNA_MIN_COVERAGE_BEFORE_ASSEMBLY and region == "wedged"):
                # Compose hit record with columns from the .psl line
                hit = {
                    "ref_name": p["q_name"],
                    "ref_size": p["q_size"],
                    "q_start": q_starts_mb,
                    "q_end": q_ends_mb,
                    "match_len": p["matches"] + p["rep_matches"] + p["mismatches"],
                    "hit_id": f"{marker_type}{hit_num}",
                    "hit_contig": p["t_name"],
                    "t_start": t_starts_mb,
                    "t_end": t_ends_mb,
                    "strand": p["strand"],
                    "matches": p["matches"] + p["rep_matches"],
                    "mismatches": p["mismatches"],
                    "coverage": pct_coverage,
                    "identity": pct_identity,
                    "score": score,
                    "wscore": wscore,
                    "region": region,
                    "gapped": gapped,
                }
                hit_num += 1
                if p["q_name"] not in raw_dna_hits:
                    raw_dna_hits[p["q_name"]] = [dict(hit)]
                else:
                    raw_dna_hits[p["q_name"]].append(dict(hit))

    # Use the name of the last contig to extract kmer size if the assembly was done within Captus
    # and determine the maximum overlap tolerated between adjacent contigs for assembly
    if hit:
        max_overlap_bp = determine_max_overlap(hit["hit_contig"])
    else:
        return False

    # Assemble hits, organized by reference
    dna_hits = {}
    for dna_ref in raw_dna_hits:
        assembled_hits = greedy_assembly_partial_hits(
            raw_dna_hits[dna_ref], max_overlap_bp, max_paralogs
        )
        if assembled_hits:
            dna_hits[dna_ref] = assembled_hits

    # Find if reference has been formatted like Angiosperms353.FAA in order to accomodate more than
    # a single reference of the same type, we can also use 'settings.REFERENCE_CLUSTER_SEPARATOR' to
    # recognize the name of the cluster of references of the same kind
    refs_have_separators = True
    for dna_ref in dna_hits:
        if settings.REF_CLUSTER_SEP not in dna_ref:
            refs_have_separators = False
            break

    if dna_hits:
        # First obtain the length of the longest match across reference sequences in reference locus
        # 'match_len' is holding temporarily the matched length of each hit
        max_len_nt_recov = {}
        for dna_ref in dna_hits:
            if refs_have_separators:
                ref_cluster = dna_ref.split(settings.REF_CLUSTER_SEP)[-1]
            else:
                ref_cluster = dna_ref
            for hit in dna_hits[dna_ref]:
                if ref_cluster not in max_len_nt_recov:
                    max_len_nt_recov[ref_cluster] = hit["match_len"]
                else:
                    if hit["match_len"] > max_len_nt_recov[ref_cluster]:
                        max_len_nt_recov[ref_cluster] = hit["match_len"]

        # Loop the object again to calculate the actual wscore
        for dna_ref in dna_hits:
            if refs_have_separators:
                ref_cluster = dna_ref.split(settings.REF_CLUSTER_SEP)[-1]
            else:
                ref_cluster = dna_ref
            for hit in dna_hits[dna_ref]:
                contigs = len(set(hit["hit_contigs"].split("\n")))
                hit["wscore"] = (
                    hit["score"]
                    * (hit["match_len"] / max_len_nt_recov[ref_cluster])
                    * (settings.EXTRA_CONTIG_PENALTY ** (contigs - 1))
                )
            # Sort hits from largest to smallest 'wscore'
            dna_hits[dna_ref] = sorted(dna_hits[dna_ref], key=lambda i: i["wscore"], reverse=True)
            # Filter by 'paralog_tolerance'
            max_wscore = dna_hits[dna_ref][0]["wscore"]
            dna_hits[dna_ref] = [
                hit for hit in dna_hits[dna_ref] if hit["wscore"] >= max_wscore / paralog_tolerance
            ]

        # If multiple references of the same kind exist in the reference, then choose the one with
        # the best hit that has the highest 'wscore', break ties by averaging the 'wscore' of all
        # secondary hits
        best_dna_hits = {}
        for dna_ref in dna_hits:
            if refs_have_separators:
                ref_cluster = dna_ref.split(settings.REF_CLUSTER_SEP)[-1]
            else:
                ref_cluster = dna_ref
            if ref_cluster not in best_dna_hits:
                best_dna_hits[ref_cluster] = dna_hits[dna_ref]
            else:
                if dna_hits[dna_ref][0]["wscore"] > best_dna_hits[ref_cluster][0]["wscore"]:
                    best_dna_hits[ref_cluster] = dna_hits[dna_ref]
                elif dna_hits[dna_ref][0]["wscore"] == best_dna_hits[ref_cluster][0]["wscore"]:
                    top_so_far = statistics.mean(
                        [
                            best_dna_hits[ref_cluster][i]["wscore"]
                            for i in range(len(best_dna_hits[ref_cluster]))
                        ]
                    )
                    tied = statistics.mean(
                        [dna_hits[dna_ref][i]["wscore"] for i in range(len(dna_hits[dna_ref]))]
                    )
                    if tied > top_so_far:
                        best_dna_hits[ref_cluster] = dna_hits[dna_ref]
        dna_hits = None
        return best_dna_hits
    else:
        return False


def taper_correction(
    fasta_in_path: Path,
    fasta_out_path: Path,
    ambig: str,
    cutoff_threshold: float,
    conservative: bool,
):
    """TAPER algorithm published in https://doi.org/10.1111/2041-210X.13696
       The code translation from Julia to Python was assisted by Anthropic's AI assistant "Claude"
       on May 2025 and then adapted to fit in the trimming stage from the 'align' module of Captus

    Args:
        fasta_in_path (Path): Path to input FASTA file
        fasta_out_path (Path): Path to output FASTA file
        mask (str): Masking character to indicate where the modifications by TAPER
        ambig (str): Ambiguity character, N for nucleotides, X for aminoacids
        cutoff_threshold (float): 3.0 as recommended in the TAPER, option available in command line
        tpars (list of dict): k,p,q,L stored in settings.py, option not available in command line
        conservative (bool): enable conservative mode for correction, option available in command line
    """

    def correct_sequences(
        sequences: list,
        k: int,
        ambig: str,
        mask: str,
        pvalue: float,
        qvalue: float,
        cutoff_threshold: float,
        conservative: bool,
    ):
        ambig_upper = ambig.upper()
        seqs_upper = [seq.upper() for seq in sequences]
        seqs_array = [list(seq) for seq in sequences]
        n = len(sequences[0])
        m = len(sequences)

        w0 = [[None for _ in range(m)] for _ in range(n)]
        for i in range(n):
            rcounts = {}
            for j in range(m):
                r = seqs_upper[j][i]
                if r not in ["-", ambig_upper]:
                    if r in rcounts:
                        rcounts[r] += 1
                    else:
                        rcounts[r] = 1
            unq = len(rcounts)
            tot = sum(rcounts.values())
            for j in range(m):
                if tot > 0:
                    if seqs_array[j][i] not in ["-", ambig]:
                        w0[i][j] = tot / (unq * rcounts[seqs_upper[j][i]])

        w1 = [[w0[i][j] for i in range(n) if w0[i][j] is not None] for j in range(m)]

        w = []
        ws = []
        if conservative:
            for arr in w1:
                if len(arr) >= k:
                    w.append([statistics.median(arr[i : i + k]) for i in range(len(arr) - k + 1)])
                    ws.append([sum(arr[i : i + k]) for i in range(len(arr) - k + 1)])
                else:
                    w.append([])
                    ws.append([])
        else:
            # The index P-1 is absolutely necessary!, in Julia it was just P
            P = int(round(k * 2 / 3)) - 1
            for arr in w1:
                if len(arr) >= k:
                    w.append([sorted(arr[i : i + k])[P] for i in range(len(arr) - k + 1)])
                    ws.append([sum(arr[i : i + k]) for i in range(len(arr) - k + 1)])
                else:
                    w.append([])
                    ws.append([])

        wsorted = [sorted(arr) for arr in w]
        wsum = []
        for arr in wsorted:
            accu = []
            s = 0
            for x in arr:
                s += x
                accu.append(s)
            wsum.append(accu)

        def f(x, y, n, m):
            return x**2 / n + (0 if n == m else (y - x) ** 2 / (m - n))

        var = []
        for i, arr in enumerate(wsum):
            if len(arr) > 0:
                var.append([f(arr[j], arr[-1], j + 1, len(arr)) for j in range(len(arr))])
            else:
                var.append([])
        wCutoff = []
        for j in range(m):
            if len(var[j]) > 0:
                max_index = var[j].index(max(var[j]))
                wCutoff.append(wsorted[j][max_index])
            else:
                wCutoff.append(0)
        cutoffSorted = sorted([wCutoff[j] for j in range(m) if len(var[j]) > 0])
        cutoffFloor = 0 if pvalue >= 1 else cutoffSorted[-int(len(cutoffSorted) * pvalue) - 1]

        sequences_corrected = []
        for j in range(m):
            wj = w[j]
            wsj = ws[j]
            L = len(wj)
            if L <= k:
                sequences_corrected.append(sequences[j])
                continue

            s = [[0, 0] for _ in range(L)]
            tiebreaker = [[0, 0] for _ in range(L)]
            bt = [[0, 0] for _ in range(L)]

            cutoff = max(
                wCutoff[j],
                cutoffFloor,
                0 if qvalue >= 1 else wsorted[j][-int(len(wsorted[j]) * qvalue) - 1],
                cutoff_threshold,
            )

            for i in range(L):
                v = 0 if wj[i] > cutoff else 1
                if i == 0:
                    s[i][0] = v
                    s[i][1] = 1 - v
                    tiebreaker[i][0] = 0
                    tiebreaker[i][1] = wsj[i]
                else:
                    s[i][0] = s[i - 1][0] + v
                    s[i][1] = s[i - 1][1] + 1 - v
                    tiebreaker[i][0] = tiebreaker[i - 1][0]
                    tiebreaker[i][1] = tiebreaker[i - 1][1] + wsj[i]
                    bt[i][0] = 1
                    bt[i][1] = 2

                if i >= k and (s[i][0], tiebreaker[i][0]) < (s[i - k][1] + v, tiebreaker[i - k][1]):
                    s[i][0] = s[i - k][1] + v
                    tiebreaker[i][0] = tiebreaker[i - k][1]
                    bt[i][0] = 2

                if i >= k and (s[i][1], tiebreaker[i][1]) < (
                    s[i - k][0] + 1 - v,
                    tiebreaker[i - k][0] + wsj[i],
                ):
                    s[i][1] = s[i - k][0] + 1 - v
                    tiebreaker[i][1] = tiebreaker[i - k][0] + wsj[i]
                    bt[i][1] = 1

            # Extract characters that are not '-' and not X
            str_chars = [
                seqs_array[j][i] for i in range(len(seqs_array[j])) if seqs_array[j][i] not in ["-", ambig]
            ]

            icur = L - 1  # Python is 0-indexed, so L-1 instead of L
            if s[L - 1][0] < s[L - 1][1]:
                for k_pos in range(L - 1, L - 1 + k):
                    if k_pos < len(str_chars):
                        str_chars[k_pos] = mask
                bcur = 2
            else:
                bcur = 1

            while True:
                if bcur == 1 and bt[icur][bcur - 1] == 1:
                    icur -= 1
                    bcur = 1
                elif bcur == 1 and bt[icur][bcur - 1] == 2:
                    icur -= k
                    bcur = 2
                    for k_pos in range(icur, min(icur + k, len(str_chars))):
                        str_chars[k_pos] = mask
                elif bcur == 2 and bt[icur][bcur - 1] == 1:
                    icur -= k
                    bcur = 1
                elif bcur == 2 and bt[icur][bcur - 1] == 2:
                    icur -= 1
                    bcur = 2
                    if icur < len(str_chars):
                        str_chars[icur] = mask
                elif bcur == 1:
                    break
                else:
                    for i in range(icur):
                        if i < len(str_chars):
                            str_chars[i] = mask
                    break

            # Reconstruct output string
            strout = ""
            i = 0
            for t in range(len(sequences[j])):
                if sequences[j][t] in [ambig, "-"]:
                    strout += sequences[j][t]
                else:
                    if i < len(str_chars):
                        strout += str_chars[i]
                        i += 1
                    else:
                        strout += sequences[j][t]  # Fallback in case of index errors

            sequences_corrected.append(strout)

        return sequences_corrected

    fasta_in = fasta_to_dict(fasta_in_path)
    seq_names = []
    seq_descs = []
    sequences = []
    for seq_name in fasta_in:
        seq_names.append(seq_name)
        seq_descs.append(fasta_in[seq_name]["description"])
        sequences.append(fasta_in[seq_name]["sequence"].upper())
    seq_array = [list(seq) for seq in sequences]
    n = len(sequences[0]) if sequences else 0
    m = len(sequences)

    for t in settings.TAPER_PARAMS:
        tmp_mask = "@"
        output = correct_sequences(
            sequences, t["k"], ambig, tmp_mask, t["p"], t["q"], cutoff_threshold, conservative
        )
        L = t["L"]
        for j in range(m):
            start = None
            cnt = 0
            for i in range(n):
                if start is None and output[j][i] == tmp_mask:
                    start = i
                    cnt = 1
                elif start is not None and output[j][i] == tmp_mask:
                    cnt += 1
                elif start is not None and output[j][i] not in [tmp_mask, "-"]:
                    if cnt < L:
                        for k in range(start, i):
                            seq_array[j][k] = "-" if output[j][k] == "-" else settings.TAPER_MASK
                    start = None
                    cnt = 0

            if start is not None and cnt < L:
                for k in range(start, n):
                    seq_array[j][k] = "-" if output[j][k] == "-" else settings.TAPER_MASK

    sequences_corrected = ["".join(seq_array[j]) for j in range(m)]
    fasta_out = {}
    for j in range(len(sequences_corrected)):
        fasta_out[seq_names[j]] = {
            "sequence": sequences_corrected[j],
            "description": seq_descs[j],
        }
    dict_to_fasta(fasta_out, fasta_out_path)

    return fasta_out_path


def write_gff3(hits, marker_type, disable_stitching, tsv_comment, out_gff_path):
    def split_coords(coords, as_strings=False):
        """
        Split coordinate pairs by segment and change coordinates to 1-based, closed system for GFF3
        """
        try:
            changed_coords = []
            for fragment in coords.split("\n"):
                fragment_coords_strings = []
                fragment_coords_tuples = []
                for coord in fragment.split(","):
                    start, end = int(coord.split("-")[0]), int(coord.split("-")[1])
                    if start > end:
                        fragment_coords_strings.append(f"{end + 1}-{start}")
                        fragment_coords_tuples.append((f"{end + 1}", f"{start}"))
                    elif start == end:
                        fragment_coords_strings.append(f"{start}-{end}")
                        fragment_coords_tuples.append((f"{start}", f"{end}"))
                    else:
                        fragment_coords_strings.append(f"{start + 1}-{end}")
                        fragment_coords_tuples.append((f"{start + 1}", f"{end}"))
                if as_strings:
                    changed_coords.append(",".join(fragment_coords_strings))
                else:
                    changed_coords.append(fragment_coords_tuples)
            return changed_coords
        except ValueError:
            print(coords)

    def split_coords_min_max(coords, as_strings=False):
        """
        Split coordinate pairs by segment and change coordinates to 1-based, closed system for GFF3
        Get only start and end coordinate for each contig to annotate mRNA if 'disable_stitching' for
        example
        """
        try:
            changed_coords = []
            for fragment in coords.split("\n"):
                min_start, max_end = math.inf, 0
                fragment_coords_strings = []
                fragment_coords_tuples = []
                for coord in fragment.split(","):
                    start, end = int(coord.split("-")[0]), int(coord.split("-")[1])
                    if start > end:
                        if end + 1 < min_start:
                            min_start = end + 1
                        if start > max_end:
                            max_end = start
                    elif start == end:
                        if start < min_start:
                            min_start = start
                        if end > max_end:
                            max_end = end
                    else:
                        if start + 1 < min_start:
                            min_start = start + 1
                        if end > max_end:
                            max_end = end
                fragment_coords_strings.append(f"{min_start}-{max_end}")
                fragment_coords_tuples.append((f"{min_start}", f"{max_end}"))
                if as_strings:
                    changed_coords.append(",".join(fragment_coords_strings))
                else:
                    changed_coords.append(fragment_coords_tuples)
            return changed_coords
        except ValueError:
            print(coords)

    if marker_type in ["NUC", "PTD", "MIT"]:
        source = urllib.parse.quote("Captus (Scipio)")
        feature_type = urllib.parse.quote(f"protein_match:{marker_type}")
    else:
        source = urllib.parse.quote("Captus (BLAT)")
        feature_type = urllib.parse.quote(f"nucleotide_match:{marker_type}")

    phase = "."
    tsv_comment_lines = tsv_comment.split("\n")
    gff = ["##gff-version 3", f"#{tsv_comment_lines[0]}", f"#{tsv_comment_lines[1]}"]
    for ref in sorted(hits):
        gff.append(f"\n# {urllib.parse.quote(ref)}")
        for h in range(len(hits[ref])):
            if len(hits[ref]) == 1:
                h_name = urllib.parse.quote(ref)
            else:
                h_name = urllib.parse.quote(f"{ref}{settings.SEQ_NAME_SEP}{h:02}")
                gff.append(f"# {h_name}")
            strands = hits[ref][h]["strand"].split("\n")
            score = f"{hits[ref][h]['score']:.3f}"
            wscore = f"""WScore={urllib.parse.quote(f"{hits[ref][h]['wscore']:.3f}")}"""
            cover_pct = f"""Coverage={urllib.parse.quote(f"{hits[ref][h]['coverage']:.2f}")}"""
            ident_pct = f"""Identity={urllib.parse.quote(f"{hits[ref][h]['identity']:.2f}")}"""
            color = f"Color={urllib.parse.quote(settings.GFF_COLORS[marker_type])}"
            ref_coords = split_coords(hits[ref][h]["ref_coords"], as_strings=True)
            hit_ids = hits[ref][h]["hit_ids"].split("\n")
            hit_contigs = hits[ref][h]["hit_contigs"].split("\n")
            if disable_stitching is True and marker_type in ["NUC", "PTD", "MIT"]:
                ref_min_max = split_coords_min_max(hits[ref][h]["ref_coords"], as_strings=True)
                hit_min_max = split_coords_min_max(hits[ref][h]["hit_coords"])
                seq_id = urllib.parse.quote(hit_contigs[0])
                strand = strands[0]
                hit_id = f"ID={hit_ids[0]}"
                name = f"Name={h_name}"
                start = str(hit_min_max[0][0][0])
                end = str(hit_min_max[0][0][1])
                query = f"""Query={urllib.parse.quote(f"{hits[ref][h]['ref_name']}:{ref_min_max[0]}")}"""
                attributes = ";".join([hit_id, name, wscore, cover_pct, ident_pct, query])
                gff.append(
                    "\t".join([seq_id, source, "mRNA", start, end, score, strand, phase, attributes])
                )
            hit_coords = split_coords(hits[ref][h]["hit_coords"])
            for c in range(len(hit_coords)):
                seq_id = urllib.parse.quote(hit_contigs[c])
                strand = strands[c]
                hit_id = f"ID={hit_ids[c]}"
                name = f"Name={h_name}"
                for p in range(len(hit_coords[c])):
                    start = str(hit_coords[c][p][0])
                    end = str(hit_coords[c][p][1])
                    query = (
                        f"""Query={urllib.parse.quote(f"{hits[ref][h]['ref_name']}:{ref_coords[c]}")}"""
                    )
                    attributes = ";".join([hit_id, name, wscore, cover_pct, ident_pct, query, color])
                    gff.append(
                        "\t".join(
                            [
                                seq_id,
                                source,
                                feature_type,
                                start,
                                end,
                                score,
                                strand,
                                phase,
                                attributes,
                            ]
                        )
                    )
    with open(out_gff_path, "w") as gff_out:
        gff_out.write("\n".join(gff) + "\n")


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


def cds_from_gff(gff_path, fasta_path, bait_length):
    """
    Extract a FASTA genome file using its corresponding GFF annotation file,
    extract and concatenate CDS sequence, classify exons as short or long
    according to `min_exon_size`, and compile a data table about exons
    and total intron lengths

    Parameters
    ----------
    gff_path : Path
        Location of GFF annotation file
    fasta_path : Path
        Location of FASTA genome file
    min_exon_size : int
        Minimum exon length to be considered as long, usually the length of the
        projected baits, e.g. 120
    """

    if f"{gff_path}".endswith(".gz"):
        opener = gzip.open
    else:
        opener = open

    cds = {}
    with opener(gff_path, "rt") as gff_in:
        for line in gff_in:
            if not line.startswith("#"):
                record = line.strip().split("\t")
                if len(record) < 8:
                    continue
                if record[2].lower() == "cds":
                    notes = {
                        n.split("=")[0].strip().lower(): n.split("=")[1].strip()
                        for n in record[8].split(";")
                        if n
                    }
                    key = "protein_id" if "protein_id" in notes else "parent"
                    shift = int(record[7]) if record[7] != "." else 0
                    start = int(record[3])
                    end = int(record[4])
                    if notes[key] not in cds:
                        cds[notes[key]] = {
                            "seq_id": record[0],
                            "strand": record[6],
                            "coords": [(start, shift, end)],
                        }
                    else:
                        cds[notes[key]]["coords"].append((start, shift, end))

    genome = fasta_to_dict(fasta_path)
    fasta_cds = {}
    fasta_long_exons = {}
    fasta_short_exons = {}
    annotations = []

    cds_data = [
        [
            "cds_id",
            "position",
            "exons",
            "exons_len",
            "long_exons",
            "long_exons_len",
            "short_exons",
            "short_exons_len",
            "introns_len",
            "gene_len",
        ]
    ]
    for cds_id in cds:
        seq_id = cds[cds_id]["seq_id"]
        strand = cds[cds_id]["strand"]
        if cds[cds_id]["strand"] == "+":
            cds[cds_id]["coords"].sort()
        elif cds[cds_id]["strand"] == "-":
            cds[cds_id]["coords"].sort(reverse=True)
        coords = cds[cds_id]["coords"]
        cds[cds_id]["seqs"] = []
        for i in range(len(coords)):
            s, e = coords[i][0], coords[i][2]
            if i == 0:
                if strand == "+":
                    s = coords[i][0] + coords[i][1]
                elif strand == "-":
                    e = coords[i][2] - coords[i][1]
            if strand == "+":
                seq = genome[seq_id]["sequence"][s - 1 : e]
            elif strand == "-":
                seq = reverse_complement(genome[seq_id]["sequence"][s - 1 : e])
            cds[cds_id]["seqs"].append(seq)
        sequence = "".join(cds[cds_id]["seqs"])
        se = ",".join([f"{t[0]}-{t[2]}" for t in coords])
        description = f"[{seq_id}:{se}:{strand}]"

        # Ignore identical annotations
        if description not in annotations:
            annotations.append(description)
            fasta_cds[cds_id] = {
                "sequence": sequence,
                "description": description,
            }
            exons = len(coords)
            exons_len = 0
            long_exons = 0
            long_exons_len = 0
            short_exons = 0
            short_exons_len = 0
            if strand == "+":
                gene_len = abs(cds[cds_id]["coords"][-1][2] - cds[cds_id]["coords"][0][0]) + 1
            elif strand == "-":
                gene_len = abs(cds[cds_id]["coords"][0][2] - cds[cds_id]["coords"][-1][0]) + 1
            for i in range(len(coords)):
                s = f"{cds[cds_id]['coords'][i][0]}"
                e = f"{cds[cds_id]['coords'][i][2]}"
                seq = cds[cds_id]["seqs"][i]
                exons_len += len(seq)
                if len(seq) >= bait_length:
                    long_exons += 1
                    long_exons_len += len(seq)
                    fasta_long_exons[f"{cds_id}_exon{i + 1}"] = {
                        "sequence": seq,
                        "description": f"[{seq_id}:{s}-{e}:{strand}]",
                    }
                else:
                    short_exons += 1
                    short_exons_len += len(seq)
                    fasta_short_exons[f"{cds_id}_exon{i + 1}"] = {
                        "sequence": seq,
                        "description": f"[{seq_id}:{s}-{e}:{strand}]",
                    }
            cds_data.append(
                [
                    cds_id,
                    f"{seq_id}:{cds[cds_id]['coords'][0][0]}-{cds[cds_id]['coords'][-1][2]}:{strand}",
                    f"{exons}",
                    f"{exons_len}",
                    f"{long_exons}",
                    f"{long_exons_len}",
                    f"{short_exons}",
                    f"{short_exons_len}",
                    f"{gene_len - exons_len}",
                    f"{gene_len}",
                ]
            )

    return fasta_cds, fasta_long_exons, fasta_short_exons, cds_data


def mmseqs_cluster(
    mmseqs_path,
    mmseqs_method,
    clustering_dir,
    clustering_input_file,
    cluster_prefix,
    clustering_tmp_dir,
    sensitivity,
    min_identity,
    seq_id_mode,
    min_coverage,
    cov_mode,
    cluster_mode,
    threads,
    ram_mb,
    proteins=False,
):
    """
    Run MMseqs2 easy-linclust but perform some parameter checking/conversion before, the FASTA input
    file has to be decompressed, we can compress it afterwards
    """
    start = time.time()
    if sensitivity < 1:
        sensitivity = 1.0
    if sensitivity > 7.5:
        sensitivity = 7.5
    if not 0 < min_identity <= 1:
        min_identity = min(1.0, round((abs(min_identity) / 100), 3))
    if not 0 < min_coverage <= 1:
        min_coverage = min(1.0, round((abs(min_coverage) / 100), 3))

    result_prefix = f"{Path(clustering_dir, cluster_prefix)}"
    mmseqs_cmd = [
        mmseqs_path,
        mmseqs_method,
        f"{clustering_input_file}",
        f"{result_prefix}",
        f"{clustering_tmp_dir}",
        "--split-memory-limit",
        f"{ram_mb * settings.MMSEQS_RAM_FRACTION:.0f}M",
        "--mask",
        f"{settings.MMSEQS_MASK_LOW_COMPLEXITY}",
        "--spaced-kmer-mode",
        f"{0}",
        "-c",
        f"{min_coverage}",
        "--cov-mode",
        f"{cov_mode}",
        "--alignment-mode",
        f"{3}",
        "--min-seq-id",
        f"{min_identity}",
        "--seq-id-mode",
        f"{seq_id_mode}",
        "--cluster-mode",
        f"{cluster_mode}",
        "--kmer-per-seq-scale",
        f"{settings.MMSEQS_KMER_PER_SEQ_SCALE}",
        "--threads",
        f"{threads}",
    ]
    if proteins is True:
        mmseqs_cmd += [
            "--gap-open",
            f"{max(1, settings.MMSEQS_GAP_OPEN_AA)}",
            "--gap-extend",
            f"{max(1, settings.MMSEQS_GAP_EXTEND_AA)}",
        ]
    else:
        mmseqs_cmd += [
            "--gap-open",
            f"{max(1, settings.MMSEQS_GAP_OPEN_NT)}",
            "--gap-extend",
            f"{max(1, settings.MMSEQS_GAP_EXTEND_NT)}",
        ]
    if mmseqs_method == "easy-cluster":
        mmseqs_cmd += ["-s", f"{sensitivity}"]
        if cluster_mode != 0:
            mmseqs_cmd += ["--cluster-reassign", f"{1}"]
    mmseqs_log_file = Path(clustering_dir, f"{cluster_prefix}_mmseqs.log")
    mmseqs_thread = ElapsedTimeThread()
    mmseqs_thread.start()
    with open(mmseqs_log_file, "w") as mmseqs_log:
        mmseqs_log.write(f"Captus' MMseqs2 Command:\n  {' '.join(mmseqs_cmd)}\n\n")
    with open(mmseqs_log_file, "a") as mmseqs_log:
        subprocess.run(mmseqs_cmd, stdout=mmseqs_log, stdin=mmseqs_log)
    mmseqs_thread.stop()
    mmseqs_thread.join()
    print()

    message = bold(f" \u2514\u2500\u2192 Clustering completed: [{elapsed_time(time.time() - start)}]")
    return message


def vsearch_cluster(
    vsearch_path,
    vsearch_method,
    clustering_dir,
    clustering_input_file,
    cluster_prefix,
    min_identity,
    seq_id_mode,
    min_coverage,
    strand,
    threads,
):
    start = time.time()
    if not 0 < min_identity <= 1:
        min_identity = min(1.0, round((abs(min_identity) / 100), 3))
    if not 0 < min_coverage <= 1:
        min_coverage = min(1.0, round((abs(min_coverage) / 100), 3))

    vsearch_cmd = [
        vsearch_path,
        vsearch_method,
        f"{Path(clustering_input_file)}",
        "--strand",
        strand,
        "--id",
        f"{min_identity}",
        "--iddef",
        f"{seq_id_mode}",
        "--query_cov",
        f"{min_coverage}",
        "--userout",
        f"{Path(clustering_dir, f'{cluster_prefix}_cluster.tsv')}",
        "--userfields",
        "target+query",
        "--maxrejects",
        "0",
        "--threads",
        f"{threads}",
    ]
    vsearch_log_file = Path(clustering_dir, f"{cluster_prefix}_vsearch.log")
    vsearch_thread = ElapsedTimeThread()
    vsearch_thread.start()
    with open(vsearch_log_file, "w") as vsearch_log:
        vsearch_log.write(f"Captus' VSEARCH Command:\n  {' '.join(vsearch_cmd)}\n\n")
    with open(vsearch_log_file, "a") as vsearch_log:
        subprocess.run(vsearch_cmd, stdout=vsearch_log, stdin=vsearch_log, stderr=vsearch_log)
    vsearch_thread.stop()
    vsearch_thread.join()
    print()

    message = bold(f" \u2514\u2500\u2192 Clustering completed: [{elapsed_time(time.time() - start)}]")
    return message


def resolve_iupac(seq):
    """
    Select a random nucleotide for IUPAC ambiguity codes, returning sequences with only A,T,G,C
    Leaves capitalization intact (important for low complexity masking)

    Parameters
    ----------
    seq : str
        Input sequence

    Returns
    -------
    str
        Output sequence with ambiguities replaced
    """
    seq_upper = seq.upper()
    seq_out = ""
    for i in range(len(seq)):
        r = seq_upper[i]
        if r in NT_IUPAC:
            r_resolved = random.choice(NT_IUPAC[r])
            if seq[i].islower():
                seq_out += r_resolved.lower()
            else:
                seq_out += r_resolved
    return seq_out


def bait_stats(bait_seq, hybrid_chem, sodium_conc, formamide_conc):
    def calc_melt_temp(bait_len, gc_content, hybrid_chem, sodium_conc, formamide_conc):
        """
        Calculate oligonucleotide melting temperature

        Parameters
        ----------
        seq_len : int
            Oligonucleotide length
        gc_content : float
            GC content in percentage
        hybrid_chem : str
            Bait hybridization chemistry, can be RNA-DNA, DNA-DNA, or RNA-RNA
        sodium_conc : float
            Na concentration
        formamide_conc : float
            Formamide concentration

        Returns
        -------
        float
            Oligonucleotide melting temperature
        """
        gc_content /= 100
        melt_temp = -1.0
        if hybrid_chem == "RNA-DNA":
            melt_temp = (
                79.8
                + (58.4 * gc_content)
                + (11.8 * math.pow(gc_content, 2))
                - (820.0 / bait_len)
                + (18.5 * math.log(sodium_conc))
                - (0.5 * formamide_conc)
            )
        elif hybrid_chem == "DNA-DNA":
            melt_temp = (
                81.5
                + (41.0 * gc_content)
                - (500.0 / bait_len)
                + (16.6 * math.log(sodium_conc))
                - (0.62 * formamide_conc)
            )
        elif hybrid_chem == "RNA-RNA":
            melt_temp = (
                79.8
                + (58.4 * gc_content)
                + (11.8 * math.pow(gc_content, 2))
                - (820.0 / bait_len)
                + (18.5 * math.log(sodium_conc))
                - (0.35 * formamide_conc)
            )
        return melt_temp

    gc_bp = []
    lower_bp = 0
    homopolymer_len = 0
    max_homopolymer_len = 1
    bait_len = len(bait_seq)
    bait_seq_upper = bait_seq.upper()

    for i in range(bait_len):
        r = bait_seq_upper[i]

        if r in NT_IUPAC:
            g_prop = NT_IUPAC[r].count("G") / len(NT_IUPAC[r])
            c_prop = NT_IUPAC[r].count("C") / len(NT_IUPAC[r])
            gc_bp.append(g_prop + c_prop)

        if i == 0:
            homopolymer_len = 1
        else:
            if r == bait_seq_upper[i - 1]:
                homopolymer_len += 1
            else:
                if homopolymer_len > max_homopolymer_len:
                    max_homopolymer_len = homopolymer_len
                homopolymer_len = 1

        if bait_seq[i].islower():
            lower_bp += 1

    gc_content = round(sum(gc_bp) / len(gc_bp) * 100, 2)

    stats = {
        "gc": gc_content,
        "melt_temp": calc_melt_temp(bait_len, gc_content, hybrid_chem, sodium_conc, formamide_conc),
        "low_complexity": round(lower_bp / bait_len * 100, 2),
        "max_homopolymer_len": max_homopolymer_len,
        "Ns": bait_seq_upper.count("N"),
    }

    return stats


def import_busco_odb1x(odb1x_tar_path: Path):
    """
    Import BUSCO databases to be used as reference target file for marker extraction, the BUSCO
    databse file has to be a .tar.gz from the odb10, containing the file 'ancestral_variants' which
    contains the aminoacid sequences representing the ortholog groups

    Parameters
    ----------
    odb10_tar_path : Path
        Any of the files from https://busco-data.ezlab.org/v5/data/lineages/

    Returns
    -------
    busco_targets
        A fasta_dict with the sequence names formatted as Captus needs them for marker extraction
    """
    busco_targets = {}
    with tarfile.open(odb1x_tar_path) as tf:
        for member in tf.getmembers():
            if Path(member.name).name == "ancestral_variants":
                seq = ""
                name = ""
                desc = ""
                for line in tf.extractfile(member).readlines():
                    line = line.decode().strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if seq:
                            busco_targets[name] = {
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
                        name = f"{name}{settings.REF_CLUSTER_SEP}{name.split('_')[0]}"
                    else:
                        seq += line
                if seq:
                    busco_targets[name] = {
                        "description": desc,
                        "sequence": seq,
                    }
    if busco_targets:
        return busco_targets
    else:
        return None


def rehead_root_msa(fasta_in: Path, fasta_out: Path, outgroup: list, remove_R_=False):
    if outgroup:
        outgroup = outgroup.split(",")
    else:
        outgroup = []

    unaligned = fasta_to_dict(fasta_in)

    if remove_R_:
        aligned = {}
        aligned_w_R_ = fasta_to_dict(fasta_out)
        for seq_name in aligned_w_R_:
            if seq_name.startswith("_R_"):
                aligned[seq_name[3:]] = aligned_w_R_[seq_name]
            else:
                aligned[seq_name] = aligned_w_R_[seq_name]
    else:
        aligned = fasta_to_dict(fasta_out)

    reheaded = {}
    for sample_name in outgroup:
        for seq_name in sorted(aligned):
            if seq_name.split(settings.SEQ_NAME_SEP)[0] == sample_name:
                reheaded[seq_name] = {
                    "description": unaligned[seq_name]["description"],
                    "sequence": aligned[seq_name]["sequence"],
                }
    for seq_name in aligned:
        if seq_name.split(settings.SEQ_NAME_SEP)[0] not in outgroup:
            reheaded[seq_name] = {
                "description": unaligned[seq_name]["description"],
                "sequence": aligned[seq_name]["sequence"],
            }
    dict_to_fasta(reheaded, fasta_out)

    return
