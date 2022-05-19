#!/usr/bin/env python3
"""
Copyright 2020 Edgardo M. Ortiz (e.ortiz.v@gmail.com)
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
import re
import statistics
import sys
import urllib

from . import settings_assembly as set_a

# Regular expression to match MEGAHIT headers assembled within Captus
CAPTUS_MEGAHIT_HEADER_REGEX = r'^NODE_\d+_length_\d+_cov_\d+\.\d+_k_\d{2,3}_flag_\d'

# Regular expression to match MEGAHIT headers
MEGAHIT_HEADER_REGEX = r'^k\d{2,3}_\d+\sflag=\d\smulti=\d+\.\d+\slen=\d+$'

# Regular expression to match SKESA headers
SKESA_HEADER_REGEX = r'^Contig_\d+_\d+\.\d+$|^Contig_\d+_\d+\.\d+_Circ$'

# Regular expression to match Spades headers
SPADES_HEADER_REGEX = r'^NODE_\d+_length_\d+_cov_\d+\.\d+|^EDGE_\d+_length_\d+_cov_\d+\.\d+'

# Set of valid nucleotides, including IUPAC ambiguities
NT_IUPAC = {
    "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"],
    "M": ["A", "C"], "R": ["G", "A"], "W": ["A", "T"],
    "S": ["G", "C"], "Y": ["C", "T"], "K": ["G", "T"],
    "V": ["A", "G", "C"], "H": ["A", "C", "T"],
    "D": ["G", "A", "T"], "B": ["G", "C", "T"],
    "N": ["A", "C", "G", "T"], "-": ["-"],
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
        NT_PIDS[f'{a}{b}'] = matches / combos
    i += 1
NT_PIDS["--"] = 0.0

# Set of valid aminoacids, including IUPAC ambiguities
AA_IUPAC = {
    "A": ["A"], "B": ["D", "N"], "C": ["C"], "D": ["D"], "E": ["E"], "F": ["F"],
    "G": ["G"], "H": ["H"], "I": ["I"], "J": ["I", "L"], "K": ["K"], "L": ["L"],
    "M": ["M"], "N": ["N"], "O": ["O"], "P": ["P"], "Q": ["Q"], "R": ["R"], "S": ["S"],
    "T": ["T"], "U": ["U"], "V": ["V"], "W": ["W"], "Y": ["Y"], "Z": ["E", "Q"],
    "X": ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "O", "P",
          "Q", "R", "S", "T", "U", "V", "W", "Y"], "*": ["*"], "-": ["-"],
}

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
        AA_PIDS[a+b] = matches / combos
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
    "A": "T", "T": "A", "G": "C", "C": "G", "a": "t", "t": "a", "g": "c", "c": "g",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "r": "y", "y": "r", "s": "s", "w": "w",
    "K": "M", "M": "K", "B": "V", "V": "B", "k": "m", "m": "k", "b": "v", "v": "b",
    "D": "H", "H": "D", "N": "N", "d": "h", "h": "d", "n": "n",
    ".": ".", "-": "-", "?": "?"
}

# Codons corresponding to aminoacid strings in GENETIC_CODES
CODONS = {
    "base1":   "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
    "base2":   "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
    "base3":   "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG",
}

# NCBI Genetic Codes (modified from ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt)
GENETIC_CODES = {
    1:  {"name": "Standard" ,
         "aa": "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "ss": "---M------**--*----M---------------M----------------------------"},
    2:  {"name": "Vertebrate Mitochondrial" ,
         "aa": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
         "ss": "----------**--------------------MMMM----------**---M------------"},
    3:  {"name": "Yeast Mitochondrial" ,
         "aa": "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "ss": "----------**----------------------MM---------------M------------"},
    4:  {"name": "Mold Mitochondrial; Protozoan Mitochondrial;"
                 " Coelenterate Mitochondrial; Mycoplasma; Spiroplasma" ,
         "aa": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "ss": "--MM------**-------M------------MMMM---------------M------------"},
    5:  {"name": "Invertebrate Mitochondrial" ,
         "aa": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
         "ss": "---M------**--------------------MMMM---------------M------------"},
    6:  {"name": "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear" ,
         "aa": "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "ss": "--------------*--------------------M----------------------------"},
    9:  {"name": "Echinoderm Mitochondrial; Flatworm Mitochondrial" ,
         "aa": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
         "ss": "----------**-----------------------M---------------M------------"},
    10: {"name": "Euplotid Nuclear" ,
         "aa": "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "ss": "----------**-----------------------M----------------------------"},
    11: {"name": "Bacterial, Archaeal and Plant Plastid" ,
         "aa": "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "ss": "---M------**--*----M------------MMMM---------------M------------"},
    12: {"name": "Alternative Yeast Nuclear" ,
         "aa": "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "ss": "----------**--*----M---------------M----------------------------"},
    13: {"name": "Ascidian Mitochondrial" ,
         "aa": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
         "ss": "---M------**----------------------MM---------------M------------"},
    14: {"name": "Alternative Flatworm Mitochondrial" ,
         "aa": "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
         "ss": "-----------*-----------------------M----------------------------"},
    15: {"name": "Blepharisma Macronuclear" ,
         "aa": "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "ss": "----------*---*--------------------M----------------------------"},
    16: {"name": "Chlorophycean Mitochondrial" ,
         "aa": "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "ss": "----------*---*--------------------M----------------------------"},
    21: {"name": "Trematode Mitochondrial" ,
         "aa": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
         "ss": "----------**-----------------------M---------------M------------"},
    22: {"name": "Scenedesmus obliquus Mitochondrial" ,
         "aa": "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "ss": "------*---*---*--------------------M----------------------------"},
    23: {"name": "Thraustochytrium Mitochondrial" ,
         "aa": "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "ss": "--*-------**--*-----------------M--M---------------M------------"},
    24: {"name": "Rhabdopleuridae Mitochondrial" ,
         "aa": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
         "ss": "---M------**-------M---------------M---------------M------------"},
    25: {"name": "Candidate Division SR1 and Gracilibacteria" ,
         "aa": "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "ss": "---M------**-----------------------M---------------M------------"},
    26: {"name": "Pachysolen tannophilus Nuclear" ,
         "aa": "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "ss": "----------**--*----M---------------M----------------------------"},
    27: {"name": "Karyorelict Nuclear" ,
         "aa": "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "ss": "--------------*--------------------M----------------------------"},
    28: {"name": "Condylostoma Nuclear" ,
         "aa": "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "ss": "----------**--*--------------------M----------------------------"},
    29: {"name": "Mesodinium Nuclear" ,
         "aa": "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "ss": "--------------*--------------------M----------------------------"},
    30: {"name": "Peritrich Nuclear" ,
         "aa": "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "ss": "--------------*--------------------M----------------------------"},
    31: {"name": "Blastocrithidia Nuclear" ,
         "aa": "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "ss": "----------**-----------------------M----------------------------"},
    32: {"name": "Balanophoraceae Plastid" ,
         "aa": "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "ss": "---M------*---*----M------------MMMM---------------M------------"},
    33: {"name": "Cephalodiscidae Mitochondrial" ,
         "aa": "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
         "ss": "---M-------*-------M---------------M---------------M------------"},
}

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

# PAM250 matrix for scoring protein alignments
PAM250 = score_matrix_to_dict([
    [".","A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","B","J","Z","X","*"],
    ["A", 2 ,-2 , 0 , 0 ,-2 , 0 , 0 , 1 ,-1 ,-1 ,-2 ,-1 ,-1 ,-3 , 1 , 1 , 1 ,-6 ,-3 , 0 , 0 ,-1 , 0 ,-1 ,-8 ],
    ["R",-2 , 6 , 0 ,-1 ,-4 , 1 ,-1 ,-3 , 2 ,-2 ,-3 , 3 , 0 ,-4 , 0 , 0 ,-1 , 2 ,-4 ,-2 ,-1 ,-3 , 0 ,-1 ,-8 ],
    ["N", 0 , 0 , 2 , 2 ,-4 , 1 , 1 , 0 , 2 ,-2 ,-3 , 1 ,-2 ,-3 , 0 , 1 , 0 ,-4 ,-2 ,-2 , 2 ,-3 , 1 ,-1 ,-8 ],
    ["D", 0 ,-1 , 2 , 4 ,-5 , 2 , 3 , 1 , 1 ,-2 ,-4 , 0 ,-3 ,-6 ,-1 , 0 , 0 ,-7 ,-4 ,-2 , 3 ,-3 , 3 ,-1 ,-8 ],
    ["C",-2 ,-4 ,-4 ,-5 ,12 ,-5 ,-5 ,-3 ,-3 ,-2 ,-6 ,-5 ,-5 ,-4 ,-3 , 0 ,-2 ,-8 , 0 ,-2 ,-4 ,-5 ,-5 ,-1 ,-8 ],
    ["Q", 0 , 1 , 1 , 2 ,-5 , 4 , 2 ,-1 , 3 ,-2 ,-2 , 1 ,-1 ,-5 , 0 ,-1 ,-1 ,-5 ,-4 ,-2 , 1 ,-2 , 3 ,-1 ,-8 ],
    ["E", 0 ,-1 , 1 , 3 ,-5 , 2 , 4 , 0 , 1 ,-2 ,-3 , 0 ,-2 ,-5 ,-1 , 0 , 0 ,-7 ,-4 ,-2 , 3 ,-3 , 3 ,-1 ,-8 ],
    ["G", 1 ,-3 , 0 , 1 ,-3 ,-1 , 0 , 5 ,-2 ,-3 ,-4 ,-2 ,-3 ,-5 , 0 , 1 , 0 ,-7 ,-5 ,-1 , 0 ,-4 , 0 ,-1 ,-8 ],
    ["H",-1 , 2 , 2 , 1 ,-3 , 3 , 1 ,-2 , 6 ,-2 ,-2 , 0 ,-2 ,-2 , 0 ,-1 ,-1 ,-3 , 0 ,-2 , 1 ,-2 , 2 ,-1 ,-8 ],
    ["I",-1 ,-2 ,-2 ,-2 ,-2 ,-2 ,-2 ,-3 ,-2 , 5 , 2 ,-2 , 2 , 1 ,-2 ,-1 , 0 ,-5 ,-1 , 4 ,-2 , 3 ,-2 ,-1 ,-8 ],
    ["L",-2 ,-3 ,-3 ,-4 ,-6 ,-2 ,-3 ,-4 ,-2 , 2 , 6 ,-3 , 4 , 2 ,-3 ,-3 ,-2 ,-2 ,-1 , 2 ,-3 , 5 ,-3 ,-1 ,-8 ],
    ["K",-1 , 3 , 1 , 0 ,-5 , 1 , 0 ,-2 , 0 ,-2 ,-3 , 5 , 0 ,-5 ,-1 , 0 , 0 ,-3 ,-4 ,-2 , 1 ,-3 , 0 ,-1 ,-8 ],
    ["M",-1 , 0 ,-2 ,-3 ,-5 ,-1 ,-2 ,-3 ,-2 , 2 , 4 , 0 , 6 , 0 ,-2 ,-2 ,-1 ,-4 ,-2 , 2 ,-2 , 3 ,-2 ,-1 ,-8 ],
    ["F",-3 ,-4 ,-3 ,-6 ,-4 ,-5 ,-5 ,-5 ,-2 , 1 , 2 ,-5 , 0 , 9 ,-5 ,-3 ,-3 , 0 , 7 ,-1 ,-4 , 2 ,-5 ,-1 ,-8 ],
    ["P", 1 , 0 , 0 ,-1 ,-3 , 0 ,-1 , 0 , 0 ,-2 ,-3 ,-1 ,-2 ,-5 , 6 , 1 , 0 ,-6 ,-5 ,-1 ,-1 ,-2 , 0 ,-1 ,-8 ],
    ["S", 1 , 0 , 1 , 0 , 0 ,-1 , 0 , 1 ,-1 ,-1 ,-3 , 0 ,-2 ,-3 , 1 , 2 , 1 ,-2 ,-3 ,-1 , 0 ,-2 , 0 ,-1 ,-8 ],
    ["T", 1 ,-1 , 0 , 0 ,-2 ,-1 , 0 , 0 ,-1 , 0 ,-2 , 0 ,-1 ,-3 , 0 , 1 , 3 ,-5 ,-3 , 0 , 0 ,-1 ,-1 ,-1 ,-8 ],
    ["W",-6 , 2 ,-4 ,-7 ,-8 ,-5 ,-7 ,-7 ,-3 ,-5 ,-2 ,-3 ,-4 , 0 ,-6 ,-2 ,-5 ,17 , 0 ,-6 ,-5 ,-3 ,-6 ,-1 ,-8 ],
    ["Y",-3 ,-4 ,-2 ,-4 , 0 ,-4 ,-4 ,-5 , 0 ,-1 ,-1 ,-4 ,-2 , 7 ,-5 ,-3 ,-3 , 0 ,10 ,-2 ,-3 ,-1 ,-4 ,-1 ,-8 ],
    ["V", 0 ,-2 ,-2 ,-2 ,-2 ,-2 ,-2 ,-1 ,-2 , 4 , 2 ,-2 , 2 ,-1 ,-1 ,-1 , 0 ,-6 ,-2 , 4 ,-2 , 2 ,-2 ,-1 ,-8 ],
    ["B", 0 ,-1 , 2 , 3 ,-4 , 1 , 3 , 0 , 1 ,-2 ,-3 , 1 ,-2 ,-4 ,-1 , 0 , 0 ,-5 ,-3 ,-2 , 3 ,-3 , 2 ,-1 ,-8 ],
    ["J",-1 ,-3 ,-3 ,-3 ,-5 ,-2 ,-3 ,-4 ,-2 , 3 , 5 ,-3 , 3 , 2 ,-2 ,-2 ,-1 ,-3 ,-1 , 2 ,-3 , 5 ,-2 ,-1 ,-8 ],
    ["Z", 0 , 0 , 1 , 3 ,-5 , 3 , 3 , 0 , 2 ,-2 ,-3 , 0 ,-2 ,-5 , 0 , 0 ,-1 ,-6 ,-4 ,-2 , 2 ,-2 , 3 ,-1 ,-8 ],
    ["X",-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-8 ],
    ["*",-8 ,-8 ,-8 ,-8 ,-8 ,-8 ,-8 ,-8 ,-8 ,-8 ,-8 ,-8 ,-8 ,-8 ,-8 ,-8 ,-8 ,-8 ,-8 ,-8 ,-8 ,-8 ,-8 ,-8 , 1 ],
])


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
    for b1, b2, b3, a, s in zip(CODONS["base1"],
                                CODONS["base2"],
                                CODONS["base3"],
                                GENETIC_CODES[genetic_code_id]["aa"],
                                GENETIC_CODES[genetic_code_id]["ss"]):
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
            return [f"{p1}{p2}{p3}" for p1 in NT_IUPAC[codon[0]]
                                    for p2 in NT_IUPAC[codon[1]]
                                    for p3 in NT_IUPAC[codon[2]]]


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
    if frame < 0: seq = reverse_complement(seq)

    codons = [seq[p:p + 3] for p in range(abs(frame) - 1, len(seq), 3)
              if len(seq[p:p + 3]) == 3]

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

    if len(codons) == 1: return seq_AA

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
            counts_Xs.append(translation.count("X"))
            counts_stops.append(translation.count("*"))
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
                "sequence": guess_frame(in_fasta_dict[seq_name]["sequence"], gc, start_as_M)}
    else:
        for seq_name in in_fasta_dict:
            translated_fasta[seq_name] = {
                "description": in_fasta_dict[seq_name]["description"],
                "sequence": translate(in_fasta_dict[seq_name]["sequence"], gc, frame, start_as_M)}

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
            f'{in_fasta_dict[seq_name]["sequence"][:-1].replace("*", "X")}'
            f'{in_fasta_dict[seq_name]["sequence"][-1]}'
        )
        seq_out = seq_out.replace(" ", "").replace("-", "")
        if in_fasta_dict[seq_name]["sequence"] != seq_out:
            has_premature_stops = True
        out_fasta_dict[seq_name] = {
            "sequence": seq_out,
            "description": in_fasta_dict[seq_name]["description"]
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


sys.setrecursionlimit(set_a.RECURSION_LIMIT)
def align_prots(s1, s2, method, scoring_matrix=PAM250):

    def gapless(s1, s2, aln):
        len_s1  = len(s1)
        len_s2  = len(s2)
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
                    pid = AA_PIDS["".join(sorted(f"{s1[i+j]}{s2[j]}"))]
                    matches += pid
                    if pid < 0.5: mismatches.append(i+j)
                    s2_aln[i+j] = s2[j]
                else:
                    pid = AA_PIDS["".join(sorted(f"{s1[j]}{s2[i+j]}"))]
                    matches += pid
                    if pid < 0.5: mismatches.append(i+j)
                    s1_aln[i+j] = s1[j]
            try:
                match_rate = matches / min_len
            except ZeroDivisionError:
                match_rate = 0.0
            if match_rate > aln["match_rate"]:
                aln["matches"]    = matches
                aln["mismatches"] = mismatches
                aln["match_rate"] = match_rate
                aln["s1_aln"]     = "".join(s1_aln)
                aln["s1_start"]   = s1_start
                aln["s1_end"]     = s1_end
                aln["s2_aln"]     = "".join(s2_aln)
                aln["s2_start"]   = s2_start
                aln["s2_end"]     = s2_end

        return aln

    def best_score(val_matrix, row, col, nt_left, nt_top, method, scoring_matrix):
        left_score = val_matrix[row  ][col-1] + scoring_matrix["gap"]
        diag_score = val_matrix[row-1][col-1] + scoring_matrix[f"{nt_left}{nt_top}"]
        up_score   = val_matrix[row-1][col  ] + scoring_matrix["gap"]
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
        seq_top  = f" {s2}"
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
            for l in range(1, len(seq_left)):
                for t in range(1, len(seq_top)):
                    val_matrix[l][t], dir_matrix[l][t] = best_score(
                        val_matrix, l, t, seq_left[l], seq_top[t], "nw", scoring_matrix
                    )
            return dir_matrix
        if method == "sw":
            top_score = 0
            opt_l = 0
            opt_t = 0
            for l in range(1, len(seq_left)):
                for t in range(1, len(seq_top)):
                    val_matrix[l][t], dir_matrix[l][t] = best_score(
                        val_matrix, l, t, seq_left[l], seq_top[t], "sw", scoring_matrix
                    )
                    if val_matrix[l][t] > top_score:
                        top_score = val_matrix[l][t]
                        opt_l = l
                        opt_t = t
            return dir_matrix, opt_l, opt_t

    def traceback_nw(dir_matrix, seq_left, seq_top):
        seq_left = f" {seq_left}"
        seq_top  = f" {seq_top}"
        traceback_nw.align_left = ""
        traceback_nw.align_top = ""
        def stepback(dir_matrix, seq_left, seq_top, l, t):
            if l < 0 or t < 0:
                return
            if dir_matrix[l][t] == 3:
                stepback(dir_matrix, seq_left, seq_top, l-1, t-1)
                traceback_nw.align_top  += seq_top[t]
                traceback_nw.align_left += seq_left[l]
            elif dir_matrix[l][t] == 2:
                stepback(dir_matrix, seq_left, seq_top, l, t-1)
                traceback_nw.align_top  += seq_top[t]
                traceback_nw.align_left += "-"
            elif dir_matrix[l][t] == 1:
                stepback(dir_matrix, seq_left, seq_top, l-1, t)
                traceback_nw.align_top  += "-"
                traceback_nw.align_left += seq_left[l]
            else:
                return
        stepback(dir_matrix, seq_left, seq_top, len(seq_left)-1, len(seq_top)-1)
        return traceback_nw.align_left, traceback_nw.align_top

    def traceback_sw(dir_matrix, opt_l, opt_t, seq_left, seq_top):
        seq_left = f" {seq_left}"
        seq_top  = f" {seq_top}"
        traceback_sw.align_left = ""
        traceback_sw.align_top = ""
        def stepback(dir_matrix, seq_left, seq_top, l, t):
            if l < 0 or t < 0:
                return
            if dir_matrix[l][t] == 3:
                stepback(dir_matrix, seq_left, seq_top, l-1, t-1)
                traceback_sw.align_top  += seq_top[t]
                traceback_sw.align_left += seq_left[l]
            elif dir_matrix[l][t] == 2:
                stepback(dir_matrix, seq_left, seq_top, l, t-1)
                traceback_sw.align_top  += seq_top[t]
                traceback_sw.align_left += "-"
            elif dir_matrix[l][t] == 1:
                stepback(dir_matrix, seq_left, seq_top, l-1, t)
                traceback_sw.align_top  += "-"
                traceback_sw.align_left += seq_left[l]
            else:
                return
        stepback(dir_matrix, seq_left, seq_top, opt_l, opt_t)
        add_top = ""
        add_left = ""
        if len(seq_left) > len(seq_top):
            add_left = seq_left[opt_l+1:]
            add_top = seq_top[opt_t+1:]
            add_top += (len(add_left) - len(add_top)) * "-"
        else:
            add_top = seq_top[opt_t+1:]
            add_left = seq_left[opt_l+1:]
            add_left += (len(add_top) - len(add_left)) * "-"
        return f"{traceback_sw.align_left}{add_left}", f"{traceback_sw.align_top}{add_top}"

    def needleman_wunsch(s1, s2, aln, scoring_matrix):
        s1_aln, s2_aln = traceback_nw(direction_matrix(s1, s2, "nw", scoring_matrix), s1, s2)
        aln["s1_start"] = len(s1_aln) - len(s1_aln.lstrip("-"))
        aln["s1_end"]   = len(s1_aln.rstrip("-"))
        aln["s1_aln"]   = s1_aln
        aln["s2_start"] = len(s2_aln) - len(s2_aln.lstrip("-"))
        aln["s2_end"]   = len(s2_aln.rstrip("-"))
        aln["s2_aln"]   = s2_aln
        aln_start = max(aln["s1_start"], aln["s2_start"])
        aln_end = min(aln["s1_end"], aln["s2_end"])
        for i in range(aln_start, aln_end):
            if "-" not in [s1_aln[i], s2_aln[i]]:
                pid = AA_PIDS["".join(sorted(f"{s1_aln[i]}{s2_aln[i]}"))]
                aln["matches"] += pid
                if pid < 0.5: aln["mismatches"].append(i)
        try:
            aln["match_rate"] = aln["matches"] / (aln_end - aln_start)
        except ZeroDivisionError:
            aln["match_rate"] = 0.0

        return aln

    def smith_waterman(s1, s2, aln, scoring_matrix):
        dir_mat, opt_l, opt_t = direction_matrix(s1, s2, "sw", scoring_matrix)
        s1_aln, s2_aln = traceback_sw(dir_mat, opt_l, opt_t, s1, s2)
        aln["s1_start"] = len(s1_aln) - len(s1_aln.lstrip("-"))
        aln["s1_end"]   = len(s1_aln.rstrip("-"))
        aln["s1_aln"]   = s1_aln
        aln["s2_start"] = len(s2_aln) - len(s2_aln.lstrip("-"))
        aln["s2_end"]   = len(s2_aln.rstrip("-"))
        aln["s2_aln"]   = s2_aln
        aln_start = max(aln["s1_start"], aln["s2_start"])
        aln_end = min(aln["s1_end"], aln["s2_end"])
        for i in range(aln_start, aln_end):
            if "-" not in [s1_aln[i], s2_aln[i]]:
                pid = AA_PIDS["".join(sorted(f"{s1_aln[i]}{s2_aln[i]}"))]
                aln["matches"] += pid
                if pid < 0.5: aln["mismatches"].append(i)
        try:
            aln["match_rate"] = aln["matches"] / (aln_end - aln_start)
        except ZeroDivisionError:
            aln["match_rate"] = 0.0

        return aln


    valid_methods = ["gapless", "nw", "sw"]
    if method not in valid_methods: return False

    # Alignment data template:
    aln = {
        "matches":    0,
        "mismatches": [],
        "match_rate": 0.0,
        "s1_start":   0,
        "s1_end":     0,
        "s1_aln":     "",
        "s2_start":   0,
        "s2_end":     0,
        "s2_aln":     "",
    }

    # Capitalize sequences
    s1 = s1.upper()
    s2 = s2.upper()

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
            if ignore_internal_gaps:
                if not "-" in pair:
                    matches += PIDS[pair]
                    aligned_length += 1
            else:
                if pair != "--":
                    matches += PIDS[pair]
                    aligned_length += 1
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
        counts = {r:site.count(r) for r in ss_site}
        combos = ((len_site * (len_site - 1)) / 2)
        if "-" in counts:
            combos -= ((counts["-"] * (counts["-"] - 1)) / 2)
        if combos == 0:
            return [0,0]
        matches = 0
        for ra in counts:
            matches += ((counts[ra] * (counts[ra] - 1)) / 2) * PIDS[ra+ra]
        i = 0
        for ra in ss_site[i:]:
            for rb in ss_site[i+1:]:
                matches += counts[ra] * counts[rb] * PIDS[ra+rb]
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
        seq  = ""
        name = ""
        desc = ""
        for line in fasta_in:
            line = line.strip("\n")
            if not line: continue
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
    in_fasta_dict, out_fasta_path, wrap=0, sort=False, append=False, write_if_empty=False
):
    """
    Saves a `in_fasta_dict` from function `fasta_to_dict()` as a FASTA file to `out_fasta_path`
    """
    if f"{out_fasta_path}".endswith(".gz"):
        opener = gzip.open
    else:
        opener = open
    if append is True:
        action = "at"
    else:
        action = "wt"
    if in_fasta_dict:
        if sort: in_fasta_dict = dict(sorted(in_fasta_dict.items(), key=lambda x: x[0]))
        with opener(out_fasta_path, action) as fasta_out:
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
            with opener(out_fasta_path, action) as fasta_out:
                fasta_out.write("")
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


def alignment_stats(fasta_path):
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
                    return False
        if length > 0:
            return length
        else:
            return False

    def num_samples(fasta_dict):
        sample_names = []
        for seq in fasta_dict:
            sample_names.append(seq.split(set_a.SEQ_NAME_SEP)[0])

        return len(set(sample_names))

    def mark_terminal_gaps(fasta_dict):
        fasta_list = []
        for seq_name in fasta_dict:
            seq = fasta_dict[seq_name]["sequence"]
            left_gaps  = "#" * (len(seq) - len(seq.lstrip("-")))
            right_gaps = "#" * (len(seq) - len(seq.rstrip("-")))
            fasta_list.append(f'{left_gaps}{seq.strip("-")}{right_gaps}')

        return fasta_list

    def transpose_aln(fasta_list:list, num_sites):
        sites = [""] * num_sites
        for seq in fasta_list:
            for pos in range(num_sites):
                sites[pos] += seq[pos]

        return sites

    def clean_patterns(sites: list, seq_type):
        """
        Replace missing data with '-' for determination of pattern type. Lists of missing data
        symbols take from:
        http://www.iqtree.org/doc/Frequently-Asked-Questions#how-does-iq-tree-treat-gapmissingambiguous-characters
        We added '#" to represent terminal gaps.

        Parameters
        ----------
        sites : list
            Transposed alignment sites (columns)
        seq_type : [type]
            'AA' if aminoacid, 'NT' if nucleotide

        Returns
        -------
        list
            Transposed alignment sites with missing data replaced by '-'
        """
        missing = []
        if seq_type == "NT":
            missing = ["?", ".", "~", "O", "N", "X", "#"]
        elif seq_type == "AA":
            missing = ["?", ".", "~", "*", "X", "#"]
        clean = []
        for site in sites:
            site_out = ""
            for r in site:
                if r in missing:
                    site_out += "-"
                else:
                    site_out += r
            clean.append(site_out)

        return clean

    def pattern_type(pattern, seq_type):
        pattern = pattern.replace("-", "").upper()
        if not pattern: return "constant"
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
    }

    aln_type = fasta_type(fasta_path)
    if aln_type == "invalid": return stats

    fasta_dict = fasta_to_dict(fasta_path)
    stats["sequences"] = len(fasta_dict)
    stats["samples"] = num_samples(fasta_dict)
    stats["avg_copies"] = stats["sequences"] / stats["samples"]

    num_sites = aligned_length(fasta_dict)
    if not num_sites: return stats
    stats["sites"] = num_sites

    sites = transpose_aln(mark_terminal_gaps(fasta_dict), num_sites)
    ids_combos = [site_pairwise_identity(site, aln_type) for site in sites]
    ids_combos = [idc for idc in ids_combos if idc is not None]
    identities = 0
    combos = 0
    for idc in ids_combos:
        identities += idc[0] * idc[1]
        combos += idc[1]
    stats["avg_pid"] = round(identities / combos, 5)

    sites = clean_patterns(sites, aln_type)
    for site in sites:
        stats[pattern_type(site, aln_type)] += 1
    stats["informativeness"] = round(stats["informative"] / stats["sites"] * 100, 5)
    stats["uninformative"] = stats["constant"] + stats["singleton"]

    stats["patterns"] = len(set(sites))

    sites_concat = "".join(sites)
    if len(sites_concat):
        stats["missingness"] = round(sites_concat.count("-") / len(sites_concat) * 100, 5)

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
        header = f'{name} {fasta_dict[name]["description"]}'.strip()
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
            header = f'{name} {fasta_dict[name]["description"]}'.strip().split()

            # MEGAHIT header example: 'k157_0 flag=1 multi=2.0000 len=274'
            new_header = "_".join([
                "NODE", header[0].split("_")[1],
                "length", header[3].split("=")[1],
                "cov", header[2].split("=")[1],
                "k", header[0].split("_")[0][1:],
                "flag", header[1].split("=")[1]
            ])
            fasta_dict_spades_headers[new_header] = {
                "description": "",
                "sequence": fasta_dict[name]["sequence"]
            }
        return fasta_dict_spades_headers, assembler

    elif assembler == "skesa":
        for name in fasta_dict:
            header = name.strip().split("_")

            # SKESA header example: 'Contig_5_7.08923', if circular: 'Contig_5_7.08923_Circ'
            new_header = [
                "NODE", header[1],
                "length", f'{len(fasta_dict[name]["sequence"])}',
                "cov", header[2]
            ]
            if "_Circ" in header:
                new_header += ["circular"]
            new_header = "_".join(new_header)
            fasta_dict_spades_headers[new_header] = {
                "description": "",
                "sequence": fasta_dict[name]["sequence"]
            }
        return fasta_dict_spades_headers, assembler

    else:
        return fasta_dict, assembler


def scipio_yaml_to_dict(
        yaml_path, min_score, min_identity, min_coverage,
        marker_type, transtable, max_paralogs, predict
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
            "    prot_len",   "    prot_start", "    prot_end",
            "    target_len", "    dna_start",  "    dna_end",
            "    matches",    "    mismatches", "    undetermined",
            "    unmatched",  "    additional", "    score",
            "        nucl_start", "        nucl_end",
            "        dna_start",  "        dna_end",
            "        prot_start", "        prot_end", "        overlap",
            "          - nucl_start", "            nucl_end",
            "            dna_start",  "            dna_end",
            "            prot_start", "            prot_end",
        ]

        yaml = {}
        with open(yaml_path, "rt") as yaml_in:
            protein = ""
            hit_num = ""
            for line in yaml_in:
                line = line.strip("\n")
                if all([not line.startswith("#"), line !=  "---", line]):
                    if not line.startswith(" "):
                        protein = line.split("_(")[0] if "_(" in line else line.split(":")[0]
                        hit_num = int(line.split("_(")[1].replace("):", "")) if "_(" in line else 0
                    elif ": " in line:
                        k, v = line.split(": ")[0], ": ".join(line.split(": ")[1:])
                        if v in nones: continue
                        if k in ints:
                            try:
                                v = int(v)
                            except ValueError:
                                pass
                        else:
                            if k == "    target" and " " in v: v = v.strip("'").split()[0]
                            if v.startswith("'"): v = v.strip("'").replace(" ", "")
                            if v.startswith("["): v = [int(n) for n in v[1:-1].split(", ")]
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
                            yaml[protein][hit_num][-1][
                                "matchings"].append(copy.deepcopy(matching))
                            yaml[protein][hit_num][-1][
                                "matchings"][-1][mkeys[k]] = v
                        elif k in mkeys:
                            yaml[protein][hit_num][-1][
                                "matchings"][-1][mkeys[k]] = v
                        elif k == "          - nucl_start":
                            yaml[protein][hit_num][-1][
                                "matchings"][-1][
                                "seqshifts"].append(copy.deepcopy(seqshift))
                            yaml[protein][hit_num][-1][
                                "matchings"][-1][
                                "seqshifts"][-1][skeys[k]] = v
                        elif k in skeys:
                            yaml[protein][hit_num][-1][
                                "matchings"][-1][
                                "seqshifts"][-1][skeys[k]] = v
                    elif line.startswith("          - "):
                        mm = int(line.replace("          - ", ""))
                        yaml[protein][hit_num][-1][
                            "matchings"][-1][
                            "mismatchlist"].append(mm)

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
                if pos in inserts: seq_out += inserts[pos]
        else:
            seq_out = mat["seq"]
        if seq_out == "": translation_out = ""

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
        rfs  = {i: translate(seq_nt, gencode, frame=i, start_as_M=False) for i in [1,2,3]}
        alns = {i: align_prots(rfs[i], seq_aa, "gapless") for i in rfs}
        rf, aln = max(alns.items(), key=(lambda x: x[1]["match_rate"])) # By GS
        aln["rf"]    = rf
        aln["lead"]  = rf - 1
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
            if (i > 0 and mod["mat_types"][i-1] in non_coding
                and mod["hit_ids"][i] == mod["hit_ids"][i-1]):
                mod["mat_nt"][i-1] += mod["mat_nt"][i][:lead]
        # Remove initial extra aminoacid from translation if needed
        if aln["s1_start"] > 0:
            mod["ref_starts"][i] += aln["s1_start"]
            mod["mat_aa"][i] = mod["mat_aa"][i][aln["s1_start"]:]

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
            if (i < mod_len - 1 and mod["mat_types"][i+1] in non_coding
                and mod["hit_ids"][i] == mod["hit_ids"][i+1]):
                mod["mat_nt"][i+1] = f'{mod["mat_nt"][i][-trail:]}{mod["mat_nt"][i+1]}'
        # Remove final extra aminoacid from translation if needed
        if aln["s1_end"] < aln["s2_end"]:
            mod["ref_ends"][i] -= aln["s2_end"] - aln["s1_end"]
            mod["mat_aa"][i] = mod["mat_aa"][i][:aln["s1_end"]]

        return mod

    def fix_model(mod: dict, gencode: dict):
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

        Returns
        -------
        dict
            Modified gene model
        """
        mod_len = len(mod["mat_types"])

        # 1. Best case scenario: gene model translates exactly as Scipio and doesn't need fixing
        cds_nt, cds_aa = concat_cds(mod)
        if translate(cds_nt, gencode, frame=1, start_as_M=False) == cds_aa:
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
        if translate(cds_nt, gencode, frame=1, start_as_M=False) == cds_aa:
            return mod

        # 3. Fix dangling nucleotides in exons when gaps are found, fixes coordinates as well
        for i in range(mod_len):
            if mod["mat_types"][i] == "exon":
                if "gap before" in mod["mat_notes"][i]:
                    mod = remove_leading(mod, i, mat_alns[i]["lead"], mat_alns[i])
                    mat_alns[i] = align_exon_nt_to_scipio_aa(
                        mod["mat_nt"][i], gencode, mod["mat_aa"][i]
                    )
                if "gap after" in mod["mat_notes"][i]:
                    mod = remove_trailing(mod, i, mat_alns[i]["trail"], mat_alns[i], mod_len)
                    mat_alns[i] = align_exon_nt_to_scipio_aa(
                        mod["mat_nt"][i], gencode, mod["mat_aa"][i]
                    )
        cds_nt, cds_aa = concat_cds(mod)
        if translate(cds_nt, gencode, frame=1, start_as_M=False) == cds_aa:
            return mod

        # 4. Remove extra dangling bases that break reading frame
        last_i = None
        for i in range(mod_len):
            if mod["mat_types"][i] == "exon":
                # Remove extra leading nucleotides/aminoacids from first exon
                if last_i is None:
                    if mat_alns[i]["lead"] > 0:
                        mod = remove_leading(mod, i, mat_alns[i]["lead"], mat_alns[i])
                # Continue with remaining exons
                else:
                    if mat_alns[last_i]["trail"] + mat_alns[i]["lead"] in [0, 3]:
                        pass
                    elif mat_alns[last_i]["trail"] == 0 and mat_alns[i]["lead"] in [1, 2]:
                        mod = remove_leading(mod, i, mat_alns[i]["lead"], mat_alns[i])
                    elif mat_alns[last_i]["trail"] == 1 and mat_alns[i]["lead"] in [0, 1]:
                        mod = remove_trailing(mod, last_i, 1, mat_alns[last_i], mod_len)
                        mod = remove_leading(mod, i, mat_alns[i]["lead"], mat_alns[i])
                    elif mat_alns[last_i]["trail"] == 2:
                        if mat_alns[i]["lead"] == 0:
                            mod = remove_trailing(mod, last_i, 2, mat_alns[last_i], mod_len)
                        elif mat_alns[i]["lead"] == 2:
                            # If last exon has an extra terminal aminoacid
                            if (mat_alns[last_i]["s1_end"]
                                < mat_alns[last_i]["s2_end"]):
                                # Trim a single base from the start of current exon
                                mod = remove_leading(mod, i, 1, mat_alns[i])
                            # If last exon has two extra bp but not extra amino
                            else:
                                # Trim a single base from the end of last exon
                                mod = remove_trailing(mod, last_i, 1, mat_alns[last_i], mod_len)
                mat_alns[i] = align_exon_nt_to_scipio_aa(mod["mat_nt"][i], gencode, mod["mat_aa"][i])
                last_i = i

        cds_nt, cds_aa = concat_cds(mod)
        cds_nt_translated = translate(cds_nt, gencode, frame=1, start_as_M=False)
        if cds_nt_translated == cds_aa:
            return mod
        else:
            aln = align_prots(cds_nt_translated, cds_aa, "gapless")
            mismatches = (len(aln["mismatches"])
                          - sum(aln["s2_aln"][i] == "X" for i in aln["mismatches"]))
            if mismatches <= set_a.SCIPIO_MAX_MISMATCHES:
                return mod
            else:
                return None

    def translate_between_exons(mod: dict, i: int, gencode: dict, predict: bool, min_identity: float):
        if (len(mod["mat_nt"][i]) > 2
            and mod["mat_types"][i-1] == mod["mat_types"][i+1] == "exon"
            and mod["hit_ids"][i-1] == mod["hit_ids"][i] == mod["hit_ids"][i+1]):

            # Align previous exon to its translation to find out its trailing bases
            prev_aln = align_exon_nt_to_scipio_aa(mod["mat_nt"][i-1], gencode, mod["mat_aa"][i-1])
            # Align next exon to to its translation to find out its leading bases
            next_aln = align_exon_nt_to_scipio_aa(mod["mat_nt"][i+1], gencode, mod["mat_aa"][i+1])
            # Take trailing bases from previous and leading bases from next and add to current chunk
            current_chunk = ""
            if prev_aln["trail"] > 0:
                current_chunk += mod["mat_nt"][i-1][-prev_aln["trail"]:]
            current_chunk += mod["mat_nt"][i]
            current_chunk += mod["mat_nt"][i+1][:next_aln["lead"]]

            # 'seq_chunks', 'lead', and 'trail' only get modified when translation is acceptable
            seq_chunks = None
            lead = None
            trail = None

            if mod["mat_types"][i] == "gap":
                # Attempt translation of current chunk and then alignment to unmatched segment of
                # reference protein ('prot_gap'):
                # ref_seq[0] = protein sequence, ref_seq[1] = starting coordinate
                ref_seq = mod["ref_seqs"][mod["hit_ids"][i]]
                prot_gap = ref_seq[0][max(mod["ref_ends"][i-1] - ref_seq[1], 0) :
                                      min(mod["ref_starts"][i+1] - ref_seq[1], len(ref_seq[0]))]
                len_gap = len(prot_gap)
                max_gap_size = mod["ref_size"] * set_a.SCIPIO_MAX_GAP_AS_REF_PROP
                min_gap_identity = (min_identity / 100) - set_a.SCIPIO_MAX_GAP_DELTA_IDENTITY

                if len_gap > 0:
                    len_chunk = len(current_chunk) // 3
                    overlap = min(len_gap, len_chunk) / max(len_gap, len_chunk)
                    # Don't align if 'current_chunk' is too long compared to the unmatched protein
                    if len_chunk <= max_gap_size and len_gap * len_chunk <= set_a.RECURSION_LIMIT:
                        rfs  = {i: translate(current_chunk, gencode, frame=i, start_as_M=False)
                                for i in [1,2,3]}
                        # Global alignment if they overlap at least 80%
                        if overlap >= 0.8:
                            alns = {i: align_prots(rfs[i], prot_gap, "nw") for i in rfs}
                        # Local alignment if the overlap is smaller
                        else:
                            alns = {i: align_prots(rfs[i], prot_gap, "sw") for i in rfs}
                        # Prioritize filling gaps without skipping leading or trailing bases

                        if (len(current_chunk) % 3 == 0 and not "*" in rfs[1]
                            and alns[1]["match_rate"] >= min_gap_identity):
                            lead, trail = 0, 0
                            seq_chunks = ["", mod["mat_nt"][i].upper()]
                            mod["mismatches"] += [pos+mod["ref_ends"][i-1]
                                                  for pos in alns[1]["mismatches"]]
                        # Rank reading frames by their match rate penalized by the number of stops
                        else:
                            for rf in alns:
                                alns[rf]["match_rate"] *= (set_a.SCIPIO_STOP_PENALTY
                                                           ** alns[rf]["s1_aln"].count("*"))
                            rf, aln = max(alns.items(), key=(lambda x: x[1]["match_rate"]))
                            if (aln["match_rate"] >= min_gap_identity
                                and not "*" in aln["s1_aln"]):
                                lead  = rf - 1
                                trail = (len(current_chunk) - lead) % 3
                                if trail > 0:
                                    seq_chunks = [f'{mod["mat_nt"][i][:lead].lower()}',
                                                  f'{mod["mat_nt"][i][lead:-trail].upper()}',
                                                  f'{mod["mat_nt"][i][-trail:].lower()}']
                                else:
                                    seq_chunks = [f'{mod["mat_nt"][i][:lead].lower()}',
                                                  f'{mod["mat_nt"][i][lead:].upper()}']
                                mod["mismatches"] += [pos+mod["ref_ends"][i-1]
                                                      for pos in aln["mismatches"]]

            if predict and mod["mat_types"][i] == "intron?":
                rf1 = translate(current_chunk, gencode, frame=1, start_as_M=False)
                if len(current_chunk) % 3 == 0 and not "*" in rf1:
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
                mod["ref_starts"][i] = mod["ref_ends"][i-1]
                mod["ref_ends"][i] = mod["ref_starts"][i+1]
                mod["hit_starts"][i] = mod["hit_ends"][i-1] + lead
                mod["hit_ends"][i] = mod["hit_starts"][i+1] - trail
                mod["mat_types"][i] += "/exon"
                if mod["mat_types"][i] == "gap":
                    mod["mat_notes"][i] += "/filled"
                elif mod["mat_types"][i] == "intron?":
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
                    if (mod["hit_ids"][last_exon] == mod["hit_ids"][i]
                        and mod["hit_starts"][i] <= mod["hit_ends"][last_exon]):
                        mod["ref_ends"][last_exon] = mod["ref_ends"][i]
                        mod["hit_ends"][last_exon] = mod["hit_ends"][i]
                        mod["ref_starts"][i] = None
                        mod["ref_ends"][i]   = None
                        mod["hit_starts"][i] = None
                        mod["hit_ends"][i]   = None
                    else:
                        last_exon = i

        return mod

    def check_gaps_and_concat_seqs(mod: dict, gencode: dict, predict: bool, min_identity: float):
        # If the model doesn't end with 'downstream' or 'exon' trim it until last piece is an 'exon'
        if mod["mat_types"][-1] not in ["exon", "downstream"]:
            for i in reversed(range(len(mod["mat_types"]))):
                if mod["mat_types"][i] == "exon":
                    for data in ["ref_starts", "ref_ends", "hit_ids", "hit_starts", "hit_ends",
                                 "mat_types", "mat_notes", "mat_nt", "mat_aa"]:
                        mod[data] = mod[data][:i+1]
                    break

        # Proceed with the checkup and concatenation
        last_exon = None
        for i in range(len(mod["mat_types"])):
            if mod["mat_types"][i] == "upstream":
                mod["seq_flanked"] += mod["mat_nt"][i].lower()
            if mod["mat_types"][i] == "exon":
                if set_a.FILL_GAP_WITH_X:
                    if last_exon is not None: # if this is not the first exon
                        gap_len = mod["ref_starts"][i] - mod["ref_ends"][last_exon]
                        if (("gap before" in mod["mat_notes"][i]
                            and "gap after" in mod["mat_notes"][last_exon])
                            or gap_len >= set_a.SCIPIO_MIN_GAP_LEN_TO_X):
                            if "/exon" not in mod["mat_types"][i-1]:
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
                if (mod["mat_types"][i-1] == "stopcodon"
                    and mod["mat_nt"][i-1] == mod["mat_nt"][i][:3]):
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
                coord = f'{abs(starts[i])}-{abs(ends[i])}'
                if last_exon is None:
                    coords += coord
                else:
                    if ids[i] == ids[last_exon]:
                        coords += f',{coord}'
                    else:
                        coords += f'\n{coord}'
                last_exon = i

        return coords

    def parse_model(
            yaml_mod: dict, protein_name: str, marker_type: str, gencode: dict, predict: bool,
            min_identity: float
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
            "ref_name":    "",    # full name of protein sequence used as reference
            "ref_seqs":    {},    # segments of reference protein matched, with start and end coords
            "ref_size":    0,     # length of protein sequence used as reference
            "ref_coords":  "",    # reference protein matching intervals
            "ref_starts":  [],    # 'prot_start' per accepted exon matching
            "ref_ends":    [],    # 'prot_end' per accepted exon matching
            "hit_ids":     [],    # NUC, PTD, or MIT followed by Scipio's hit 'ID'(s)
            "hit_contigs": [],    # accepted 'target' names (= contig names in assembly)
            "hit_coords":  "",    # matching intervals per 'target'
            "hit_starts":  [],    # 'dna_start' for each accepted exon
            "hit_ends":    [],    # 'dna_end' for each accepted exon
            "mat_types":   [],    # 'type' of each 'matching' across 'target'(s)
            "mat_notes":   [],    # info about the 'matching', e.g.: exon followed by gap
            "mat_nt":      [],    # 'seq' in nucleotides for each accepted 'matching'
            "mat_aa":      [],    # 'translation' in aminoacids for each accepted exon
            "strand":      [],    # 'strand' matched per 'target'
            "matches":     [],    # number of aminoacid matches per 'target'
            "mismatches":  [],    # number of aminoacid mismatches per 'target'
            "coverage":    0.0,   # (matches + mismatches) / ref_size * 100
            "identity":    0.0,   # matches / (matches + mismatches) * 100
            "score":       0.0,   # (matches - mismatches) / ref_size
            "lwscore":     0.0,   # score * (length AA / max length AA across refs)
            "gapped":      False, # recovered protein has gaps with respect to the reference
            "seq_flanked": "",    # concatenation of 'mat_nt'
            "seq_gene":    "",    # 'seq_flanked' without 'upstream' or 'downstream' nucleotides
            "seq_nt":      "",    # concatenation of 'mat_nt' only for exons
            "seq_aa":      "",    # concatenation of 'mat_aa'
        }
        add_separator = False

        for ctg in yaml_mod: # ctg = contig (gene models can span several contigs)
            if ctg["prot_start"] < ctg["prot_end"]: # only accept protein stretches longer than 0 bp
                ref_starts = []
                ref_ends =   []
                hit_ids =    []
                hit_starts = []
                hit_ends =   []
                mat_types =  []
                mat_notes =  []
                mat_nt =     []
                mat_aa =     []
                mismatches = []
                has_overlap =  False
                hit_id = f'{marker_type}{ctg["ID"]}'
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
                    mat_nt.append(set_a.SCIPIO_CONTIG_SEPARATOR)
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
                        if (current_exon == 1
                            and "gap to querystart" in ctg["reason"]):
                            note += "gap before/"
                        if (current_exon == 1
                            and "gap to previous hit" in ctg["reason"]):
                            note += "gap before/"
                        if (current_exon == exons_in_ctg
                            and "gap to next hit" in ctg["reason"]):
                            note += "gap after/"
                        if (current_exon == exons_in_ctg
                            and "gap to queryend" in ctg["reason"]):
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
                            if mat["overlap"]: has_overlap = True
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
                    mod["ref_name"]         = protein_name
                    mod["ref_seqs"][hit_id] = (ctg["prot_seq"],
                                               ctg["prot_start"],
                                               ctg["prot_end"])
                    mod["ref_size"]         = ctg["prot_len"]
                    mod["ref_starts"]      += ref_starts
                    mod["ref_ends"]        += ref_ends
                    mod["hit_ids"]         += hit_ids
                    mod["hit_contigs"]     += [ctg["target"]]
                    mod["hit_starts"]      += hit_starts
                    mod["hit_ends"]        += hit_ends
                    mod["mat_types"]       += mat_types
                    mod["mat_notes"]       += mat_notes
                    mod["mat_nt"]          += mat_nt
                    mod["mat_aa"]          += mat_aa
                    mod["strand"]          += [ctg["strand"]]
                    mod["matches"]         += [ctg["matches"]]
                    mod["mismatches"]      += mismatches
                    add_separator = False if has_overlap else True

        mod = fix_model(mod, gencode)
        if not mod: return None

        mod = check_gaps_and_concat_seqs(mod, gencode, predict, min_identity)
        prot_len_predicted = sum(mod["hit_ends"][i] - mod["hit_starts"][i]
                                 for i in range(len(mod["mat_types"]))
                                 if "predicted" in mod["mat_notes"][i])
        mod = merge_adjoining_exons(mod)
        mod["ref_coords"]  = concat_coords(mod["hit_ids"], mod["ref_starts"], mod["ref_ends"])
        mod["hit_coords"]  = concat_coords(mod["hit_ids"], mod["hit_starts"], mod["hit_ends"])
        mod["hit_ids"]     = "\n".join(list(dict.fromkeys(list(filter(None, mod["hit_ids"])))))
        mod["hit_contigs"] = "\n".join(mod["hit_contigs"])
        mod["strand"]      = "\n".join(mod["strand"])
        prot_len           = (len(mod["seq_nt"].replace("n", "")) - prot_len_predicted) // 3
        mismatch_rate      = len(set(mod["mismatches"])) / prot_len
        mismatches         = mismatch_rate * prot_len
        matches            = prot_len - mismatches
        mod["coverage"]    = (matches + mismatches) / mod["ref_size"] * 100
        mod["identity"]    = matches / (matches + mismatches) * 100
        mod["score"]       = (matches - mismatches) / mod["ref_size"]
        # mod["lwscore"]     = mod["score"] * (matches + mismatches) / mod["ref_size"]
        # actual lwscore is computed in the filtered_models loop, for now just store length
        mod["lwscore"]     = min(prot_len, mod["ref_size"])
        mod["gapped"]      = bool("gap" in "".join(mod["mat_notes"]))

        return mod


    gencode = genetic_code(transtable)
    yaml = load_scipio_yaml(yaml_path)
    unfiltered_models = {}
    for protein in yaml: # prot = protein (reference protein)
        for yaml_model in yaml[protein]: # yaml_mod = gene model (hit, or paralog)
            model = parse_model(yaml[protein][yaml_model], protein, marker_type, gencode, predict)
            if model:
                if protein in unfiltered_models:
                    unfiltered_models[protein][yaml_model] = model
                else:
                    unfiltered_models[protein] = {yaml_model: model}

    # Separate reference protein names formatted like in the Angiosperms353.FAA file to get
    # the name of the protein cluster only when EVERY reference protein has the
    # 'set_a.REFERENCE_CLUSTER_SEPARATOR'
    refs_have_separators = True
    for protein in unfiltered_models:
        if set_a.REFERENCE_CLUSTER_SEPARATOR not in protein:
            refs_have_separators = False
            break

    # Keep track of longest recovered AA length per locus
    max_len_aa_recov = {}
    for protein in unfiltered_models:
        for model in unfiltered_models[protein]:
            if refs_have_separators:
                prot_cluster = protein.split(set_a.REFERENCE_CLUSTER_SEPARATOR)[-1]
            else:
                prot_cluster = protein
            # At this point ["lwscore"] contains only the AA length recovered (later calculated)
            if prot_cluster in max_len_aa_recov:
                if unfiltered_models[protein][model]["lwscore"] > max_len_aa_recov[prot_cluster]:
                    max_len_aa_recov[prot_cluster] = unfiltered_models[protein][model]["lwscore"]
            else:
                max_len_aa_recov[prot_cluster] = unfiltered_models[protein][model]["lwscore"]

    # Filter models by 'min_identity' and 'min_coverage'
    filtered_models = {}
    for protein in unfiltered_models:
        accepted_models = []
        for model in unfiltered_models[protein]:
            if refs_have_separators:
                prot_cluster = protein.split(set_a.REFERENCE_CLUSTER_SEPARATOR)[-1]
            else:
                prot_cluster = protein
            unfiltered_models[protein][model]["lwscore"] = (
                unfiltered_models[protein][model]["score"]
                * (unfiltered_models[protein][model]["lwscore"]
                   / max_len_aa_recov[prot_cluster])
            )
            if (unfiltered_models[protein][model]["score"] >= min_score
                and unfiltered_models[protein][model]["identity"] >= min_identity
                and unfiltered_models[protein][model]["coverage"] >= min_coverage):
                accepted_models.append(unfiltered_models[protein][model])
        if accepted_models:
            accepted_models = sorted(accepted_models, key=lambda i: i["lwscore"], reverse=True)
            if max_paralogs > -1:
                accepted_models = accepted_models[:max_paralogs+1]
            filtered_models[protein] = accepted_models
    unfiltered_models = None

    # Keep only the best hit 'models[protein][0]' with highest score and its paralogs for each
    # protein cluster. Check 'settings_assembly.py' for a more detailed description of how to
    # format your protein reference files under 'REFERENCE_CLUSTER_SEPARATOR'
    if filtered_models:
        best_models = {}
        for protein in filtered_models:
            if refs_have_separators:
                protein_cluster = protein.split(set_a.REFERENCE_CLUSTER_SEPARATOR)[-1]
            else:
                protein_cluster = protein
            if protein_cluster not in best_models:
                best_models[protein_cluster] = filtered_models[protein]
            else:
                if filtered_models[protein][0]["lwscore"] > best_models[protein_cluster][0]["lwscore"]:
                    best_models[protein_cluster] = filtered_models[protein]
        filtered_models = None
        return best_models
    else:
        return None


def blat_misc_dna_psl_to_dict(
        psl_path, target_dict, min_identity, min_coverage, marker_type, max_paralogs
):
    """
    Parse .psl from BLAT, assemble greedily the partial hits, and return the best set of hits if
    the reference contains more than a single sequence of the same type (analogous to the reference
    file Angiosperms353.FAA)
    """

    def split_blat_hit(q_size, block_sizes, q_starts, t_starts, q_end_pos, t_end_pos):
        """
        Split aligned block so no block contains an insertion as long as the query size
        """
        t_inserts = 0
        q_start, t_start = [q_starts[0]], [t_starts[0]]
        q_end, t_end = [], []
        for i in range(len(block_sizes) - 1):
            if (t_starts[i + 1] - (t_starts[i] + block_sizes[i])
                >= min(q_size * set_a.DNA_MAX_INSERT_PROP, set_a.DNA_MAX_INSERT_SIZE)):
                t_inserts += t_starts[i + 1] - (t_starts[i] + block_sizes[i])
                q_start.append(q_starts[i + 1])
                t_start.append(t_starts[i + 1])
                q_end.append(q_starts[i] + block_sizes[i])
                t_end.append(t_starts[i] + block_sizes[i])
        q_end.append(q_end_pos)
        t_end.append(t_end_pos)

        return q_start, q_end, t_start, t_end, t_inserts

    def calculate_psl_identity(
            q_end, q_start, t_end, t_start, q_inserts, matches, rep_matches, mismatches
    ):
        """
        Adapted for the case of DNA vs DNA search only from:
        https://genome-source.gi.ucsc.edu/gitlist/kent.git/raw/master/src/utils/pslScore/pslScore.pl
        """
        millibad = 0
        q_ali_size = q_end - q_start
        t_ali_size = t_end - t_start
        ali_size = q_ali_size
        if t_ali_size < q_ali_size:
            ali_size = t_ali_size
        if ali_size <= 0:
            return millibad
        size_dif = abs(q_ali_size - t_ali_size)
        insert_factor = q_inserts
        total = matches + rep_matches + mismatches
        if total != 0:
            round_away_from_zero = 3 * math.log(1 + size_dif)
            if round_away_from_zero < 0:
                round_away_from_zero = int(round_away_from_zero - 0.5)
            else:
                round_away_from_zero = int(round_away_from_zero + 0.5)
            millibad = (1000 * (mismatches + insert_factor + round_away_from_zero)) / total

        return 100.0 - millibad * 0.1

    def determine_matching_region(q_size, q_start, q_end, t_size, t_start, t_end, t_strand):
        """
        Determine if a contig matches entirely the query, or if it is a partial hit. In cases of
        partial hits determine if it is a split hit (proximal, middle, or distal) or if the hit
        is partial and subsumed within a larger stretch of sequence unrelated to the query
        """
        region = ""
        if q_end - q_start >= q_size * (1 - set_a.DNA_TOLERANCE_PROP):
            region = "full"
        elif (q_start <= q_size * set_a.DNA_TOLERANCE_PROP
              and (q_size - q_end) <= q_size * set_a.DNA_TOLERANCE_PROP):
            region = "full"
        elif (q_start >= q_size * set_a.DNA_TOLERANCE_PROP
              and (q_size - q_end) >= q_size * set_a.DNA_TOLERANCE_PROP):
            region = "middle"
        elif q_start <= q_size * set_a.DNA_TOLERANCE_PROP:
            region = "proximal"
        elif (q_size - q_end) <= q_size * set_a.DNA_TOLERANCE_PROP:
            region = "distal"
        else:  # hit is only partial and surrounded by a large proportion of unmatched sequence
            region = "wedged"
        # Now check if the flanks in the contig are too long to be part of a multi-part hit
        left_flank_too_long = t_start / (t_end - t_start) < set_a.DNA_TOLERANCE_PROP
        right_flank_too_long = (t_size - t_end) / (t_end - t_start) < set_a.DNA_TOLERANCE_PROP
        if region == "middle":
            if left_flank_too_long or right_flank_too_long:
                region = "wedged"
        elif ((region == "proximal" and t_strand == "+")
              or (region == "distal" and t_strand == "-")):
            if right_flank_too_long: region = "wedged"
        elif ((region == "proximal" and t_strand == "-")
              or (region == "distal" and t_strand == "+")):
            if left_flank_too_long: region = "wedged"
        return region

    def determine_max_overlap(contig_name):
        template = re.compile(CAPTUS_MEGAHIT_HEADER_REGEX)
        if template.match(contig_name):
            max_overlap_bp = int(contig_name.split("_")[7])
        else:
            max_overlap_bp = set_a.DNA_MAX_OVERLAP_BP
        return max_overlap_bp

    def greedy_assembly_partial_hits(hits_list, max_overlap_bp, max_paralogs: int):
        """
        Partial hits are assembled greedily starting by the hits with the highest 'lwscore',
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
            # Sort partial hits by their 'lwscore'
            part_hits = sorted(part_hits, key=lambda i: i["lwscore"], reverse=True)
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

        assembly_paths = ([[full_hits[f]] for f in range(len(full_hits))] + sorted_paths)

        # Search FASTA assembly input ('target_dict') and extract/stitch needed sequence fragments
        assembly = extract_and_stitch_edges(assembly_paths, max_overlap_bp, max_paralogs)
        assembly_paths = None
        return assembly

    def pair_is_compatible(h1, h2, max_overlap_bp):
        if (max(h1["identity"], h2["identity"]) * (1 - set_a.DNA_TOLERANCE_PROP)
            > min(h1["identity"], h2["identity"])):
            return False
        overlap = min(h1["q_end"][-1], h2["q_end"][-1]) - max(h1["q_start"][0], h2["q_start"][0])
        if overlap < 1:
            return True
        if (h1["match_len"] * (1 - set_a.DNA_TOLERANCE_PROP) <= overlap
            or h2["match_len"] * (1 - set_a.DNA_TOLERANCE_PROP) <= overlap):
            return False
        return bool(overlap <= max_overlap_bp)

    def extract_and_stitch_edges(assembly_paths, max_overlap_bp: int, max_paralogs: int):
        """
        Returns a list of stitched sequences with metadata, sorted by relevance: 'full' hits first
        sorted by 'lwscore', followed by assembled hits sorted by 'lwscore', and finally
        partial unassembled hits sorted by 'lwscore'.
        Assembly: For overlaps follow the coordinate system of the query and remove the overlap from
        the partial hit with the lower 'identity', for non-overlapped partial hits, concatenate hits
        intercalating them with as many Ns as indicated by gap in query
        """
        raw_assembly = []
        for path in assembly_paths:
            len_path = len(path)
            if len_path > 1:
                for i in range(1, len_path):
                    if path[i]["region"] == "proximal": path[i]["region"] = "middle"
                for i in range(len_path - 1):
                    if path[i]["region"] == "distal": path[i]["region"] = "middle"
            asm_hit = {
                "ref_name": path[0]["ref_name"],        # full name of non-coding reference sequence
                "ref_size": path[0]["ref_size"],           # length of non-coding reference sequence
                "ref_coords": join_coords(path[0]["q_start"], path[0]["q_end"]),     # coords in ref
                "hit_ids": path[0]["hit_id"],                                   # 'DNA##' or 'CLR##'
                "hit_contigs": path[0]["hit_contig"],              # contig name(s) used in assembly
                "hit_coords": join_coords(path[0]["t_start"], path[0]["t_end"]),   # coords in match
                "strand": path[0]["strand"],                                     # contigs strand(s)
                "matches": path[0]["matches"],                  # accumulated matches across targets
                "mismatches": path[0]["mismatches"],         # accumulated mismatches across targets
                "coverage": path[0]["coverage"],         # ((matches + mismatches) / ref_size) * 100
                "identity": path[0]["identity"],          # (matches / (matches + mismatches)) * 100
                "score": path[0]["score"],  # Scipio-like score as (matches - mismatches) / ref_size
                "lwscore": path[0]["lwscore"], # Scipio-like * (len matched / locus max len matched)
                "gapped": path[0]["gapped"],          # set to True when is assembly of partial hits
                "region": path[0]["region"],    # full, proximal, middle, distal with respect to ref
                # assembled sequence match plus upstream and downstream buffer
                "seq_flanked": extract_psl_sequence(target_dict, path[0],
                                                    set_a.DNA_UP_DOWN_STREAM_BP,
                                                    flanked=True),
                # assembled sequence match
                "seq_gene": extract_psl_sequence(target_dict, path[0], 0),
                "seq_nt": "",  # not used
                "seq_aa": "",  # not used
                "match_len": path[0]["match_len"],  # 'matches' + 'rep_matches' + 'mismatches'
            }

            if len(path) > 1:
                ref_starts = [list(path[0]["q_start"])]
                ref_ends = [list(path[0]["q_end"])]
                hit_ids = [path[0]["identity"]]
                match_props = [path[0]["matches"] / len(asm_hit["seq_gene"])]
                mismatch_props = [path[0]["mismatches"] / len(asm_hit["seq_gene"])]
                overlap_sum = 0
                for h in range(len(path) - 1):
                    asm_hit["hit_ids"] += f'\n{path[h+1]["hit_id"]}'
                    asm_hit["hit_contigs"] += f'\n{path[h+1]["hit_contig"]}'
                    asm_hit["hit_coords"] += (
                        f'\n{join_coords(path[h+1]["t_start"], path[h+1]["t_end"])}'
                    )
                    asm_hit["strand"] += f'\n{path[h+1]["strand"]}'
                    hit_ids.append(path[h+1]["identity"])
                    asm_hit["region"] += f',{path[h+1]["region"]}'
                    next_seq_flanked = extract_psl_sequence(target_dict, path[h+1],
                                                            set_a.DNA_UP_DOWN_STREAM_BP,
                                                            flanked=True)
                    next_seq_gene = extract_psl_sequence(target_dict, path[h+1], 0)
                    match_props.append(path[h+1]["matches"] / len(next_seq_gene))
                    mismatch_props.append(path[h+1]["mismatches"] / len(next_seq_gene))

                    overlap = path[h]["q_end"][-1] - path[h+1]["q_start"][0]
                    # Negative 'overlap' is a gap that has to be filled with 'n's
                    gap_len = abs(overlap) if overlap < 0 else 0
                    if gap_len > 0:
                        asm_hit["seq_flanked"] = stitch_contigs(
                            asm_hit["seq_flanked"], asm_hit["hit_contigs"].split("\n")[-1],
                            next_seq_flanked, path[h+1]["hit_contig"], gap_len
                        )
                        asm_hit["seq_gene"] += f'{"n" * gap_len}{next_seq_gene}'
                        ref_starts.append(list(path[h+1]["q_start"]))
                        ref_ends.append(list(path[h+1]["q_end"]))
                    else:
                        overlap_sum += overlap
                        # Ignore overlapped portion from the hit with smaller 'identity' to the ref
                        if path[h]["identity"] >= path[h+1]["identity"]:
                            asm_hit["seq_flanked"] += next_seq_flanked[overlap:]
                            asm_hit["seq_gene"] += next_seq_gene[overlap:]
                            ref_starts.append(list(path[h+1]["q_start"]))
                            ref_starts[-1][0] += overlap
                            ref_ends.append(list(path[h+1]["q_end"]))
                        else:
                            asm_hit["seq_flanked"] = (
                                f'{asm_hit["seq_flanked"][:-overlap]}{next_seq_flanked}'
                            )
                            asm_hit["seq_gene"] = f'{asm_hit["seq_gene"][:-overlap]}{next_seq_gene}'
                            ref_ends[-1][-1] -= overlap
                            ref_starts.append(list(path[h+1]["q_start"]))
                            ref_ends.append(list(path[h+1]["q_end"]))

                # Reformat the matched query coordinates for decorating the FASTA description line
                asm_hit["ref_coords"] = "\n".join([join_coords(s, e) for s, e
                                                   in zip(ref_starts, ref_ends)])
                # To avoid inflating the 'score' and 'coverage' of hits with insertions we need to
                # find out first the 'matched_len' discounting gaps and insertions to ref
                match_len = ref_ends[-1][-1] - ref_starts[0][0] - asm_hit["seq_gene"].count("n")
                asm_hit["coverage"] = match_len / asm_hit["ref_size"] * 100.0
                # Calculate the mean 'identity' of all the partial hits used in the assembled path
                asm_hit["identity"] = statistics.mean(hit_ids)
                # Recalculate the 'score' and 'lwscore' using sum of matches/mismatches from all
                # partial hits used in the assemble path
                ave_match_prop = statistics.mean(match_props)
                ave_mismatch_prop = 1 - ave_match_prop
                asm_hit["score"] = (((ave_match_prop * (match_len + overlap_sum))
                                     - (ave_mismatch_prop * (match_len + overlap_sum)))
                                    / (asm_hit["ref_size"] + overlap_sum))
                asm_hit["lwscore"] = asm_hit["score"] * (match_len / asm_hit["ref_size"])
                asm_hit["gapped"] = bool("n" in asm_hit["seq_gene"])
                asm_hit["match_len"] = match_len

            # Append hits to the global assembly 'raw_assembly'
            raw_assembly.append(dict(asm_hit))

        # Sorting hits by 'lwscore'
        raw_assembly = sorted(raw_assembly, key=lambda i: i["lwscore"], reverse=True)
        # Filter hit by 'min_coverage' and 'min_identity'
        assembly = []
        for hit in raw_assembly:
            if hit["identity"] >= min_identity and hit["coverage"] >= min_coverage:
                assembly.append(hit)
        # Filter by 'max_paralogs'
        if max_paralogs > -1:
            assembly = assembly[:max_paralogs+1]
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
        ends   = list(hit_dict["t_end"])
        region = hit_dict["region"]
        sequence = ""
        if flanked:
            for i in range(len(starts)):
            # When we extract the flanked matches and there is intervening blocks of unmatched
            # sequence we must extract those as well, for example matching CDS against a genome will
            # match the exons, in the flanked version we recover the introns as well
                sequence += fasta_dict[contig]["sequence"][starts[i]:ends[i]].upper()
                try:
                    sequence += fasta_dict[contig]["sequence"][ends[i]:starts[i+1]].lower()
                except IndexError:
                    continue
            # Now we extract the flanks, being aware of strand
            seq_len = len(fasta_dict[contig]["sequence"])
            start_flank = max((starts[0] - up_down_stream_bp), 0)
            end_flank = min((ends[-1] + up_down_stream_bp), seq_len)
            left_flank = fasta_dict[contig]["sequence"][start_flank:starts[0]].lower()
            right_flank = fasta_dict[contig]["sequence"][ends[-1]:end_flank].lower()
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
                current_trail = current_ctg_seq[i+1:]
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
        if check_len > set_a.DNA_MIN_OVERLAP_BP:
            while check_len >= set_a.DNA_MIN_OVERLAP_BP:
                if overlap_matches(current_trail[-check_len:], next_lead[0:check_len]):
                    scaffold_overlap = check_len
                    break
                check_len -= 1
        scaffold = ""
        sep_len = gap_len - (len(current_trail) + len(next_lead))
        separator = "n" * sep_len if sep_len > 0 else set_a.DNA_CONTIG_SEPARATOR

        if scaffold_overlap == 0:
            scaffold = f'{current_ctg_seq}{separator}{next_ctg_seq}'
        else:
            if current_ctg_cov >= next_ctg_cov:
                scaffold = f'{current_ctg_seq}{next_ctg_seq[scaffold_overlap:]}'
            else:
                scaffold = f'{current_ctg_seq[:-scaffold_overlap]}{next_ctg_seq}'

        return scaffold

    def overlap_matches(seq1, seq2):
        matches = 0
        seq_len = len(seq1)
        for n1, n2 in zip(seq1.upper(), seq2.upper()):
            matches += NT_PIDS["".join(sorted([n1, n2]))]
        if matches == seq_len:
            return True
        if matches / seq_len >= (1 - set_a.DNA_TOLERANCE_PROP):
            return True
        if seq_len - matches <= 1.0 and seq_len <= 10:
            return True
        return False

    raw_dna_hits = {}
    hit = {}
    hit_num = 1
    with open(psl_path) as psl_in:
        for line in psl_in:
            cols = line.split()
            matches, mismatches, rep_matches = int(cols[0]), int(cols[1]), int(cols[2])
            q_inserts, t_strand, q_name, q_size = int(cols[4]), cols[8], cols[9], int(cols[10])
            t_base_inserts, t_name, t_size = int(cols[7]), cols[13], int(cols[14])
            # When an insertion as long as the query size is detected we must attempt splitting hit
            if t_base_inserts >= min(q_size * set_a.DNA_MAX_INSERT_PROP,
                                     set_a.DNA_MAX_INSERT_SIZE):
                block_sizes = [int(s) for s in cols[18].strip(",").split(",")]
                q_starts = [int(p) for p in cols[19].strip(",").split(",")]
                t_starts = [int(p) for p in cols[20].strip(",").split(",")]
                q_start, q_end, t_start, t_end, _ = split_blat_hit(
                    q_size, block_sizes, q_starts, t_starts, int(cols[12]), int(cols[16])
                )
            else:
                q_start, q_end = [int(cols[11])], [int(cols[12])]
                t_start, t_end = [int(cols[15])], [int(cols[16])]

            coverage = (matches + rep_matches + mismatches) / q_size * 100
            identity = calculate_psl_identity(
                q_end[-1], q_start[0], t_end[-1], t_start[0],
                q_inserts, matches, rep_matches, mismatches
            )
            score = (matches + rep_matches - mismatches) / q_size
            lwscore = score * ((matches + rep_matches + mismatches) / q_size)
            region = determine_matching_region(q_size, q_start[0], q_end[-1],
                                               t_size, t_start[0], t_end[-1], t_strand)
            gapped = not bool(region == "full")

            # Ignore hits with not enough coverage of the reference, or partial hits immersed in
            # larger segments of unrelated sequence
            if not (coverage < set_a.DNA_MIN_COVERAGE_BEFORE_ASSEMBLY and region == "wedged"):

                # Compose hit record with columns from the .psl line
                hit = {
                    "ref_name":   q_name,
                    "ref_size":   q_size,
                    "q_start":    q_start,                      # list of starts
                    "q_end":      q_end,                          # list of ends
                    "match_len":  matches + rep_matches + mismatches,
                    "hit_id":     f"{marker_type}{hit_num}",
                    "hit_contig": t_name,
                    "t_start":    t_start,                      # list of starts
                    "t_end":      t_end,                          # list of ends
                    "strand":     t_strand,
                    "matches":    matches + rep_matches,
                    "mismatches": mismatches,
                    "coverage":   coverage,
                    "identity":   identity,
                    "score":      score,
                    "lwscore":    lwscore,
                    "region":     region,
                    "gapped":     gapped,
                }
                hit_num += 1
                if q_name not in raw_dna_hits:
                    raw_dna_hits[q_name] = [dict(hit)]
                else:
                    raw_dna_hits[q_name].append(dict(hit))

    # Use the name of the last contig to extract kmer size if the assembly was done within Captus
    # and determine the maximum overlap tolerated between adjacent contigs for assembly
    if hit:
        max_overlap_bp = determine_max_overlap(hit["hit_contig"])
    else:
        return False

    # Assemble hits, organized by reference
    dna_hits = {}
    for dna_ref in raw_dna_hits:
        assembled_hits = greedy_assembly_partial_hits(raw_dna_hits[dna_ref],
                                                      max_overlap_bp,
                                                      max_paralogs)
        if assembled_hits:
            dna_hits[dna_ref] = assembled_hits

    # Find if reference has been formatted like Angiosperms353.FAA in order to accomodate more than
    # a single reference of the same type, we can also use 'set_a.REFERENCE_CLUSTER_SEPARATOR' to
    # recognize the name of the cluster of references of the same kind
    refs_have_separators = True
    for dna_ref in dna_hits:
        if set_a.REFERENCE_CLUSTER_SEPARATOR not in dna_ref:
            refs_have_separators = False
            break

    if dna_hits:
        # First obtain the length of the longest match across reference sequences in reference locus
        # 'match_len' is holding temporarily the matched length of each hit
        max_len_nt_recov = {}
        for dna_ref in dna_hits:
            if refs_have_separators:
                dna_ref_cluster = dna_ref.split(set_a.REFERENCE_CLUSTER_SEPARATOR)[-1]
            else:
                dna_ref_cluster = dna_ref
            for hit in dna_hits[dna_ref]:
                if dna_ref_cluster not in max_len_nt_recov:
                    max_len_nt_recov[dna_ref_cluster] = hit["match_len"]
                else:
                    if hit["match_len"] > max_len_nt_recov[dna_ref_cluster]:
                        max_len_nt_recov[dna_ref_cluster] = hit["match_len"]

        # Loop the object again to calculate the actual lwscore
        for dna_ref in dna_hits:
            if refs_have_separators:
                dna_ref_cluster = dna_ref.split(set_a.REFERENCE_CLUSTER_SEPARATOR)[-1]
            else:
                dna_ref_cluster = dna_ref
            for hit in dna_hits[dna_ref]:
                hit["lwscore"] = hit["score"] * (hit["match_len"] / max_len_nt_recov[dna_ref_cluster])
            # Sort hits from largest to smallest 'lwscore'
            dna_hits[dna_ref] = sorted(dna_hits[dna_ref], key=lambda i: i["lwscore"], reverse=True)

        # If multiple references of the same kind exist in the reference, then choose the one with
        # the best hit that has the highest 'lwscore'
        best_dna_hits = {}
        for dna_ref in dna_hits:
            if refs_have_separators:
                dna_ref_cluster = dna_ref.split(set_a.REFERENCE_CLUSTER_SEPARATOR)[-1]
            else:
                dna_ref_cluster = dna_ref
            if dna_ref_cluster not in best_dna_hits:
                best_dna_hits[dna_ref_cluster] = dna_hits[dna_ref]
            else:
                if dna_hits[dna_ref][0]["lwscore"] > best_dna_hits[dna_ref_cluster][0]["lwscore"]:
                    best_dna_hits[dna_ref_cluster] = dna_hits[dna_ref]
        dna_hits = None
        return best_dna_hits
    else:
        return False


def write_gff3(hits, marker_type, out_gff_path):

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

    if marker_type in ["NUC", "PTD", "MIT"]:
        source = urllib.parse.quote("Captus (Scipio)")
        feature_type = urllib.parse.quote(f"protein_match:{marker_type}")
    else:
        source = urllib.parse.quote("Captus (BLAT)")
        feature_type = urllib.parse.quote(f"nucleotide_match:{marker_type}")

    phase = "."

    gff = ["##gff-version 3"]
    for ref in hits:
        gff.append(f"\n# {urllib.parse.quote(ref)}")
        for h in range(len(hits[ref])):
            if len(hits[ref]) == 1:
                h_name = urllib.parse.quote(ref)
            else:
                h_name = urllib.parse.quote(f"{ref}{set_a.SEQ_NAME_SEP}{h:02}")
                gff.append(f"# {h_name}")
            ref_coords = split_coords(hits[ref][h]["ref_coords"], as_strings=True)
            hit_ids = hits[ref][h]["hit_ids"].split("\n")
            hit_contigs = hits[ref][h]["hit_contigs"].split("\n")
            hit_coords = split_coords(hits[ref][h]["hit_coords"])
            strands = hits[ref][h]["strand"].split("\n")
            score = f'{hits[ref][h]["score"]:.3f}'
            lwscore = f"""LWScore={urllib.parse.quote(f'{hits[ref][h]["lwscore"]:.3f}')}"""
            cover_pct = f"""Coverage={urllib.parse.quote(f'{hits[ref][h]["coverage"]:.2f}')}"""
            ident_pct = f"""Identity={urllib.parse.quote(f'{hits[ref][h]["identity"]:.2f}')}"""
            color = f"Color={urllib.parse.quote(set_a.GFF_COLORS[marker_type])}"
            for c in range(len(hit_coords)):
                seq_id = urllib.parse.quote(hit_contigs[c])
                strand = strands[c]
                hit_id = f"ID={hit_ids[c]}"
                name = f"Name={h_name}"
                for p in range(len(hit_coords[c])):
                    start = str(hit_coords[c][p][0])
                    end = str(hit_coords[c][p][1])
                    query = (
                        f"""Query={urllib.parse.quote(f'{hits[ref][h]["ref_name"]}:{ref_coords[c]}')}"""
                    )
                    attributes = ";".join([
                        hit_id, name, lwscore, query, cover_pct, ident_pct, color
                    ])
                    gff.append("\t".join([
                        seq_id, source, feature_type, start, end, score, strand, phase, attributes
                    ]))
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

