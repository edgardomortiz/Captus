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


# import ast
import gzip
import math
import re
import statistics
import urllib
from collections import OrderedDict
from pathlib import Path

from . import settings_assembly as set_a
from .misc import bold, dim, red

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
    "N": ["A", "C", "G", "T"], "-": [],
}

# Set of valid aminoacids, including IUPAC ambiguities
AA_IUPAC = {
    "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N",
    "O", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z", "X", "*", "-"
}

# Ambiguous aminoacids used in fuzzy translation
AA_AMBIGS = {
    ("D", "N"): "N",  # "B"
    ("E", "Q"): "Q",  # "Z"
    ("I", "L"): "L",  # "J"
}

# Since the list of nucleotides is a subset of the list of aminoacids we need a list of only the
# aminoacid letters that are not also nucleotide letters
AA_NOT_IN_NT = AA_IUPAC - set(NT_IUPAC)

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
    "base1": "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
    "base2": "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
    "base3": "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG",
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


def translate(seq, genetic_code: dict, frame=1):
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

    if frame < 0:
        seq = reverse_complement(seq)

    codons = [seq[p:p + 3].upper() for p in range(abs(frame) - 1, len(seq), 3)
              if len(seq[p:p + 3]) == 3]

    if len(codons) < 2:
        return ""

    seq_AA = ""

    starts_fuzzy = []
    for triplet in resolve_codon(codons[0]):
        if triplet in genetic_code["starts"]:
            starts_fuzzy.append(genetic_code["starts"][triplet])
        elif triplet in genetic_code["aminos"]:
            starts_fuzzy.append(genetic_code["aminos"][triplet])
        else:
            starts_fuzzy.append("X")
    seq_AA += unresolve_aminoacid(starts_fuzzy)

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


def translate_fasta_dict(in_fasta_dict: dict, genetic_code_id: int, frame="guess"):
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
    def guess_frame(seq, genetic_code):
        translations, counts_Xs, counts_stops = [], [], []
        for f in [1, 2, 3, -1, -2, -3]:
            translation = translate(seq, genetic_code, frame=f)
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
    translated_fasta = OrderedDict()
    if frame == "guess":
        for seq_name in in_fasta_dict:
            translated_fasta[seq_name] = {
                "description": in_fasta_dict[seq_name]["description"],
                "sequence": guess_frame(in_fasta_dict[seq_name]["sequence"], gc)}
    else:
        for seq_name in in_fasta_dict:
            translated_fasta[seq_name] = {
                "description": in_fasta_dict[seq_name]["description"],
                "sequence": translate(in_fasta_dict[seq_name]["sequence"], gc)}
    return translated_fasta


def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence
    """
    return "".join([REV_COMP_DICT[n] for n in seq[::-1]])


def fasta_to_dict(fasta_path, ordered=False):
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
    if ordered is True:
        fasta_out = OrderedDict()
    else:
        fasta_out = {}
    with opener(fasta_path, "rt") as fasta_in:
        for line in fasta_in:
            line = line.strip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if len(line.split()) > 1:
                    name = line[1:].split()[0]
                    desc = " ".join(line[1:].split()[1:])
                else:
                    name = line[1:].rstrip()
                    desc = ""
                fasta_out[name] = {"description": desc, "sequence": ""}
            else:
                fasta_out[name]["sequence"] += line
    return fasta_out


def dict_to_fasta(in_fasta_dict, out_fasta_path, wrap=0, append=False):
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
    if bool(in_fasta_dict):
        with opener(out_fasta_path, action) as fasta_out:
            if wrap == 0:
                for name in in_fasta_dict:
                    if in_fasta_dict[name]["description"]:
                        header = f'>{name} {in_fasta_dict[name]["description"]}'.strip()
                    else:
                        header = f">{name}".strip()
                    fasta_out.write(f'{header}\n{in_fasta_dict[name]["sequence"]}\n')
            elif wrap > 0:
                for name in in_fasta_dict:
                    if in_fasta_dict[name]["description"]:
                        header = f'>{name} {in_fasta_dict[name]["description"]}'.strip()
                    else:
                        header = f">{name}".strip()
                    sequence = in_fasta_dict[name]["sequence"]
                    wrapped = "\n".join([sequence[i:i + wrap] for i in range(0, len(sequence), wrap)])
                    fasta_out.write(f'{header}\n{wrapped}\n')
    return out_fasta_path


def is_fasta_nt(fasta_path):
    """
    Verify FASTA format and that sequence is only nucleotides, returns 'True' when the file contains
    only nucleotides or nucleotide IUPAC ambiguity codes and 'False' when it contains aminoacids.
    When other characters are found in the sequence or the files doesn't have headers it returns
    'not a FASTA'
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
                    line_set = set(line.upper().replace("U", "T"))
                    if line_set - set(NT_IUPAC):
                        if not line_set - set(AA_IUPAC):
                            has_aminos = True
                            break
                else:
                    return "not a FASTA"
    if has_headers and has_aminos:
        return False
    elif has_headers and not has_aminos:
        return True
    else:
        return "not a FASTA"


def fasta_headers_to_spades(fasta_dict, ordered=True):
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

    if ordered is True:
        fasta_dict_spades_headers = OrderedDict()
    else:
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


def scipio_yaml_to_dict(yaml_path, min_identity, min_coverage, marker_type):
    """
    Takes a Scipio's YAML output file 'yaml_path' and returns a dictionary of all the hits per
    protein and their curated sequence
    The function creates a new score based on Scipio's own score and the length of the hit
    """

    def insert_shifts(seq_in, shifts_dict):
        """
        Given a dictionary of the form {pos: "nn"}, loop the sequence and insert the values AFTER
        the index coordinate
        """
        seq_out = ""
        for pos in range(len(seq_in)):
            seq_out += seq_in[pos]
            if pos in shifts_dict:
                seq_out += shifts_dict[pos]
        return seq_out

    # Initialize variables before starting the loop along the YAML file
    raw = {}  # raw models dictionary, need to be filtered by 'min_coverage' and 'min_identity' after
    protein = ""
    hit = 0
    match_type = ""
    overlapped = False
    include = False
    # Elements within target separated by ',' between targets by '\n'
    hit_data = {
        "ref_name": "",     # full name of reference protein sequence
        "ref_size": 0,      # length of reference protein sequence
        "ref_coords": "",   # interval(s) of protein matched
        "hit_id": "",       # 'NUC', 'PTD', or 'MIT' followed by Scipio's ID(s)
        "hit_contig": "",   # contig name(s) used in assembly
        "hit_coords": "",   # matching interval(s) in hit(s)
        "strand": "",       # strand(s)
        "matches": 0,       # accumulated number of matches across targets
        "mismatches": 0,    # accumulated number of mismatches across targets
        "coverage": 0.0,    # reference coverage as ((matches + mismatches) / ref_size) * 100
        "identity": 0.0,    # identity percentage as (matches / (matches + mismatches)) * 100
        "score": 0.0,       # Scipio-like score as (matches - mismatches) / ref_size
        "lwscore": 0.0,     # Scipio-like score multiplied by ((matches + mismatches) / ref_size)
        "region": "",       # included here to match BLAT's data for non-coding extraction
        "gapped": False,    # hit contains gaps with respect to reference
        "seq_flanked": "",  # concatenation of contigs, overlaps merged or 50 'n's in between
        "seq_gene": "",     # gene sequence excluding upstream and downstream regions
        "seq_nt": "",       # concatenation of CDS in nucleotide (includes STOP codon)
        "seq_aa": "",       # concatenation of CDS in aminoacid (STOP codon excluded)
    }

    # The processing loop for the YAML file starts here
    with open(yaml_path) as yaml_in:
        for line in yaml_in:
            line = line.strip("\n")

            # Skip empty and comment lines
            if line and not line.startswith("---") and not line.startswith("#"):

                # Lines with no indentation mark the start of each model
                if not line.startswith(" "):

                    # Only paralogs end in '_(1)', '_(2)', etc.  Since the best model doesn't have
                    # a subindex we call it '0'. Extract 'protein' ignoring paralog subindex, and
                    # set 'hit' to paralog subindex or to '0' for best model
                    if "_(" in line:
                        protein = line.split("_(")[0]
                    else:
                        protein = line.split(":")[0]

                    # Add the protein hit to 'raw' dictionary
                    if protein not in raw:
                        raw[protein] = [dict(hit_data)] # first hit in list has index 0
                    else:
                        raw[protein].append(dict(hit_data)) # paralogs will have index > 0
                    raw[protein][-1]["ref_name"] = protein

                # Parse model info
                if line.startswith("  "):
                    if len(line.split(": ")) > 1:
                        value = line.split(": ")[1]
                        if value in ["~", "''"]:
                            value = ""

                    if line.startswith("  - ID:"):
                        hit_id = value

                    if line.startswith("    reason:"):
                        if "gap" in line:
                            gapped = True
                        else:
                            gapped = raw[protein][-1]["gapped"]

                    if line.startswith("    prot_len:"):
                        prot_len = int(value)

                    if line.startswith("    prot_start:"):
                        prot_start = abs(int(value))

                    if line.startswith("    prot_end:"):
                        include = False if abs(int(value)) < prot_start else True

                    if include:
                        if raw[protein][-1]["hit_id"] == "":
                            raw[protein][-1]["hit_id"] = f"{marker_type}{hit_id}"
                        else:
                            raw[protein][-1]["hit_id"] += f"\n{marker_type}{hit_id}"

                    if include:
                        raw[protein][-1]["gapped"] = gapped

                    if include:
                        raw[protein][-1]["ref_size"] = prot_len

                    if include and line.startswith("    target:"):
                        target = value.split()[0].strip("'")
                        if raw[protein][-1]["hit_contig"] == "":
                            raw[protein][-1]["hit_contig"] = target
                        else:
                            raw[protein][-1]["hit_contig"] += f"\n{target}"
                            raw[protein][-1]["ref_coords"] += "\n"
                            raw[protein][-1]["hit_coords"] += "\n"
                            if not overlapped:
                                raw[protein][-1]["seq_flanked"] += set_a.SCIPIO_CONTIG_SEPARATOR
                                raw[protein][-1]["seq_gene"] += set_a.SCIPIO_CONTIG_SEPARATOR

                    if include and line.startswith("    strand:"):
                        if raw[protein][-1]["strand"] == "":
                            raw[protein][-1]["strand"] = value.strip("'")
                        else:
                            raw[protein][-1]["strand"] += f"""\n{value.strip("'")}"""

                    if include and line.startswith("    matches:"):
                        raw[protein][-1]["matches"] += int(value)

                    if include and line.startswith("    mismatches:"):
                        raw[protein][-1]["mismatches"] += int(value)
                        raw[protein][-1]["coverage"] = (
                            ((raw[protein][-1]["matches"] + raw[protein][-1]["mismatches"])
                            / raw[protein][-1]["ref_size"]) * 100
                        )
                        try:
                            raw[protein][-1]["identity"] = (
                                (raw[protein][-1]["matches"]
                                / (raw[protein][-1]["matches"]
                                    + raw[protein][-1]["mismatches"])) * 100
                            )
                        except ZeroDivisionError:
                            raw[protein][-1]["identity"] = 0.0
                        raw[protein][-1]["score"] = (
                            (raw[protein][-1]["matches"] - raw[protein][-1]["mismatches"])
                            / raw[protein][-1]["ref_size"]
                        )
                        raw[protein][-1]["lwscore"] = (
                            raw[protein][-1]["score"] * (
                                (raw[protein][-1]["matches"] + raw[protein][-1]["mismatches"])
                                / raw[protein][-1]["ref_size"]
                            )
                        )

                    if include and line.startswith("    upstream:"):
                        raw[protein][-1]["seq_flanked"] += value

                    if include and line.startswith("    upstream_gap:"):
                        raw[protein][-1]["seq_flanked"] += value
                        # raw[protein][-1]["seq_gene"] += value

                    if include and line.startswith("      - type:"):
                        match_type = value.strip("?")
                        if match_type == "exon":
                            shifts_dict = {}

                    if include and match_type == "exon" and line.startswith("        dna_start:"):
                        exon_start = int(value)
                        match_interval = value.strip("-")

                    if include and match_type == "exon" and line.startswith("        dna_end:"):
                        match_interval += f'-{value.strip("-")}'
                        if (raw[protein][-1]["hit_coords"] == "" or
                            raw[protein][-1]["hit_coords"].endswith("\n")):
                            raw[protein][-1]["hit_coords"] += match_interval
                        else:
                            raw[protein][-1]["hit_coords"] += f",{match_interval}"

                    if include and match_type == "exon" and line.startswith("        prot_start:"):
                        prot_interval = value

                    if include and match_type == "exon" and line.startswith("        prot_end:"):
                        prot_interval += f"-{value}"
                        if (raw[protein][-1]["ref_coords"] == "" or
                            raw[protein][-1]["ref_coords"].endswith("\n")):
                            raw[protein][-1]["ref_coords"] += prot_interval
                        else:
                            raw[protein][-1]["ref_coords"] += f",{prot_interval}"

                    if include and line.startswith("        seq:"):
                        if match_type == "intron" or match_type == "gap":
                            raw[protein][-1]["seq_flanked"] += value
                            raw[protein][-1]["seq_gene"] += value
                        if match_type == "exon":
                            exon_raw_seq = value

                    if include and match_type == "exon" and line.startswith("            dna_start:"):
                        shift_start = int(value)
                    if include and match_type == "exon" and line.startswith("            dna_end:"):
                        shift_end = int(value)
                        if (shift_end - shift_start) % 3 > 0:
                            shifts_dict[shift_start - exon_start - 1] = (
                                (3 - (shift_end - shift_start)) * "N"
                            )

                    if include and line.startswith("        translation:"):
                        raw[protein][-1]["seq_aa"] += value.strip("'")
                        if shifts_dict:
                            exon_seq = insert_shifts(exon_raw_seq, shifts_dict).upper()
                        else:
                            exon_seq = exon_raw_seq.upper()
                        raw[protein][-1]["seq_flanked"] += exon_seq
                        raw[protein][-1]["seq_gene"] += exon_seq
                        raw[protein][-1]["seq_nt"] += exon_seq
                        shifts_dict = {}

                    if include and line.startswith("        overlap:"):
                        overlapped = bool(value)

                    if include and line.startswith("    stopcodon:"):
                        raw[protein][-1]["seq_nt"] += value.upper()

                    if include and line.startswith("    downstream:"):
                        raw[protein][-1]["seq_flanked"] += value

                    if include and line.startswith("    downstream_gap:"):
                        raw[protein][-1]["seq_flanked"] += value
                        # raw[protein][-1]["seq_gene"] += value

    models = {}  # for models filtered by 'min_coverage' and 'min_identity'
    for protein in raw:
        accepted_hits = []
        for hit in raw[protein]:
            if hit["identity"] >= min_identity and hit["coverage"] >= min_coverage:
                accepted_hits.append(hit)
        if accepted_hits:
            models[protein] = accepted_hits
    raw = None

    # Separate reference protein names formatted like in the Angiosperms353.FAA file to get
    # the name of the protein cluster only when EVERY reference protein has the
    # 'set_a.REFERENCE_CLUSTER_SEPARATOR'
    refs_have_separators = True
    for protein in models:
        if set_a.REFERENCE_CLUSTER_SEPARATOR not in protein:
            refs_have_separators = False
            break

    # Keep only the best hit 'models[protein][0]' with highest score and its paralogs for each
    # protein cluster. Check 'settings_assembly.py' for a more detailed description of how to
    # format your protein reference files under 'REFERENCE_CLUSTER_SEPARATOR'
    if models:
        best_models = {}
        for protein in models:
            if refs_have_separators:
                protein_cluster = protein.split(set_a.REFERENCE_CLUSTER_SEPARATOR)[-1]
            else:
                protein_cluster = protein
            if protein_cluster not in best_models:
                best_models[protein_cluster] = models[protein]
            else:
                if models[protein][0]["lwscore"] > best_models[protein_cluster][0]["lwscore"]:
                    best_models[protein_cluster] = models[protein]
        models = None
        return best_models
    else:
        return None


def blat_misc_dna_psl_to_dict(psl_path, target_dict, min_identity, min_coverage, marker_type):
    """
    Parse .psl from BLAT, assemble greedily the partial hits, and return the best set of hits if
    the reference contains more than a single sequence of the same type (analogous to the reference
    file Angiosperms353.FAA)
    """

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

    def determine_matching_region(q_size, q_start, q_end, t_size, t_start, t_end):
        """
        Determine if a contig matches entirely the query, or if it is a partial hit. In cases of
        partial hits determine if it is a split hit (proximal, middle, or distal) or if the hit
        is partial and subsumed within a larger stretch of sequence unrelated to the query
        """
        if q_end - q_start >= q_size * (1 - set_a.DNA_TOLERANCE_PROP):
            return "full"
        elif t_end - t_start >= t_size * (1 - (set_a.DNA_TOLERANCE_PROP * 2)):
            return "middle"
        elif (q_start <= q_size * set_a.DNA_TOLERANCE_PROP
              and t_end >= t_size * (1 - set_a.DNA_TOLERANCE_PROP)):
            return "proximal"
        elif (q_end >= q_size * (1 - set_a.DNA_TOLERANCE_PROP)
              and t_start <= t_size * set_a.DNA_TOLERANCE_PROP):
            return "distal"
        else:  # hit is only partial and surrounded by a large proportion of unknown sequence
            return "wedged"

    def extract_psl_sequence(
            fasta_dict, seq_name, strand, start_pos, end_pos, match_part, up_down_stream_bp
    ):
        """
        Extract sequence from a 'fasta_dict' object using BLAT's PSL coordinate style
        """
        if up_down_stream_bp > 0:
            seq_len = len(fasta_dict[seq_name]["sequence"])
            if match_part == "full" or match_part == "proximal":
                start_pos = max((start_pos - up_down_stream_bp), 0)
            if match_part == "full" or match_part == "distal":
                end_pos = min((end_pos + up_down_stream_bp), seq_len)
        if strand == "+":
            return fasta_dict[seq_name]["sequence"][start_pos:end_pos]
        elif strand == "-":
            return reverse_complement(fasta_dict[seq_name]["sequence"][start_pos:end_pos])

    def extract_and_stitch_edges(assembly_paths):
        """
        Returns a list of stitched sequences with metadata, sorted by relevance: 'full' hits first
        sorted by 'lwscore', followed by assembled hits sorted by 'lwscore', and finally
        partial unassembled hits sorted by 'lwscore'.
        Assembly: For overlaps follow the coordinate system of the query and remove the overlap from
        the partial hit with the lower 'match_id', for non-overlapped partial hits, concatenate hits
        intercalating them with as many Ns as indicated by gap in query
        """
        full_hits = []
        sewn_hits = []
        part_hits = []
        asm_num = 1
        for path in assembly_paths:
            identities = [path[0]["identity"]]
            asm_hit = {
                "ref_name": path[0]["ref_name"],        # full name of non-coding reference sequence
                "ref_size": path[0]["ref_size"],           # length of non-coding reference sequence
                "ref_coords": f'{path[0]["q_start"]}-{path[0]["q_end"]}',     # AA coords in protein
                "hit_id": path[0]["hit_id"],                                    # 'DNA##' or 'CLR##'
                "hit_contig": path[0]["hit_contig"],               # contig name(s) used in assembly
                "hit_coords": f'{path[0]["t_start"]}-{path[0]["t_end"]}',       # NT coords in match
                "strand": path[0]["strand"],                                     # contigs strand(s)
                "matches": path[0]["matches"],                  # accumulated matches across targets
                "mismatches": path[0]["mismatches"],         # accumulated mismatches across targets
                "coverage": path[0]["coverage"],         # ((matches + mismatches) / ref_size) * 100
                "identity": path[0]["identity"],          # (matches / (matches + mismatches)) * 100
                "score": path[0]["score"],  # Scipio-like score as (matches - mismatches) / ref_size
                "lwscore": path[0]["lwscore"],   # Scipio-like * ((matches + mismatches) / ref_size)
                "gapped": path[0]["gapped"],          # set to True when is assembly of partial hits
                "region": path[0]["region"],    # full, proximal, middle, distal with respect to ref
                # assembled sequence match plus upstream and downstream buffer
                "seq_flanked": extract_psl_sequence(
                    target_dict,
                    path[0]["hit_contig"],
                    path[0]["strand"],
                    path[0]["t_start"],
                    path[0]["t_end"],
                    path[0]["region"],
                    set_a.DNA_UP_DOWN_STREAM_BP
                ),
                # assembled sequence match
                "seq_gene": extract_psl_sequence(
                    target_dict,
                    path[0]["hit_contig"],
                    path[0]["strand"],
                    path[0]["t_start"],
                    path[0]["t_end"],
                    path[0]["region"],
                    0
                ),
                "seq_nt": "",  # not used
                "seq_aa": "",  # not used
            }

            if len(path) > 1:
                for h in range(len(path) - 1):
                    asm_hit["ref_coords"] += f'\n{path[h + 1]["q_start"]}-{path[h + 1]["q_end"]}'
                    asm_hit["hit_id"] += f'\n{path[h + 1]["hit_id"]}'
                    asm_hit["hit_contig"] += "\n" + path[h + 1]["hit_contig"]
                    asm_hit["hit_coords"] += f'\n{path[h + 1]["t_start"]}-{path[h + 1]["t_end"]}'
                    asm_hit["strand"] += f'\n{path[h + 1]["strand"]}'
                    asm_hit["matches"] += path[h + 1]["matches"]
                    asm_hit["mismatches"] += path[h + 1]["mismatches"]
                    asm_hit["coverage"] = min((asm_hit["coverage"] + path[h + 1]["coverage"]), 100.0)
                    identities.append(path[h + 1]["identity"])
                    asm_hit["region"] += f',{path[h + 1]["region"]}'
                    next_seq_flanked = extract_psl_sequence(
                        target_dict,
                        path[h + 1]["hit_contig"],
                        path[h + 1]["strand"],
                        path[h + 1]["t_start"],
                        path[h + 1]["t_end"],
                        path[h + 1]["region"],
                        set_a.DNA_UP_DOWN_STREAM_BP
                    )
                    next_seq_gene = extract_psl_sequence(
                        target_dict,
                        path[h + 1]["hit_contig"],
                        path[h + 1]["strand"],
                        path[h + 1]["t_start"],
                        path[h + 1]["t_end"],
                        path[h + 1]["region"],
                        0
                    )
                    overlap = path[h]["q_end"] - path[h + 1]["q_start"]

                    # Negative 'overlap' is a gap that has to be filled with 'n's
                    if overlap < 0:
                        asm_hit["seq_flanked"] += f'{"n" * abs(overlap)}{next_seq_flanked}'
                        asm_hit["seq_gene"] += f'{"n" * abs(overlap)}{next_seq_gene}'
                    else:

                        # Ignore overlapped portion from the hit with smaller protein 'identity'
                        if path[h]["identity"] >= path[h + 1]["identity"]:
                            asm_hit["seq_flanked"] += next_seq_flanked[overlap:]
                            asm_hit["seq_gene"] += next_seq_gene[overlap:]
                        else:
                            asm_hit["seq_flanked"] = (
                                f'{asm_hit["seq_flanked"][:-overlap]}{next_seq_flanked}'
                            )
                            asm_hit["seq_gene"] = f'{asm_hit["seq_gene"][:-overlap]}{next_seq_gene}'

                # Calculate the mean 'match_id' of all the partial hits used in the assembled path
                asm_hit["identity"] = statistics.mean(identities)

                # Likewise, recalculate the 'score' and 'lwscore' using sum of matches/mismatches
                # from all partial hits used in the assemble path
                asm_hit["score"] = (
                    (asm_hit["matches"] - asm_hit["mismatches"]) / asm_hit["ref_size"]
                )
                asm_hit["lwscore"] = (
                    asm_hit["score"] * (
                        (asm_hit["matches"] + asm_hit["mismatches"]) / asm_hit["ref_size"]
                    )
                )
                asm_hit["gapped"] = True if "n" in asm_hit["seq_gene"] else False
                asm_num += 1

            # Separate hits by type
            if asm_hit["region"] == "full":
                full_hits.append(dict(asm_hit))
            elif "," in asm_hit["region"]:
                sewn_hits.append(dict(asm_hit))
            else:
                part_hits.append(dict(asm_hit))

        # Sort sets of hits by 'lwscore' and concatenate in a single 'assembly' list
        if full_hits:
            full_hits = sorted(full_hits, key=lambda i: i["lwscore"], reverse=True)
        if sewn_hits:
            sewn_hits = sorted(sewn_hits, key=lambda i: i["lwscore"], reverse=True)
        if part_hits:
            part_hits = sorted(part_hits, key=lambda i: i["lwscore"], reverse=True)
        raw_assembly = full_hits + sewn_hits + part_hits

        # Final filtering by 'min_coverage' and 'min_identity'
        assembly = []
        for hit in raw_assembly:
            if hit["identity"] >= min_identity and hit["coverage"] >= min_coverage:
                assembly.append(hit)
        raw_assembly = None
        return assembly

    def greedy_assembly_partial_hits(hits_list, max_overlap_bp):
        """
        TODO: check this
        Start with the most 'proximal' hit and look for 'middle' or 'distal' hits that are at most
        overlapped by min('tolerance_prop'*'q_size', max_overlap_bp) and with 'match_id' >= proximal
        hit's 'match_id'*'tolerance_prop'. Contruct a list of paths to proceed with the function
        'extract_and_stitch_edges' to finally return a list of ranked assembled sequences to
        repopulate the dictionary of 'dna_hits'
        """
        full_hits, part_hits = [], []
        compatible_pairs = []
        all_paths, last_path = [], []
        contigs_used, valid_paths = [], []
        hit_paths = []
        for hit in hits_list:
            if hit["region"] == "full":
                full_hits.append(hit)
            else:
                part_hits.append(hit)

        # Sort partial hits by their start coordinate of the match to the reference
        part_hits = sorted(part_hits, key=lambda i: i["q_start"])

        # Find pair of campatible hits, not too much overlap, not too different id. to the reference
        for h1 in range(len(part_hits)):
            for h2 in range(h1 + 1, len(part_hits)):
                if (part_hits[h1]["q_end"] - max_overlap_bp <= part_hits[h2]["q_start"] and
                    max(part_hits[h1]["identity"], part_hits[h2]["identity"])
                    * set_a.DNA_TOLERANCE_PROP
                    <= min(part_hits[h1]["identity"], part_hits[h2]["identity"])):
                    compatible_pairs.append((part_hits[h1]["hit_id"], part_hits[h2]["hit_id"]))

        # Find concatenation paths greedily, starting by the hits closest to the first query start
        if compatible_pairs:
            for p1 in range(len(compatible_pairs)):
                if last_path == []:
                    last_path = [compatible_pairs[p1][0], compatible_pairs[p1][1]]
                for p2 in range(p1 + 1, len(compatible_pairs)):
                    if last_path[-1] == compatible_pairs[p2][0]:
                        last_path.append(compatible_pairs[p2][1])
                all_paths.append(last_path)
                last_path = []

        # Filter redundant paths so every hit is used only once
        if all_paths:
            contigs_used, valid_paths = list(all_paths[0]), list([all_paths[0]])
            for path in all_paths:
                if set(path) - set(contigs_used) == set(path):
                    valid_paths.append(path)
                    contigs_used += path

        # Recover all information associated with the hit names
        if valid_paths:
            for path in valid_paths:
                path_out = []
                for step in path:
                    for hit in part_hits:
                        if step == hit["hit_id"]:
                            path_out.append(hit)
                hit_paths.append(path_out)

        # Put together all the paths to proceed with the assembly, simple full hits first, then
        # the assembled partial hits, and finally the orphaned partial hits
        assembly_paths = (
            [[full_hits[f]] for f in range(len(full_hits))] + hit_paths
            + [[part_hits[p]] for p in range(len(part_hits))
               if part_hits[p]["hit_id"] not in contigs_used]
        )

        # Search FASTA assembly input ('target_dict') and extract/stitch needed sequence fragments
        assembly = extract_and_stitch_edges(assembly_paths)
        assembly_paths = None
        return assembly

    def determine_max_overlap(contig_name):
        template = re.compile(CAPTUS_MEGAHIT_HEADER_REGEX)
        if template.match(contig_name):
            max_overlap_bp = int(contig_name.split("_")[7])
        else:
            max_overlap_bp = set_a.DNA_MAX_OVERLAP_BP
        return max_overlap_bp

    raw_dna_hits = {}
    hit_num = 1
    with open(psl_path) as psl_in:
        for line in psl_in:
            cols = line.split()
            matches, mismatches, rep_matches = int(cols[0]), int(cols[1]), int(cols[2])
            q_inserts, t_strand = int(cols[4]), cols[8]
            q_base_inserts, t_base_inserts = int(cols[5]), int(cols[7])
            q_name, q_size, q_start, q_end = cols[9], int(cols[10]), int(cols[11]), int(cols[12])
            t_name, t_size, t_start, t_end = cols[13], int(cols[14]), int(cols[15]), int(cols[16])

            coverage = ((q_end - q_start - q_base_inserts - t_base_inserts) / q_size) * 100.0
            identity = calculate_psl_identity(q_end, q_start, t_end, t_start, q_inserts,
                                              matches, rep_matches, mismatches)
            score = ((matches + rep_matches) - mismatches) / q_size
            lwscore = score * ((matches + rep_matches + mismatches) / q_size)
            region = determine_matching_region(q_size, q_start, q_end, t_size, t_start, t_end)
            gapped = False if region == "full" else True

            # Ignore hits with not enough coverage of the reference, or partial hits immersed in
            # larger segments of unrelated sequence
            if coverage >= set_a.DNA_MIN_COVERAGE_BEFORE_ASSEMBLY and region != "wedged":

                # Compose hit record with columns from the .psl line
                hit = {
                    "ref_name": q_name,
                    "ref_size": q_size,
                    "q_start": q_start,
                    "q_end": q_end,
                    "hit_id": f"{marker_type}{hit_num}",
                    "hit_contig": t_name,
                    "t_start": t_start,
                    "t_end": t_end,
                    "strand": t_strand,
                    "matches": matches + rep_matches,
                    "mismatches": mismatches,
                    "coverage": coverage,
                    "identity": identity,
                    "score": score,
                    "lwscore": lwscore,
                    "region": region,
                    "gapped": gapped
                }
                hit_num += 1
                if q_name not in raw_dna_hits:
                    raw_dna_hits[q_name] = [dict(hit)]
                else:
                    raw_dna_hits[q_name].append(dict(hit))

    # Use the name of the last contig to extract kmer size if the assembly was done within Captus
    # and determine the maximum overlap tolerated between adjacent contigs for assembly
    max_overlap_bp = determine_max_overlap(hit["hit_contig"])

    # Assemble hits, organized by reference
    dna_hits = {}
    for dna_ref in raw_dna_hits:
        assembled_hits = greedy_assembly_partial_hits(raw_dna_hits[dna_ref], max_overlap_bp)
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

    # If multiple references of the same kind exist in the reference, then choose the one with the
    # best hit that has the highest 'lwscore'
    if dna_hits:
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
                h_name = urllib.parse.quote(f"{ref}|{h:02}")
                gff.append(f"# {h_name}")
            ref_coords = split_coords(hits[ref][h]["ref_coords"], as_strings=True)
            hit_ids = hits[ref][h]["hit_id"].split("\n")
            hit_contigs = hits[ref][h]["hit_contig"].split("\n")
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

