#!/usr/bin/env python3
"""
Copyright 2020 Edgardo M. Ortiz (e.ortiz.v@gmail.com)
https://github.com/edgardomortiz/Captus

Captus' version is stored here in a separate file so it can exist in only one place.

This file is part of Captus. Captus is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Captus is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Captus. If
not, see <http://www.gnu.org/licenses/>.
"""

import math
import statistics

REV_COMP_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                 'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
                 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                 'D': 'H', 'H': 'D', 'N': 'N',
                 'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
                 'k': 'm', 'm': 'k', 'b': 'v', 'v': 'b',
                 'd': 'h', 'h': 'd', 'n': 'n',
                 '.': '.', '-': '-', '?': '?'}

def reverse_complement(seq):
    return "".join([REV_COMP_DICT[n] for n in seq[::-1]])


def fasta_to_dict(fasta_path):
    """
    Turns a FASTA file given with 'fasta_path' into a dictionary. For example, for the sequence:
        >k157_0 flag=1 multi=2.0000 len=274
        ATATTGATATTTCATAATAATAGTTTTTGAACTAAAAAGAAATTTTTCCTCCAATTATGTGGG
    Returns the dictionary:
        {'k157_0' : {'description': 'flag=1 multi=2.0000 len=274', 
                     'sequence': 'ATATTGATATTTCATAATAATAGTTTTTGAACTAAAAAGAAATTTTTCCTCCAATTATGTGGG'}}
    """
    fasta_out = {}
    if fasta_path.endswith(".gz"):
        opener = gzip.open
    else:
        opener = open
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


def blat_dna_psl_to_dict(psl_path, target_dict, query_dict, min_coverage, best_hits_only):
    """
    Parse .psl from BLAT to extract best hits, greedily completing coverage of the reference, 
    scoring length, then identity to reference, only considering partial hits when they are at the 
    ends of their respective contigs
    """

    def calc_psl_identity(q_end, q_start, t_end, t_start, q_inserts, matches, rep_matches, 
                          mismatches):
        """
        Adapted from:
        https://genome-source.gi.ucsc.edu/gitlist/kent.git/raw/master/src/utils/pslScore/pslScore.pl
        for the case of DNA vs DNA search only
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

    def get_matching_part(q_size, q_start, q_end, t_size, t_start, t_end, tolerance_prop=0.05):
        """
        Determine if a contig matches entirely the query, or if it is a partial hit. In cases of 
        partial hits determine if it is a split hit (proximal, middle, or distal) or if the hit
        is partial and subsumed within a larger stretch of sequence unrelated to the query
        """
        if q_end - q_start >= q_size * (1 - tolerance_prop):
            return "full"
        elif t_end - t_start >= t_size * (1 - (tolerance_prop * 2)):
            return "middle"
        elif q_start <= q_size * tolerance_prop and t_end >= t_size * (1 - tolerance_prop):
            return "proximal"
        elif q_end >= q_size * (1 - tolerance_prop) and t_start <= t_size * tolerance_prop:
            return "distal"
        else: # hit is only partial and surrounded by a large proportion of unknown sequence
            return "wedged"

    def extract_psl_sequence(fasta_dict, seq_name, strand, start_pos, end_pos, match_part, 
                             up_down_stream):
        """
        Extract sequence from a 'fasta_dict' object using BLAT's PSL coordinate style
        """
        if up_down_stream > 0:
            seq_len = len(fasta_dict[seq_name]["sequence"])
            if match_part == "full" or match_part == "proximal":
                start_pos = max((start_pos - up_down_stream), 0)
            if match_part == "full" or match_part == "distal":
                end_pos = min((end_pos + up_down_stream), seq_len)
        if strand == "+":
            return fasta_dict[seq_name]["sequence"][start_pos:end_pos]
        elif strand == "-":
            return reverse_complement(fasta_dict[seq_name]["sequence"])[start_pos:end_pos]

    def extract_and_stitch_edges(assembly_paths, target_dict, up_down_stream=1000):
        """
        Returns a list of stitched sequences with metadata, sorted by relevance: 'full' hits first
        sorted by 'match_score', followed by assembled hits sorted by 'match_score', and finally
        partial unassembled hits sorted by 'match_score'.
        Assembly: For overlaps follow the coordinate system of the query and remove the overlap from 
        the partial hit with the lower 'match_id', for non-overlapped partial hits, concatenate hits
        intercalating them with as many Ns as indicated by gap in query
        """
        full_hits = []
        sewn_hits = []
        part_hits = []
        asm_num = 1
        for path in assembly_paths:
            identities = [path[0]["match_id"]]
            asm_hit = {"hit": path[0]["hit"], 
                       "matches": path[0]["matches"],
                       "mismatches": path[0]["mismatches"],
                       "q_name": path[0]["q_name"], 
                       "q_size": path[0]["q_size"], 
                       "q_coords": str(path[0]["q_start"]) + "-" + str(path[0]["q_end"]), 
                       "q_cov": path[0]["q_cov"], 
                       "t_name": path[0]["t_name"], 
                       "t_coords": str(path[0]["t_start"]) + "-" + str(path[0]["t_end"]), 
                       "t_strand": path[0]["t_strand"],
                       "match_id": path[0]["match_id"], 
                       "match_score": path[0]["match_score"], 
                       "match_part": path[0]["match_part"], 
                       "seq_match": extract_psl_sequence(target_dict, 
                                                   path[0]["t_name"], 
                                                   path[0]["t_strand"], 
                                                   path[0]["t_start"], 
                                                   path[0]["t_end"],
                                                   path[0]["match_part"],
                                                   0),
                       "seq_full": extract_psl_sequence(target_dict, 
                                                   path[0]["t_name"], 
                                                   path[0]["t_strand"], 
                                                   path[0]["t_start"], 
                                                   path[0]["t_end"],
                                                   path[0]["match_part"],
                                                   up_down_stream)}
            if len(path) > 1:
                for h in range(len(path) -1):
                    asm_hit["matches"] += path[h+1]["matches"]
                    asm_hit["mismatches"] += path[h+1]["mismatches"]
                    asm_hit["q_coords"] += ","+str(path[h+1]["q_start"])+"-"+str(path[h+1]["q_end"])
                    asm_hit["q_cov"] = min((asm_hit["q_cov"] + path[h+1]["q_cov"]), 100.0)
                    asm_hit["t_name"] += ","+path[h+1]["t_name"]
                    asm_hit["t_coords"] += ","+str(path[h+1]["t_start"])+"-"+str(path[h+1]["t_end"])
                    asm_hit["t_strand"] += ","+path[h+1]["t_strand"]
                    asm_hit["match_part"] += ","+path[h+1]["match_part"]
                    identities.append(path[h+1]["match_id"])
                    next_seq_match = extract_psl_sequence(target_dict, 
                                                          path[h+1]["t_name"], 
                                                          path[h+1]["t_strand"], 
                                                          path[h+1]["t_start"], 
                                                          path[h+1]["t_end"],
                                                          path[h+1]["match_part"],
                                                          0)
                    next_seq_full = extract_psl_sequence(target_dict, 
                                                          path[h+1]["t_name"], 
                                                          path[h+1]["t_strand"], 
                                                          path[h+1]["t_start"], 
                                                          path[h+1]["t_end"],
                                                          path[h+1]["match_part"],
                                                          up_down_stream)
                    overlap = path[h]["q_end"] - path[h+1]["q_start"]
                    if overlap < 0:
                        asm_hit["seq_match"] += ("N" * abs(overlap)) + next_seq_match
                        asm_hit["seq_full"] += ("N" * abs(overlap)) + next_seq_full
                    else:
                        if path[h]["match_id"] >= path[h+1]["match_id"]:
                            asm_hit["seq_match"] += next_seq_match[overlap:]
                            asm_hit["seq_full"] += next_seq_full[overlap:]
                        else:
                            asm_hit["seq_match"] = asm_hit["seq_match"][:-overlap] + next_seq_match
                            asm_hit["seq_full"] = asm_hit["seq_full"][:-overlap] + next_seq_full
                asm_hit["hit"] = asm_hit["q_name"] + "|" + "asm_"+str(asm_num).zfill(2)
                asm_hit["match_id"] = statistics.mean(identities)
                asm_hit["match_score"] = ((asm_hit["matches"] + asm_hit["mismatches"]) 
                                          * ((asm_hit["matches"] - asm_hit["mismatches"]) 
                                             / asm_hit["q_size"]))
                asm_num += 1
            if asm_hit["match_part"] == "full":
                full_hits.append(dict(asm_hit))
            elif "," in asm_hit["match_part"]:
                sewn_hits.append(dict(asm_hit))
            else:
                part_hits.append(dict(asm_hit))
        if full_hits:
            full_hits = sorted(full_hits, key=lambda i: i["match_score"], reverse=True)
        if sewn_hits:
            sewn_hits = sorted(sewn_hits, key=lambda i: i["match_score"], reverse=True)
        if part_hits:
            part_hits = sorted(part_hits, key=lambda i: i["match_score"], reverse=True)
        assembly = full_hits + sewn_hits + part_hits
        for a in range(len(assembly)):
            assembly[a]["hit"] = a
        return assembly

    def greedy_assembly_partial_hits(hits_list, tolerance_prop=0.05, max_overlap_bp=50):
        """
        Start with the most 'proximal' hit and look for 'middle' or 'distal' hits that are at most
        overlapped by min('tolerance_prop'*'q_size', max_overlap_bp) and with 'match_id' >= proximal 
        hit's 'match_id'*'tolerance_prop'. Contruct a list of paths to proceed with the function
        'extract_and_stitch_edges' to finally return a list of ranked assembled sequences to 
        repopulate the dictionary of 'nc_hits'
        """
        full_hits, part_hits = [], []
        compatible_pairs = []
        all_paths, last_path = [], []
        contigs_used, valid_paths = [], []
        hit_paths = []
        for hit in hits_list:
            if hit["match_part"] == "full":
                full_hits.append(hit)
            else:
                part_hits.append(hit)
        part_hits = sorted(part_hits, key=lambda i: i["q_start"])
        for h1 in range(len(part_hits)):
            for h2 in range(h1+1, len(part_hits)):
                tolerance_bp = min(int(tolerance_prop * part_hits[h1]["q_size"]), max_overlap_bp)
                if (part_hits[h1]["q_end"] - tolerance_bp <= part_hits[h2]["q_start"] and 
                    max(part_hits[h1]["match_id"], part_hits[h2]["match_id"]) * tolerance_prop <= 
                    min(part_hits[h1]["match_id"], part_hits[h2]["match_id"])):
                    compatible_pairs.append((part_hits[h1]["hit"], part_hits[h2]["hit"]))
        if compatible_pairs:
            for p1 in range(len(compatible_pairs)):
                if last_path == []:
                    last_path = [compatible_pairs[p1][0], compatible_pairs[p1][1]]
                for p2 in range(p1 + 1, len(compatible_pairs)):
                    if last_path[-1] == compatible_pairs[p2][0]:
                        last_path.append(compatible_pairs[p2][1])
                all_paths.append(last_path)
                last_path = []
        if all_paths:
            contigs_used, valid_paths = list(all_paths[0]), list([all_paths[0]])
            for path in all_paths:
                if set(path) - set(contigs_used) == set(path):
                    valid_paths.append(path)
                    contigs_used += path
        if valid_paths:
            for path in valid_paths:
                path_out = []
                for step in path:
                    for hit in part_hits:
                        if step == hit["hit"]:
                            path_out.append(hit)
                hit_paths.append(path_out)
        assembly_paths = ([[full_hits[f]] for f in range(len(full_hits))] 
                          + hit_paths
                          + [[part_hits[h]] for h in range(len(part_hits)) 
                             if part_hits[h]["hit"] not in contigs_used])
        assembly = extract_and_stitch_edges(assembly_paths, target_dict)
        return assembly

    nc_hits = {}
    with open(psl_path) as psl_in:
        for line in psl_in:
            cols = line.split()
            matches, mismatches, rep_matches = int(cols[0]), int(cols[1]), int(cols[2])
            q_inserts, t_strand = int(cols[4]), cols[8]
            q_name, q_size, q_start, q_end = cols[9], int(cols[10]), int(cols[11]), int(cols[12])
            t_name, t_size, t_start, t_end = cols[13], int(cols[14]), int(cols[15]), int(cols[16])
            q_cov = ((q_end - q_start) / q_size) * 100.0
            match_score = ((matches + rep_matches + mismatches) 
                           * (((matches + rep_matches) - mismatches) / q_size))
            match_id = calc_psl_identity(q_end, q_start, t_end, t_start, q_inserts, matches, 
                                         rep_matches, mismatches)
            match_part = get_matching_part(q_size, q_start, q_end, t_size, t_start, t_end)
            hit = {"hit": t_name+"|"+t_strand+"|"+str(q_start)+"-"+str(q_end)+"|"+str(match_id),
                   "matches": matches + rep_matches,
                   "mismatches": mismatches,
                   "q_name": q_name,
                   "q_size": q_size,
                   "q_start": q_start,
                   "q_end": q_end,
                   "q_cov": q_cov,
                   "t_name": t_name,
                   "t_size": t_size,
                   "t_start": t_start,
                   "t_end": t_end,
                   "t_strand": t_strand,
                   "match_id": match_id,
                   "match_score": match_score,
                   "match_part": match_part}
            if q_cov >= min_coverage and match_part != "wedged":
                if q_name not in nc_hits:
                    nc_hits[q_name] = [dict(hit)]
                else:
                    nc_hits[q_name].append(dict(hit))
    for nc_ref in nc_hits:
        nc_hits[nc_ref] = greedy_assembly_partial_hits(nc_hits[nc_ref])

    # Find if reference have been formatted like Angiosperms353.FAA in order to accomodate more than 
    # a single reference of the same type, we can also use 'set_a.REFERENCE_CLUSTER_SEPARATOR'
    refs_have_separators = True
    for nc_ref in nc_hits:
        if "-" not in nc_ref:
            refs_have_separators = False
            break

    # If multiple references of the same type exist in the reference, then choos the one with the
    # highest 'match_score'
    if nc_hits:
        best_nc_hits = {}
        for nc_ref in nc_hits:
            if refs_have_separators:
                nc_ref_cluster = nc_ref.split("-")[-1]
            else:
                nc_ref_cluster = nc_ref
            if nc_ref_cluster not in best_nc_hits:
                best_nc_hits[nc_ref_cluster] = nc_hits[nc_ref]
            else:
                if (nc_hits[nc_ref][0]["match_score"] 
                    > best_nc_hits[nc_ref_cluster][0]["match_score"]):
                    best_nc_hits[nc_ref_cluster] = nc_hits[nc_ref]
        print(best_nc_hits)
        return best_nc_hits
    else:
        return False



psl_path = "/Users/emortiz/Desktop/captus_tests/cucumis/blat_id75ms15.psl"
target_dict = fasta_to_dict("/Users/emortiz/Desktop/captus_tests/cucumis/cu009.fasta")
query_dict = fasta_to_dict("/Users/emortiz/Desktop/captus_tests/cucumis/non_coding_test.fasta")
min_coverage = 10
best_hits_only = False

blat_dna_psl_to_dict(psl_path, target_dict, query_dict, min_coverage, best_hits_only)
