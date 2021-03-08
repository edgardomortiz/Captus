
def scipio_yaml_to_dict(yaml_path, min_identity, min_coverage, marker_type):
    """
    Takes a Scipio's YAML output file 'yaml_path' and returns a dictionary of all the hits per
    protein and their curated sequence
    The function creates a new score based on Scipio's own score and the length of the hit
    """

    # Initialize variables before starting the loop along the YAML file
    raw = {}  # raw models dictionary, need to be filtered by 'min_coverage' and 'min_identity' after
    protein = ""
    hit = 0
    match_type = ""
    overlapped = False
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
        "region": "",       # included here just to match BLAT's data for non-coding extraction
        "gapped": False,    # hit contains gaps with respect to reference
        "seq_flanked": "",  # concatenation of contigs, overlaps merged or 50 'n's in between
        "seq_gene": "",     # gene sequence excluding upstream and downstream regions
        "seq_nt": "",       # concatenation of CDS in nucleotide (includes STOP codon)
        "seq_aa": "",       # concatenation of CDS in aminoacid (STOP codon excluded)
    }

