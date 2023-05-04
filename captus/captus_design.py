#!/usr/bin/env python3
"""
Copyright 2020-2023 Edgardo M. Ortiz (e.ortiz.v@gmail.com)
https://github.com/edgardomortiz/Captus

This is the control program for the bait design pipeline of Captus.

Multi-level argparse modified from:
https://chase-seibert.github.io/blog/2014/03/21/python-multilevel-argparse.html

This file is part of Captus. Captus is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Captus is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Captus. If
not, see <http://www.gnu.org/licenses/>.
"""


import argparse
import sys

from . import settings as settings
from .cluster import cluster
from .misc import MyHelpFormatter, bold, red
from .version import __version__


class CaptusDesign(object):

    ################################################################################################
    ########################################################################## CAPTUS-DESIGN SECTION
    def __init__(self):

        description = bold(f'Captus {__version__}:'
                           " Assembly of Phylogenomic Datasets from High-Throughput Sequencing data")
        parser=argparse.ArgumentParser(
            usage="captus_assembly command [options]",
            description=description,
            formatter_class=MyHelpFormatter,
            epilog="For help on a particular command: captus_assembly command -h",
            add_help=False
        )

        required_group = parser.add_argument_group("Captus-design commands")
        required_group.add_argument(
            "command",
            help="B|Program commands (in typical order of execution)\n"
                 "cluster = Find clusters in CDS or transcripts across several samples. Genomes"
                 " in FASTA + annotation track in GFF, or simply transcripts or CDS in FASTA are"
                 " allowed\n"
                 "select = Filter the clusters according to several selection parameters\n"
                 "bait = Create baits based on the selected clusters. It can process the clusters"
                 " selected in the previous step or any set of input alignments"
        )

        help_group = parser.add_argument_group("Help")
        help_group.add_argument(
            "-h", "--help",
            action="help",
            help="Show this help message and exit"
        )
        help_group.add_argument(
            "--version",
            action="version",
            version=f"Captus v{__version__}",
            help="Show Captus' version number"
        )

        if len(sys.argv) == 1:
            parser.print_help()
            exit(red("\nERROR: Missing command\n"))

        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            parser.print_help()
            exit(red(f"\nERROR: Unrecognized command: {sys.argv[1]}\n"))
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()


    ################################################################################################
    ################################################################################ CLUSTER SECTION
    def cluster(self):
        description = bold(
            "Captus-design: Cluster; find clusters of markers across genomes or transcripts\n"
        )
        parser=argparse.ArgumentParser(
            usage="captus_design cluster -m MARKERS_TO_CLUSTER [MARKERS_TO_CLUSTER ...] [options]",
            description=description,
            formatter_class=MyHelpFormatter,
            epilog="For more information, please see https://github.com/edgardomortiz/Captus",
            add_help=False
        )

        input_group = parser.add_argument_group("Input")
        input_group.add_argument(
            "-m", "--markers_to_cluster",
            action="store",
            nargs="+",
            type=str,
            required=True,
            dest="markers_to_cluster",
            help="B|Input directory containing data to cluster. Valid extensions are: .fasta,"
                 " .fasta.gz, .fna, .fna.gz, .gff, .gff.gz, .gff3, .gff3.gz, .gtf, .gtf.gz."
                 " Everything before he extension will be used as sample name. If a FASTA file is"
                 " found with its corresponding annotation track in GFF, then the CDS features will"
                 " be used for clustering, otherwise, every sequence in the FASTA file will be used"
                 " for clustering. There are a few ways to provide the input filenames:\n"
                 "A directory = path to directory containing the files (e.g.: -m ./genomic_data)\n"
                 "A list = filenames separated by space (e.g.: -m A.fasta A.gff B.fna.gz C.fna)\n"
                 "A pattern = UNIX matching expression (e.g.: -m ./genomic_data/A.*)"
        )
        input_group.add_argument(
            "-b", "--bait_length",
            action="store",
            default=120,
            type=int,
            dest="bait_length",
            help="Projected length for the baits (a.k.a. probes), only exons of at least this length"
                 " are used for the creation of baits"
        )

        output_group = parser.add_argument_group("Output")
        output_group.add_argument(
            "-o", "--out",
            action="store",
            default="./01_clustered_markers",
            type=str,
            dest="out",
            help="Output directory name"
        )
        output_group.add_argument(
            "--keep_all",
            action="store_true",
            dest="keep_all",
            help="Do not delete any intermediate files"
        )
        output_group.add_argument(
            "--overwrite",
            action="store_true",
            dest="overwrite",
            help="Overwrite previous results"
        )

        clustering_group = parser.add_argument_group("Clustering of markers")
        clustering_group.add_argument(
            "--clust_program",
            action="store",
            default="mmseqs",
            type=str,
            dest="clust_program",
            choices=["mmseqs", "vsearch"],
            help="Clustering software to use for deduplication and clustering"
        )
        clustering_group.add_argument(
            "--max_seq_len",
            action="store",
            default="auto",
            type=str,
            dest="max_seq_len",
            help="When left to 'auto' Captus will filter sequences longer than 65535 bp if MMseqs2"
                 " is used or longer than 5000 bp if VSEARCH is used. Otherwise provide a length in"
                 " bp to remove longer sequences prior to deduplication and clustering"
        )
        clustering_group.add_argument(
            "-s", "--strand",
            action="store",
            default="both",
            type=str,
            dest="strand",
            choices=["both", "plus"],
            help="Strand to compare for deduplication or clustering when using VSEARCH (MMSeqs2"
                 " processes both strands by default). Using both strands increases"
                 " processing time in VSEARCH"
        )
        clustering_group.add_argument(
            "-d", "--dedup_threshold",
            action="store",
            default=99.5,
            type=float,
            dest="dedup_threshold",
            help="Percent identity threshold for within-sample marker deduplication"
        )
        clustering_group.add_argument(
            "-c", "--clust_threshold",
            action="store",
            default=75,
            type=float,
            dest="clust_threshold",
            help="Percent identity threshold for across-samples marker clustering"
        )

        align_group = parser.add_argument_group("Alignment of clusters")
        align_group.add_argument(
            "--align_singletons",
            action="store_true",
            dest="align_singletons",
            help="Align clusters containing sequences from a single species"
        )
        align_group.add_argument(
            "--timeout",
            action="store",
            default=21600,
            type=int,
            dest="timeout",
            help="Maximum allowed time in seconds for a single alignment"
        )

        curation_group = parser.add_argument_group("Curation of aligned clusters")
        curation_group.add_argument(
            "--gaps",
            action="store",
            default=0.2,
            type=float,
            dest="gaps",
            help="Maximum proportion of missing data allowed from leading and trailing columns on"
                 " aligned clusters"
        )
        curation_group.add_argument(
            "--focal_species",
            action="store",
            type=str,
            dest="focal_species",
            help="Comma-separated list (no spaces) of species whose presence must be tracked in the"
                 " alignments"
        )
        curation_group.add_argument(
            "--outgroup_species",
            action="store",
            type=str,
            dest="outgroup_species",
            help="Comma-separated list (no spaces) of outgroup species names"
        )
        curation_group.add_argument(
            "--addon_samples",
            action="store",
            type=str,
            dest="addon_samples",
            help="Comma-separated list (no spaces) of samples to exclude from the maximum proportion"
                 " of missing data, these can be poor quality assemblies for example"
        )

        other_group = parser.add_argument_group("Other")
        other_group.add_argument(
            "--redo_from",
            action="store",
            choices=[
                "importation",
                "dereplication",
                "clustering",
                "alignment",
                "curation",
            ],
            type=str,
            dest="redo_from",
            help="B|Repeat analysis from a particular stage:\n"
                 "import = Delete all subdirectories and import input data again\n"
                 "deduplication = Delete all deduplicated FASTA files and restart\n"
                 "clustering = Delete clustering subdirectory and restart\n"
                 "alignment = Delete alignments subdirectory and restart\n"
                 "curation = Delete curated alignments subdirectory and restart"
        )
        other_group.add_argument(
            "--mmseqs_path",
            action="store",
            default="mmseqs",
            type=str,
            dest="mmseqs_path",
            help="Path to MMseqs2"
        )
        other_group.add_argument(
            "--vsearch_path",
            action="store",
            default="vsearch",
            type=str,
            dest="vsearch_path",
            help="Path to VSEARCH"
        )
        other_group.add_argument(
            "--mafft_path",
            action="store",
            default="mafft",
            type=str,
            dest="mafft_path",
            help="Path to MAFFT (default: mafft/mafft.bat)"
        )
        other_group.add_argument(
            "--ram",
            action="store",
            default="auto",
            type=str,
            dest="ram",
            help="Maximum RAM in GB (e.g.: 4.5) dedicated to Captus, 'auto' uses"
                 f" {settings.RAM_FRACTION:.0%}% of available RAM"
        )
        other_group.add_argument(
            "--threads",
            action="store",
            default="auto",
            type=str,
            dest="threads",
            help="Maximum number of CPUs dedicated to Captus, 'auto' uses all available CPUs"
        )
        other_group.add_argument(
            "--concurrent",
            action="store",
            default='auto',
            type=str,
            dest="concurrent",
            help="Captus will attempt to execute this many alignments concurrently. CPUs will be"
                 " divided by this value for each individual process. If set to 'auto', Captus will"
                 " set as many processes as to at least have 2 threads available for each MAFFT"
                 " process. Clustering is done sample by sample using as many threads as possible."
        )
        other_group.add_argument(
            "--debug",
            action="store_true",
            dest="debug",
            help="Enable debugging mode, parallelization is disabled so errors are logged to screen"
        )
        other_group.add_argument(
            "--show_more",
            action="store_true",
            dest="show_more",
            help="Show individual alignment information during the run. Detailed information is"
                 " written regardless to the log"
        )

        help_group = parser.add_argument_group("Help")
        help_group.add_argument(
            "-h", "--help",
            action="help",
            help="Show this help message and exit"
        )
        help_group.add_argument(
            "--version",
            action="version",
            version=f"Captus v{__version__}",
            help="Show Captus' version number"
        )

        if len(sys.argv) == 2:
            parser.print_help()
            exit(red("\nERROR: Missing required argument -m/--markers_to_cluster\n"))

        full_command = " ".join(sys.argv)
        args = parser.parse_args(sys.argv[2:])
        cluster(full_command, args)
        exit(1)


    ################################################################################################
    ################################################################################# SELECT SECTION
    def select(self):
        description = bold("Captus-design: Select; perform de novo assembly using MEGAHIT\n")
        parser=argparse.ArgumentParser(
            usage="captus_design select -c CAPTUS_CLUSTERS_DIR [options]",
            description=description,
            formatter_class=MyHelpFormatter,
            epilog="For more information, please see https://github.com/edgardomortiz/Captus",
            add_help=False
        )

        input_group = parser.add_argument_group("Input")
        input_group.add_argument(
            "-c", "--captus_clusters_dir",
            action="store",
            default="./01_clustered_markers",
            type=str,
            required=True,
            dest="captus_clusters_dir",
            help="Path to the output directory that contains the curated aligned clusters found"
                 " during the previous step of Captus-design. This directory is called"
                 " '01_clustered_markers' if you did not specify a different name during the"
                 " 'cluster' step"
        )
        input_group.add_argument(
            "-d", "--dry_run",
            action="store_true",
            dest="dry_run",
            help="Just preview the number of markers resulting from the selecting parameters as well"
                 " as the total projected capture footprint"
        )

        output_group = parser.add_argument_group("Output")
        output_group.add_argument(
            "-o", "--out",
            action="store",
            default="./02_selected_markers",
            type=str,
            dest="out",
            help="Output directory name. This directory will contain a copy of only the selected"
                 " markers"
        )
        output_group.add_argument(
            "--keep_all",
            action="store_true",
            dest="keep_all",
            help="Do not delete any intermediate files"
        )
        output_group.add_argument(
            "--overwrite",
            action="store_true",
            dest="overwrite",
            help="Overwrite previous results"
        )

        filter_group = parser.add_argument_group("Filtering of markers")
        filter_group.add_argument(
            "--par",
            action="store_true",
            dest="include_paralogs",
            help="Include paralog clusters in the selection"
        )
        filter_group.add_argument(
            "--len",
            action="store",
            type=str,
            dest="length",
            default="720,20000",
            help="Alignment range length in bp (min,max)"
        )
        filter_group.add_argument(
            "--pid",
            action="store",
            type=str,
            dest="pairwise_identity",
            default="75.0,99.99",
            help="Average pairwise percentage identity range (min,max)"
        )
        filter_group.add_argument(
            "--gaps",
            action="store",
            type=str,
            dest="gappiness",
            default="0.0,5.0",
            help="Range of percentage of internal gaps allowed (min,max)"
        )
        filter_group.add_argument(
            "--sequences",
            action="store",
            type=str,
            dest="num_sequences",
            default="1,1000",
            help="Range of number of sequences per alignment (min,max)"
        )
        filter_group.add_argument(
            "--samples",
            action="store",
            type=str,
            dest="num_samples",
            default="1,1000",
            help="Range of number of samples per alignment (min,max)"
        )
        filter_group.add_argument(
            "--focal",
            action="store",
            type=str,
            dest="num_focal_species",
            default="0,1000",
            help="Range of number of focal species per alignment (min,max)"
        )
        filter_group.add_argument(
            "--outgroup",
            action="store",
            type=str,
            dest="num_outgroup_species",
            default="0,1000",
            help="Range of number of outgroup species per alignment (min,max)"
        )
        filter_group.add_argument(
            "--addon",
            action="store",
            type=str,
            dest="num_add_on_samples",
            default="0,1000",
            help="Range of number of add-on samples per alignment (min,max)"
        )
        filter_group.add_argument(
            "--species",
            action="store",
            type=str,
            dest="num_species",
            default="0,1000",
            help="Range of number of species per alignment (min,max)"
        )
        filter_group.add_argument(
            "--genera",
            action="store",
            type=str,
            dest="num_genera",
            default="0,1000",
            help="Range of number of genera per alignment (min,max)"
        )
        filter_group.add_argument(
            "--cdl",
            action="store",
            type=str,
            dest="cds_length",
            default="0,20000",
            help="Original length range of entire CDS in bp before trimming (min,max)"
        )
        filter_group.add_argument(
            "--llr",
            action="store",
            type=str,
            dest="length_long_exon_retained",
            default="0,20000",
            help="Length range of sequence retained (in bp) corresponding to long exons (min,max)"
        )
        filter_group.add_argument(
            "--lsr",
            action="store",
            type=str,
            dest="length_short_exon_retained",
            default="0,20000",
            help="Length range of sequence retained (in bp) corresponding to short exons (min,max)"
        )
        filter_group.add_argument(
            "--ptr",
            action="store",
            type=str,
            dest="percentage_total_cds_retained",
            default="0.0,100.0",
            help="Percentage range of the CDS original length retained (min,max)"
        )
        filter_group.add_argument(
            "--plr",
            action="store",
            type=str,
            dest="percentage_long_exon_retained",
            default="0.0,100.0",
            help="Percentage range of sequence retained corresponding to long exons (min,max)"
        )
        filter_group.add_argument(
            "--psr",
            action="store",
            type=str,
            dest="percentage_short_exon_retained",
            default="0.0,15.0",
            help="Percentage range of sequence retained corresponding to short exons (min,max)"
        )

        help_group = parser.add_argument_group("Help")
        help_group.add_argument(
            "-h", "--help",
            action="help",
            help="Show this help message and exit"
        )
        help_group.add_argument(
            "--version",
            action="version",
            version=f"Captus v{__version__}",
            help="Show Captus' version number"
        )

        if len(sys.argv) == 2:
            parser.print_help()
            exit(red("\nERROR: Missing required argument -c/--captus_clusters_dir\n"))

        full_command = " ".join(sys.argv)
        args = parser.parse_args(sys.argv[2:])
        assemble(full_command, args)
        exit(1)


    ################################################################################################
    ################################################################################### BAIT SECTION
    def bait(self):
        description = bold("Captus-design: Bait; create baits from selected alignments\n")
        parser=argparse.ArgumentParser(
            usage="captus_design bait -s CAPTUS_SELECTED_DIR [options]",
            description=description,
            formatter_class=MyHelpFormatter,
            epilog="For more information, please see https://github.com/edgardomortiz/Captus",
            add_help=False
        )

        input_group = parser.add_argument_group("Input")
        input_group.add_argument(
            "-s", "--captus_selected_dir",
            action="store",
            default="./02_selected_markers",
            type=str,
            required=True,
            dest="captus_selected_dir",
            help="Path to an output directory from the 'select' step of Captus-design which is"
                 " tipically called '02_selected_markers'. Captus will search for FASTA files"
                 " inside this directory and every subdirectory"
        )
        input_group.add_argument(
            "-c", "--captus_clusters_dir",
            action="store",
            default="./01_clustered_markers",
            type=str,
            dest="captus_clusters_dir",
            help="Path to the output directory that contains the curated aligned clusters found"
                 " during the previous 'cluster' step of Captus-design. This directory is called"
                 " '01_clustered_markers' if you did not specify a different name during the"
                 " 'cluster' step. Captus uses the exon data table to create baits that do not span"
                 " exonic boundaries"
        )

        output_group = parser.add_argument_group("Output")
        output_group.add_argument(
            "-o", "--out",
            action="store",
            default="./03_baits",
            type=str,
            dest="out",
            help="Output directory name"
        )
        output_group.add_argument(
            "--keep_all",
            action="store_true",
            dest="keep_all",
            help="Do not delete any intermediate files"
        )
        output_group.add_argument(
            "--overwrite",
            action="store_true",
            dest="overwrite",
            help="Overwrite previous results"
        )

        bait_group = parser.add_argument_group("Bait creation and filtering")
        bait_group.add_argument(
            "-b", "--bait-length",
            action="store",
            type=int,
            dest="bait_length",
            default=120,
            help="Length in bp for bait design"
        )
        bait_group.add_argument(
            "--iN", "--include-Ns",
            action="store_true",
            dest="include_Ns",
            help="If enabled, baits with Ns will also be considered for design. When omitted, baits"
                 " with Ns are discarded"
        )
        bait_group.add_argument(
            "--ki", "--keep-iupac",
            action="store_true",
            dest="keep_iupac",
            help="If enabled, IUPAC ambiguity codes will be kept in the final baits. When omitted,"
                 " the IUPAC codes are disambiguated choosing one nucleotide in the ambiguitiy at"
                 " random"
        )
        bait_group.add_argument(
            "--er", "--exclude-reference",
            type=str,
            dest="exclude_reference",
            help="If you want to exclude baits that map to a reference sequence, you can specify the"
                 " path to a FASTA file"
        )
        bait_group.add_argument(
            "--gc", "--gc-content",
            action="store",
            type=str,
            dest="gc_content",
            default="30.0,50.0",
            help="GC content allowed in bait"
        )
        bait_group.add_argument(
            "--mt", "--melting-temperature",
            action="store",
            type=str,
            dest="melting_temperature",
            default="0.0,120.0",
            help="Melting temperature allowed for bait"
        )
        bait_group.add_argument(
            "--hc", "--hybridization-chemistry",
            action="store",
            type=str,
            dest="hybridization_chemistry",
            choices={"RNA-DNA", "RNA-RNA", "DNA-DNA"},
            default="RNA-DNA",
            help="Probe-Target chemistry, used in melting temperature calculation"
        )
        bait_group.add_argument(
            "-s", "--sodium",
            action="store",
            type=float,
            dest="sodium",
            default=0.9,
            help="Na concentration for melting temperature calculation"
        )
        bait_group.add_argument(
            "-f", "--formamide",
            action="store",
            type=float,
            dest="formamide",
            default=0.0,
            help="Formamide concentration for melting temperature calculation"
        )
        bait_group.add_argument(
            "--mmp", "--max-masked-percentage",
            action="store",
            type=float,
            dest="max_masked_percentage",
            default=25.0,
            help="Low complexity regions of baits are masked using VSEARCH, define maximum"
                 " percentage of sequence allowed to be masked"
        )
        bait_group.add_argument(
            "--mhl", "--max-homopolymer-length",
            type=int,
            dest="max_homopolymer_length",
            default=6,
            help="Maximum length in bp for homopolymers within a bait"
        )

        clustering_group = parser.add_argument_group("Bait clustering and tiling")
        clustering_group.add_argument(
            "--tpo", "--tiling-percentage-overlap",
            type=float,
            dest="tiling_percentage_overlap",
            default=50.0,
            help="Maximum percentage of the length of the bait allowed to overlap with the adjacent"
                 " bait in the tiling design, e.g. if you use 25%% of overlap for 120 bp baits, they"
                 " will overlap by at most 30 bp"
        )
        clustering_group.add_argument(
            "--ct", "--clustering-threshold",
            action="store",
            type=int,
            dest="clustering_threshold",
            default=86,
            help="The entire set of potential baits are finally clustered at this identity, keeping"
                 " a single centroid bait sequence per cluster. It has been shown than a bait can"
                 " hybridize with a target sequence that is at least 75%% identical, therefore do"
                 " not use clustering thresholds lower than 0.75"
        )

        other_group = parser.add_argument_group("Other")
        other_group.add_argument(
            "--vsearch_path",
            action="store",
            default="vsearch",
            type=str,
            dest="vsearch_path",
            help="Path to VSEARCH"
        )
        other_group.add_argument(
            "--threads",
            action="store",
            default="auto",
            type=str,
            dest="threads",
            help="Maximum number of CPUs dedicated to Captus, 'auto' uses all available CPUs"
        )
        other_group.add_argument(
            "--concurrent",
            action="store",
            default="auto",
            type=str,
            dest="concurrent",
            help="Captus will attempt to execute this many extractions concurrently. RAM and CPUs"
                 " will be divided by this value for each individual process. If set to 'auto',"
                 " Captus will set as many processes as to at least have 2GB of RAM available for"
                 " each process due to the RAM requirements of BLAT"
        )
        other_group.add_argument(
            "--debug",
            action="store_true",
            dest="debug",
            help="Enable debugging mode, parallelization is disabled so errors are logged to screen"
        )
        other_group.add_argument(
            "--show_less",
            action="store_true",
            dest="show_less",
            help="Do not show individual sample information during the run, the information is still"
                 " written to the log"
        )

        help_group = parser.add_argument_group("Help")
        help_group.add_argument(
            "-h", "--help",
            action="help",
            help="Show this help message and exit"
        )
        help_group.add_argument(
            "--version",
            action="version",
            version=f"Captus v{__version__}",
            help="Show Captus' version number"
        )

        if len(sys.argv) == 2:
            parser.print_help()
            exit(red("\nERROR: Missing required argument -a/--captus_assemblies_dir\n"))

        full_command = " ".join(sys.argv)
        args = parser.parse_args(sys.argv[2:])
        extract(full_command, args)
        exit(1)


def main():
    CaptusDesign()


if __name__ == "__main__":
    main()

