#!/usr/bin/env python3
"""
Copyright 2020-2025 Edgardo M. Ortiz (e.ortiz.v@gmail.com)
https://github.com/edgardomortiz/Captus

This is the control program for the assembly pipeline of Captus.

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

from . import settings
from .align import align
from .assemble import assemble
from .clean import clean
from .extract import extract
from .misc import MyHelpFormatter, bold, red
from .version import __version__


class CaptusAssembly(object):
    ################################################################################################
    ######################################################################## CAPTUS-ASSEMBLY SECTION
    def __init__(self):
        description = bold(
            f"Captus {__version__}:"
            " Assembly of Phylogenomic Datasets from High-Throughput Sequencing data"
        )
        parser = argparse.ArgumentParser(
            usage="captus command [options]",
            description=description,
            formatter_class=MyHelpFormatter,
            epilog="For help on a particular command: captus command -h",
            add_help=False,
        )

        required_group = parser.add_argument_group("Captus-assembly commands")
        required_group.add_argument(
            "command",
            help="B|Program commands (in typical order of execution)\n"
            "clean = Trim adaptors and quality filter reads with BBTools, run FastQC on the"
            " raw and cleaned reads\n"
            "assemble = Perform de novo assembly with MEGAHIT and estimate contig depth of"
            " coverage with Salmon. Assembling reads that were cleaned with the 'clean' command"
            " is recommended, but reads cleaned elsewhere are also allowed\n"
            "extract = Recover targeted markers with BLAT and Scipio: Extracting markers"
            " from the assembly obtained with the 'assemble' command is recommended, but any"
            " other assemblies in FASTA format are also allowed\n"
            "align = Align extracted markers across samples with MAFFT or MUSCLE: Marker"
            " alignment depends on the directory structure created by the 'extract' command."
            " This step also performs paralog filtering and alignment trimming using TAPER and"
            " ClipKIT",
        )

        help_group = parser.add_argument_group("Help")
        help_group.add_argument("-h", "--help", action="help", help="Show this help message and exit")
        help_group.add_argument(
            "--version",
            action="version",
            version=f"Captus v{__version__}",
            help="Show Captus' version number",
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
    ################################################################################## CLEAN SECTION
    def clean(self):
        description = bold(
            "Captus-assembly: Clean; remove adaptors and quality-filter reads with BBTools\n"
        )
        parser = argparse.ArgumentParser(
            usage="captus clean -r READS [READS ...] [options]",
            description=description,
            formatter_class=MyHelpFormatter,
            epilog="For more information, please see https://github.com/edgardomortiz/Captus",
            add_help=False,
        )

        input_group = parser.add_argument_group("Input")
        input_group.add_argument(
            "-r",
            "--reads",
            action="store",
            nargs="+",
            type=str,
            required=True,
            dest="reads",
            help="B|FASTQ files. Valid filename extensions are: .fq, .fastq, .fq.gz, and .fastq.gz."
            " The names must include the string '_R1' (and '_R2' when pairs are provided)."
            " Everything before the string '_R1' will be used as sample name. There are a few"
            " ways to provide the FASTQ files:\n"
            "A directory = path to directory containing FASTQ files (e.g.: -r ./raw_reads)\n"
            "A list = filenames separated by space (e.g.: -r A_R1.fq A_R2.fq B_R1.fq C_R1.fq)\n"
            "A pattern = UNIX matching expression (e.g.: -r ./raw_reads/*.fastq.gz)",
        )

        output_group = parser.add_argument_group("Output")
        output_group.add_argument(
            "-o",
            "--out",
            action="store",
            default="./01_clean_reads",
            type=str,
            dest="out",
            help="Output directory name",
        )
        output_group.add_argument(
            "--keep_all",
            action="store_true",
            dest="keep_all",
            help="Do not delete any intermediate files",
        )
        output_group.add_argument(
            "--overwrite", action="store_true", dest="overwrite", help="Overwrite previous results"
        )

        adaptor_group = parser.add_argument_group("Adaptor trimming")
        adaptor_group.add_argument(
            "--adaptor_set",
            action="store",
            default="ALL",
            type=str,
            dest="adaptor_set",
            help="B|Set of adaptors to remove\n"
            "Illumina = Illumina adaptors included in BBTools\n"
            "BGI = BGISEQ, DNBSEG, or MGISEQ adaptors\n"
            "ALL = Illumina + BGI\n"
            "Alternatively, you can provide a path to a FASTA file containing your own adaptors",
        )
        adaptor_group.add_argument(
            "--rna", action="store_true", dest="rna", help="Trim ploy-A tails from RNA-Seq reads"
        )

        quality_group = parser.add_argument_group("Quality trimming and filtering")
        quality_group.add_argument(
            "--trimq",
            action="store",
            default=13,
            type=int,
            dest="trimq",
            help="Leading and trailing read regions with average PHRED quality score below this"
            " value will be trimmed",
        )
        quality_group.add_argument(
            "--maq",
            action="store",
            default=16,
            type=int,
            dest="maq",
            help="After quality trimming, reads with average PHRED quality score below this value"
            " will be removed",
        )
        quality_group.add_argument(
            "--ftl",
            action="store",
            default=0,
            type=int,
            dest="ftl",
            help="Trim any base to the left of this position. For example, if you want to remove 4"
            " bases from the left of the reads set this number to 5",
        )
        quality_group.add_argument(
            "--ftr",
            action="store",
            default=0,
            type=int,
            dest="ftr",
            help="Trim any base to the right of this position. For example, if you want to truncate"
            " your reads length to 100 bp set this number to 100",
        )

        qcstats_group = parser.add_argument_group("QC Statistics")
        qcstats_group.add_argument(
            "--qc_program",
            action="store",
            choices=[
                "fastqc",
                "falco",
            ],
            default="fastqc",
            type=str,
            dest="qc_program",
            help="Which program to use to obtain the statistics from the raw and cleaned FASTQ files."
            " Falco produces identical results to FastQC while being much faster",
        )
        qcstats_group.add_argument(
            "--skip_qc_stats",
            action="store_true",
            dest="skip_qc_stats",
            help="Skip FastQC/Falco analysis on raw and cleaned reads",
        )

        other_group = parser.add_argument_group("Other")
        other_group.add_argument(
            "--bbduk_path",
            action="store",
            default="bbduk.sh",
            type=str,
            dest="bbduk_path",
            help="Path to bbduk.sh",
        )
        other_group.add_argument(
            "--falco_path",
            action="store",
            default="falco",
            type=str,
            dest="falco_path",
            help="Path to Falco",
        )
        other_group.add_argument(
            "--fastqc_path",
            action="store",
            default="fastqc",
            type=str,
            dest="fastqc_path",
            help="Path to FastQC",
        )
        other_group.add_argument(
            "--ram",
            action="store",
            default="auto",
            type=str,
            dest="ram",
            help="Maximum RAM in GB (e.g.: 4.5) dedicated to Captus, 'auto' uses"
            f" {settings.RAM_FRACTION:.0%}% of available RAM",
        )
        other_group.add_argument(
            "--threads",
            action="store",
            default="auto",
            type=str,
            dest="threads",
            help="Maximum number of CPUs dedicated to Captus, 'auto' uses all available CPUs",
        )
        other_group.add_argument(
            "--concurrent",
            action="store",
            default="auto",
            type=str,
            dest="concurrent",
            help="Captus will attempt to run FastQC concurrently on this many samples. If set to"
            f" 'auto', Captus will run at most {settings.MAX_HDD_READ_INSTANCES} instances of"
            " FastQC or as many CPU cores are available, whatever number is lower",
        )
        other_group.add_argument(
            "--debug",
            action="store_true",
            dest="debug",
            help="Enable debugging mode, parallelization is disabled so errors are logged to screen",
        )
        other_group.add_argument(
            "--show_less",
            action="store_true",
            dest="show_less",
            help="Do not show individual sample information during the run, the information is still"
            " written to the log",
        )

        help_group = parser.add_argument_group("Help")
        help_group.add_argument("-h", "--help", action="help", help="Show this help message and exit")
        help_group.add_argument(
            "--version",
            action="version",
            version=f"Captus v{__version__}",
            help="Show Captus' version number",
        )

        if len(sys.argv) == 2:
            parser.print_help()
            exit(red("\nERROR: Missing required argument -r/--reads\n"))

        full_command = " ".join(sys.argv)
        args = parser.parse_args(sys.argv[2:])
        clean(full_command, args)
        exit(1)

    ################################################################################################
    ############################################################################### ASSEMBLE SECTION
    def assemble(self):
        description = bold("Captus-assembly: Assemble; perform de novo assembly using MEGAHIT\n")
        parser = argparse.ArgumentParser(
            usage="captus assemble -r READS [READS ...] [options]",
            description=description,
            formatter_class=MyHelpFormatter,
            epilog="For more information, please see https://github.com/edgardomortiz/Captus",
            add_help=False,
        )

        input_group = parser.add_argument_group("Input")
        input_group.add_argument(
            "-r",
            "--reads",
            action="store",
            default="./01_clean_reads",
            nargs="+",
            type=str,
            required=True,
            dest="reads",
            help="B|FASTQ files. Valid file name extensions are: .fq, .fastq, .fq.gz, and .fastq.gz."
            " The names must include the string '_R1' (and '_R2' when pairs are provided). "
            " Everything before the string '_R1' will be used as sample name. There are a few"
            " ways to provide the FASTQ files:\n"
            "A directory = path to directory containing FASTQ files (e.g.: -r ./raw_reads)\n"
            "A list = file names separated by space (e.g.: -r A_R1.fq A_R2.fq B_R1.fq C_R1.fq)\n"
            "A pattern = UNIX matching expression (e.g.: -r ./raw_reads/*.fastq.gz)",
        )
        input_group.add_argument(
            "--sample_reads_target",
            action="store",
            default=0,
            type=int,
            dest="sample_reads_target",
            help="Use this number of read pairs (or reads if single-end) for assembly. Reads are"
            " randomly subsampled with 'reformat.sh' from BBTools (option:"
            " srt/samplereadstarget). Useful for limiting the amount of data of samples with"
            " very high sequencing depth. To use all the reads, set this value to 0",
        )

        output_group = parser.add_argument_group("Output")
        output_group.add_argument(
            "-o",
            "--out",
            action="store",
            default="./02_assemblies",
            type=str,
            dest="out",
            help="Output directory name. Inside this directory the output for each sample will be"
            " stored in a subdirectory named as 'Sample_name__captus-asm'",
        )
        output_group.add_argument(
            "--keep_all",
            action="store_true",
            dest="keep_all",
            help="Do not delete any intermediate files",
        )
        output_group.add_argument(
            "--overwrite", action="store_true", dest="overwrite", help="Overwrite previous results"
        )

        megahit_group = parser.add_argument_group("MEGAHIT")
        megahit_group.add_argument(
            "--k_list",
            action="store",
            type=str,
            dest="k_list",
            help="Comma-separated list of kmer sizes, all must be odd values in the range 15-255, in"
            " increments of at most 28. If not provided, a list optimized for hybridization"
            " capture and/or genome skimming data will be used. The final kmer size will be"
            " adjusted automatically so it never exceeds the mean read length of the sample by"
            " more than 31",
        )
        megahit_group.add_argument(
            "--min_count",
            action="store",
            type=int,
            dest="min_count",
            help="Minimum contig depth (a.k.a. multiplicity in MEGAHIT), accepted values are"
            " integers >= 1. Reducing it to 1 may increase the amount of low-depth contigs"
            " likely produced from reads with errors. Increase above 2 if the data has high"
            " and even sequencing depth",
        )
        megahit_group.add_argument(
            "--prune_level",
            action="store",
            type=int,
            dest="prune_level",
            help="Prunning strength for low-coverage edges during graph cleaning. Increasing the"
            " value beyond 2 can speed up the assembly at the cost of losing low-depth contigs."
            " Accepted values are integers between 0 and 3",
        )
        megahit_group.add_argument(
            "--merge_level",
            action="store",
            default="20,0.95",
            type=str,
            dest="merge_level",
            help="Merge complex bubbles, the first number multiplied by the kmer size represents the"
            " maximum bubble length to merge, the second number represents the minimum"
            " similarity required to merge bubbles",
        )
        megahit_group.add_argument(
            "--preset",
            action="store",
            choices=[
                "CAPSKIM",
                "RNA",
                "WGS",
            ],
            default="CAPSKIM",
            type=str,
            dest="preset",
            help="B|The default preset is 'CAPSKIM', these settings work well with either"
            " hybridization capture or genome skimming data (or a combination of both). You can"
            " assemble RNA-Seq reads with the 'RNA' preset or high-coverage Whole Genome"
            " Sequencing reads with the 'WGS' preset, however, both presets require a minimum"
            " of 8GB of RAM to work well (default: CAPSKIM)\n"
            f"CAPSKIM = --k-list {settings.MEGAHIT_PRESETS['CAPSKIM']['k_list']}"
            f" --min-count {settings.MEGAHIT_PRESETS['CAPSKIM']['min_count']}"
            f" --prune-level {settings.MEGAHIT_PRESETS['CAPSKIM']['prune_level']}\n"
            f"RNA = --k-list {settings.MEGAHIT_PRESETS['RNA']['k_list']}"
            f" --min-count {settings.MEGAHIT_PRESETS['RNA']['min_count']}"
            f" --prune-level {settings.MEGAHIT_PRESETS['RNA']['prune_level']}\n"
            f"WGS = --k-list {settings.MEGAHIT_PRESETS['WGS']['k_list']}"
            f" --min-count {settings.MEGAHIT_PRESETS['WGS']['min_count']}"
            f" --prune-level {settings.MEGAHIT_PRESETS['WGS']['prune_level']}\n",
        )
        megahit_group.add_argument(
            "--min_contig_len",
            action="store",
            default="auto",
            type=str,
            dest="min_contig_len",
            help="Minimum contig length in output assembly, 'auto' is mean input read length +"
            " smallest kmer size in '--k_list'",
        )
        megahit_group.add_argument(
            "--tmp_dir",
            action="store",
            default="$HOME",
            type=str,
            dest="tmp_dir",
            help="Location to create the temporary directory 'captus_assembly_tmp' for MEGAHIT"
            " assembly. Sometimes, when working on external hard drives MEGAHIT will refuse to"
            " run unless this directory is created in an internal hard drive",
        )
        megahit_group.add_argument(
            "--max_contig_gc",
            action="store",
            default=100.0,
            type=float,
            dest="max_contig_gc",
            help="Maximum GC percentage allowed per contig. Useful to filter contamination. For"
            " example, bacteria usually exceed 60%% GC content while eukaryotes rarely exceed"
            " that limit. 100.0 disables the GC filter",
        )
        megahit_group.add_argument(
            "--disable_mapping",
            action="store_true",
            dest="disable_mapping",
            help="Disable mapping the reads back to the contigs using Salmon for accurate depth"
            " estimation. If disabled, the approximate depth estimation given by MEGAHIT will"
            " be used instead",
        )
        megahit_group.add_argument(
            "--min_contig_depth",
            action="store",
            default="auto",
            type=str,
            dest="min_contig_depth",
            help="Minimum contig depth of coverage in output assembly; 'auto' will retain contigs"
            " with depth of coverage greater than 1.0x when '--disable_mapping' is chosen,"
            " otherwise it will retain only contigs of at least 1.5x. Accepted values are"
            " decimals greater or equal to 0. Use 0 to disable the filter",
        )
        megahit_group.add_argument(
            "--redo_filtering",
            action="store_true",
            dest="redo_filtering",
            help="Enable if you want to try different values for `--max_contig_gc` or"
            " `--min_contig_depth`. Only the filtering step will be repeated",
        )

        other_group = parser.add_argument_group("Other")
        other_group.add_argument(
            "--reformat_path",
            action="store",
            default="reformat.sh",
            type=str,
            dest="reformat_path",
            help="Path to reformat.sh",
        )
        other_group.add_argument(
            "--megahit_path",
            action="store",
            default="megahit",
            type=str,
            dest="megahit_path",
            help="Path to MEGAHIT",
        )
        other_group.add_argument(
            "--megahit_toolkit_path",
            action="store",
            default="megahit_toolkit",
            type=str,
            dest="megahit_toolkit_path",
            help="Path to MEGAHIT's toolkit",
        )
        other_group.add_argument(
            "--salmon_path",
            action="store",
            default="salmon",
            type=str,
            dest="salmon_path",
            help="Path to Salmon",
        )
        other_group.add_argument(
            "--ram",
            action="store",
            default="auto",
            type=str,
            dest="ram",
            help="Maximum RAM in GB (e.g.: 4.5) dedicated to Captus, 'auto' uses"
            f" {settings.RAM_FRACTION:.0%}% of available RAM",
        )
        other_group.add_argument(
            "--threads",
            action="store",
            default="auto",
            type=str,
            dest="threads",
            help="Maximum number of CPUs dedicated to Captus, 'auto' uses all available CPUs",
        )
        other_group.add_argument(
            "--concurrent",
            action="store",
            default="auto",
            type=str,
            dest="concurrent",
            help="Captus will attempt to assemble this many samples concurrently. RAM and CPUs"
            " will be divided by this value for each individual MEGAHIT process. For example"
            " if you set --threads to 12 and --concurrent to 3, then each MEGAHIT assembly will"
            " be done using --threads/--concurrent = 4 CPUs. If set to 'auto', Captus will run"
            " as many concurrent assemblies as possible with a minimum of 4 CPUs and 4 GB of"
            " RAM per assembly (8 GB if presets RNA or WGS are used)",
        )
        other_group.add_argument(
            "--debug",
            action="store_true",
            dest="debug",
            help="Enable debugging mode, parallelization is disabled so errors are logged to screen",
        )
        other_group.add_argument(
            "--show_less",
            action="store_true",
            dest="show_less",
            help="Do not show individual sample information during the run, the information is still"
            " written to the log",
        )

        help_group = parser.add_argument_group("Help")
        help_group.add_argument("-h", "--help", action="help", help="Show this help message and exit")
        help_group.add_argument(
            "--version",
            action="version",
            version=f"Captus v{__version__}",
            help="Show Captus' version number",
        )

        if len(sys.argv) == 2:
            parser.print_help()
            exit(red("\nERROR: Missing required argument -r/--reads\n"))

        full_command = " ".join(sys.argv)
        args = parser.parse_args(sys.argv[2:])
        assemble(full_command, args)
        exit(1)

    ################################################################################################
    ################################################################################ EXTRACT SECTION
    def extract(self):
        description = bold("Captus-assembly: Extract; recover markers from FASTA assemblies\n")
        parser = argparse.ArgumentParser(
            usage="captus extract -a CAPTUS_ASSEMBLIES_DIR [options]",
            description=description,
            formatter_class=MyHelpFormatter,
            epilog="For more information, please see https://github.com/edgardomortiz/Captus",
            add_help=False,
        )

        input_group = parser.add_argument_group("Input")
        input_group.add_argument(
            "-a",
            "--captus_assemblies_dir",
            action="store",
            default="./02_assemblies",
            type=str,
            required=True,
            dest="captus_assemblies_dir",
            help="Path to an output directory from the 'assemble' step of Captus-assembly which is"
            " tipically called '02_assemblies'. If you DID NOT assemble any sample within"
            " Captus and want to start exclusivey with FASTAs assembled elsewhere, the path"
            " provided here will be created in order to contain your assemblies provided with"
            " '-f' into a proper directory structure needed by Captus",
        )
        input_group.add_argument(
            "-f",
            "--fastas",
            action="store",
            nargs="+",
            type=str,
            dest="fastas",
            help="B|FASTA assembly file(s) that were not assembled with Captus. Valid file name"
            " extensions are: .fa, .fna, .fasta, .fa.gz, .fna.gz, .fasta.gz. These FASTA files"
            " must contain only nucleotides (no aminoacids). All the text before the extension"
            " of the filename will be used as sample name. These FASTAs will be automatically"
            " copied to the path provided with '-a'/'--captus_assemblies_dir' using the correct"
            " directory structure needed by Captus. There are a few ways to provide the FASTA"
            " files:\n"
            "A directory = path to directory containing FASTA files (e.g.: -f ./my_fastas)\n"
            "A list = filenames separated by space (e.g.: -f speciesA.fa speciesB.fasta.gz)\n"
            "A pattern = UNIX matching expression (e.g.: -f ./my_fastas/*.fasta)",
        )

        output_group = parser.add_argument_group("Output")
        output_group.add_argument(
            "-o",
            "--out",
            action="store",
            default="./03_extractions",
            type=str,
            dest="out",
            help="Output directory name",
        )
        output_group.add_argument(
            "--ignore_depth",
            action="store_true",
            dest="ignore_depth",
            help="Do not filter contigs based on their depth of coverage",
        )
        output_group.add_argument(
            "--disable_stitching",
            action="store_true",
            dest="disable_stitching",
            help="Use only if you are sure your target loci will be found in a single contig (for"
            " example if you have a chromosome-level assembly), otherwise Captus tries to join"
            " partial matches to a target that are scattered across multiple contigs if their"
            " structure and overlaps are compatible",
        )
        output_group.add_argument(
            "--max_paralogs",
            action="store",
            default=-1,
            type=int,
            dest="max_paralogs",
            help="Maximum number of secondary hits (copies) of any particular reference target"
            " marker allowed in the output. We recommend disabling the removal of paralogs"
            " (secondary hits/copies) during the 'extract' step because the 'align' step uses a"
            " more sophisticated filter for paralogs. -1 disables the removal of paralogs",
        )
        output_group.add_argument(
            "--paralog_tolerance",
            action="store",
            default=5.0,
            type=float,
            dest="paralog_tolerance",
            help="Only paralogs with a wscore >= (locus best hit wscore / paralog_tolerance) will be"
            " retained",
        )
        output_group.add_argument(
            "--max_loci_files",
            action="store",
            default=0,
            type=int,
            dest="max_loci_files",
            help="When the number of markers in a reference target file exceeds this number, Captus"
            " will not write a separate FASTA file per sample per marker to not overload I/O."
            " The single FASTA file containing all recovered markers per sample needed by the"
            " 'align' step is still produced as are the rest of output files",
        )
        output_group.add_argument(
            "--keep_all",
            action="store_true",
            dest="keep_all",
            help="Do not delete any intermediate files",
        )
        output_group.add_argument(
            "--overwrite", action="store_true", dest="overwrite", help="Overwrite previous results"
        )

        protein_group = parser.add_argument_group("Proteins extraction global options (Scipio)")
        protein_group.add_argument(
            "--max_loci_scipio_x2",
            action="store",
            default=2000,
            type=int,
            dest="max_loci_scipio_x2",
            help="When the number of loci in a protein reference target file exceeds this number,"
            " Captus will not run a second, more exhaustive round of Scipio. Usually the"
            " results from the first round are extremely similar and sufficient, the second"
            " round can become extremely slow as the number of reference target proteins grows",
        )
        protein_group.add_argument(
            "--predict",
            action="store_true",
            dest="predict",
            help="Scipio flags introns as dubious when the splice signals are not found at the"
            " exon edges, this may indicate that there are additional aminoacids in the"
            " recovered protein that are not present in the reference target protein. Enable"
            " this flag to attempt translation of these dubious introns, if the translation"
            " does not introduce premature stop codons they will be added to the recovered"
            " protein",
        )

        scipio_nuc_group = parser.add_argument_group("Nuclear proteins extraction (Scipio)")
        scipio_nuc_group.add_argument(
            "-n",
            "--nuc_refs",
            action="store",
            type=str,
            dest="nuc_refs",
            help="B|Set of nuclear protein reference target sequences, options are:\n"
            "Angiosperms353 = The original set of target proteins from Angiosperms353\n"
            "Mega353 = The improved set of target proteins from Angiosperms353\n"
            "Alternatively, provide a path to a FASTA file containing your reference target"
            " protein sequences in either nucleotide or aminoacid. When the FASTA file is in"
            " nucleotides, '--nuc_transtable' will be used to translate it to aminoacids",
        )
        scipio_nuc_group.add_argument(
            "--nuc_transtable",
            action="store",
            default=settings.DEFAULT_GENETIC_CODES["NUC"]["id"],
            type=int,
            dest="nuc_transtable",
            help="Genetic code table to translate your nuclear proteins. Complete list of tables at:"
            " https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi (default:"
            f" {settings.DEFAULT_GENETIC_CODES['NUC']['id']}:"
            f" {settings.DEFAULT_GENETIC_CODES['NUC']['name']})",
        )
        scipio_nuc_group.add_argument(
            "--nuc_min_score",
            action="store",
            default=0.13,
            type=float,
            dest="nuc_min_score",
            help="Minimum Scipio score to retain hits to reference target proteins.",
        )
        scipio_nuc_group.add_argument(
            "--nuc_min_identity",
            action="store",
            default=65,
            type=float,
            dest="nuc_min_identity",
            help="Minimum identity percentage to retain hits to reference target proteins",
        )
        scipio_nuc_group.add_argument(
            "--nuc_min_coverage",
            action="store",
            default=20,
            type=float,
            dest="nuc_min_coverage",
            help="Minimum coverage percentage of reference target protein to consider a hit by a contig",
        )
        scipio_nuc_group.add_argument(
            "--nuc_depth_tolerance",
            action="store",
            default=2.0,
            type=float,
            dest="nuc_depth_tolerance",
            help="Minimum depth = 10^(log(depth of contig with best hit in locus) / "
            "nuc_depth_tolerance), values must be >= 1.0 (1 is the most strict)",
        )

        scipio_ptd_group = parser.add_argument_group("Plastidial proteins extraction (Scipio)")
        scipio_ptd_group.add_argument(
            "-p",
            "--ptd_refs",
            action="store",
            type=str,
            dest="ptd_refs",
            help="B|Set of plastidial protein reference target sequences, options are:\n"
            "SeedPlantsPTD = A set of plastidial proteins for Seed Plants, curated by us\n"
            "Alternatively, provide a path to a FASTA file containing your reference target"
            " protein sequences in either nucleotide or aminoacid. When the FASTA file is in"
            " nucleotides, '--ptd_transtable' will be used to translate it to aminoacids",
        )
        scipio_ptd_group.add_argument(
            "--ptd_transtable",
            action="store",
            default=settings.DEFAULT_GENETIC_CODES["PTD"]["id"],
            type=int,
            dest="ptd_transtable",
            help="Genetic code table to translate your plastidial proteins. Complete list of tables"
            " at: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi (default:"
            f" {settings.DEFAULT_GENETIC_CODES['PTD']['id']}:"
            f" {settings.DEFAULT_GENETIC_CODES['PTD']['name']})",
        )
        scipio_ptd_group.add_argument(
            "--ptd_min_score",
            action="store",
            default=0.2,
            type=float,
            dest="ptd_min_score",
            help="Minimum Scipio score to retain hits to reference target proteins",
        )
        scipio_ptd_group.add_argument(
            "--ptd_min_identity",
            action="store",
            default=65,
            type=float,
            dest="ptd_min_identity",
            help="Minimum identity percentage to retain hits to reference target proteins",
        )
        scipio_ptd_group.add_argument(
            "--ptd_min_coverage",
            action="store",
            default=20,
            type=float,
            dest="ptd_min_coverage",
            help="Minimum coverage percentage of reference target protein to consider a hit by a contig",
        )
        scipio_ptd_group.add_argument(
            "--ptd_depth_tolerance",
            action="store",
            default=2.0,
            type=float,
            dest="ptd_depth_tolerance",
            help="Minimum depth = 10^(log(depth of contig with best hit in locus) / "
            "ptd_depth_tolerance), values must be >= 1.0 (1 is the most strict)",
        )

        scipio_mit_group = parser.add_argument_group("Mitochondrial proteins extraction (Scipio)")
        scipio_mit_group.add_argument(
            "-m",
            "--mit_refs",
            action="store",
            type=str,
            dest="mit_refs",
            help="B|Set of mitochondrial protein reference target sequences, options are:\n"
            "SeedPlantsMIT = A set of mitochondrial proteins for Seed Plants, curated by us\n"
            "Alternatively, provide a path to a FASTA file containing your reference target"
            " protein sequences in either nucleotide or aminoacid. When the FASTA file is in"
            " nucleotides, '--mit_transtable' will be used to translate it to aminoacids",
        )
        scipio_mit_group.add_argument(
            "--mit_transtable",
            action="store",
            default=settings.DEFAULT_GENETIC_CODES["MIT"]["id"],
            type=int,
            dest="mit_transtable",
            help="Genetic code table to translate your mitochondrial proteins. Complete list of"
            " tables at: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi (default:"
            f" {settings.DEFAULT_GENETIC_CODES['MIT']['id']}:"
            f" {settings.DEFAULT_GENETIC_CODES['MIT']['name']})",
        )
        scipio_mit_group.add_argument(
            "--mit_min_score",
            action="store",
            default=0.2,
            type=float,
            dest="mit_min_score",
            help="Minimum Scipio score to retain hits to reference target proteins",
        )
        scipio_mit_group.add_argument(
            "--mit_min_identity",
            action="store",
            default=65,
            type=float,
            dest="mit_min_identity",
            help="Minimum identity percentage to retain hits to reference target proteins",
        )
        scipio_mit_group.add_argument(
            "--mit_min_coverage",
            action="store",
            default=20,
            type=float,
            dest="mit_min_coverage",
            help="Minimum coverage percentage of reference target protein to consider a hit by a contig",
        )
        scipio_mit_group.add_argument(
            "--mit_depth_tolerance",
            action="store",
            default=2.0,
            type=float,
            dest="mit_depth_tolerance",
            help="Minimum depth = 10^(log(depth of contig with best hit in locus) / "
            "mit_depth_tolerance), values must be >= 1.0 (1 is the most strict)",
        )

        non_coding_group = parser.add_argument_group("Miscellaneous DNA extraction (BLAT)")
        non_coding_group.add_argument(
            "-d",
            "--dna_refs",
            action="store",
            type=str,
            dest="dna_refs",
            help="Path to a FASTA nucleotide file of miscellaneous DNA reference targets",
        )
        non_coding_group.add_argument(
            "--dna_min_identity",
            action="store",
            default=80,
            type=float,
            dest="dna_min_identity",
            help="Minimum identity percentage to reference target sequences to retain matches",
        )
        non_coding_group.add_argument(
            "--dna_min_coverage",
            action="store",
            default=20,
            type=float,
            dest="dna_min_coverage",
            help="Minimum coverage percentage of reference target sequence to retain matches",
        )
        non_coding_group.add_argument(
            "--dna_depth_tolerance",
            action="store",
            default=2.0,
            type=float,
            dest="dna_depth_tolerance",
            help="Minimum depth = 10^(log(depth of contig with best hit in locus) / "
            "dna_depth_tolerance), values must be >= 1.0 (1 is the most strict)",
        )

        mmseqs2_group = parser.add_argument_group("Assemblies clustering (MMseqs2)")
        mmseqs2_group.add_argument(
            "-c",
            "--cluster_leftovers",
            action="store_true",
            dest="cluster_leftovers",
            help="Enable MMseqs2 clustering across samples of the contigs that had no hits to the"
            " reference target markers. A new miscellaneous DNA reference is built from the"
            " best representative of each cluster in order to perform a miscellaneous DNA"
            " marker extraction.",
        )
        mmseqs2_group.add_argument(
            "--mmseqs_method",
            action="store",
            choices=[
                "easy-linclust",
                "easy-cluster",
            ],
            default="easy-linclust",
            type=str,
            dest="mmseqs_method",
            help="B|MMseqs2 clustering algorithm, options are:\n"
            "easy-linclust = Fast linear time (for huge datasets), less sensitive clustering\n"
            "easy-cluster = Sensitive homology search (recommended but slower)",
        )
        mmseqs2_group.add_argument(
            "--cl_mode",
            action="store",
            default=2,
            type=int,
            dest="cl_mode",
            choices=[0, 1, 2],
            help="B|MMseqs2 clustering mode (https://github.com/soedinglab/mmseqs2/wiki#clustering-"
            "modes), options are:\n"
            "0 = Greedy set cover\n"
            "1 = Connected component\n"
            "2 = Greedy incremental (analogous to CD-HIT)",
        )
        mmseqs2_group.add_argument(
            "--cl_sensitivity",
            action="store",
            default=settings.MMSEQS_SENSITIVITY,
            type=float,
            dest="cl_sensitivity",
            help="MMseqs2 sensitivity, from 1 to 7.5, only applicable when using 'easy-cluster'."
            " Common reference points are: 1 (faster), 4 (fast), 7.5 (sens)",
        )
        mmseqs2_group.add_argument(
            "--cl_min_identity",
            action="store",
            default="auto",
            type=str,
            dest="cl_min_identity",
            help="Minimum identity percentage between sequences in a cluster, when set to 'auto'"
            f" it becomes {settings.MMSEQS_BLAT_DNA_IDENTITY_FACTOR:.0%}% of the '--dna_min"
            f"_identity' value but never less than {settings.MMSEQS_MIN_AUTO_MIN_IDENTITY}%%",
        )
        mmseqs2_group.add_argument(
            "--cl_seq_id_mode",
            action="store",
            default=1,
            type=int,
            dest="cl_seq_id_mode",
            help="B|MMseqs2 sequence identity mode, options are:\n"
            "0 = Alignment length\n"
            "1 = Shorter sequence\n"
            "2 = Longer sequence",
        )
        mmseqs2_group.add_argument(
            "--cl_min_coverage",
            action="store",
            default=80,
            type=float,
            dest="cl_min_coverage",
            help="Any sequence in a cluster has to be at least this percent included in the length"
            " of the longest sequence in the cluster",
        )
        mmseqs2_group.add_argument(
            "--cl_cov_mode",
            action="store",
            default=settings.MMSEQS_COV_MODE,
            type=int,
            dest="cl_cov_mode",
            help="MMseqs2 sequence coverage mode (https://github.com/soedinglab/mmseqs2/wiki#how-"
            "to-set-the-right-alignment-coverage-to-cluster)",
        )
        mmseqs2_group.add_argument(
            "--cl_min_samples",
            action="store",
            default="auto",
            type=str,
            dest="cl_min_samples",
            help="Minimum number of samples per cluster, if set to 'auto' the number is adjusted to"
            f" {settings.CLR_MIN_SAMPLE_PROP:.0%}% of the total number of samples or at least 4",
        )
        mmseqs2_group.add_argument(
            "--cl_rep_min_len",
            action="store",
            default=500,
            type=int,
            dest="cl_rep_min_len",
            help="After clustering is finished, only accept cluster representatives of at least this"
            " length to be part of the new miscellaneous DNA reference targets. Use 0 to"
            " disable this filter",
        )
        mmseqs2_group.add_argument(
            "--cl_max_seq_len",
            action="store",
            default=5000,
            type=int,
            dest="cl_max_seq_len",
            help="Do not cluster sequences longer than this length in bp, the maximum allowed by"
            " MMseqs2 is 65535. Use 0 to disable this filter",
        )
        mmseqs2_group.add_argument(
            "--cl_max_copies",
            action="store",
            default=3,
            type=float,
            dest="cl_max_copies",
            help="Maximum average number of sequences per sample in a cluster. This can exclude loci"
            " that are extremely paralogous",
        )
        mmseqs2_group.add_argument(
            "--cl_tmp_dir",
            action="store",
            default="$HOME",
            type=str,
            dest="cl_tmp_dir",
            help="Where to create the temporary directory 'captus_mmseqs_tmp' for MMseqs2."
            " Clustering can become slow when done on external drives, set this location to a"
            " fast, preferably local, drive",
        )

        other_group = parser.add_argument_group("Other")
        other_group.add_argument(
            "--scipio_path",
            action="store",
            default="bundled",
            type=str,
            dest="scipio_path",
            help="Path to Scipio",
        )
        other_group.add_argument(
            "--blat_path",
            action="store",
            default="bundled",
            type=str,
            dest="blat_path",
            help="Path to BLAT >= 36x7 (this version is the first one that guarantees the same"
            " result both in Mac and Linux)",
        )
        other_group.add_argument(
            "--mmseqs_path",
            action="store",
            default="mmseqs",
            type=str,
            dest="mmseqs_path",
            help="Path to MMseqs2",
        )
        other_group.add_argument(
            "--mafft_path",
            action="store",
            default="mafft",
            type=str,
            dest="mafft_path",
            help="Path to MAFFT (default: mafft/mafft.bat)",
        )
        other_group.add_argument(
            "--ram",
            action="store",
            default="auto",
            type=str,
            dest="ram",
            help="Maximum RAM in GB (e.g.: 4.5) dedicated to Captus, 'auto' uses"
            f" {settings.RAM_FRACTION:.0%}% of available RAM",
        )
        other_group.add_argument(
            "--threads",
            action="store",
            default="auto",
            type=str,
            dest="threads",
            help="Maximum number of CPUs dedicated to Captus, 'auto' uses all available CPUs",
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
            " each process due to the RAM requirements of BLAT",
        )
        other_group.add_argument(
            "--debug",
            action="store_true",
            dest="debug",
            help="Enable debugging mode, parallelization is disabled so errors are logged to screen",
        )
        other_group.add_argument(
            "--show_less",
            action="store_true",
            dest="show_less",
            help="Do not show individual sample information during the run, the information is still"
            " written to the log",
        )

        help_group = parser.add_argument_group("Help")
        help_group.add_argument("-h", "--help", action="help", help="Show this help message and exit")
        help_group.add_argument(
            "--version",
            action="version",
            version=f"Captus v{__version__}",
            help="Show Captus' version number",
        )

        if len(sys.argv) == 2:
            parser.print_help()
            exit(red("\nERROR: Missing required argument -a/--captus_assemblies_dir\n"))

        full_command = " ".join(sys.argv)
        args = parser.parse_args(sys.argv[2:])
        extract(full_command, args)
        exit(1)

    ################################################################################################
    ################################################################################## ALIGN SECTION
    def align(self):
        description = bold("Captus-assembly: Align; collect, align, and curate aligned markers\n")
        parser = argparse.ArgumentParser(
            usage="captus align -e CAPTUS_EXTRACTIONS_DIR [options]",
            description=description,
            formatter_class=MyHelpFormatter,
            epilog="For more information, please see https://github.com/edgardomortiz/Captus",
            add_help=False,
        )

        input_group = parser.add_argument_group("Input")
        input_group.add_argument(
            "-e",
            "--captus_extractions_dir",
            action="store",
            default="./03_extractions",
            type=str,
            required=True,
            dest="captus_extractions_dir",
            help="Path to the output directory that contains the assemblies and extractions from"
            " previous steps of Captus-assembly. This directory is called '02_assemblies' if"
            " you did not specify a different name during the 'assemble' or 'extract' steps",
        )
        input_group.add_argument(
            "-m",
            "--markers",
            action="store",
            default="all",
            type=str,
            dest="markers",
            help="B|Which markers to align, you can provide a comma-separated list, no spaces"
            " (default: all)\n"
            "NUC = Nuclear proteins inside directories '01_coding_NUC'\n"
            "PTD = Plastidial proteins inside directories '02_coding_PTD'\n"
            "MIT = Mitochondrial proteins inside directories '03_coding_MIT'\n"
            "DNA = Miscellaneous DNA markers inside directories '04_misc_DNA'\n"
            "CLR = Cluster-derived DNA markers inside directories '05_clusters'\n"
            "ALL = Shortcut for NUC,PTD,MIT,DNA,CLR",
        )
        input_group.add_argument(
            "-f",
            "--formats",
            action="store",
            default="AA,NT,GE,MA",
            type=str,
            dest="formats",
            help="B|Which alignment format(s) to prepare for each marker category, you can provide a"
            " comma-separated list, no spaces (default: AA,NT,GE,MA)\n"
            "Valid types for NUC, PTD, and MIT markers:\n"
            "AA = Coding sequences in aminoacids\n"
            "NT = Coding sequences in nucleotides\n"
            "GE = Complete gene sequences (exons + introns) without flanking upstream or"
            " downstream basepairs\n"
            "GF = Complete gene sequences with flanking upstream and downstream basepairs\n"
            "Valid types for miscellaneous DNA and CLusteR-derived markers:\n"
            "MA = Matched sequences without flanking upstream or downstream basepairs\n"
            "MF = Matched sequences with flanking upstream and downstream basepairs\n"
            "ALL = Shortcut for AA,NT,GE,GF,MA,MF",
        )
        input_group.add_argument(
            "--max_paralogs",
            action="store",
            default=5,
            type=int,
            dest="max_paralogs",
            help="Maximum number of secondary hits (copies) per sample to import from the"
            " extraction step. Large numbers of marker copies per sample can increase"
            " alignment times. Hits (copies) are ranked from best to worst during the"
            " 'extract' step. -1 disables the removal of paralogs and aligns them all, which"
            " might be useful if you expect very high ploidy levels for example",
        )
        input_group.add_argument(
            "-s",
            "--min_samples",
            action="store",
            default=4,
            type=int,
            dest="min_samples",
            help="Minimum number of samples in a marker to proceed with alignment. Markers with"
            " fewer samples will be skipped",
        )

        output_group = parser.add_argument_group("Output")
        output_group.add_argument(
            "-o",
            "--out",
            action="store",
            default="./04_alignments",
            type=str,
            dest="out",
            help="Output directory name",
        )
        output_group.add_argument(
            "--collect_only",
            action="store_true",
            dest="collect_only",
            help="Only collect the markers from the extraction folder and exit (skips addition of"
            " reference target sequences and subsequent steps)",
        )
        output_group.add_argument(
            "--keep_w_refs",
            action="store_true",
            dest="keep_w_refs",
            help="Keep the directories containing the alignments that include target reference"
            " sequences. The deletion of these directories is performed after alignments are"
            " filtered for paralogs and before trimming. Enable if you plan to repeat the paralog"
            " filtering in the future with '--redo_from filtering'",
        )
        output_group.add_argument(
            "--keep_all",
            action="store_true",
            dest="keep_all",
            help="Do not delete any intermediate files",
        )
        output_group.add_argument(
            "--overwrite", action="store_true", dest="overwrite", help="Overwrite previous results"
        )

        align_group = parser.add_argument_group("Alignment")
        align_group.add_argument(
            "--align_method",
            action="store",
            choices=[
                "mafft_auto",
                "mafft_genafpair",
                "mafft_localpair",
                "mafft_globalpair",
                "mafft_retree1",
                "mafft_retree2",
                "muscle_align",
                "muscle_super5",
            ],
            default="mafft_auto",
            type=str,
            dest="align_method",
            help="For MAFFT's algorithms see:"
            " https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html"
            " For MUSCLE's algorithms see:"
            " https://drive5.com/muscle5/manual/commands.html",
        )
        align_group.add_argument(
            "--timeout",
            action="store",
            default=21600,
            type=int,
            dest="timeout",
            help="Maximum allowed time in seconds for a single alignment",
        )
        align_group.add_argument(
            "--disable_codon_align",
            action="store_true",
            dest="disable_codon_align",
            help="Do not align nucleotide coding sequences based on their protein alignment",
        )
        align_group.add_argument(
            "--outgroup",
            action="store",
            default=None,
            type=str,
            dest="outgroup",
            help="Outgroup sample names, separated by commas, no spaces. Since phylogenetic programs"
            " usually root the resulting trees at the first taxon in the alignment, Captus will"
            " place these samples at the beggining of every alignment in the given order",
        )

        paralog_group = parser.add_argument_group("Paralog filtering")
        paralog_group.add_argument(
            "--filter_method",
            action="store",
            choices=["naive", "informed", "both", "none"],
            default="both",
            type=str,
            dest="filter_method",
            help="B|Methods for filtering paralogous sequences:\n"
            "naive = Only the best hit for each sample (marked as hit=00) is retained\n"
            "informed = Only keep the copy (regardless of hit ranking) that is most similar to"
            " the reference target sequence that was chosen most frequently among all other"
            " samples in the alignment"
            "both = Two separate folders will be created, each containing the results from each"
            " filtering method"
            "none = Skip paralog removal, just remove reference target sequences from the"
            " alignments. Useful for phylogenetic methods that allow paralogs like ASTRAL-Pro",
        )
        paralog_group.add_argument(
            "--tolerance",
            action="store",
            default=4.0,
            type=float,
            dest="tolerance",
            help="Only applicable to the 'informed' filter. If the selected copy's identity to the"
            " most commonly chosen reference target  is below this number of Standard"
            " Deviations from the mean, it will also be removed (the lower the number the"
            " stricter the filter). Use -1 to disable the filter",
        )

        trimming_group = parser.add_argument_group("Trimming")
        trimming_group.add_argument(
            "-t",
            "--taper_cutoff",
            action="store",
            default=3.0,
            type=float,
            dest="taper_cutoff",
            help="TAPER cutoff threshold, values greater than 1.0 are recommended, the lower the"
            " value the more aggressive the correction, 3.0 recommended by TAPER's authors",
        )
        trimming_group.add_argument(
            "--taper_conservative",
            action="store_true",
            dest="taper_conservative",
            help="Enable the more conservative mode of TAPER. Captus uses the aggressive mode by"
            " default, see 'correction_multi_aggressive.jl' at https://github.com/chaoszhang/TAPER",
        )
        trimming_group.add_argument(
            "--taper_unfiltered",
            action="store_true",
            dest="taper_unfiltered",
            help="Enable TAPER correction even for alignments than have not been paralog-filtered,"
            " TAPER is only able to distinguish error when an unfiltered alignment contains copies"
            " of the locus that are not extremely divergent",
        )
        trimming_group.add_argument(
            "--disable_taper",
            action="store_true",
            dest="disable_taper",
            help="Disable the TAPER algorithm for masking for erroneous regions in alignments, see"
            " https://doi.org/10.1111/2041-210X.13696",
        )
        trimming_group.add_argument(
            "--clipkit_method",
            action="store",
            choices=[
                "smart-gap",
                "gappy",
                "kpic",
                "kpic-smart-gap",
                "kpic-gappy",
                "kpi",
                "kpi-smart-gap",
                "kpi-gappy",
            ],
            default="gappy",
            type=str,
            dest="clipkit_method",
            help="ClipKIT's algorithm, see https://jlsteenwyk.com/ClipKIT/advanced/index.html#modes",
        )
        trimming_group.add_argument(
            "-g",
            "--clipkit_gaps",
            action="store",
            default=0.90,
            type=float,
            dest="clipkit_gaps",
            help="Gappyness threshold per position. Accepted values between 0 and 1. This argument"
            " is ignored when using the 'kpi' and 'kpic' algorithms or intermediate steps that"
            " use 'smart-gap', the higher the value the more columns are preserved",
        )
        trimming_group.add_argument(
            "-d",
            "--min_data_per_column",
            action="store",
            default=0,
            type=int,
            dest="min_data_per_column",
            help="Minimum number of non-missing sites per column. When this parameter is > 0, Captus"
            " will dynamically calculate a '--clipkit_gaps' threshold per alignment to keep"
            " this minimum amount of data per column",
        )
        trimming_group.add_argument(
            "--ends_only",
            action="store_true",
            dest="ends_only",
            help="Trim only the ends of the alignments (do not trim internal gaps)",
        )
        trimming_group.add_argument(
            "-c",
            "--min_coverage",
            action="store",
            default=0.40,
            type=float,
            dest="min_coverage",
            help="Minimum coverage of sequence as proportion of the mean of sequence lengths in the"
            " alignment, ignoring gaps. Accepted values between 0 and 1. After ClipKIT finishes"
            " trimming columns, Captus will also remove short sequences below this threshold",
        )

        other_group = parser.add_argument_group("Other")
        other_group.add_argument(
            "--redo_from",
            action="store",
            choices=[
                "alignment",
                "filtering",
                "removal",
                "trimming",
            ],
            type=str,
            dest="redo_from",
            help="B|Repeat analysis from a particular stage:\n"
            "alignment = Delete all subdirectories with alignments and restart\n"
            "filtering = Delete all subdirectories with filtered alignments and restart\n"
            "removal = Delete all subdirectories with alignments with reference targets removed"
            " and restart\n"
            "trimming = Delete all subdirectories with trimmed alignments and restart",
        )
        other_group.add_argument(
            "--mafft_path",
            action="store",
            default="mafft",
            type=str,
            dest="mafft_path",
            help="Path to MAFFT (default: mafft/mafft.bat)",
        )
        other_group.add_argument(
            "--muscle_path",
            action="store",
            default="muscle",
            type=str,
            dest="muscle_path",
            help="Path to MUSCLE (default: muscle)",
        )
        other_group.add_argument(
            "--clipkit_path",
            action="store",
            default="clipkit",
            type=str,
            dest="clipkit_path",
            help="Path to ClipKIT",
        )
        other_group.add_argument(
            "--ram",
            action="store",
            default="auto",
            type=str,
            dest="ram",
            help="Maximum RAM in GB (e.g.: 4.5) dedicated to Captus, 'auto' uses"
            f" {settings.RAM_FRACTION:.0%}% of available RAM",
        )
        other_group.add_argument(
            "--threads",
            action="store",
            default="auto",
            type=str,
            dest="threads",
            help="Maximum number of CPUs dedicated to Captus, 'auto' uses all available CPUs",
        )
        other_group.add_argument(
            "--concurrent",
            action="store",
            default="auto",
            type=str,
            dest="concurrent",
            help="Captus will attempt to execute this many alignments concurrently. CPUs will be"
            " divided by this value for each individual process. If set to 'auto', Captus will"
            " set as many processes as to at least have 2 threads available for each MAFFT"
            " or MUSCLE process",
        )
        other_group.add_argument(
            "--debug",
            action="store_true",
            dest="debug",
            help="Enable debugging mode, parallelization is disabled so errors are logged to screen",
        )
        other_group.add_argument(
            "--show_more",
            action="store_true",
            dest="show_more",
            help="Show individual alignment information during the run. Detailed information is"
            " written regardless to the log",
        )

        help_group = parser.add_argument_group("Help")
        help_group.add_argument("-h", "--help", action="help", help="Show this help message and exit")
        help_group.add_argument(
            "--version",
            action="version",
            version=f"Captus v{__version__}",
            help="Show Captus' version number",
        )

        if len(sys.argv) == 2:
            parser.print_help()
            exit(red("\nERROR: Missing required argument -a/--captus_assemblies_dir\n"))

        full_command = " ".join(sys.argv)
        args = parser.parse_args(sys.argv[2:])
        align(full_command, args)
        exit(1)

def main():
    CaptusAssembly()

if __name__ == "__main__":
    main()
