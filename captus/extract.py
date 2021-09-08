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

import json
import math
import shutil
import statistics
import subprocess
import time
from collections import OrderedDict
from pathlib import Path

from tqdm import tqdm

from . import log
from . import settings_assembly as settings
from .bioformats import (blat_misc_dna_psl_to_dict, dict_to_fasta, fasta_headers_to_spades,
                         fasta_to_dict, fasta_type, fix_premature_stops, scipio_yaml_to_dict,
                         split_mmseqs_clusters_file, translate_fasta_dict, write_gff3)
from .misc import (ElapsedTimeThread, bioperl_get_version, blat_path_version, bold, bold_green,
                   bold_yellow, compress_list_files, dim, elapsed_time, format_dep_msg,
                   has_valid_ext, is_dir_empty, make_output_dir, make_tmp_dir_within,
                   mmseqs_path_version, python_library_check, quit_with_error, red,
                   remove_formatting, scipio_path_version, set_ram, set_threads,
                   tqdm_parallel_async_run, tqdm_serial_run, yaml_perl_get_version)
from .version import __version__


def extract(full_command, args):

    captus_start = time.time()
    out_dir, out_dir_msg = make_output_dir(args.out)
    log.logger = log.Log(Path(args.out, "captus-assembly_extract.log"), stdout_verbosity_level=1)
    mar = 25  # Margin for aligning parameters and values

    ################################################################################################
    ############################################################################### STARTING SECTION
    log.log_section_header("Starting Captus-assembly: Extract", single_newline=False)
    log.log_explanation(
        "Welcome to the marker extraction step of Captus-assembly. In this step, Captus will use"
        " Scipio to search within your FASTA assemblies and recover any set of reference proteins"
        " provided through '--nuc_refs', '--ptd_refs', and/or '--mit_refs'. Captus includes some"
        " reference protein sets to work with Plants like the 'Angiosperms353' protein set for"
        " nuclear genes, and two organellar protein reference sets 'SeedPlantsPTD' for plastids,"
        " and 'SeedPlantsMIT' for mitochondria.", extra_empty_lines_after=0
    )
    log.log_explanation(
        "If you have references that are not proteins (e.g. non-coding regions, gene sequences"
        " including introns, full mRNAs, etc.) You can provide these miscellaneous DNA references"
        " with '--dna_refs' and Captus will use BLAT to find and extract matches in your assemblies.",
        extra_empty_lines_after=0
    )
    log.log_explanation(
        "Finally, if you want to explore the usefulness of the contigs that were not hit by any"
        " protein or other DNA references after the Captus' extraction process, or simply if you do"
        " not have reference sets to test, you can try the option '--cluster_leftovers' to attempt"
        " sequence clustering across samples in order to discover homologous markers.",
        extra_empty_lines_after=0
    )
    log.log_explanation("For more information, please see https://github.com/edgardomortiz/Captus")

    log.log(f'{"Captus version":>{mar}}: {bold(f"v{__version__}")}')
    log.log(f'{"Command":>{mar}}: {bold(full_command)}')
    ram_B, ram_MB, ram_GB, ram_GB_total = set_ram(args.ram)
    log.log(f'{"Max. RAM":>{mar}}: {bold(f"{ram_GB:.1f}GB")} {dim(f"(out of {ram_GB_total:.1f}GB)")}')
    threads_max, threads_total = set_threads(args.threads)
    log.log(f'{"Max. Threads":>{mar}}: {bold(threads_max)} {dim(f"(out of {threads_total})")}')
    log.log("")

    # Set up concurrency now to be used either at extraction or during clustering/extraction
    concurrent, ram_B_per_extraction = adjust_blat_concurrency(args.concurrent, threads_max, ram_B)

    # Set up software usage and actions based on given arguments
    if any([args.nuc_refs, args.ptd_refs, args.mit_refs, args.dna_refs]):
        skip_extraction = False
        _, scipio_version, scipio_status = scipio_path_version(args.scipio_path)
        _, bioperl_version, bioperl_status = bioperl_get_version()
        _, yaml_version, yaml_status = yaml_perl_get_version()
        _, blat_version, blat_status = blat_path_version(args.blat_path)
    else:
        skip_extraction = True
        scipio_version, bioperl_version, yaml_version, blat_version = [""] * 4
        scipio_status, bioperl_status, yaml_status, blat_status = ["not used"] * 4

    if args.cluster_leftovers:
        skip_clustering = False
        _, mmseqs_version, mmseqs_status = mmseqs_path_version(args.mmseqs2_path)
        _, blat_version, blat_status = blat_path_version(args.blat_path)
    else:
        skip_clustering = True
        mmseqs_version, mmseqs_status = "", "not used"

    log.log(f'{"Dependencies":>{mar}}:')
    log.log(format_dep_msg(f'{"Scipio   ":>{mar}}: ', scipio_version, scipio_status))
    branch_bioperl = "\u251C\u2500BioPerl"
    log.log(format_dep_msg(f'{branch_bioperl:>{mar}}: ', bioperl_version, bioperl_status))
    branch_yaml = "\u2514\u2500YAML   "
    log.log(format_dep_msg(f'{branch_yaml:>{mar}}: ', yaml_version, yaml_status))
    log.log(format_dep_msg(f'{"BLAT     ":>{mar}}: ', blat_version, blat_status))
    log.log(format_dep_msg(f'{"MMseqs2  ":>{mar}}: ', mmseqs_version, mmseqs_status))
    log.log("")

    log.log(f'{"Python libraries":>{mar}}:')
    numpy_found, numpy_version, numpy_status = python_library_check("numpy")
    pandas_found, pandas_version, pandas_status = python_library_check("pandas")
    plotly_found, plotly_version, plotly_status = python_library_check("plotly")
    log.log(format_dep_msg(f'{"numpy":>{mar}}: ', numpy_version, numpy_status))
    log.log(format_dep_msg(f'{"pandas":>{mar}}: ', pandas_version, pandas_status))
    log.log(format_dep_msg(f'{"plotly":>{mar}}: ', plotly_version, plotly_status))
    log.log("")

    captus_assemblies_dir, _ = make_output_dir(args.captus_assemblies_dir)
    _, asms_dir_status_msg = status_captus_assemblies_dir(captus_assemblies_dir, args.fastas, mar)
    log.log(f'{"Captus assemblies dir":>{mar}}: {bold(captus_assemblies_dir)}')
    log.log(asms_dir_status_msg)
    log.log("")

    # Last dependency verification before modifying Captus' assembly directory
    if any(dep_status == "not found" for dep_status in [
        scipio_status, bioperl_status, yaml_status, blat_status
    ]):
        skip_extraction = True
        if skip_clustering:
            quit_with_error(
                "At least one of Scipio's dependencies could not be found, please check your"
                " '--scipio_path', '--blat_path', and that your perl installation includes the"
                " modules 'BioPerl' and 'YAML'. Additionally, '--cluster_leftovers' was not enabled."
            )
        else:
            log.log(
                f"{bold('WARNING:')} At least one of Scipio's dependencies could not be found,"
                " please check your '--scipio_path', '--blat_path', and that your perl installation"
                " includes the modules 'BioPerl' and 'YAML'. Clustering with MMseqs2 will be"
                " attempted."
            )
    if mmseqs_status == "not found":
        skip_clustering = True
        if skip_extraction:
            quit_with_error(
                "MMseqs2 could not be found, please check your '--mseqs_path'. Additionally, no"
                " reference protein or nucleotide sets were provided for extraction."
            )
        else:
            log.log(
                f"{bold('WARNING:')} MMseqs2 could not be found, please check your '--mseqs_path'."
                " Extraction of reference protein sets will be attempted."
            )

    # Check and import extra FASTA assemblies
    fastas_to_extract = find_fasta_assemblies(captus_assemblies_dir, out_dir)
    if args.fastas:
        log.log_explanation(
            "Now Captus will verify the additional assembly files provided with '--fastas'. A new"
            f" sample directory will be created within '{out_dir}' for each of the new FASTA"
            " assemblies to contain the imported file"
        )
        asms_before_import = len(fastas_to_extract)
        find_and_copy_fastas(args.fastas, captus_assemblies_dir, args.overwrite,
                             threads_max, args.show_less)
        fastas_to_extract = find_fasta_assemblies(captus_assemblies_dir, out_dir)
        log.log("")
        asms_imported = len(fastas_to_extract) - asms_before_import
        log.log(f'{"Assemblies imported":>{mar}}: {bold(asms_imported)}')
    log.log(f'{"Total assemblies found":>{mar}}: {bold(len(fastas_to_extract))}')
    log.log("")
    log.log(f'{"Output directory":>{mar}}: {bold(out_dir)}')
    log.log(f'{"":>{mar}}  {dim(out_dir_msg)}')
    log.log("")


    ################################################################################################
    ############################################################################# EXTRACTION SECTION
    log.log_section_header("Protein and DNA Markers Extraction from FASTA assemblies")
    log.log_explanation(
        "Now Captus will search within your FASTA assemblies for the set of reference proteins"
        " provided by '--nuc_refs', '--ptd_refs', and/or '--mit_refs' and extract the sequences both"
        " in nucleotide and translated to aminoacid using the specified '--nuc_transtable',"
        " '--ptd_transtable' and '--mit_transtable' respectively.", extra_empty_lines_after=0
    )
    log.log_explanation(
        "Additionally, miscellaneous DNA references provided with '--dna_refs' will also be"
        " extracted."
    )
    if skip_extraction:
        log.log(red(
            "Skipping extraction step... (to enable provide any/all of these options: '--nuc_refs',"
            " '--ptd_refs', '--mit_refs', '--dna_refs', and verify that all the software"
            " dependencies are correctly installed)"
        ))
        log.log("")
    else:
        if not fastas_to_extract:
            quit_with_error(
                "Captus could not find FASTA assemblies to process, please verify your "
                "'--captus_assemblies_dir' and/or '--fastas' argument"
            )

        log.log(f'{"Concurrent extractions":>{mar}}: {bold(concurrent)}')
        log.log(f'{"RAM per extraction":>{mar}}: {bold(f"{ram_B_per_extraction / 1024 ** 3:.1f}GB")}')
        log.log(f'{"Threads per extraction":>{mar}}: {bold("1")}'
                f' {dim("(Scipio and BLAT are single-threaded)")}')
        log.log("")
        protein_refs = prepare_protein_refs(args.nuc_refs,
                                            args.ptd_refs,
                                            args.mit_refs,
                                            args.nuc_transtable,
                                            args.ptd_transtable,
                                            args.mit_transtable)
        log.log("")
        log.log(bold(f'{"Nuclear proteins":>{mar}}:'))
        log.log(f'{"reference":>{mar}}: {protein_refs["NUC"]["AA_msg"]}')
        if protein_refs["NUC"]["AA_path"]:
            nuc_query = fasta_to_dict(protein_refs["NUC"]["AA_path"], ordered=True)
            nuc_query_info = reference_info(nuc_query)
            log.log(f'{"reference info":>{mar}}: {nuc_query_info["info_msg"]}')
            log.log(f'{"min_score":>{mar}}: {bold(args.nuc_min_score)}')
            nuc_min_identity = adjust_min_identity(args.nuc_min_identity, args.nuc_transtable)
            log.log(f'{"min_identity":>{mar}}: {bold(nuc_min_identity)}')
            nuc_min_coverage = adjust_min_coverage(args.nuc_min_coverage)
            log.log(f'{"min_coverage":>{mar}}: {bold(nuc_min_coverage)}')
            log.log(f'{"translation table":>{mar}}: {bold(args.nuc_transtable)}')
        log.log("")

        log.log(bold(f'{"Plastidial proteins":>{mar}}:'))
        log.log(f'{"reference":>{mar}}: {protein_refs["PTD"]["AA_msg"]}')
        if protein_refs["PTD"]["AA_path"]:
            ptd_query = fasta_to_dict(protein_refs["PTD"]["AA_path"], ordered=True)
            ptd_query_info = reference_info(ptd_query)
            log.log(f'{"reference info":>{mar}}: {ptd_query_info["info_msg"]}')
            log.log(f'{"min_score":>{mar}}: {bold(args.ptd_min_score)}')
            ptd_min_identity = adjust_min_identity(args.ptd_min_identity, args.ptd_transtable)
            log.log(f'{"min_identity":>{mar}}: {bold(ptd_min_identity)}')
            ptd_min_coverage = adjust_min_coverage(args.ptd_min_coverage)
            log.log(f'{"min_coverage":>{mar}}: {bold(ptd_min_coverage)}')
            log.log(f'{"translation table":>{mar}}: {bold(args.ptd_transtable)}')
        log.log("")

        log.log(bold(f'{"Mitochondrial proteins":>{mar}}:'))
        log.log(f'{"reference":>{mar}}: {protein_refs["MIT"]["AA_msg"]}')
        if protein_refs["MIT"]["AA_path"]:
            mit_query = fasta_to_dict(protein_refs["MIT"]["AA_path"], ordered=True)
            mit_query_info = reference_info(mit_query)
            log.log(f'{"reference info":>{mar}}: {mit_query_info["info_msg"]}')
            log.log(f'{"min_score":>{mar}}: {bold(args.mit_min_score)}')
            mit_min_identity = adjust_min_identity(args.mit_min_identity, args.mit_transtable)
            log.log(f'{"min_identity":>{mar}}: {bold(mit_min_identity)}')
            mit_min_coverage = adjust_min_coverage(args.mit_min_coverage)
            log.log(f'{"min_coverage":>{mar}}: {bold(mit_min_coverage)}')
            log.log(f'{"translation table":>{mar}}: {bold(args.mit_transtable)}')
        log.log("")

        dna_ref = prepare_dna_refs(args.dna_refs)
        log.log(bold(f'{"Miscellaneous DNA":>{mar}}:'))
        log.log(f'{"reference":>{mar}}: {dna_ref["DNA"]["NT_msg"]}')
        if dna_ref["DNA"]["NT_path"]:
            dna_query = fasta_to_dict(dna_ref["DNA"]["NT_path"], ordered=True)
            dna_query_info = reference_info(dna_query)
            log.log(f'{"reference info":>{mar}}: {dna_query_info["info_msg"]}')
            log.log(f'{"min_identity":>{mar}}: {bold(args.dna_min_identity)}')
            log.log(f'{"min_coverage":>{mar}}: {bold(args.dna_min_coverage)}')
        log.log("")
        log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
        log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
        log.log(f'{"Samples to process":>{mar}}: {bold(len(fastas_to_extract))}')

        extract_coding = bool(any([protein_refs["NUC"]["AA_path"],
                                   protein_refs["PTD"]["AA_path"],
                                   protein_refs["MIT"]["AA_path"]]))

        if extract_coding or bool(dna_ref["DNA"]["NT_path"]):
            scipio_params = []
            blat_params = []
            cleanup_params = []
            for sample in fastas_to_extract:
                if protein_refs["NUC"]["AA_path"]:
                    scipio_params.append((
                        args.scipio_path,
                        args.nuc_min_score,
                        nuc_min_identity,
                        args.nuc_min_coverage,
                        args.blat_path,
                        args.overwrite,
                        args.keep_all,
                        fastas_to_extract[sample]["assembly_path"],
                        fastas_to_extract[sample]["sample_dir"],
                        sample,
                        protein_refs["NUC"]["AA_path"],
                        nuc_query,
                        nuc_query_info,
                        "NUC",
                        args.nuc_transtable,
                        args.max_loci_files,
                        args.max_loci_scipio2x
                    ))
                if protein_refs["PTD"]["AA_path"]:
                    scipio_params.append((
                        args.scipio_path,
                        args.ptd_min_score,
                        ptd_min_identity,
                        args.ptd_min_coverage,
                        args.blat_path,
                        args.overwrite,
                        args.keep_all,
                        fastas_to_extract[sample]["assembly_path"],
                        fastas_to_extract[sample]["sample_dir"],
                        sample,
                        protein_refs["PTD"]["AA_path"],
                        ptd_query,
                        ptd_query_info,
                        "PTD",
                        args.ptd_transtable,
                        args.max_loci_files,
                        args.max_loci_scipio2x
                    ))
                if protein_refs["MIT"]["AA_path"]:
                    scipio_params.append((
                        args.scipio_path,
                        args.mit_min_score,
                        mit_min_identity,
                        args.mit_min_coverage,
                        args.blat_path,
                        args.overwrite,
                        args.keep_all,
                        fastas_to_extract[sample]["assembly_path"],
                        fastas_to_extract[sample]["sample_dir"],
                        sample,
                        protein_refs["MIT"]["AA_path"],
                        mit_query,
                        mit_query_info,
                        "MIT",
                        args.mit_transtable,
                        args.max_loci_files,
                        args.max_loci_scipio2x
                    ))
                if dna_ref["DNA"]["NT_path"]:
                    blat_params.append((
                        args.blat_path,
                        args.dna_min_identity,
                        args.dna_min_coverage,
                        args.overwrite,
                        args.keep_all,
                        fastas_to_extract[sample]["assembly_path"],
                        fastas_to_extract[sample]["sample_dir"],
                        sample,
                        dna_ref["DNA"]["NT_path"],
                        dna_query,
                        dna_query_info,
                        "DNA",
                        args.max_loci_files
                    ))
                cleanup_params.append((
                    sample,
                    fastas_to_extract[sample]["sample_dir"],
                    fastas_to_extract[sample]["assembly_path"],
                    args.keep_all,
                    args.overwrite,
                    skip_clustering
                ))
            log.log(f'{"Extractions to process":>{mar}}:'
                    f' {bold(f"{len(scipio_params)} protein")} and'
                    f' {bold(f"{len(blat_params)} nucleotide")}')
            log.log("")
            if protein_refs["NUC"]["AA_path"]:
                log.log(f'{"Nuclear proteins":>{mar}}:'
                        f' {bold(f"{out_dir}/[Sample_name]__captus-ext/")}'
                        f'{bold(settings.MARKER_DIRS["NUC"])}')
            if protein_refs["PTD"]["AA_path"]:
                log.log(f'{"Plastidial proteins":>{mar}}:'
                        f' {bold(f"{out_dir}/[Sample_name]__captus-ext/")}'
                        f'{bold(settings.MARKER_DIRS["PTD"])}')
            if protein_refs["MIT"]["AA_path"]:
                log.log(f'{"Mitochondrial proteins":>{mar}}:'
                        f' {bold(f"{out_dir}/[Sample_name]__captus-ext/")}'
                        f'{bold(settings.MARKER_DIRS["MIT"])}')
            if dna_ref["DNA"]["NT_path"]:
                log.log(f'{"Miscellaneous DNA markers":>{mar}}:'
                        f' {bold(f"{out_dir}/[Sample_name]__captus-ext/")}'
                        f'{bold(settings.MARKER_DIRS["DNA"])}')
            log.log(f'{"Annotated assemblies":>{mar}}: '
                    f'{bold(f"{out_dir}/[Sample_name]__captus-ext/06_assembly_annotated")}')
            log.log(f'{"":>{mar}}  '
                    f'{dim("(All these output directories will be created for each sample)")}')
            log.log("")

            json_path = update_refs_json({**protein_refs, **dna_ref}, out_dir) # concat two dict
            log.log(f'{"Paths to references used":>{mar}}: {bold(f"{json_path}")}')
            log.log("")

            if scipio_params:
                d_msg = "Extracting protein-coding markers with Scipio"
                f_msg = "Protein-coding markers: finished processing"
                if args.debug:
                    tqdm_serial_run(scipio_coding, scipio_params, d_msg, f_msg,
                                    "extraction", args.show_less)
                else:
                    tqdm_parallel_async_run(scipio_coding, scipio_params, d_msg, f_msg,
                                            "extraction", concurrent, args.show_less)
                log.log("")

            if blat_params:
                d_msg = "Extracting miscellaneous DNA markers with BLAT"
                f_msg = "Miscellaneous DNA markers: finished processing"
                if args.debug:
                    tqdm_serial_run(blat_misc_dna, blat_params, d_msg, f_msg,
                                    "extraction", args.show_less)
                else:
                    tqdm_parallel_async_run(blat_misc_dna, blat_params, d_msg, f_msg,
                                            "extraction", concurrent, args.show_less)
                log.log("")

            d_msg = "Merging GFFs, summarizing recovery stats, and cleaning up"
            f_msg = "Merged GFF and summary recovery stats created"
            if args.debug:
                tqdm_serial_run(cleanup_post_extraction, cleanup_params, d_msg, f_msg,
                                "sample", args.show_less)
            else:
                tqdm_parallel_async_run(cleanup_post_extraction, cleanup_params, d_msg, f_msg,
                                        "sample", concurrent, args.show_less)
            log.log("")

        # Nothing to extract
        else:
            log.log(red(
                "Skipping extraction step... (no valid references found, verify the paths provided"
                " through '--nuc_refs', '--ptd_refs', '--mit_refs', '--dna_refs')"
            ))
            log.log("")


    ################################################################################################
    ####################################################################### MEGAHIT ASSEMBLY SECTION
    log.log_section_header("Clustering and Extracting Markers Across Samples with MMseqs2/BLAT")
    log.log_explanation(
        "Now Captus will cluster the assembly sequences across samples to discover homologous"
        " markers. If protein and/or miscellaneous DNA references sequences were extracted in the"
        " previous step, the '06_assembly_annotated/leftover_contigs.fasta.gz' files from each"
        " sample will be used for clustering, otherwise the '01_assembly/assembly.fasta' files are"
        " used. USE WITH CAUTION: Clustering time and/or success will depend on the size of your"
        " assembly files, it should be safe to use if you assembled captured data, RNAseq reads, or"
        " even genome skimming data.", extra_empty_lines_after=0
    )
    log.log_explanation(
        "After clustering has been completed, only the clusters with at least '--cl_min_samples'"
        " will be considered and a single representative sequence from each cluster will be selected"
        " to create a new reference file to perform a final DNA marker extraction across all samples."
        " The extraction of these new markers is done with the same settings chosen for"
        " miscellaneous DNA markers."
    )
    if skip_clustering:
        log.log(red(
            "Skipping clustering step... (to enable use '--cluster_leftovers', and verify that all"
            " the software dependencies are correctly installed)"
        ))
        log.log("")
    else:
        if not fastas_to_extract:
            quit_with_error(
                "Captus could not find FASTA assemblies to process, please verify your "
                "'--captus_assemblies_dir' and/or '--fastas' argument"
            )

        clustering_dir, clustering_dir_msg = make_output_dir(Path(out_dir, "01_clustering_data"))
        clustering_input_file = Path(clustering_dir, "clustering_input.fasta")

        # When 'cl_min_identity' is set to 'auto' it becomes 90% of 'dna_min_identity' never
        # becoming less than 75%
        if args.cl_min_identity == 'auto':
            cl_min_identity = max(settings.MMSEQS2_BLAT_DNA_IDENTITY_FACTOR * args.dna_min_identity,
                                  settings.MMSEQS_MIN_AUTO_MIN_IDENTITY)
        else:
            cl_min_identity = float(args.cl_min_identity)
        clust_tmp_dir = make_tmp_dir_within(args.cl_tmp_dir, "captus_mmseqs2_tmp")
        fastas_to_cluster, num_leftovers = find_fasta_leftovers(fastas_to_extract)
        if args.cl_min_samples == "auto":
            cl_min_samples = int(max(4, len(fastas_to_cluster) * settings.CLR_MIN_SAMPLE_PROP))
        else:
            cl_min_samples = min(int(args.cl_min_samples), len(fastas_to_cluster))

        log.log(bold_yellow("  \u25BA STEP 1 OF 3: Clustering contigs across samples with MMseqs2"))
        log.log("")
        log.log(f'{"min_seq_id":>{mar}}: {bold(cl_min_identity)}')
        log.log(f'{"cov":>{mar}}: {bold(args.cl_min_coverage)}')
        log.log(f'{"gap_open":>{mar}}: {bold(args.cl_gap_open)}')
        log.log(f'{"gap_extend":>{mar}}: {bold(args.cl_gap_extend)}')
        log.log(f'{"max_seq_len":>{mar}}: {bold(args.cl_max_seq_len)}')
        log.log(f'{"tmp_dir":>{mar}}: {bold(clust_tmp_dir)}')
        log.log("")
        log.log(f'{"Min. samples per cluster":>{mar}}: {bold(cl_min_samples)}')
        log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
        log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
        log.log("")
        log.log(f'{"Using leftover contigs":>{mar}}: {bold(num_leftovers)}')
        log.log(f'{"Using entire assembly":>{mar}}: {bold(len(fastas_to_cluster) - num_leftovers)}')
        log.log(bold(f'{"Total samples to cluster":>{mar}}: {len(fastas_to_cluster)}'))
        log.log("")
        log.log(f'{"Clustering directory":>{mar}}: {bold(clustering_dir)}')
        log.log(f'{"":>{mar}}  {dim(clustering_dir_msg)}')
        log.log("")

        if is_dir_empty(clustering_dir) or clustering_input_file.is_file() or args.overwrite:
            if clustering_input_file.is_file() and args.overwrite is False:
                log.log(
                    f"{bold('WARNING:')} The input FASTA file for clustering was found in"
                    f" '{bold(clustering_input_file)}' and it will be used, to recreate it enable"
                    " '--overwrite'"
                )
                log.log("")
            else:
                rehead_and_concatenate_fastas(fastas_to_cluster, clustering_dir,
                                              clustering_input_file,
                                              min(settings.MAX_WRITING_INSTANCES, threads_max),
                                              args.show_less)
            log.log("")
            log.log(bold("Clustering contigs across samples (this may take a while, please wait):"))
            clust_prefix, num_clusters = mmseqs2_cluster(args.mmseqs2_path, clustering_input_file,
                                                         clustering_dir, cl_min_identity,
                                                         args.cl_min_coverage, args.cl_gap_open,
                                                         args.cl_gap_extend, args.cl_max_seq_len,
                                                         clust_tmp_dir, threads_max, args.keep_all)
            log.log("")
            log.log(bold("Filtering clusters:"))
            captus_clusters_ref = filter_clusters(clust_prefix, num_clusters, cl_min_samples)
            log.log("")
            log.log("")
            log.log(bold_yellow(
                "  \u25BA STEP 2 OF 3: Extracting cluster-derived DNA markers with BLAT"
            ))
            log.log("")
            clust_ref = prepare_dna_refs(captus_clusters_ref, cluster=True)
            log.log(f'{"reference":>{mar}}: {clust_ref["CLR"]["NT_msg"]}')
            if clust_ref["CLR"]["NT_path"]:
                clust_query = fasta_to_dict(clust_ref["CLR"]["NT_path"], ordered=True)
                clust_query_info = reference_info(clust_query)
                log.log(f'{"reference info":>{mar}}: {clust_query_info["info_msg"]}')
                log.log(f'{"dna_min_identity":>{mar}}: {bold(args.dna_min_identity)}')
                log.log(f'{"dna_min_coverage":>{mar}}: {bold(args.dna_min_coverage)}')
            log.log("")
            log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
            log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
            log.log(f'{"Samples to process":>{mar}}: {bold(len(fastas_to_extract))}')
            if clust_ref["CLR"]["NT_path"]:
                blat_clusters_params = []
                for sample in fastas_to_extract:
                    blat_clusters_params.append((
                        args.blat_path,
                        args.dna_min_identity,
                        args.dna_min_coverage,
                        args.overwrite,
                        args.keep_all,
                        fastas_to_extract[sample]["assembly_path"],
                        fastas_to_extract[sample]["sample_dir"],
                        sample,
                        clust_ref["CLR"]["NT_path"],
                        clust_query,
                        clust_query_info,
                        "CLR",
                        args.max_loci_files
                    ))
                log.log(f'{"Extractions to process":>{mar}}: {bold(len(blat_clusters_params))}')
                log.log("")
                log.log(f'{"Clustering output":>{mar}}:'
                        f' {bold(f"{out_dir}/[Sample_name]__captus-ext/")}'
                        f'{bold(settings.MARKER_DIRS["CLR"])}')
                log.log(f'{"":>{mar}}  {dim("A directory will be created for each sample")}')
                log.log("")

                json_path = update_refs_json(clust_ref, out_dir)
                log.log(f'{"Paths to references used":>{mar}}: {bold(f"{json_path}")}')
                log.log("")

                d_msg = "Extracting cluster-derived DNA markers with BLAT"
                f_msg = "Cluster-derived DNA markers: finished processing"
                if args.debug:
                    tqdm_serial_run(blat_misc_dna, blat_clusters_params, d_msg, f_msg,
                                    "extraction", args.show_less)
                else:
                    tqdm_parallel_async_run(blat_misc_dna, blat_clusters_params, d_msg, f_msg,
                                            "extraction", concurrent, args.show_less)
                log.log("")
                log.log("")
                log.log(bold_yellow(
                    "  \u25BA STEP 3 OF 3: Final Merging of GFF Annotations and File Cleanup"
                ))
                log.log("")
                cleanup_params = []
                for sample in fastas_to_extract:
                    cleanup_params.append((
                        sample,
                        fastas_to_extract[sample]["sample_dir"],
                        fastas_to_extract[sample]["assembly_path"],
                        args.keep_all,
                        True, # overwrite
                        skip_clustering,
                        True, # cluster
                    ))
                d_msg = "Merging GFFs, summarizing recovery stats, and cleaning up"
                f_msg = "Merged GFF and summary recovery stats created"

                if args.debug:
                    tqdm_serial_run(cleanup_post_extraction, cleanup_params, d_msg, f_msg,
                                    "sample", args.show_less)
                else:
                    tqdm_parallel_async_run(cleanup_post_extraction, cleanup_params, d_msg, f_msg,
                                            "sample", concurrent, args.show_less)
                log.log("")
                cleanup_clustering_dir(clustering_dir, clust_prefix, clust_tmp_dir,
                                       threads_max, args.keep_all)
                log.log("")

        # Nothing to cluster, just skip already processed
        else:
            log.log(red(
                "Skipping clustering step... Captus found output files in: '{}', to replace them"
                " enable --overwrite".format(clustering_dir)
            ))
            log.log("")


    ################################################################################################
    ############################################################################## SUMMARIZE SECTION
    log.log_section_header("Statistics Summarization")
    log.log_explanation(
        "Now Captus will collect the extraction statistics from each sample to compile a"
        " comprehensive table and a HTML report for visualization of extraction statistics."
    )
    ext_stats_tsv = collect_ext_stats(out_dir)
    if ext_stats_tsv:
        log.log(f'{"Extraction statistics":>{mar}}: {bold(ext_stats_tsv)}')
        log.log("")
        if all([numpy_found, pandas_found, plotly_found]):

            from .report import build_extraction_report

            log.log_explanation(
                "Generating Marker Extraction report..."
            )
            ext_html_report, ext_html_msg = build_extraction_report(out_dir, ext_stats_tsv)
            log.log(f'{"Extraction report":>{mar}}: {bold(ext_html_report)}')
            log.log(f'{"":>{mar}}  {dim(ext_html_msg)}')
        else:
            log.log(
                f"{bold('WARNING:')} Captus uses 'numpy', 'pandas', and 'plotly' to generate  an HTML"
                " report based on the marker recovery statistics. At least one of these libraries could"
                " not be found, please verify these libraries are installed and available."
            )
    else:
        log.log(red("Skipping summarization step... (no extraction statistics files were found)"))
    log.log("")


    ################################################################################################
    ################################################################################# ENDING SECTION
    log.log_section_header(
        "Captus-assembly: Extract -> successfully completed"
        f" [{elapsed_time(time.time() - captus_start)}]"
    )
    log.log("")


def adjust_blat_concurrency(concurrent, threads_max, ram_B):
    """
    Adjust the proposed number of 'concurrent' BLAT processes so 'RAM_per_assembly' is never smaller
    than 'settings.BLAT_MIN_RAM_B'. Once the right 'concurrent' has been found 'ram_per_assembly' is
    ireadjusted
    """
    if concurrent == "auto":
        concurrent = threads_max
    else:
        try:
            concurrent = int(concurrent)
        except ValueError:
            quit_with_error("Invalid value for '--concurrent', set it to 'auto' or use a number")

    if concurrent > 1:
        while ram_B // concurrent < settings.BLAT_MIN_RAM_B:
            concurrent -= 1
    ram_B_per_extraction = ram_B // concurrent
    return concurrent, ram_B_per_extraction


def status_captus_assemblies_dir(captus_assemblies_dir, extra_fastas, margin):
    """
    Verify that the internal directory structure of `captus_assemblies_dir` matches the structure
    produced by the 'assemble' step of 'captus_assembly'
    """
    indent = (margin + 2) * " "
    num_assemblies = 0
    sample_dirs_lacking_assembly = []
    sample_dirs = list(Path(captus_assemblies_dir).resolve().rglob("*__captus-asm"))
    for sample_dir in sample_dirs:
        if not list(Path(sample_dir).rglob("*assembly.fasta")):
            sample_dirs_lacking_assembly.append(Path(sample_dir).name)
        else:
            num_assemblies += 1
    if sample_dirs:
        if len(sample_dirs) == num_assemblies and not sample_dirs_lacking_assembly:
            message = bold_green(
                f"{indent}VALID Captus assemblies directory with {len(sample_dirs)} samples"
                f" and {num_assemblies} 'assembly.fasta' files"
            )
            status = "valid"
            return status, message
        elif num_assemblies < len(sample_dirs):
            message = (
                f"{indent}WARNING: Captus assembly directory with {len(sample_dirs)} samples"
                f" and only {num_assemblies} 'assembly.fasta' files\n"
            )
            message += (
                f"{indent}{len(sample_dirs_lacking_assembly)} sample(s) lacking FASTA assemblies:\n"
                f"{indent}- "
            )
            message += f"\n{indent}- ".join(
                [name.replace("__captus-asm", "") for name in sample_dirs_lacking_assembly]
            )
            status = "incomplete"
            return status, message
    else:
        if extra_fastas:
            message = dim(f"{indent}Captus assembly directory sucessfully created")
            status = "ready"
            return status, message
        else:
            quit_with_error(
                "No assemblies to process, please provide a valid Captus assemblies directory with"
                " '--captus_assemblies_dir' and/or provide other pre-assembled FASTA assemblies with"
                " '--fastas'"
            )


def find_fasta_assemblies(captus_assemblies_dir, out_dir):
    captus_assemblies_dir = Path(captus_assemblies_dir)
    fastas_to_extract = {}
    if captus_assemblies_dir.exists():
        sample_dirs = list(Path(captus_assemblies_dir).resolve().rglob("*__captus-asm"))
        for sample_dir in sample_dirs:
            fastas = list(Path(sample_dir.resolve()).rglob("*/assembly.fasta"))
            for fasta in fastas:
                if f"{fasta.parent.parent}".endswith("__captus-asm"):
                    sample_name = fasta.parent.parent.parts[-1].replace("__captus-asm", "")
                    sample_dir = Path(out_dir, f"{sample_name}__captus-ext")
                    fastas_to_extract[sample_name] = {
                        "assembly_path": fasta,
                        "sample_dir": sample_dir
                    }
    return fastas_to_extract


def find_and_copy_fastas(fastas, captus_assemblies_dir, overwrite, threads_max, show_less):
    """
    Receives a list of files or a folder name. Only FASTA nucleotide files are accepted. They can
    have the following extensions: .fa, .fna, .fasta, .fa.gz, .fna.gz, .fasta.gz
    Sample names are derived from the filenames by removing the extensions.
    Search within a provided directory is recursive
    """
    valid_exts = settings.FASTA_VALID_EXTENSIONS
    if type(fastas) is not list:
        fastas = [fastas]
    if len(fastas) == 1 and Path(fastas[0]).is_dir():
        fastas = [file for file in Path(fastas[0]).resolve().rglob("*")
                  if has_valid_ext(file, valid_exts)]
    else:
        fastas = [Path(any_file).resolve() for any_file in fastas
                  if Path(any_file).is_file() and has_valid_ext(any_file, valid_exts)]
    check_and_copy_found_fasta_params = []
    for fasta in fastas:
        check_and_copy_found_fasta_params.append((
                fasta,
                valid_exts,
                captus_assemblies_dir,
                overwrite
        ))
    tqdm_parallel_async_run(check_and_copy_found_fasta ,check_and_copy_found_fasta_params,
                            "Verifying and importing provided FASTA files",
                            "Verification and copy completed", "assembly", threads_max, show_less)


def check_and_copy_found_fasta(fasta_path, valid_exts, captus_assemblies_dir, overwrite):
    start = time.time()
    fasta_name = fasta_path.name
    for ext in valid_exts:
        if fasta_name.lower().endswith(ext.lower()):
            ext_idx = fasta_name.lower().find(ext.lower())
            break
    sample_name = fasta_name[:ext_idx]
    sample_assembly_dir = Path(captus_assemblies_dir, f"{sample_name}__captus-asm", "01_assembly")
    sample_assembly_file = Path(sample_assembly_dir, "assembly.fasta")
    if overwrite is True or not sample_assembly_file.exists():
        if sample_assembly_dir.exists():
            shutil.rmtree(sample_assembly_dir, ignore_errors=True)
        sample_assembly_dir.mkdir(parents=True)
        if fasta_type(fasta_path) == "NT":
            fasta_out, _ = fasta_headers_to_spades(fasta_to_dict(fasta_path, ordered=True),
                                                   ordered=True)
            dict_to_fasta(fasta_out, Path(sample_assembly_dir, "assembly.fasta"), wrap=80)

            # Write a log indicating original file location of the 'assembly.fasta'
            with open(Path(sample_assembly_dir, "assembly.log"), "wt") as copy_log:
                copy_log.write(
                    "The following 'assembly.fasta' file was not assembled with Captus:\n"
                    f"  {sample_assembly_file}\nIt was successfully imported from the original file:\n"
                    f"  {fasta_path}\n"
                )
            message = (
                f"'{fasta_name}': copied to '[captus_assemblies_dir]/"
                f'{Path(f"{sample_name}__captus-asm","01_assembly","assembly.fasta")}'
                f"' [{elapsed_time(time.time() - start)}]"
            )
        else:
            message = dim(f"'{fasta_name}': skipped, this FASTA contains aminoacids")
    else:
        message = (
            f"'{fasta_name}': skipped, '[captus_assemblies_dir]/"
            f'{Path(f"{sample_name}__captus-asm","01_assembly","assembly.fasta")}'
            "' already exists"
        )
    return message


def prepare_protein_refs(
    nuc_refset, ptd_refset, mit_refset, nuc_transtable, ptd_transtable, mit_transtable
):
    """
    Looks for bundled sets of reference proteins or verify given paths, translates FASTA file
    if it only contains nucleotides
    """
    def check_refset(refset, transtable, marker: str):
        """
        Verify that `refset` is either a bundled reference set or a valid FASTA file. If the FASTA
        file contains only nucleotides then this function assumes it is only coding sequence and
        translates it using `transtable`

        Parameters
        ----------
        refset : str
            Name of bundled reference set or path to reference FASTA file
        marker : str
            Genomic compartment of the protein set, allowed values are "NUC", "PTD", and "MIT"
        """
        aa_path, nt_path, aa_msg, nt_msg = None, None, dim("not used"), dim("not used")
        if refset is None:
            return {"AA_path": aa_path, "AA_msg": aa_msg, "NT_path": nt_path, "NT_msg": nt_msg}
        elif f"{refset}".lower() in settings.PROT_REFS[marker]:
            aa_path = settings.PROT_REFS[marker][f"{refset}".lower()]["AA"]
            nt_path = settings.PROT_REFS[marker][f"{refset}".lower()]["NT"]
            aa_msg = f"{bold(refset)} {dim(aa_path)}"
            nt_msg = f"{bold(refset)} {dim(nt_path)}"
        elif Path(refset).is_file() and fasta_type(refset) == "NT":
            start = time.time()
            suffix = settings.TRANSLATED_REF_SUFFIX
            if refset.endswith(".gz"):
                aa_path = Path(Path(refset).resolve().parent,
                               f'{Path(refset.replace(".gz", "")).stem}{suffix}')
            else:
                aa_path = Path(Path(refset).resolve().parent, f"{Path(refset).stem}{suffix}")
            nt_path = Path(refset).resolve()
            amino_refset = translate_fasta_dict(fasta_to_dict(refset, ordered=True), transtable)
            amino_refset_fixed = fix_premature_stops(amino_refset)
            if amino_refset_fixed is None:
                dict_to_fasta(amino_refset, aa_path, wrap=80)
                log.log(
                    f"Translated '{bold(refset)}' using Genetic Code: {bold(transtable)}"
                    f" [{elapsed_time(time.time() - start)}]"
                )
            else:
                dict_to_fasta(amino_refset_fixed, aa_path, wrap=80)
                log.log(
                    f"Translated '{bold(refset)}' using Genetic Code: {bold(transtable)}"
                    f" [{elapsed_time(time.time() - start)}]"
                )
                log.log(
                    f"WARNING: Premature stop codons were found in {refset}"
                    " and automatically converted to X"
                )
            aa_msg = f'{bold(aa_path)} {dim("(Translated from")} {dim(refset)}{dim(")")}'
            nt_msg = f"{nt_path}"
        elif Path(refset).is_file() and fasta_type(refset) == "AA":
            amino_refset = fasta_to_dict(Path(refset).resolve())
            amino_refset_fixed = fix_premature_stops(amino_refset)
            if amino_refset_fixed is None:
                aa_path = Path(Path(refset).resolve().parent, f"{Path(refset).stem}.faa")
                dict_to_fasta(amino_refset, aa_path, wrap=80)
            else:
                aa_path = Path(Path(refset).resolve().parent, f"{Path(refset).stem}_fixed.faa")
                dict_to_fasta(amino_refset_fixed, aa_path, wrap=80)
                log.log(
                    f"WARNING: Premature stop codons were found in {refset}"
                    " and automatically converted to X"
                )
            aa_msg = bold(aa_path)
        elif Path(refset).is_file() and fasta_type(refset) == "invalid":
            aa_msg = red("not a valid FASTA")
        else:
            aa_msg = red("file not found")
        return {"AA_path": aa_path, "AA_msg": aa_msg, "NT_path": nt_path, "NT_msg": nt_msg}

    log.log("Verifying protein reference sets...")
    protein_refs = {
        "NUC": check_refset(nuc_refset, nuc_transtable, "NUC"),
        "PTD": check_refset(ptd_refset, ptd_transtable, "PTD"),
        "MIT": check_refset(mit_refset, mit_transtable, "MIT"),
    }

    return protein_refs


def prepare_dna_refs(dna_refs, cluster=False):
    nt_path, nt_msg = None, dim("not used")
    if dna_refs is None:
        return {"DNA": {"AA_path": None, "AA_msg": None, "NT_path": nt_path, "NT_msg": nt_msg}}
    elif f"{dna_refs}".lower() in settings.DNA_REFS:
        nt_path = settings.DNA_REFS[dna_refs.lower()].resolve()
        nt_msg = f'{bold(dna_refs)} {dim(nt_path)}'
    elif Path(dna_refs).is_file() and fasta_type(dna_refs) == "NT":
        nt_path = Path(dna_refs).resolve()
        nt_msg = bold(nt_path)
    else:
        nt_msg = red("not a valid FASTA")

    if cluster:
        return {"CLR": {"AA_path": None, "AA_msg": None, "NT_path": nt_path, "NT_msg": nt_msg}}
    else:
        return {"DNA": {"AA_path": None, "AA_msg": None, "NT_path": nt_path, "NT_msg": nt_msg}}


def update_refs_json(refs_paths: dict, out_dir):
    """
    Attempts to load a JSON file with the paths to the references last used and update it with the
    ones used currently. Creates a new JSON if the file is not found.

    Parameters
    ----------
    refs_paths : dict
        A dictionary with the reference paths for each marker type except CLR
    """
    json_path = Path(out_dir, settings.JSON_REFS)
    try:
        with open(json_path, "rt") as jin:
            refs_json = json.load(jin)
    except FileNotFoundError:
        refs_json = {
            "NUC" : {"AA_path": None, "AA_msg": "not used", "NT_path": None, "NT_msg": "not used"},
            "PTD" : {"AA_path": None, "AA_msg": "not used", "NT_path": None, "NT_msg": "not used"},
            "MIT" : {"AA_path": None, "AA_msg": "not used", "NT_path": None, "NT_msg": "not used"},
            "DNA" : {"AA_path": None, "AA_msg": None,       "NT_path": None, "NT_msg": "not used"},
            "CLR" : {"AA_path": None, "AA_msg": None,       "NT_path": None, "NT_msg": "not used"},
        }

    if refs_json:
        for marker in refs_paths:
            for info in refs_paths[marker]:
                if refs_paths[marker][info]:
                    if info.endswith("_path"):
                        refs_json[marker][info] = f"{refs_paths[marker][info]}"
                    else:
                        msg_clean = remove_formatting(f"{refs_paths[marker][info]}")
                        if msg_clean != "not used":
                            refs_json[marker][info] = msg_clean

    with open(json_path, "wt") as jout:
        json.dump(refs_json, jout, indent=4)

    return json_path


def adjust_min_identity(min_identity, transtable):
    """
    Reduce 'min_identity' if 'transtable' represents a divergent genetic code that potentially will
    find matches with many mismatches because BLAT only uses the Standard Code (Table 1).
    Detect if by mistake the number was provided as a decimal between 0 and 1 and fix it.
    """
    if 1 >= min_identity > 0:
        min_identity = round(min_identity * 100, 2)
    if transtable in settings.DIVERGENT_GENETIC_CODES:
        return min(settings.SCIPIO_MAX_IDENTITY_DIV_CODE, min_identity)
    else:
        return min_identity


def adjust_min_coverage(min_coverage):
    if 1 >= min_coverage > 0:
        return round(min_coverage * 100, 2)
    else:
        return min_coverage


def scipio_coding(
        scipio_path, min_score, min_identity, min_coverage, blat_path, overwrite, keep_all,
        target_path, sample_dir, sample_name, query_path, query_dict, query_info, marker_type,
        transtable, max_loci_files, max_loci_scipio2x
):
    """
    Perform two consecutive rounds of Scipio, the first run with mostly default settings and with a
    'min_score' multiplied by 'settings.SCIPIO_SCORE_FACTOR'. This round narrows down the set of
    target contigs (your assembly) to only those with hits and the set of reference proteins to
    those that produced the best hits (just when you provided multiple reference proteins of the
    same cluster, for more details check the section 'REFERENCE_CLUSTER_SEPARATOR' in the file
    'settings_assembly.py')
    """

    start = time.time()

    genes = {"NUC": "nuclear genes", "PTD": "plastidial genes", "MIT": "mitochondrial genes"}

    initial_contigs = fasta_to_dict(target_path, ordered=True)

    # Set programs' paths in case of using the bundled versions
    if scipio_path == "bundled":
        scipio_path = settings.BUNDLED_SCIPIO
    if blat_path == "bundled":
        blat_path = settings.BUNDLED_BLAT

    # Group Scipio's basic parameters, the 'sample_dir', and the 'marker_type' into a list
    scipio_params = [
        scipio_path,    # 0
        min_score,      # 1
        min_identity,   # 2
        min_coverage,   # 3
        blat_path,      # 4
        transtable,     # 5
        sample_dir,     # 6
        marker_type,    # 7
        keep_all,       # 8
    ]

    # Run Scipio twice if 'query_info["num_loci"]' does not exceed 'max_loci_scipio2x'
    if query_info["num_loci"] <= max_loci_scipio2x:

        # Use the function run_scipio_command() that to run Scipio and also get the name of the YAML
        # output 'yaml_out_file'
        yaml_initial_file = run_scipio_command(scipio_params, target_path, query_path,
                                               overwrite, stage="initial")
        if yaml_initial_file is None:
            message = dim(
                f"'{sample_name}': extraction of {genes[marker_type]}"
                " skipped (output files already exist)"
            )
            return message
        else:
            yaml_initial_dir = yaml_initial_file.parent
            initial_models = scipio_yaml_to_dict(yaml_initial_file, min_identity,
                                                 min_coverage, marker_type)

        # Parse YAML to subselect only the best proteins from the 'query' and the contigs with hits
        # from the assembly ('target')
        if initial_models is None:
            message = red(f"'{sample_name}': FAILED extraction of {genes[marker_type]}")
            return message
        else:
            final_target, final_query = filter_query_and_target(query_dict, initial_contigs,
                                                                yaml_initial_dir, initial_models,
                                                                marker_type)

        # Perform final Scipio's run (more exhaustive but with fewer contigs and reference proteins)
        yaml_final_file = run_scipio_command(scipio_params, final_target, final_query,
                                             overwrite, stage="final")

    # Run a single Scipio run when 'num_refs' exceeds 'max_loci_scipio2x'
    else:
        yaml_final_file = run_scipio_command(scipio_params, target_path, query_path,
                                             overwrite, stage="single")

    if yaml_final_file is None:
        message = dim(
            f"'{sample_name}': extraction of {genes[marker_type]}"
            " skipped (output files already exist)"
        )
        return message
    else:
        yaml_final_dir = yaml_final_file.parent
        final_models = scipio_yaml_to_dict(yaml_final_file, min_identity, min_coverage, marker_type)

    # Parse final YAML to produce output FASTA files and reports
    if final_models is None:
        message = red(f"'{sample_name}': FAILED extraction of {genes[marker_type]}")
        return message
    else:
        write_gff3(final_models, marker_type, Path(yaml_final_dir, f"{marker_type}_contigs.gff"))
        recovery_stats = write_fastas_and_report(final_models, sample_name, initial_contigs,
                                                 yaml_final_dir, marker_type, overwrite,
                                                 max_loci_files)
        message = (
            f"'{sample_name}': recovered {recovery_stats['num_loci']} {genes[marker_type]}"
            f' ({recovery_stats["num_loci"] / query_info["num_loci"]:.1%} of {query_info["num_loci"]}),'
            f' {recovery_stats["total_length_best_hits"] / query_info["total_length_loci"]:.1%} of total'
            f' reference length, {recovery_stats["num_paralogs"]} paralogs'
            f" found [{elapsed_time(time.time() - start)}]"
        )
        return message


def reference_info(query_dict):
    num_seqs = len(query_dict)
    separators_found = 0

    for seq_name in query_dict:
        if settings.REFERENCE_CLUSTER_SEPARATOR in seq_name:
            if len(list(filter(None, seq_name.split(settings.REFERENCE_CLUSTER_SEPARATOR)))) > 1:
                separators_found += 1

    loci_lengths = {}
    if separators_found == num_seqs:
        for seq_name in query_dict:
            locus_name = seq_name.split(settings.REFERENCE_CLUSTER_SEPARATOR)[-1]
            if locus_name not in loci_lengths:
                loci_lengths[locus_name] = [len(query_dict[seq_name]["sequence"])]
            else:
                loci_lengths[locus_name].append(len(query_dict[seq_name]["sequence"]))
        for locus_name in loci_lengths:
            loci_lengths[locus_name] = statistics.mean(loci_lengths[locus_name])
    else:
        for seq_name in query_dict:
            loci_lengths[seq_name] = len(query_dict[seq_name]["sequence"])

    num_loci = len(loci_lengths)
    total_length_loci = round(sum(list(loci_lengths.values())))
    info_msg = bold(f'{num_loci} loci, {num_seqs} sequences ')
    if separators_found == num_seqs:
        if num_loci == num_seqs:
            info_msg += dim(f'(loci names found, detected a single sequence per locus)')
        elif num_loci < num_seqs:
            info_msg += dim(f'(loci names found, detected multiple sequences per locus)')
    else:
        info_msg += dim(
            f'(locus name separator "{settings.REFERENCE_CLUSTER_SEPARATOR}" missing in '
            f'{num_seqs - separators_found} sequences, each sequence taken as a different locus)')

    ref_info = {
        "num_seqs": num_seqs,
        "separators_found": separators_found,
        "num_loci": num_loci,
        "total_length_loci": total_length_loci,
        "info_msg": info_msg
    }

    return ref_info


def run_scipio_command(scipio_params, target, query, overwrite, stage):
    # Set output directory and files according to 'sample_dir', 'genome', and 'stage'
    marker_type = scipio_params[7]
    scipio_out_dir = Path(scipio_params[6], settings.MARKER_DIRS[marker_type])
    if stage == "initial":
        scipio_out_dir   = Path(scipio_out_dir, f"00_initial_scipio_{marker_type}")
        blat_out_file    = Path(scipio_out_dir, f"{marker_type}_scipio_initial.psl")
        scipio_out_file  = Path(scipio_out_dir, f"{marker_type}_scipio_initial.yaml")
        scipio_log_file  = Path(scipio_out_dir, f"{marker_type}_scipio_initial.log")

        # Adjust Scipio's 'min_score' for initial round allowing more lenient matching
        scipio_min_score = round(scipio_params[1] * settings.SCIPIO_SCORE_INITIAL_FACTOR, 5)
    else:
        blat_out_file   = Path(scipio_out_dir, f"{marker_type}_scipio_final.psl")
        scipio_out_file = Path(scipio_out_dir, f"{marker_type}_scipio_final.yaml")
        scipio_log_file = Path(scipio_out_dir, f"{marker_type}_scipio_final.log")
        scipio_min_score = scipio_params[1]

    scipio_out_dir, _ = make_output_dir(scipio_out_dir)
    if is_dir_empty(scipio_out_dir) is True or overwrite is True:

        # Set minimum BLAT score according to number of references (too many hits created when using
        # a large number of reference proteins), higher score, fewer hits
        if stage == "single":
            blat_score = f"{settings.SCIPIO_1X_BLAT_MIN_SCORE}"
        else:
            blat_score = f"{settings.SCIPIO_2X_BLAT_MIN_SCORE}"

        # Build the basic part of the command
        basic = [
            "perl", f"{scipio_params[0]}",
            f"--blat_output={blat_out_file}",
            "--overwrite",
            "--verbose",
            "--keep_blat_output",
            "--show_blatline",
            f"--min_score={scipio_min_score}",
            f"--min_identity={scipio_params[2]}",
            f"--min_coverage={scipio_params[3]}",
            "--max_mismatch=0",  # 0 means infinite
            "--multiple_results",
            f"--blat_bin={scipio_params[4]}",
            f"--blat_score={blat_score}",
            f"--blat_identity={scipio_params[2] * settings.SCIPIO_BLAT_IDENTITY_FACTOR:.0f}",
            f"--transtable={scipio_params[5]}",
            f"--accepted_intron_penalty={settings.SCIPIO_ACCEPTED_INTRON_PENALTY}"
        ]

        # Initial settings, genome-specific
        initial_settings = settings.SCIPIO_GENOME_BASIC_SETTINGS[marker_type]

        # Get extra parameters for final Scipio's run according to the 'genome'
        final_settings = settings.SCIPIO_GENOME_EXTRA_SETTINGS[marker_type]

        # Concatenate pieces of the command
        if stage == "final":
            scipio_cmd = basic + final_settings + [f"{target}", f"{query}"]
        else:
            scipio_cmd = basic + initial_settings + [f"{target}", f"{query}"]

        # Execute Scipio's command
        with open(scipio_out_file, "w") as yaml_out:
            with open(scipio_log_file, "w") as yaml_log:
                subprocess.run(scipio_cmd, stdout=yaml_out, stderr=yaml_log)

        # Erase BLAT .psl files if 'keep_all' is disabled
        if not scipio_params[8]:
            blat_out_file.unlink()

        return scipio_out_file

    else:

        return None


def filter_query_and_target(query_dict, target_dict, yaml_initial_dir, initial_models, marker_type):
    """
    Retain only best model proteins and contigs that were hit after initial Scipio run
    """
    best_proteins, hit_contigs = {}, {}
    for protein in initial_models:
        best_protein_name = initial_models[protein][0]["ref_name"]
        best_proteins[best_protein_name] = query_dict[best_protein_name]
        for h in range(len(initial_models[protein])):
            for contig in initial_models[protein][h]["hit_contig"].split("\n"):
                if contig not in hit_contigs:
                    hit_contigs[contig] = target_dict[contig]
    final_query = Path(yaml_initial_dir, f"{marker_type}_best_proteins.faa")
    final_target = Path(yaml_initial_dir, f"{marker_type}_hit_contigs.fna")
    dict_to_fasta(best_proteins, final_query, wrap=80)
    dict_to_fasta(hit_contigs, final_target, wrap=80)
    return final_target, final_query


def write_fastas_and_report(
        hits, sample_name, target_dict, out_dir, marker_type, overwrite, max_loci_files
):
    """
    Take a 'hits' dictionary from coding or miscellaneous DNA extraction and produce FASTA outputs
    according to type of marker and type of extraction
    """

    def format_coords(ref_coords):
        """
        Turn 'ref_coords' from '32-56,123-983\n23122-22685' to '33-56,124-983;22686-23122'
        Making starts always the smallest values and adding 1 to start (1-based, closed coords)
        """
        formatted_coords = []
        for segment in ref_coords.split("\n"):
            segment_coords = []
            for coord in segment.split(","):
                start, end = int(coord.split("-")[0]), int(coord.split("-")[1])
                if start == end:
                    segment_coords.append(f"{start}-{end}")
                else:
                    segment_coords.append(f"{min(start, end) + 1}-{max(start, end)}")
            formatted_coords.append(",".join(segment_coords))
        return ";".join(formatted_coords)

    num_loci, num_paralogs = 0, 0
    lengths_best_hits, coverages_best_hits = [], []
    flanked_seqs, gene_seqs, cds_aa_seqs, cds_nt_seqs, hit_contigs = {}, {}, {}, {}, {}
    stats_header = "\t".join(["sample_name", "marker_type",
                              "locus", "ref_name", "ref_coords", "ref_type", "ref_len_matched",
                              "hit", "pct_recovered", "pct_identity", "score", "lwscore",
                              "hit_len", "cds_len", "intron_len", "flanks_len", "frameshifts",
                              "ctg_names", "ctg_strands", "ctg_coords"])
    stats = []
    for ref in sorted(hits):
        num_loci += 1
        for h in range(len(hits[ref])):

            length = hits[ref][h]["matches"] + hits[ref][h]["mismatches"]
            if h == 0:
                lengths_best_hits.append(length)
                coverages_best_hits.append(hits[ref][h]["coverage"])
            else:
                num_paralogs += 1

            description = (
                f'hit={h:02}|lwscore={hits[ref][h]["lwscore"]:.3f}|'
                f'cover={hits[ref][h]["coverage"]:.2f}|ident={hits[ref][h]["identity"]:.2f}|'
                f'score={hits[ref][h]["score"]:.3f}'
            )

            seq_flanked = hits[ref][h]["seq_flanked"]
            len_flanked = f"|length={len(seq_flanked)}"
            seq_gene = hits[ref][h]["seq_gene"]
            len_gene = f"|length={len(seq_gene)}"
            seq_nt = hits[ref][h]["seq_nt"]
            len_nt = f"|length={len(seq_nt)}"
            seq_aa = hits[ref][h]["seq_aa"]
            len_aa = f"|length={len(seq_aa)}"

            ref_coords = format_coords(hits[ref][h]["ref_coords"])
            query = (
                f'|query={hits[ref][h]["ref_name"]}:{ref_coords}'
                f'|contigs={hits[ref][h]["hit_contig"]}'.replace("\n", ";")
            )

            shifts_flanked, shifts_gene, shifts_nt, shifts_aa = "", "", "", ""
            if marker_type in ["NUC", "PTD", "MIT"]:
                shifts_flanked = [str(p + 1) for p in range(len(seq_flanked))
                                  if seq_flanked[p] == "N"]
                if shifts_flanked:
                    shifts_flanked = f'|frameshifts={",".join(shifts_flanked)}'
                else:
                    shifts_flanked = ""
                shifts_gene = [str(p + 1) for p in range(len(seq_gene)) if seq_gene[p] == "N"]
                if shifts_gene:
                    shifts_gene = f'|frameshifts={",".join(shifts_gene)}'
                else:
                    shifts_gene = ""
                shifts_nt = [str(p + 1) for p in range(len(seq_nt)) if seq_nt[p] == "N"]
                shifts_aa = [str(math.ceil(int(p) / 3)) for p in shifts_nt]
                if shifts_nt:
                    shifts_aa = f'|frameshifts={",".join(shifts_aa)}'
                    shifts_nt = f'|frameshifts={",".join(shifts_nt)}'
                else:
                    shifts_aa, shifts_nt = "", ""

            if len(hits[ref]) == 1:
                seq_name = "|".join([sample_name, ref])
            else:
                seq_name = "|".join([sample_name, ref, f"{h:02}"])

            flanked_seqs[seq_name] = {
                "description": f"{description}{len_flanked}{shifts_flanked}{query}",
                "sequence": seq_flanked,
                "ref_name": ref
            }
            gene_seqs[seq_name] = {
                "description": f"{description}{len_gene}{shifts_gene}{query}",
                "sequence": seq_gene,
                "ref_name": ref
            }

            if marker_type in ["NUC", "PTD", "MIT"]:
                cds_nt_seqs[seq_name] = {
                    "description": f"{description}{len_nt}{shifts_nt}{query}",
                    "sequence": seq_nt,
                    "ref_name": ref
                }
                cds_aa_seqs[seq_name] = {
                    "description": f"{description}{len_aa}{shifts_aa}{query}",
                    "sequence": seq_aa,
                    "ref_name": ref
                }

            # Format data for summary table
            hit_len = len(seq_flanked.replace("n", ""))
            flanks_len = len(seq_flanked.replace("n", "")) - len(seq_gene.replace("n", ""))
            if marker_type in ["NUC", "PTD", "MIT"]:
                intron_len = max(len(seq_gene.replace("n", "")) - len(seq_nt), 0)
                stats_row = {"ref_type": "prot",
                             "cds_len": f'{len(seq_nt)}',
                             "intron_len": f'{intron_len}',
                             "frameshifts": shifts_nt.replace("|frameshifts=", "")}
            else:
                stats_row = {"ref_type": "nucl",
                             "cds_len": "NA",
                             "intron_len": "NA",
                             "frameshifts": "NA"}
            stats.append("\t".join([sample_name,
                                    marker_type,
                                    ref,
                                    hits[ref][h]["ref_name"],
                                    ref_coords,
                                    stats_row["ref_type"],
                                    f"{length}",
                                    f"{h:02}",
                                    f'{hits[ref][h]["coverage"]:.2f}',
                                    f'{hits[ref][h]["identity"]:.2f}',
                                    f'{hits[ref][h]["score"]:.3f}',
                                    f'{hits[ref][h]["lwscore"]:.3f}',
                                    f'{hit_len}',
                                    stats_row["cds_len"],
                                    stats_row["intron_len"],
                                    f'{flanks_len}',
                                    stats_row["frameshifts"],
                                    f'{hits[ref][h]["hit_contig"]}'.replace("\n", ";"),
                                    f'{hits[ref][h]["strand"]}'.replace("\n", ";"),
                                    format_coords(hits[ref][h]["hit_coords"]),
                                    ]))

            for contig in hits[ref][h]["hit_contig"].split("\n"):
                if contig not in hit_contigs:
                    hit_contigs[contig] = dict(target_dict[contig])

    # Write statistics table
    with open(Path(out_dir, f"{marker_type}_recovery_stats.tsv"), "w") as stats_out:
        stats_out.write(stats_header + "\n" + "\n".join(stats) + "\n")

    # Write multi-sequence FASTAs and setup directories for locus-wise files, only prepare a
    # separate file per marker if the number of references is not greater than args.max_loci_files
    with open(Path(out_dir, f"{marker_type}_contigs.list.txt"), "w") as hit_list_out:
        for contig in hit_contigs:
            hit_list_out.write(f"{contig}\n")

    if marker_type in ["NUC", "PTD", "MIT"]:
        dict_to_fasta(cds_aa_seqs, Path(out_dir, f'{marker_type}{settings.FORMAT_SUFFIXES["AA"]}'))
        dict_to_fasta(cds_nt_seqs, Path(out_dir, f'{marker_type}{settings.FORMAT_SUFFIXES["NT"]}'))
        dict_to_fasta(gene_seqs, Path(out_dir, f'{marker_type}{settings.FORMAT_SUFFIXES["GE"]}'))
        dict_to_fasta(flanked_seqs, Path(out_dir, f'{marker_type}{settings.FORMAT_SUFFIXES["GF"]}'))
    elif marker_type in ["DNA", "CLR"]:
        dict_to_fasta(gene_seqs, Path(out_dir, f'{marker_type}{settings.FORMAT_SUFFIXES["MA"]}'))
        dict_to_fasta(flanked_seqs, Path(out_dir, f'{marker_type}{settings.FORMAT_SUFFIXES["MF"]}'))

    if len(gene_seqs) <= max_loci_files:
        if marker_type in ["NUC", "PTD", "MIT"]:
            cds_aa_seqs_dir, _ = make_output_dir(Path(out_dir, settings.FORMAT_DIRS["AA"]))
            cds_nt_seqs_dir, _ = make_output_dir(Path(out_dir, settings.FORMAT_DIRS["NT"]))
            gene_seqs_dir, _ = make_output_dir(Path(out_dir, settings.FORMAT_DIRS["GE"]))
            flanked_seqs_dir, _ = make_output_dir(Path(out_dir, settings.FORMAT_DIRS["GF"]))
            if is_dir_empty(cds_aa_seqs_dir) is True or overwrite is True:
                for faa_file in cds_aa_seqs_dir.glob("*.faa"):
                    faa_file.unlink()
            if is_dir_empty(cds_nt_seqs_dir) is True or overwrite is True:
                for fna_file in cds_nt_seqs_dir.glob("*.fna"):
                    fna_file.unlink()
        elif marker_type in ["DNA", "CLR"]:
            gene_seqs_dir, _ = make_output_dir(Path(out_dir, settings.FORMAT_DIRS["MA"]))
            flanked_seqs_dir, _ = make_output_dir(Path(out_dir, settings.FORMAT_DIRS["MF"]))

        if is_dir_empty(gene_seqs_dir) is True or overwrite is True:
            for fna_file in gene_seqs_dir.glob("*.fna"):
                fna_file.unlink()
        if is_dir_empty(flanked_seqs_dir) is True or overwrite is True:
            for fna_file in flanked_seqs_dir.glob("*.fna"):
                fna_file.unlink()

        for seq_name_full in sorted(flanked_seqs):
            if len(seq_name_full.split("|")) == 3:
                seq_name_short = "|".join(seq_name_full.split("|")[::2])
            else:
                seq_name_short = seq_name_full.split("|")[0]
            dict_to_fasta({seq_name_short: dict(flanked_seqs[seq_name_full])},
                          Path(flanked_seqs_dir, f'{flanked_seqs[seq_name_full]["ref_name"]}.fna'),
                          append=True)
            dict_to_fasta({seq_name_short: dict(gene_seqs[seq_name_full])},
                          Path(gene_seqs_dir, f'{gene_seqs[seq_name_full]["ref_name"]}.fna'),
                          append=True)
            if marker_type in ["NUC", "PTD", "MIT"]:
                dict_to_fasta({seq_name_short: dict(cds_aa_seqs[seq_name_full])},
                              Path(cds_aa_seqs_dir, f'{cds_aa_seqs[seq_name_full]["ref_name"]}.faa'),
                              append=True)
                dict_to_fasta({seq_name_short: dict(cds_nt_seqs[seq_name_full])},
                              Path(cds_nt_seqs_dir, f'{cds_nt_seqs[seq_name_full]["ref_name"]}.fna'),
                              append=True)

    recovery_stats = {
        "num_loci": num_loci,
        "num_paralogs": num_paralogs,
        "mean_length_best_hits": round(statistics.mean(lengths_best_hits), 2),
        "total_length_best_hits": sum(lengths_best_hits),
        "mean_coverage_best_hits": round(statistics.mean(coverages_best_hits), 2)
    }

    return recovery_stats


def blat_misc_dna(
        blat_path, min_identity, min_coverage, overwrite, keep_all, target, sample_dir, sample_name,
        query_path, query_dict, query_info, marker_type, max_loci_files
):
    """
    Extract matches of miscellaneous DNA sequences by comparing the assemblies to a set of
    references, these can be formatted as the proteins references '>sample-locus_name'
    """

    start = time.time()

    dna_target = fasta_to_dict(target)

    # Set BLAT path in case of using the bundled version
    if blat_path == "bundled":
        blat_path = settings.BUNDLED_BLAT

    # Create output directory
    blat_dna_out_dir  = Path(sample_dir, settings.MARKER_DIRS[marker_type])
    blat_dna_out_file = Path(blat_dna_out_dir, f"{marker_type}_blat_search.psl")
    blat_dna_log_file = Path(blat_dna_out_dir, f"{marker_type}_blat_search.log")
    dna_gff_file = Path(blat_dna_out_dir, f"{marker_type}_contigs.gff")
    blat_dna_out_dir, _ = make_output_dir(blat_dna_out_dir)

    if is_dir_empty(blat_dna_out_dir) is True or overwrite is True:
        blat_cmd = [
            f"{blat_path}",
            "-t=dna",
            "-q=dna",
            "-noHead",
            f"-minIdentity={min_identity}",
            f"{target}",
            f"{query_path}",
            f"{blat_dna_out_file}"
        ]
        with open(blat_dna_log_file, "w") as blat_log:
            blat_log.write(f"Captus' BLAT command:\n  {' '.join(blat_cmd)}\n\n")
        with open(blat_dna_log_file, "a") as blat_log:
            subprocess.run(blat_cmd, stdout=blat_log, stderr=blat_log)

        dna_hits = blat_misc_dna_psl_to_dict(blat_dna_out_file, dna_target, min_identity,
                                             min_coverage, marker_type)
        if not dna_hits:
            message = red(f"'{sample_name}': FAILED extraction of miscellaneous DNA markers")
            return message
        else:
            if not keep_all:
                Path(blat_dna_out_file).unlink()
            write_gff3(dna_hits, marker_type, dna_gff_file)
            recovery_stats = write_fastas_and_report(dna_hits, sample_name, dna_target,
                                                     blat_dna_out_dir, marker_type, overwrite,
                                                     max_loci_files)
            message = (
                f"'{sample_name}': recovered {recovery_stats['num_loci']} DNA markers"
                f' ({recovery_stats["num_loci"] / query_info["num_loci"]:.1%} of'
                f' {query_info["num_loci"]}),'
                f' {recovery_stats["total_length_best_hits"] / query_info["total_length_loci"]:.1%}'
                f' of total reference length, {recovery_stats["num_paralogs"]} paralogs found'
                f" [{elapsed_time(time.time() - start)}]"
            )
            return message

    else:
        message = dim(
            f"'{sample_name}': extraction of miscellaneous DNA markers"
            " skipped (output files already exist)"
        )
        return message


def cleanup_post_extraction(
        sample_name, sample_dir, assembly_path, keep_all, overwrite, skip_clustering, cluster=False
):
    """
    Concatenate all '.gff' from extraction folders to make a master annotation file and a single
    FASTA with all hit contigs as well as a FASTA with the leftover contigs that can be used for
    later clustering. Also concatenate all '_recovery_stats.tsv' tables without repeating the header.
    Switch clust to True in order to create a 'leftover_contigs_after_clustering.fasta.gz'
    """

    start = time.time()

    annotated_assembly_dir, _ = make_output_dir(Path(sample_dir, "06_assembly_annotated"))

    gff_file = Path(annotated_assembly_dir, f"{sample_name}_hit_contigs.gff")
    stats_file = Path(annotated_assembly_dir, f"{sample_name}_recovery_stats.tsv")
    hit_contigs_file = Path(annotated_assembly_dir, f"{sample_name}_hit_contigs.fasta")
    leftovers_clust_file = Path(annotated_assembly_dir, "leftover_contigs_after_clustering.fasta.gz")
    leftovers_file = Path(annotated_assembly_dir, "leftover_contigs.fasta.gz")
    out_files = [gff_file, stats_file, hit_contigs_file, leftovers_clust_file, leftovers_file]

    if is_dir_empty(annotated_assembly_dir) is True or overwrite is True:
        # Concatenate all GFF from different genomes
        sample_gffs = list(sample_dir.resolve().rglob("[A-Z]*_contigs.gff"))
        sample_gffs = [gff for gff in sample_gffs
                       if gff.parts[-2] != "06_assembly_annotated"]
        if not cluster and not skip_clustering:
            sample_gffs = [gff for gff in sample_gffs
                           if gff.parts[-2] != settings.MARKER_DIRS["CLR"]]
        gff_lines = 0
        with open(gff_file, "wt") as gff_out:
            gff_out.write("##gff-version 3\n")
            if sample_gffs:
                for gff in sorted(sample_gffs):
                    with open(gff, "rt") as gff_in:
                        for line in gff_in:
                            if line != "##gff-version 3\n":
                                gff_out.write(line)
                                gff_lines += 1

        if gff_lines == 0:
            for file in out_files:
                if file.exists():
                    file.unlink()
            message = red(
                f"'{sample_name}': GFF was empty, marker extraction FAILED"
                f" [{elapsed_time(time.time() - start)}]"
            )
            return message

        # Concatenate all recovery statistics tables
        stats_header = "\t".join(["sample_name", "marker_type",
                                  "locus", "ref_name", "ref_coords", "ref_type", "ref_len_matched",
                                  "hit", "pct_recovered", "pct_identity", "score", "lwscore",
                                  "hit_len", "cds_len", "intron_len", "flanks_len", "frameshifts",
                                  "ctg_names", "ctg_strands", "ctg_coords"]) + "\n"
        sample_stats = list(sample_dir.resolve().rglob("[A-Z]*_recovery_stats.tsv"))
        sample_stats = [tsv for tsv in sample_stats
                        if tsv.parts[-2] != "06_assembly_annotated"]
        if not cluster and not skip_clustering:
            sample_stats = [tsv for tsv in sample_stats
                            if tsv.parts[-2] != settings.MARKER_DIRS["CLR"]]
        tsv_lines = 0
        with open(stats_file, "wt") as tsv_out:
            tsv_out.write(stats_header)
            if sample_stats:
                for table in sorted(sample_stats):
                    with open(table, "rt") as stats_in:
                        for line in stats_in:
                            if line != stats_header:
                                tsv_out.write(line)
                                tsv_lines += 1
        if tsv_lines == 0:
            tsv_out.unlink()

        # Write FASTAs of contigs with and without hits
        names_hit_contigs = []
        contig_lists = list(sample_dir.resolve().rglob("[A-Z]*_contigs.list.txt"))
        if not cluster and not skip_clustering:
            contig_lists = [cl for cl in contig_lists
                            if cl.parts[-2] != settings.MARKER_DIRS["CLR"]]
        if contig_lists:
            for cl in contig_lists:
                with open(cl, "rt") as clin:
                    for line in clin:
                        names_hit_contigs.append(line.strip("\n"))
        all_contigs = fasta_to_dict(assembly_path, ordered=True)
        hit_contigs = OrderedDict()
        leftover_contigs = OrderedDict()
        for contig_name in all_contigs:
            if contig_name in names_hit_contigs:
                hit_contigs[contig_name] = dict(all_contigs[contig_name])
            else:
                leftover_contigs[contig_name] = dict(all_contigs[contig_name])
        dict_to_fasta(hit_contigs, hit_contigs_file, wrap=80)
        if cluster:
            dict_to_fasta(leftover_contigs, leftovers_clust_file, wrap=80)
        else:
            dict_to_fasta(leftover_contigs, leftovers_file, wrap=80)

        # Erase unnecessary files if '--keep_all' is disabled
        if not keep_all:
            files_to_delete = []
            files_to_delete += list(sample_dir.resolve().rglob("*.psl"))
            files_to_delete += list(sample_dir.resolve().rglob("*.yaml"))
            files_to_delete += list(sample_dir.resolve().rglob("*_hit_contigs.fna"))
            for del_file in files_to_delete:
                del_file.unlink()

        message = (
            f"'{sample_name}': GFF annotations and recovery stats merged, unnecessary files removed"
            f" [{elapsed_time(time.time() - start)}]"
        )
        return message
    else:
        message = dim(
            f"'{sample_name}': cleanup and merging of GFFs and stats skipped"
            " (output files already exist)"
        )
        return message


def find_fasta_leftovers(fastas_to_extract):
    """
    The directory '06_assembly_annotated' will be searched within each sample directory, if the
    file 'leftover_contigs.fasta.gz' is found it will be used for clustering, if not it will revert
    to the main assembly file of the sample (unfiltered contigs from MEGAHIT)
    """
    fastas_to_cluster = {}
    num_leftovers = 0
    for sample in fastas_to_extract:
        leftovers_file = Path(fastas_to_extract[sample]["sample_dir"],
                              "06_assembly_annotated",
                              "leftover_contigs.fasta.gz")
        if leftovers_file.is_file():
            fastas_to_cluster[sample] = leftovers_file.resolve()
            num_leftovers += 1
        else:
            fastas_to_cluster[sample] = fastas_to_extract[sample]["assembly_path"]
    return fastas_to_cluster, num_leftovers


def rehead_and_concatenate_fastas(
        fastas_to_cluster, clustering_dir, clustering_input_file, threads_max, show_less
):
    """
    Since all the FASTAs coming from all the samples will be joined in a single file for clustering,
    we have to include the sample names in the headers to be able to distinguish them later.
    The sample's name + '|' + original contig name will be used as new headers, so try to avoid
    the use of '|' to name the samples. The descriptions will be lost (MMseqs ignores them afaik).
    """
    start = time.time()
    rehead_params = []
    for sample in fastas_to_cluster:
        rehead_params.append((sample, fastas_to_cluster[sample], clustering_dir))
    tqdm_parallel_async_run(rehead_fasta_with_sample_name, rehead_params,
                            "Preparing FASTA assemblies for clustering",
                            "FASTA pre-processing completed", "file", threads_max, show_less)
    log.log("")
    fastas_to_concatenate = list(clustering_dir.glob("*_leftover_contigs.fasta"))
    log.log(bold("Concatenating input for clustering:"))
    tqdm_cols = min(shutil.get_terminal_size().columns, 120)
    with tqdm(total=len(fastas_to_concatenate), ncols=tqdm_cols, unit="file") as pbar:
        for cat_fasta in fastas_to_concatenate:
            dict_to_fasta(fasta_to_dict(cat_fasta), clustering_input_file, append=True)
            cat_fasta.unlink()
            pbar.update()
    log.log(bold(
        f" \u2514\u2500\u2192 File '{clustering_input_file.name}'"
        f" prepared in {elapsed_time(time.time() - start)}(s)"
    ))


def rehead_fasta_with_sample_name(sample_name, sample_fasta_path, clustering_dir):
    """
    Prepend sample name to sequence header followed by a '|'
    """
    start = time.time()
    sample_fasta = fasta_to_dict(sample_fasta_path)
    reheaded_fasta = {}
    reheaded_fasta_file = Path(clustering_dir, f"{sample_name}_leftover_contigs.fasta")
    for header in sample_fasta:
        reheaded_fasta[f"{sample_name}|{header}"] = dict(sample_fasta[header])
    dict_to_fasta(reheaded_fasta, reheaded_fasta_file)
    message = (
        f"'{sample_name}': file '{sample_fasta_path.name}'"
        f" reheaded [{elapsed_time(time.time() - start)}]"
    )
    return message


def mmseqs2_cluster(
        mmseqs2_path, clustering_input_file, clustering_dir, clust_min_identity, clust_min_coverage,
        clust_gap_open, clust_gap_extend, clust_max_seq_len, clust_tmp_dir, threads, keep_all
):
    """
    Run MMseqs easy-linclust but perform some parameter checking/conversion before, the FASTA input
    file has to be decompressed, we can compress it afterwards
    """
    start = time.time()
    if not 0 < clust_min_identity <= 1:
        clust_min_identity = min(1.0, round((abs(clust_min_identity) / 100), 3))
    if not 0 < clust_min_coverage <= 1:
        clust_min_coverage = min(1.0, round((abs(clust_min_coverage) / 100), 3))

    files_prefix = f"clust_id{clust_min_identity * 100:.2f}_cov{clust_min_coverage * 100:.2f}"
    result_prefix = f"{Path(clustering_dir, files_prefix)}"
    mmseqs2_cmd = [
        mmseqs2_path,
        "easy-linclust",
        f"{clustering_input_file}",
        f"{result_prefix}",
        f"{clust_tmp_dir}",
        "--min-seq-id", f"{clust_min_identity}",
        "-c", f"{clust_min_coverage}",
        "--cov-mode", f"{settings.MMSEQS2_COV_MODE}",
        "--cluster-mode", "0",
        "--gap-open", f"{max(1, clust_gap_open)}",
        "--gap-extend", f"{max(1, clust_gap_extend)}",
        "--kmer-per-seq-scale", f"{settings.MMSEQS2_KMER_PER_SEQ_SCALE}",
        "--threads", f"{threads}"
    ]
    mmseqs2_log_file = Path(clustering_dir, f"{files_prefix}.log")
    mmseqs2_thread = ElapsedTimeThread()
    mmseqs2_thread.start()
    with open(mmseqs2_log_file, "w") as mmseqs2_log:
        mmseqs2_log.write(f"Captus' MMseqs2 Command:\n  {' '.join(mmseqs2_cmd)}\n\n")
    with open(mmseqs2_log_file, "a") as mmseqs2_log:
        subprocess.run(mmseqs2_cmd, stdout=mmseqs2_log, stdin=mmseqs2_log)
    num_clusters = 0
    with open(Path(clustering_dir, f"{files_prefix}_rep_seq.fasta"), "r") as clusters:
        for line in clusters:
            if line.startswith(">"):
                num_clusters += 1
    mmseqs2_thread.stop()
    mmseqs2_thread.join()
    print()

    message = bold(f" \u2514\u2500\u2192 Clustering completed: [{elapsed_time(time.time() - start)}]")
    log.log(message)
    return result_prefix, num_clusters


def filter_clusters(clust_prefix, num_clusters, clust_min_samples):
    """
    Write clusters that have at least 'clust_min_samples' to 'clustering_dir'
    """
    start = time.time()

    clusters_raw = split_mmseqs_clusters_file(Path(f"{clust_prefix}_all_seqs.fasta"))

    rep_seq_filtered = []
    passed = []
    failed = []
    singletons = []  # singletons have a single sequence, or just sequences from the same sample

    tqdm_cols = min(shutil.get_terminal_size().columns, 120)
    with tqdm(total=len(clusters_raw), ncols=tqdm_cols, unit="cluster") as pbar:
        for cluster in clusters_raw:
            sample_names = [cluster[i][1:].split("|")[0] for i in range(0, len(cluster), 2)]
            if len(set(sample_names)) >= clust_min_samples:
                passed.append(cluster)
                if len(set(sample_names)) == 1:
                    singletons.append(cluster)
                rep_seq_filtered += [cluster[0], cluster[1]]
                pbar.update()
                continue
            elif len(set(sample_names)) == 1:
                singletons.append(cluster)
                pbar.update()
                continue
            else:
                failed.append(cluster)
                pbar.update()
                continue

    if singletons:
        singletons_file = Path(f"{clust_prefix}_singletons.fasta")
        with open(singletons_file, "w") as singletons_out:
            singletons_out.write("\n//\n".join(["\n".join(s_clust) for s_clust in singletons]))
    if failed:
        failed_file = Path(f"{clust_prefix}_failed.fasta")
        with open(failed_file, "w") as failed_out:
            failed_out.write("\n//\n".join(["\n".join(f_clust) for f_clust in failed]))
    if passed:
        passed_file = Path(f"{clust_prefix}_passed.fasta")
        with open(passed_file, "w") as passed_out:
            passed_out.write("\n//\n".join(["\n".join(p_clust) for p_clust in passed]))

    p1 = bold(
        f" \u2514\u2500\u2192 Clusters retained: {len(passed)} [{elapsed_time(time.time() - start)}]"
    )
    p2 = dim(f"(Filtered out: {len(singletons)} singletons and {len(failed)} by '--cl_min_samples')")
    log.log(f"{p1} {p2}")
    log.log("")

    new_names = ["captus_cluster_name\toriginal_sequence_name"]
    num_digits = len(str(len(rep_seq_filtered) // 2))
    num_cluster = 0
    for i in range(0, len(rep_seq_filtered), 2):
        num_cluster += 1

        # Make new names for a reference file for a Captus' miscellaneous DNA extraction:
        # >sample_name-marker_name
        new_header = (
            f'{rep_seq_filtered[i].split("|")[0]}{settings.REFERENCE_CLUSTER_SEPARATOR}'
            f"captus_{num_cluster:0{num_digits}}"
        )

        # Record new and old header in 'names_equivalence_file'
        new_names.append(f"{new_header[1:]}\t{rep_seq_filtered[i][1:]}")

        # Replace the header before writing the FASTA reference
        rep_seq_filtered[i] = new_header

    # Save a list of new reference names and their respective original sequence names
    names_equivalence_file = Path(f"{clust_prefix}_original_names.tsv")
    with open(names_equivalence_file, "w") as names_out:
        names_out.write("\n".join(new_names))

    # Save the new reference file
    captus_clusters_file = Path(f"{clust_prefix}_captus_clusters_refs.fasta")
    with open(captus_clusters_file, "w") as reps_out:
        reps_out.write("\n".join(rep_seq_filtered))
    p1 = "Created new reference file for a miscellaneous DNA marker extraction:"
    p2 = bold(f" \u2514\u2500\u2192 {captus_clusters_file}")
    log.log(f"{p1}\n{p2}")
    return captus_clusters_file


def collect_ext_stats(out_dir):
    samples_stats = sorted(list(Path(out_dir).resolve().rglob("*_recovery_stats.tsv")))
    samples_stats = [tsv for tsv in samples_stats if tsv.parts[-2] == "06_assembly_annotated"]
    if not samples_stats:
        return None
    else:
        stats_file_out = Path(out_dir, "captus-assembly_extract.stats.tsv")
        header = (
            "\t".join(["sample_name", "marker_type", "locus",
                       "ref_name", "ref_coords", "ref_type", "ref_len_matched",
                       "hit", "pct_recovered", "pct_identity", "score", "lwscore",
                       "hit_len", "cds_len", "intron_len", "flanks_len", "frameshifts",
                       "ctg_names", "ctg_strands", "ctg_coords"]) + "\n"
        )
        with open(stats_file_out, "wt") as tsv_out:
            tsv_out.write(header)
            for tsv in samples_stats:
                with open(tsv, "rt") as tsv_in:
                    for line in tsv_in:
                        if line != header:
                            tsv_out.write(line)
        return stats_file_out


def cleanup_clustering_dir(clustering_dir, clust_prefix, clust_tmp_dir, threads, keep_all):
    """
    Erase junk, compress big files
    """
    shutil.rmtree(clust_tmp_dir, ignore_errors=True)
    f_all_seqs = Path(f"{clust_prefix}_all_seqs.fasta")
    f_captus_refs = Path(f"{clust_prefix}_captus_clusters_refs.fasta")
    f_cluster_tsv = Path(f"{clust_prefix}_cluster.tsv")
    f_failed = Path(f"{clust_prefix}_failed.fasta")
    f_original_names = Path(f"{clust_prefix}_original_names.tsv")
    f_passed = Path(f"{clust_prefix}_passed.fasta")
    f_rep_seq = Path(f"{clust_prefix}_rep_seq.fasta")
    f_singletons = Path(f"{clust_prefix}_singletons.fasta")
    f_cluster_input = Path(clustering_dir, "clustering_input.fasta")
    files_to_compress = [f_failed, f_passed, f_singletons]  # Don't compress 'f_captus_refs' for now
    files_to_delete = [f_all_seqs, f_cluster_tsv, f_rep_seq, f_cluster_input]
    if keep_all:
        files_to_compress += files_to_delete
    else:
        for f_del in files_to_delete:
            if f_del.exists():
                f_del.unlink()
    compress_list_files(files_to_compress, threads)

