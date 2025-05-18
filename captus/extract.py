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

import json
import math
import platform
import random
import shutil
import statistics
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from tqdm import tqdm

from . import log, settings
from .bioformats import (
    blat_misc_dna_psl_to_dict,
    dict_to_fasta,
    fasta_headers_to_spades,
    fasta_to_dict,
    fasta_type,
    fix_premature_stops,
    import_busco_odb1x,
    mmseqs_cluster,
    parse_psl_record,
    rehead_root_msa,
    scipio_yaml_to_dict,
    split_mmseqs_clusters_file,
    translate_fasta_dict,
    write_gff3,
)
from .misc import (
    bioperl_get_version,
    blat_path_version,
    bold,
    bold_green,
    bold_yellow,
    compress_list_files,
    dim,
    dir_is_empty,
    elapsed_time,
    file_is_empty,
    format_dep_msg,
    has_valid_ext,
    mafft_path_version,
    make_output_dir,
    make_tmp_dir_within,
    mmseqs_path_version,
    python_library_check,
    quit_with_error,
    red,
    remove_formatting,
    scipio_path_version,
    set_ram,
    set_threads,
    successful_exit,
    tqdm_parallel_async_run,
    tqdm_parallel_nested_run,
    tqdm_serial_run,
    yaml_perl_get_version,
)
from .version import __version__


def extract(full_command, args):
    captus_start = time.time()
    out_dir, out_dir_msg = make_output_dir(args.out)
    log.logger = log.Log(Path(args.out, "captus-extract.log"), stdout_verbosity_level=1)
    mar = 25  # Margin for aligning parameters and values

    ################################################################################################
    ############################################################################### STARTING SECTION
    log.log_section_header("Starting Captus-assembly: EXTRACT", single_newline=False)
    log.log_explanation(
        "Welcome to the marker extraction step of Captus-assembly. In this step, Captus will use"
        " Scipio to search within your FASTA assemblies and recover any set of reference proteins"
        " provided through '--nuc_refs', '--ptd_refs', and/or '--mit_refs'. Captus includes some"
        " reference protein sets to work with Plants like the 'Angiosperms353' protein set for"
        " nuclear genes, and two organellar protein reference sets 'SeedPlantsPTD' for plastids,"
        " and 'SeedPlantsMIT' for mitochondria.",
        extra_empty_lines_after=0,
    )
    log.log_explanation(
        "If you have references that are not proteins (e.g. non-coding regions, gene sequences"
        " including introns, full mRNAs, etc.) You can provide these miscellaneous DNA references"
        " with '--dna_refs' and Captus will use BLAT to find and extract matches in your assemblies.",
        extra_empty_lines_after=0,
    )
    log.log_explanation(
        "Finally, if you want to explore the usefulness of the contigs that were not hit by any"
        " protein or other DNA references after the Captus' extraction process, or simply if you do"
        " not have reference sets to test, you can try the option '--cluster_leftovers' to attempt"
        " sequence clustering across samples in order to discover homologous markers.",
        extra_empty_lines_after=0,
    )
    log.log_explanation("For more information, please see https://github.com/edgardomortiz/Captus")

    log.log(f"{'Captus version':>{mar}}: {bold(f'v{__version__}')}")
    log.log(f"{'Command':>{mar}}: {bold(full_command)}")
    log.log(f"{'OS':>{mar}}: {bold(platform.platform())}")
    log.log(f"{'Host':>{mar}}: {bold(platform.node())}")
    tsv_comment = f"#Captus v{__version__}\n#Command: {full_command}\n"
    ram_B, ram_MB, ram_GB, ram_GB_total = set_ram(args.ram)
    log.log(f"{'Max. RAM':>{mar}}: {bold(f'{ram_GB:.1f}GB')} {dim(f'(out of {ram_GB_total:.1f}GB)')}")
    threads_max, threads_total = set_threads(args.threads)
    log.log(f"{'Max. Threads':>{mar}}: {bold(threads_max)} {dim(f'(out of {threads_total})')}")
    log.log("")

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
        _, mmseqs_version, mmseqs_status = mmseqs_path_version(args.mmseqs_path)
        _, mafft_version, mafft_status = mafft_path_version(args.mafft_path)
        _, blat_version, blat_status = blat_path_version(args.blat_path)
    else:
        skip_clustering = True
        mmseqs_version, mmseqs_status = "", "not used"
        mafft_version, mafft_status = "", "not used"

    log.log(f"{'Dependencies':>{mar}}:")
    log.log(format_dep_msg(f"{'Scipio   ':>{mar}}: ", scipio_version, scipio_status))
    branch_bioperl = "\u251c\u2500BioPerl"
    log.log(format_dep_msg(f"{branch_bioperl:>{mar}}: ", bioperl_version, bioperl_status))
    branch_yaml = "\u2514\u2500YAML   "
    log.log(format_dep_msg(f"{branch_yaml:>{mar}}: ", yaml_version, yaml_status))
    log.log(format_dep_msg(f"{'BLAT     ':>{mar}}: ", blat_version, blat_status))
    log.log(format_dep_msg(f"{'MMseqs2  ':>{mar}}: ", mmseqs_version, mmseqs_status))
    log.log(format_dep_msg(f"{'MAFFT    ':>{mar}}: ", mafft_version, mafft_status))
    log.log("")

    log.log(f"{'Python libraries':>{mar}}:")
    numpy_found, numpy_version, numpy_status = python_library_check("numpy")
    pandas_found, pandas_version, pandas_status = python_library_check("pandas")
    plotly_found, plotly_version, plotly_status = python_library_check("plotly")
    log.log(format_dep_msg(f"{'numpy':>{mar}}: ", numpy_version, numpy_status))
    log.log(format_dep_msg(f"{'pandas':>{mar}}: ", pandas_version, pandas_status))
    log.log(format_dep_msg(f"{'plotly':>{mar}}: ", plotly_version, plotly_status))
    log.log("")

    captus_assemblies_dir, _ = make_output_dir(args.captus_assemblies_dir)
    _, asms_dir_status_msg = status_captus_assemblies_dir(captus_assemblies_dir, args.fastas, mar)
    log.log(f"{'Captus assemblies dir':>{mar}}: {bold(captus_assemblies_dir)}")
    log.log(asms_dir_status_msg)
    log.log("")

    # Last dependency verification before modifying Captus' assembly directory
    if any(
        dep_status == "not found"
        for dep_status in [scipio_status, bioperl_status, yaml_status, blat_status]
    ):
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
                "MMseqs2 could not be found, please check your '--mmseqs_path'. Additionally, no"
                " reference protein or nucleotide sets were provided for extraction."
            )
        else:
            log.log(
                f"{bold('WARNING:')} MMseqs2 could not be found, please check your '--mmseqs_path'."
                " Extraction of reference protein sets will be attempted."
            )
    if mafft_status == "not found":
        skip_clustering = True
        if skip_extraction:
            quit_with_error(
                "MAFFT could not be found, please check your '--mafft_path'. Additionally, no"
                " reference protein or nucleotide sets were provided for extraction."
            )
        else:
            log.log(
                f"{bold('WARNING:')} MAFFT could not be found, please check your '--mafft_path'."
                " Extraction of reference protein sets will be attempted."
            )

    # Check and import extra FASTA assemblies
    fastas_to_extract, skipped_extract = find_fasta_assemblies(captus_assemblies_dir, out_dir)
    if args.fastas:
        log.log_explanation(
            "Now Captus will verify the additional assembly files provided with '--fastas'. A new"
            f" sample directory will be created within '{out_dir}' for each of the new FASTA"
            " assemblies to contain the imported file"
        )
        asms_before_import = len(fastas_to_extract)
        find_and_copy_fastas(
            args.fastas, captus_assemblies_dir, args.overwrite, threads_max, args.show_less
        )
        fastas_to_extract, skipped_extract = find_fasta_assemblies(captus_assemblies_dir, out_dir)
        log.log("")
        asms_imported = len(fastas_to_extract) - asms_before_import
        log.log(f"{'Assemblies imported':>{mar}}: {bold(asms_imported)}")
    log.log(f"{'Total assemblies found':>{mar}}: {bold(len(fastas_to_extract))}")
    log.log("")
    log.log(f"{'Output directory':>{mar}}: {bold(out_dir)}")
    log.log(f"{'':>{mar}}  {dim(out_dir_msg)}")
    log.log("")

    if skipped_extract:
        log.log(f"{bold('WARNING:')} {len(skipped_extract)} sample(s) will be skipped")
        for msg in skipped_extract:
            log.log(msg)
        log.log("")

    ################################################################################################
    ############################################################################# EXTRACTION SECTION
    log.log_section_header("Protein and DNA Markers Extraction from FASTA assemblies")
    log.log_explanation(
        "Now Captus will search within your FASTA assemblies for the set of reference proteins"
        " provided by '--nuc_refs', '--ptd_refs', and/or '--mit_refs' and extract the sequences both"
        " in nucleotide and translated to aminoacid using the specified '--nuc_transtable',"
        " '--ptd_transtable' and '--mit_transtable' respectively.",
        extra_empty_lines_after=0,
    )
    log.log_explanation(
        "Additionally, miscellaneous DNA references provided with '--dna_refs' will also be extracted."
    )
    if skip_extraction:
        log.log(
            red(
                "Skipping extraction step... (to enable provide any/all of these options: '--nuc_refs',"
                " '--ptd_refs', '--mit_refs', '--dna_refs', and verify that all the software"
                " dependencies are correctly installed)"
            )
        )
        log.log("")
    else:
        if not fastas_to_extract:
            quit_with_error(
                "Captus could not find FASTA assemblies to process, please verify your "
                "'--captus_assemblies_dir' and/or '--fastas' argument"
            )

        num_refs = 0
        if args.nuc_refs:
            num_refs += 1
            if args.nuc_refs.lower() in settings.PROT_REF_TARGETS["NUC"]:
                args.nuc_transtable = settings.PROT_REF_TARGETS["NUC"][args.nuc_refs.lower()][
                    "transtable"
                ]
        if args.ptd_refs:
            num_refs += 1
            if args.ptd_refs.lower() in settings.PROT_REF_TARGETS["PTD"]:
                args.ptd_transtable = settings.PROT_REF_TARGETS["PTD"][args.ptd_refs.lower()][
                    "transtable"
                ]
        if args.mit_refs:
            num_refs += 1
            if args.mit_refs.lower() in settings.PROT_REF_TARGETS["MIT"]:
                args.mit_transtable = settings.PROT_REF_TARGETS["MIT"][args.mit_refs.lower()][
                    "transtable"
                ]
        num_prot_extractions = len(fastas_to_extract) * num_refs

        protein_refs = prepare_protein_refs(
            args.nuc_refs,
            args.ptd_refs,
            args.mit_refs,
            args.nuc_transtable,
            args.ptd_transtable,
            args.mit_transtable,
            out_dir,
        )

        log.log("")
        if any(
            [
                protein_refs["NUC"]["AA_path"],
                protein_refs["PTD"]["AA_path"],
                protein_refs["MIT"]["AA_path"],
            ]
        ):
            log.log(bold(f"{'PROTEIN OPTIONS':>{mar}}:"))
            prot_concurrent, prot_threads, prot_ram = adjust_concurrency(
                fastas_to_extract, num_prot_extractions, args.concurrent, threads_max, ram_B, "protein"
            )
            if protein_refs["NUC"]["AA_path"]:
                nuc_query = fasta_to_dict(protein_refs["NUC"]["AA_path"])
                nuc_query_parts_paths = split_refs(nuc_query, out_dir, "NUC", prot_threads)
                nuc_query_info = reference_info(nuc_query)
            if protein_refs["PTD"]["AA_path"]:
                ptd_query = fasta_to_dict(protein_refs["PTD"]["AA_path"])
                ptd_query_parts_paths = split_refs(ptd_query, out_dir, "PTD", prot_threads)
                ptd_query_info = reference_info(ptd_query)
            if protein_refs["MIT"]["AA_path"]:
                mit_query = fasta_to_dict(protein_refs["MIT"]["AA_path"])
                mit_query_parts_paths = split_refs(mit_query, out_dir, "MIT", prot_threads)
                mit_query_info = reference_info(mit_query)
            log.log(f"{'Max. loci for Scipio x2':>{mar}}: {bold(args.max_loci_scipio_x2)}")
            log.log(f"{'Predict dubious introns':>{mar}}: {bold(args.predict)}")
            log.log(f"{'Concurrent extractions':>{mar}}: {bold(prot_concurrent)}")
            log.log(f"{'RAM per extraction':>{mar}}: {bold(f'{prot_ram / 1024**3:.1f}GB')}")
            log.log(f"{'Threads per extraction':>{mar}}: {bold(prot_threads)}")
            log.log("")
        log.log(bold(f"{'Nuclear proteins':>{mar}}:"))
        log.log(f"{'reference':>{mar}}: {protein_refs['NUC']['AA_msg']}")
        if protein_refs["NUC"]["AA_path"]:
            log.log(f"{'reference info':>{mar}}: {nuc_query_info['info_msg']}")
            log.log(f"{'translation table':>{mar}}: {bold(args.nuc_transtable)}")
            log.log(f"{'min_score':>{mar}}: {bold(args.nuc_min_score)}")
            nuc_min_identity = adjust_min_identity(args.nuc_min_identity, args.nuc_transtable)
            log.log(f"{'min_identity':>{mar}}: {bold(nuc_min_identity)}")
            nuc_min_coverage = adjust_min_coverage(args.nuc_min_coverage)
            log.log(f"{'min_coverage':>{mar}}: {bold(nuc_min_coverage)}")
            nuc_dt, nuc_dt_msg = depth_tolerance_check(args.nuc_depth_tolerance, args.ignore_depth)
            log.log(f"{'depth_tolerance':>{mar}}: {nuc_dt_msg}")
        log.log("")
        log.log(bold(f"{'Plastidial proteins':>{mar}}:"))
        log.log(f"{'reference':>{mar}}: {protein_refs['PTD']['AA_msg']}")
        if protein_refs["PTD"]["AA_path"]:
            log.log(f"{'reference info':>{mar}}: {ptd_query_info['info_msg']}")
            log.log(f"{'translation table':>{mar}}: {bold(args.ptd_transtable)}")
            log.log(f"{'min_score':>{mar}}: {bold(args.ptd_min_score)}")
            ptd_min_identity = adjust_min_identity(args.ptd_min_identity, args.ptd_transtable)
            log.log(f"{'min_identity':>{mar}}: {bold(ptd_min_identity)}")
            ptd_min_coverage = adjust_min_coverage(args.ptd_min_coverage)
            log.log(f"{'min_coverage':>{mar}}: {bold(ptd_min_coverage)}")
            ptd_dt, ptd_dt_msg = depth_tolerance_check(args.ptd_depth_tolerance, args.ignore_depth)
            log.log(f"{'depth_tolerance':>{mar}}: {ptd_dt_msg}")
        log.log("")
        log.log(bold(f"{'Mitochondrial proteins':>{mar}}:"))
        log.log(f"{'reference':>{mar}}: {protein_refs['MIT']['AA_msg']}")
        if protein_refs["MIT"]["AA_path"]:
            log.log(f"{'reference info':>{mar}}: {mit_query_info['info_msg']}")
            log.log(f"{'translation table':>{mar}}: {bold(args.mit_transtable)}")
            log.log(f"{'min_score':>{mar}}: {bold(args.mit_min_score)}")
            mit_min_identity = adjust_min_identity(args.mit_min_identity, args.mit_transtable)
            log.log(f"{'min_identity':>{mar}}: {bold(mit_min_identity)}")
            mit_min_coverage = adjust_min_coverage(args.mit_min_coverage)
            log.log(f"{'min_coverage':>{mar}}: {bold(mit_min_coverage)}")
            mit_dt, mit_dt_msg = depth_tolerance_check(args.mit_depth_tolerance, args.ignore_depth)
            log.log(f"{'depth_tolerance':>{mar}}: {mit_dt_msg}")
        log.log("")
        log.log("")

        dna_ref_size = 0
        dna_ref = prepare_dna_refs(args.dna_refs)
        log.log(bold(f"{'MISCELLANEOUS DNA OPTIONS':>{mar}}:"))
        if dna_ref["DNA"]["NT_path"]:
            num_dna_extractions = len(fastas_to_extract)
            dna_concurrent, dna_threads, dna_ram = adjust_concurrency(
                fastas_to_extract, num_dna_extractions, args.concurrent, threads_max, ram_B, "dna"
            )
            dna_query = fasta_to_dict(dna_ref["DNA"]["NT_path"])
            dna_query_parts_paths = split_refs(dna_query, out_dir, "DNA", dna_threads)
            dna_query_info = reference_info(dna_query)
            dna_ref_size = max(dna_ref_size, dna_query_info["total_size"])
            log.log(f"{'Concurrent extractions':>{mar}}: {bold(dna_concurrent)}")
            log.log(f"{'RAM per extraction':>{mar}}: {bold(f'{dna_ram / 1024**3:.1f}GB')}")
            log.log(f"{'Threads per extraction':>{mar}}: {bold(dna_threads)}")
            log.log("")
        log.log(f"{'reference':>{mar}}: {dna_ref['DNA']['NT_msg']}")
        if dna_ref["DNA"]["NT_path"]:
            log.log(f"{'reference info':>{mar}}: {dna_query_info['info_msg']}")
            log.log(f"{'min_identity':>{mar}}: {bold(args.dna_min_identity)}")
            log.log(f"{'min_coverage':>{mar}}: {bold(args.dna_min_coverage)}")
            dna_dt, dna_dt_msg = depth_tolerance_check(args.dna_depth_tolerance, args.ignore_depth)
            log.log(f"{'depth_tolerance':>{mar}}: {dna_dt_msg}")
        log.log("")
        log.log("")

        log.log(bold(f"{'OUTPUT OPTIONS':>{mar}}:"))
        log.log(f"{'Ignore depth of coverage':>{mar}}: {bold(args.ignore_depth)}")
        log.log(f"{'Disable contig stitching':>{mar}}: {bold(args.disable_stitching)}")
        max_paralogs_msg = dim("(Keep all paralogs)") if args.max_paralogs == -1 else ""
        log.log(f"{'Max. paralogs':>{mar}}: {bold(args.max_paralogs)} {max_paralogs_msg}")
        log.log(f"{'Paralog tolerance':>{mar}}: {bold(args.paralog_tolerance)}")
        loci_files_msg = ""
        if args.max_loci_files == 0:
            loci_files_msg = dim("(Do not write separate loci files per sample)")
        log.log(f"{'Max. separate loci files':>{mar}}: {bold(args.max_loci_files)} {loci_files_msg}")
        log.log(f"{'Overwrite files':>{mar}}: {bold(args.overwrite)}")
        log.log(f"{'Keep all files':>{mar}}: {bold(args.keep_all)}")
        log.log("")
        log.log(f"{'Samples to process':>{mar}}: {bold(len(fastas_to_extract))}")

        extract_coding = bool(
            any(
                [
                    protein_refs["NUC"]["AA_path"],
                    protein_refs["PTD"]["AA_path"],
                    protein_refs["MIT"]["AA_path"],
                ]
            )
        )

        if extract_coding or bool(dna_ref["DNA"]["NT_path"]):
            scipio_params = []
            blat_params = []
            cleanup_params = []
            for sample in sorted(
                fastas_to_extract, key=lambda x: fastas_to_extract[x]["assembly_size"], reverse=True
            ):
                if protein_refs["NUC"]["AA_path"]:
                    scipio_params.append(
                        (
                            args.scipio_path,
                            args.nuc_min_score,
                            nuc_min_identity,
                            args.nuc_min_coverage,
                            args.blat_path,
                            fastas_to_extract[sample]["assembly_path"],
                            fastas_to_extract[sample]["sample_dir"],
                            sample,
                            protein_refs["NUC"]["AA_path"],
                            nuc_query,
                            nuc_query_parts_paths,
                            nuc_query_info,
                            "NUC",
                            args.nuc_transtable,
                            args.ignore_depth,
                            nuc_dt,
                            args.disable_stitching,
                            args.max_loci_files,
                            args.max_loci_scipio_x2,
                            args.max_paralogs,
                            args.paralog_tolerance,
                            args.predict,
                            prot_threads,
                            prot_ram,
                            tsv_comment,
                            args.debug,
                            args.overwrite,
                            args.keep_all,
                        )
                    )
                if protein_refs["PTD"]["AA_path"]:
                    scipio_params.append(
                        (
                            args.scipio_path,
                            args.ptd_min_score,
                            ptd_min_identity,
                            args.ptd_min_coverage,
                            args.blat_path,
                            fastas_to_extract[sample]["assembly_path"],
                            fastas_to_extract[sample]["sample_dir"],
                            sample,
                            protein_refs["PTD"]["AA_path"],
                            ptd_query,
                            ptd_query_parts_paths,
                            ptd_query_info,
                            "PTD",
                            args.ptd_transtable,
                            args.ignore_depth,
                            ptd_dt,
                            args.disable_stitching,
                            args.max_loci_files,
                            args.max_loci_scipio_x2,
                            args.max_paralogs,
                            args.paralog_tolerance,
                            args.predict,
                            prot_threads,
                            prot_ram,
                            tsv_comment,
                            args.debug,
                            args.overwrite,
                            args.keep_all,
                        )
                    )
                if protein_refs["MIT"]["AA_path"]:
                    scipio_params.append(
                        (
                            args.scipio_path,
                            args.mit_min_score,
                            mit_min_identity,
                            args.mit_min_coverage,
                            args.blat_path,
                            fastas_to_extract[sample]["assembly_path"],
                            fastas_to_extract[sample]["sample_dir"],
                            sample,
                            protein_refs["MIT"]["AA_path"],
                            mit_query,
                            mit_query_parts_paths,
                            mit_query_info,
                            "MIT",
                            args.mit_transtable,
                            args.ignore_depth,
                            mit_dt,
                            args.disable_stitching,
                            args.max_loci_files,
                            args.max_loci_scipio_x2,
                            args.max_paralogs,
                            args.paralog_tolerance,
                            args.predict,
                            prot_threads,
                            prot_ram,
                            tsv_comment,
                            args.debug,
                            args.overwrite,
                            args.keep_all,
                        )
                    )
                if dna_ref["DNA"]["NT_path"]:
                    blat_params.append(
                        (
                            args.blat_path,
                            args.dna_min_identity,
                            args.dna_min_coverage,
                            fastas_to_extract[sample]["assembly_path"],
                            fastas_to_extract[sample]["sample_dir"],
                            sample,
                            dna_ref["DNA"]["NT_path"],
                            dna_query,
                            dna_query_parts_paths,
                            dna_query_info,
                            "DNA",
                            args.ignore_depth,
                            dna_dt,
                            args.disable_stitching,
                            args.max_loci_files,
                            args.max_paralogs,
                            args.paralog_tolerance,
                            tsv_comment,
                            dna_threads,
                            dna_ram,
                            args.debug,
                            args.overwrite,
                            args.keep_all,
                        )
                    )
                cleanup_params.append(
                    (
                        sample,
                        fastas_to_extract[sample]["sample_dir"],
                        fastas_to_extract[sample]["assembly_path"],
                        tsv_comment,
                        args.keep_all,
                        args.overwrite,
                        skip_clustering,
                    )
                )
            log.log(
                f"{'Extractions to process':>{mar}}:"
                f" {bold(f'{len(scipio_params)} protein')} and"
                f" {bold(f'{len(blat_params)} nucleotide')}"
            )
            log.log("")

            if protein_refs["NUC"]["AA_path"]:
                log.log(
                    f"{'Nuclear proteins':>{mar}}:"
                    f" {bold(f'{out_dir}/[Sample_name]__captus-ext/')}"
                    f"{bold(settings.MARKER_DIRS['NUC'])}"
                )
            if protein_refs["PTD"]["AA_path"]:
                log.log(
                    f"{'Plastidial proteins':>{mar}}:"
                    f" {bold(f'{out_dir}/[Sample_name]__captus-ext/')}"
                    f"{bold(settings.MARKER_DIRS['PTD'])}"
                )
            if protein_refs["MIT"]["AA_path"]:
                log.log(
                    f"{'Mitochondrial proteins':>{mar}}:"
                    f" {bold(f'{out_dir}/[Sample_name]__captus-ext/')}"
                    f"{bold(settings.MARKER_DIRS['MIT'])}"
                )
            if dna_ref["DNA"]["NT_path"]:
                log.log(
                    f"{'Miscellaneous DNA markers':>{mar}}:"
                    f" {bold(f'{out_dir}/[Sample_name]__captus-ext/')}"
                    f"{bold(settings.MARKER_DIRS['DNA'])}"
                )
            log.log(
                f"{'Annotated assemblies':>{mar}}: "
                f"{bold(f'{out_dir}/[Sample_name]__captus-ext/06_assembly_annotated')}"
            )
            log.log(f"{'':>{mar}}  {dim('(These output directories will be created for each sample)')}")
            log.log("")

            json_path = update_refs_json({**protein_refs, **dna_ref}, out_dir)  # concat two dict
            log.log(f"{'Paths to references used':>{mar}}: {bold(f'{json_path}')}")
            log.log("")

            if scipio_params:
                d_msg = "Extracting protein-coding markers with Scipio"
                f_msg = "Protein-coding markers: finished processing"
                if args.debug:
                    tqdm_serial_run(
                        scipio_coding, scipio_params, d_msg, f_msg, "extraction", args.show_less
                    )
                else:
                    tqdm_parallel_nested_run(
                        scipio_coding,
                        scipio_params,
                        d_msg,
                        f_msg,
                        "extraction",
                        prot_concurrent,
                        args.show_less,
                    )
                log.log("")

            if blat_params:
                d_msg = "Extracting miscellaneous DNA markers with BLAT"
                f_msg = "Miscellaneous DNA markers: finished processing"
                if args.debug:
                    tqdm_serial_run(
                        blat_misc_dna, blat_params, d_msg, f_msg, "extraction", args.show_less
                    )
                else:
                    tqdm_parallel_nested_run(
                        blat_misc_dna,
                        blat_params,
                        d_msg,
                        f_msg,
                        "extraction",
                        dna_concurrent,
                        args.show_less,
                    )
                log.log("")

            d_msg = "Merging GFFs, summarizing recovery stats, and cleaning up"
            f_msg = "Merged GFF and summary recovery stats created"
            cleanup_concurrent = min(settings.MAX_HDD_WRITE_INSTANCES, threads_max)
            if args.debug:
                tqdm_serial_run(
                    cleanup_post_extraction, cleanup_params, d_msg, f_msg, "sample", args.show_less
                )
            else:
                tqdm_parallel_async_run(
                    cleanup_post_extraction,
                    cleanup_params,
                    d_msg,
                    f_msg,
                    "sample",
                    cleanup_concurrent,
                    args.show_less,
                )
            if not args.keep_all:
                shutil.rmtree(Path(out_dir, settings.REF_TARGETS_SPLIT_DIR), ignore_errors=True)
            log.log("")

        # Nothing to extract
        else:
            log.log(
                red(
                    "Skipping extraction step... (no valid references found, verify the paths provided"
                    " through '--nuc_refs', '--ptd_refs', '--mit_refs', '--dna_refs')"
                )
            )
            log.log("")

    ################################################################################################
    ############################################################################# CLUSTERING SECTION
    log.log_section_header("Clustering and Extracting Markers Across Samples with MMseqs2/BLAT")
    log.log_explanation(
        "Now Captus will cluster the assembly sequences across samples to discover homologous"
        " markers. If protein and/or miscellaneous DNA references sequences were extracted in the"
        " previous step, the '06_assembly_annotated/leftover_contigs.fasta.gz' files from each"
        " sample will be used for clustering, otherwise the '01_assembly/assembly.fasta' files are"
        " used. USE WITH CAUTION: Clustering time and/or success will depend on the size of your"
        " assembly files, it should be safe to use if you assembled captured data, RNAseq reads, or"
        " even genome skimming data.",
        extra_empty_lines_after=0,
    )
    log.log_explanation(
        "After clustering has been completed, only the clusters with at least '--cl_min_samples'"
        " will be considered and a single representative sequence from each cluster will be selected"
        " to create a new reference file to perform a final DNA marker extraction across all samples."
        " The extraction of these new markers is done with the same settings chosen for"
        " miscellaneous DNA markers."
    )
    if skip_clustering:
        log.log(
            red(
                "Skipping clustering step... (to enable use '--cluster_leftovers', and verify that all"
                " the software dependencies are correctly installed)"
            )
        )
        log.log("")
    else:
        if not fastas_to_extract:
            quit_with_error(
                "Captus could not find FASTA assemblies to process, please verify your "
                "'--captus_assemblies_dir' and/or '--fastas' argument"
            )

        clustering_dir, clustering_dir_msg = make_output_dir(Path(out_dir, "02_clustering_data"))
        clustering_input_file = Path(clustering_dir, "clustering_input.fasta")

        # When 'cl_min_identity' is set to 'auto' it becomes 90% of 'dna_min_identity' never
        # becoming less than 75%
        if args.cl_min_identity == "auto":
            dna_min_identity = max(args.dna_min_identity, settings.MMSEQS_MIN_AUTO_MIN_IDENTITY)
            cl_min_identity = max(
                settings.MMSEQS_BLAT_DNA_IDENTITY_FACTOR * dna_min_identity,
                settings.MMSEQS_MIN_AUTO_MIN_IDENTITY,
            )
        else:
            dna_min_identity = cl_min_identity = float(args.cl_min_identity)
        if args.cl_mode != 2:
            args.cl_seq_id_mode = 0
        clust_tmp_dir = make_tmp_dir_within(args.cl_tmp_dir, "captus_mmseqs_tmp")
        fastas_to_cluster, num_leftovers = find_fasta_leftovers(fastas_to_extract)
        num_samples = len(fastas_to_cluster)
        if args.cl_min_samples == "auto":
            cl_min_samples = int(max(4, num_samples * settings.CLR_MIN_SAMPLE_PROP))
        else:
            cl_min_samples = min(int(args.cl_min_samples), num_samples)

        log.log(bold_yellow("  \u25ba STEP 1 OF 3: Clustering contigs across samples with MMseqs2"))
        log.log("")
        log.log(f"{'MMseqs2 method':>{mar}}: {bold(args.mmseqs_method)}")
        log.log(f"{'cluster_mode':>{mar}}: {bold(args.cl_mode)}")
        log.log(f"{'sensitivity':>{mar}}: {bold(args.cl_sensitivity)}")
        log.log(f"{'min_seq_id':>{mar}}: {bold(cl_min_identity)}")
        log.log(f"{'seq_id_mode':>{mar}}: {bold(args.cl_seq_id_mode)}")
        log.log(f"{'cov':>{mar}}: {bold(args.cl_min_coverage)}")
        log.log(f"{'cov_mode':>{mar}}: {bold(args.cl_cov_mode)}")
        log.log(f"{'gap_open':>{mar}}: {bold(settings.MMSEQS_GAP_OPEN)}")
        log.log(f"{'gap_extend':>{mar}}: {bold(settings.MMSEQS_GAP_EXTEND)}")
        log.log(f"{'max_seq_len':>{mar}}: {bold(args.cl_max_seq_len)}")
        log.log(f"{'tmp_dir':>{mar}}: {bold(clust_tmp_dir)}")
        log.log("")
        log.log(f"{'Min. locus length':>{mar}}: {bold(args.cl_rep_min_len)}")
        log.log(f"{'Min. samples per cluster':>{mar}}: {bold(cl_min_samples)}")
        log.log(f"{'Max. copies per cluster':>{mar}}: {bold(args.cl_max_copies)}")
        log.log(f"{'Overwrite files':>{mar}}: {bold(args.overwrite)}")
        log.log(f"{'Keep all files':>{mar}}: {bold(args.keep_all)}")
        log.log("")
        log.log(f"{'Using leftover contigs':>{mar}}: {bold(num_leftovers)}")
        log.log(f"{'Using entire assembly':>{mar}}: {bold(len(fastas_to_cluster) - num_leftovers)}")
        log.log(bold(f"{'Total samples to cluster':>{mar}}: {len(fastas_to_cluster)}"))
        log.log("")
        log.log(f"{'Clustering directory':>{mar}}: {bold(clustering_dir)}")
        log.log(f"{'':>{mar}}  {dim(clustering_dir_msg)}")
        log.log("")

        if dir_is_empty(clustering_dir) or clustering_input_file.is_file() or args.overwrite:
            if clustering_input_file.is_file() and args.overwrite is False:
                log.log(
                    f"{bold('WARNING:')} The input FASTA file for clustering was found in"
                    f" '{bold(clustering_input_file)}' and it will be used, to recreate it enable"
                    " '--overwrite'"
                )
                log.log("")
            else:
                rehead_and_concatenate_fastas(
                    fastas_to_cluster,
                    clustering_dir,
                    clustering_input_file,
                    args.cl_max_seq_len,
                    min(settings.MAX_HDD_WRITE_INSTANCES, threads_max),
                    args.show_less,
                )
            captus_cluster_refs = cluster_and_select_refs(
                num_samples,
                cl_min_samples,
                args.cl_max_copies,
                args.cl_rep_min_len,
                args.mmseqs_path,
                args.mmseqs_method,
                args.cl_mode,
                args.cl_sensitivity,
                clustering_input_file,
                clustering_dir,
                cl_min_identity,
                args.cl_seq_id_mode,
                args.cl_min_coverage,
                args.cl_cov_mode,
                clust_tmp_dir,
                args.mafft_path,
                threads_max,
                ram_MB,
                args.debug,
            )
            log.log("")
            log.log("")
            log.log(
                bold_yellow("  \u25ba STEP 2 OF 3: Extracting cluster-derived DNA markers with BLAT")
            )
            log.log("")
            clust_ref_size = 0
            clust_ref = prepare_dna_refs(captus_cluster_refs, cluster=True)
            if clust_ref["CLR"]["NT_path"]:
                num_clr_extractions = len(fastas_to_extract)
                clust_concurrent, clust_threads, clust_ram = adjust_concurrency(
                    fastas_to_extract, num_clr_extractions, args.concurrent, threads_max, ram_B, "dna"
                )
                clust_query = fasta_to_dict(clust_ref["CLR"]["NT_path"])
                clust_query_parts_paths = split_refs(clust_query, out_dir, "CLR", clust_threads)
                clust_query_info = reference_info(clust_query)
                clust_ref_size = max(clust_ref_size, clust_query_info["total_size"])
                log.log(f"{'Concurrent extractions':>{mar}}: {bold(clust_concurrent)}")
                log.log(f"{'RAM per extraction':>{mar}}: {bold(f'{clust_ram / 1024**3:.1f}GB')}")
                log.log(f"{'Threads per extraction':>{mar}}: {bold(clust_threads)}")
                log.log("")
            log.log(f"{'reference':>{mar}}: {clust_ref['CLR']['NT_msg']}")
            if clust_ref["CLR"]["NT_path"]:
                log.log(f"{'reference info':>{mar}}: {clust_query_info['info_msg']}")
                log.log(f"{'dna_min_identity':>{mar}}: {bold(dna_min_identity)}")
                log.log(f"{'dna_min_coverage':>{mar}}: {bold(args.dna_min_coverage)}")
                clr_dt, clr_dt_msg = depth_tolerance_check(args.dna_depth_tolerance, args.ignore_depth)
                log.log(f"{'depth_tolerance':>{mar}}: {clr_dt_msg}")
            log.log("")
            log.log(f"{'Overwrite files':>{mar}}: {bold(args.overwrite)}")
            log.log(f"{'Keep all files':>{mar}}: {bold(args.keep_all)}")
            log.log(f"{'Samples to process':>{mar}}: {bold(num_clr_extractions)}")
            if clust_ref["CLR"]["NT_path"]:
                blat_clusters_params = []
                for sample in sorted(
                    fastas_to_extract,
                    key=lambda x: fastas_to_extract[x]["assembly_size"],
                    reverse=True,
                ):
                    blat_clusters_params.append(
                        (
                            args.blat_path,
                            dna_min_identity,
                            args.dna_min_coverage,
                            fastas_to_extract[sample]["assembly_path"],
                            fastas_to_extract[sample]["sample_dir"],
                            sample,
                            clust_ref["CLR"]["NT_path"],
                            clust_query,
                            clust_query_parts_paths,
                            clust_query_info,
                            "CLR",
                            args.ignore_depth,
                            clr_dt,
                            args.disable_stitching,
                            args.max_loci_files,
                            args.max_paralogs,
                            args.paralog_tolerance,
                            tsv_comment,
                            clust_threads,
                            clust_ram,
                            args.debug,
                            args.overwrite,
                            args.keep_all,
                        )
                    )
                log.log(f"{'Extractions to process':>{mar}}: {bold(len(blat_clusters_params))}")
                log.log("")
                log.log(
                    f"{'Clustering output':>{mar}}:"
                    f" {bold(f'{out_dir}/[Sample_name]__captus-ext/')}"
                    f"{bold(settings.MARKER_DIRS['CLR'])}"
                )
                log.log(f"{'':>{mar}}  {dim('(A directory will be created for each sample)')}")
                log.log("")

                json_path = update_refs_json(clust_ref, out_dir)
                log.log(f"{'Paths to references used':>{mar}}: {bold(f'{json_path}')}")
                log.log("")

                d_msg = "Extracting cluster-derived DNA markers with BLAT"
                f_msg = "Cluster-derived DNA markers: finished processing"
                if args.debug:
                    tqdm_serial_run(
                        blat_misc_dna,
                        blat_clusters_params,
                        d_msg,
                        f_msg,
                        "extraction",
                        args.show_less,
                    )
                else:
                    tqdm_parallel_nested_run(
                        blat_misc_dna,
                        blat_clusters_params,
                        d_msg,
                        f_msg,
                        "extraction",
                        clust_concurrent,
                        args.show_less,
                    )
                log.log("")
                log.log("")
                log.log(
                    bold_yellow(
                        "  \u25ba STEP 3 OF 3: Final Merging of GFF Annotations and File Cleanup"
                    )
                )
                log.log("")
                cleanup_params = []
                for sample in fastas_to_extract:
                    cleanup_params.append(
                        (
                            sample,
                            fastas_to_extract[sample]["sample_dir"],
                            fastas_to_extract[sample]["assembly_path"],
                            tsv_comment,
                            args.keep_all,
                            True,  # overwrite
                            skip_clustering,
                            True,  # cluster
                        )
                    )
                d_msg = "Merging GFFs, summarizing recovery stats, and cleaning up"
                f_msg = "Merged GFF and summary recovery stats created"
                cleanup_concurrent = min(settings.MAX_HDD_WRITE_INSTANCES, threads_max)
                if args.debug:
                    tqdm_serial_run(
                        cleanup_post_extraction,
                        cleanup_params,
                        d_msg,
                        f_msg,
                        "sample",
                        args.show_less,
                    )
                else:
                    tqdm_parallel_async_run(
                        cleanup_post_extraction,
                        cleanup_params,
                        d_msg,
                        f_msg,
                        "sample",
                        cleanup_concurrent,
                        args.show_less,
                    )
                log.log("")

        # Nothing to cluster, just skip already processed
        else:
            log.log(
                red(
                    f"Skipping clustering step... Captus found output files in: '{clustering_dir}',"
                    " to replace them enable --overwrite"
                )
            )
            log.log("")

    ################################################################################################
    ############################################################################## SUMMARIZE SECTION
    log.log_section_header("Statistics Summarization")
    log.log_explanation(
        "Now Captus will collect the extraction statistics from each sample to compile a"
        " comprehensive table and a HTML report for visualization of extraction statistics."
    )
    ext_stats_tsv = collect_ext_stats(out_dir, tsv_comment)
    if ext_stats_tsv:
        log.log(f"{'Extraction statistics':>{mar}}: {bold(ext_stats_tsv)}")
        log.log("")
        if all([numpy_found, pandas_found, plotly_found]):
            from .report import build_extraction_report

            log.log_explanation("Generating Marker Extraction report...")
            ext_html_report, ext_html_msg = build_extraction_report(out_dir, ext_stats_tsv)
            log.log(f"{'Extraction report':>{mar}}: {bold(ext_html_report)}")
            log.log(f"{'':>{mar}}  {dim(ext_html_msg)}")
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
    successful_exit(
        "Captus-assembly: EXTRACT -> successfully completed"
        f" [{elapsed_time(time.time() - captus_start)}]"
    )


def adjust_concurrency(fastas_to_extract, num_samples, concurrent, threads_max, ram_B, ref_type):
    """
    Adjust the proposed number of 'concurrent' BLAT/Scipio processes so 'RAM_per_extraction' is
    never smaller than 'settings.EXTRACTION_MIN_RAM_B'. Once the right 'concurrent' has been found,
    'ram_per_extraction' is readjusted
    """
    if concurrent == "auto":
        concurrent = min(num_samples, threads_max)
    else:
        try:
            concurrent = int(concurrent)
        except ValueError:
            quit_with_error("Invalid value for '--concurrent', set it to 'auto' or use a number")

    asm_sizes = [fastas_to_extract[sample]["assembly_size"] for sample in fastas_to_extract]
    i = int(round(len(asm_sizes) * 2 / 3)) - 1
    asm_size = math.ceil(sorted(asm_sizes)[i] / 512**3) * 512**3 # Round up to the 0.5GB
    if ref_type == "protein":
        min_ram_b = asm_size * settings.BLAT_PROT_RAM_FACTOR
    elif ref_type == "dna":
        min_ram_b = asm_size * settings.BLAT_DNA_RAM_FACTOR
    min_threads = settings.EXTRACT_MIN_THREADS

    if min_ram_b >= ram_B or min_threads >= threads_max:
        return 1, 1, ram_B

    threads_max = min(ram_B // min_ram_b, threads_max)

    if concurrent > 1:
        while threads_max // concurrent < min_threads:
            concurrent -= 1
    ram_B_per_extraction = ram_B // concurrent

    if ram_B_per_extraction < min_ram_b:
        while ram_B // concurrent < min_ram_b:
            concurrent -= 1

    threads_per_extraction = threads_max // concurrent
    ram_B_per_extraction = ram_B // concurrent
    return concurrent, threads_per_extraction, ram_B_per_extraction


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
    skipped = []
    if captus_assemblies_dir.exists():
        sample_dirs = list(Path(captus_assemblies_dir).resolve().rglob("*__captus-asm"))
        for sample_dir in sample_dirs:
            if settings.SEQ_NAME_SEP in f"{sample_dir}".replace("__captus-asm", ""):
                sample_name = sample_dir.parts[-1].replace("__captus-asm", "")
                skipped.append(
                    f"'{sample_dir.parts[-1]}': SKIPPED, pattern"
                    f" '{settings.SEQ_NAME_SEP}' not allowed in sample name"
                    f" '{sample_name}'"
                )
            else:
                fastas = list(Path(sample_dir.resolve()).rglob("*/assembly.fasta"))
                for fasta in fastas:
                    if f"{fasta.parent.parent}".endswith("__captus-asm"):
                        sample_name = fasta.parent.parent.parts[-1].replace("__captus-asm", "")
                        sample_dir = Path(out_dir, f"{sample_name}__captus-ext")
                        fastas_to_extract[sample_name] = {
                            "assembly_path": fasta,
                            "assembly_size": fasta.stat().st_size,
                            "sample_dir": sample_dir,
                        }
    return fastas_to_extract, skipped


def find_and_copy_fastas(fastas, captus_assemblies_dir, overwrite, threads_max, show_less):
    """
    Receives a list of files or a folder name. Only FASTA nucleotide files are accepted. They can
    have the following extensions: .fa, .fna, .fasta, .fa.gz, .fna.gz, .fasta.gz
    Sample names are derived from the filenames by removing the extensions.
    Search within a provided directory is recursive
    """
    valid_exts = settings.FASTA_VALID_EXTENSIONS
    if not isinstance(fastas, list):
        fastas = [fastas]
    if len(fastas) == 1 and Path(fastas[0]).is_dir():
        fastas = [
            file for file in Path(fastas[0]).resolve().rglob("*") if has_valid_ext(file, valid_exts)
        ]
    else:
        fastas = [
            Path(any_file).resolve()
            for any_file in fastas
            if Path(any_file).is_file() and has_valid_ext(any_file, valid_exts)
        ]
    check_and_copy_found_fasta_params = []
    for fasta in fastas:
        check_and_copy_found_fasta_params.append((fasta, valid_exts, captus_assemblies_dir, overwrite))
    tqdm_parallel_async_run(
        check_and_copy_found_fasta,
        check_and_copy_found_fasta_params,
        "Verifying and importing provided FASTA files",
        "Verification and copy completed",
        "assembly",
        threads_max,
        show_less,
    )


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
            fasta_out, _ = fasta_headers_to_spades(fasta_to_dict(fasta_path))
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
                f"{Path(f'{sample_name}__captus-asm', '01_assembly', 'assembly.fasta')}"
                f"' [{elapsed_time(time.time() - start)}]"
            )
        else:
            message = dim(f"'{fasta_name}': SKIPPED, this FASTA contains aminoacids")
    else:
        message = (
            f"'{fasta_name}': SKIPPED, '[captus_assemblies_dir]/"
            f"{Path(f'{sample_name}__captus-asm', '01_assembly', 'assembly.fasta')}"
            "' already exists"
        )
    return message


def prepare_protein_refs(
    nuc_refset, ptd_refset, mit_refset, nuc_transtable, ptd_transtable, mit_transtable, out_dir
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
        elif f"{refset}".lower() in settings.PROT_REF_TARGETS[marker]:
            aa_path = settings.PROT_REF_TARGETS[marker][f"{refset}".lower()]["AA"]
            nt_path = settings.PROT_REF_TARGETS[marker][f"{refset}".lower()]["NT"]
            aa_msg = f"{bold(refset)} {dim(aa_path)}"
            nt_msg = f"{bold(refset)} {dim(nt_path)}"
        elif Path(refset).is_file() and fasta_type(refset) == "NT":
            start = time.time()
            suffix = settings.TRANSLATED_REF_SUFFIX
            if refset.endswith(".gz"):
                aa_path = Path(
                    Path(refset).resolve().parent, f"{Path(refset.replace('.gz', '')).stem}{suffix}"
                )
            else:
                aa_path = Path(Path(refset).resolve().parent, f"{Path(refset).stem}{suffix}")
            nt_path = Path(refset).resolve()
            amino_refset = translate_fasta_dict(fasta_to_dict(refset), transtable)
            amino_refset_fixed = fix_premature_stops(amino_refset)
            if amino_refset_fixed is None:
                dict_to_fasta(amino_refset, aa_path)
                log.log(
                    f"Translated '{bold(refset)}' using Genetic Code: {bold(transtable)}"
                    f" [{elapsed_time(time.time() - start)}]"
                )
            else:
                dict_to_fasta(amino_refset_fixed, aa_path)
                log.log(
                    f"Translated '{bold(refset)}' using Genetic Code: {bold(transtable)}"
                    f" [{elapsed_time(time.time() - start)}]"
                )
                log.log(
                    f"WARNING: Premature stop codons were found in {refset}"
                    " and automatically converted to X"
                )
            aa_msg = f"{bold(aa_path)} {dim('(Translated from')} {dim(refset)}{dim(')')}"
            nt_msg = f"{nt_path}"
        elif Path(refset).is_file() and fasta_type(refset) == "AA":
            amino_refset = fasta_to_dict(Path(refset).resolve())
            amino_refset_fixed = fix_premature_stops(amino_refset)
            if amino_refset_fixed is None:
                aa_path = Path(Path(refset).resolve().parent, f"{Path(refset).stem}.faa")
                dict_to_fasta(amino_refset, aa_path)
            else:
                aa_path = Path(Path(refset).resolve().parent, f"{Path(refset).stem}_fixed.faa")
                dict_to_fasta(amino_refset_fixed, aa_path)
                log.log(
                    f"WARNING: {refset} contained gaps that were removed"
                    " and/or premature stops that were converted to X"
                )
            aa_msg = bold(aa_path)
        elif Path(refset).is_file() and "odb1" in refset and refset.endswith(".tar.gz"):
            log.log(
                f"'{Path(refset).name}' seems to be a BUSCO lineage"
                " database, Captus will attempt to import it..."
            )
            amino_refset = import_busco_odb1x(Path(refset))
            if amino_refset is None:
                aa_msg = red("BUSCO lineage database could not be imported")
            else:
                refset_dir = Path(out_dir, settings.REF_TARGETS_DIR)
                refset_stem = Path(refset).name.replace("".join(Path(refset).suffixes), "")
                make_output_dir(refset_dir)
                amino_refset_fixed = fix_premature_stops(amino_refset)
                if amino_refset_fixed is None:
                    aa_path = Path(refset_dir, f"{refset_stem}.faa")
                    dict_to_fasta(amino_refset, aa_path)
                else:
                    aa_path = Path(refset_dir, f"{refset_stem}_fixed.faa")
                    dict_to_fasta(amino_refset_fixed, aa_path)
                    log.log(
                        f"WARNING: {refset_stem} contained gaps that were removed"
                        " and/or premature stops that were converted to X"
                    )
                log.log(f"'{refset_stem}' was successfully imported")
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
    elif f"{dna_refs}".lower() in settings.DNA_REF_TARGETS:
        nt_path = settings.DNA_REF_TARGETS[dna_refs.lower()].resolve()
        nt_msg = f"{bold(dna_refs)} {dim(nt_path)}"
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
            "NUC": {"AA_path": None, "AA_msg": "not used", "NT_path": None, "NT_msg": "not used"},
            "PTD": {"AA_path": None, "AA_msg": "not used", "NT_path": None, "NT_msg": "not used"},
            "MIT": {"AA_path": None, "AA_msg": "not used", "NT_path": None, "NT_msg": "not used"},
            "DNA": {"AA_path": None, "AA_msg": None, "NT_path": None, "NT_msg": "not used"},
            "CLR": {"AA_path": None, "AA_msg": None, "NT_path": None, "NT_msg": "not used"},
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
        return min(settings.SCIPIO_MAX_IDENT_DIV_CODE, min_identity)
    else:
        return min_identity


def adjust_min_coverage(min_coverage):
    if 1 >= min_coverage > 0:
        return round(min_coverage * 100, 2)
    else:
        return min_coverage


def split_refs(query_dict, out_dir, marker_type, threads, final_round=False):
    # Split reference file in groups of roughly REFS_SPLIT_CHUNK_SIZE to run Scipio in parallel
    # on each

    ref_split_chunk_size = settings.REF_SPLIT_CHUNK_SIZE
    ref_split_dir = Path(out_dir, settings.REF_TARGETS_SPLIT_DIR, settings.MARKER_DIRS[marker_type])
    if final_round is True:
        ref_split_dir = Path(out_dir)
        ref_split_chunk_size = ref_split_chunk_size // settings.BEST_N_TARGETS_INITIAL
    seq_names = list(query_dict)
    random.Random(settings.RANDOM_SEED).shuffle(seq_names)
    num_seqs = len(seq_names)
    if num_seqs < threads:
        seqs_per_chunk = num_seqs
    else:
        cc = math.ceil(num_seqs / ref_split_chunk_size)
        cf = math.floor(num_seqs / ref_split_chunk_size)
        if cc % threads != 0:
            cc += threads - (cc % threads)
        if cf % threads != 0:
            cf += threads - (cf % threads)
        chunks_ceiling = max(threads, cc)
        chunks_floor = max(threads, cf)
        if abs(ref_split_chunk_size - (num_seqs / chunks_ceiling)) <= abs(
            ref_split_chunk_size - (num_seqs / chunks_floor)
        ):
            seqs_per_chunk = math.ceil(num_seqs / chunks_ceiling)
        else:
            seqs_per_chunk = math.ceil(num_seqs / chunks_floor)
    make_output_dir(ref_split_dir)
    ref_seqs_paths = {}
    part = 1
    for i in range(0, len(seq_names), seqs_per_chunk):
        split_fasta = {}
        for seq_name in seq_names[i : i + seqs_per_chunk]:
            split_fasta[seq_name] = query_dict[seq_name]
        split_fasta_path = Path(ref_split_dir, f"{marker_type}_part{part}.fasta")
        dict_to_fasta(split_fasta, split_fasta_path)
        ref_seqs_paths[split_fasta_path] = set(split_fasta)
        part += 1
    return ref_seqs_paths


def reference_info(query_dict):
    num_seqs = len(query_dict)
    total_size = 0
    separators_found = 0

    for seq_name in query_dict:
        total_size += len(query_dict[seq_name]["sequence"])
        if settings.REF_CLUSTER_SEP in seq_name:
            if len(list(filter(None, seq_name.split(settings.REF_CLUSTER_SEP)))) > 1:
                separators_found += 1

    loci_lengths = {}
    if separators_found == num_seqs:
        for seq_name in query_dict:
            locus_name = seq_name.split(settings.REF_CLUSTER_SEP)[-1]
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
    info_msg = bold(f"{num_loci:,} loci, {num_seqs:,} sequences ")
    if separators_found == num_seqs:
        if num_loci == num_seqs:
            info_msg += dim("(loci names found, detected a single sequence per locus)")
        elif num_loci < num_seqs:
            info_msg += dim("(loci names found, detected multiple sequences per locus)")
    else:
        info_msg += dim(
            f'(locus name separator "{settings.REF_CLUSTER_SEP}" missing in '
            f"{num_seqs - separators_found} sequences, each sequence taken as a different locus)"
        )

    ref_info = {
        "num_seqs": num_seqs,
        "total_size": total_size,
        "separators_found": separators_found,
        "num_loci": num_loci,
        "total_length_loci": total_length_loci,
        "info_msg": info_msg,
    }

    return ref_info


def depth_tolerance_check(depth_tolerance, ignore_depth):
    depth_tolerance = max(1, depth_tolerance)
    if ignore_depth is True:
        msg = dim("Filter disabled, '--ignore_depth' is being used")
    else:
        msg = bold(f"{depth_tolerance} ")
        msg += dim("(min = 10^(log(depth of the contig with best hit in locus) / ")
        msg += bold(f"{depth_tolerance}")
        msg += dim("))")
    return depth_tolerance, msg


def scipio_coding(
    scipio_path,
    min_score,
    min_identity,
    min_coverage,
    blat_path,
    target_path,
    sample_dir,
    sample_name,
    query_path,
    query_dict,
    query_parts_paths,
    query_info,
    marker_type,
    transtable,
    ignore_depth,
    depth_tolerance,
    disable_stitching,
    max_loci_files,
    max_loci_scipio_x2,
    max_paralogs,
    paralog_tolerance,
    predict,
    threads,
    ram_bytes,
    tsv_comment,
    debug,
    overwrite,
    keep_all,
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

    # Set programs' paths in case of using the bundled versions
    if scipio_path == "bundled":
        scipio_path = settings.BUNDLED_SCIPIO
    if blat_path == "bundled":
        blat_path = settings.BUNDLED_BLAT

    # Group Scipio's basic parameters, the 'sample_dir', and the 'marker_type' into a dictionary
    scipio_params = {
        "scipio_path": scipio_path,
        "min_score": min_score,
        "min_identity": min_identity,
        "min_coverage": min_coverage,
        "blat_path": blat_path,
        "transtable": transtable,
        "sample_dir": sample_dir,
        "marker_type": marker_type,
    }

    # Run Scipio twice if 'query_info["num_loci"]' does not exceed 'max_loci_scipio_x2'
    if query_info["num_loci"] <= max_loci_scipio_x2:
        # The function parallel_scipio() splits the PSL file according to groups of targets and runs
        # Scipio simultaneouslo on each part
        psl_initial_file = parallel_scipio(
            scipio_params,
            target_path,
            query_parts_paths,
            query_info,
            ignore_depth,
            depth_tolerance,
            disable_stitching,
            overwrite,
            threads,
            ram_bytes,
            debug,
            keep_all,
            stage="initial",
        )
        if psl_initial_file is None:
            message = dim(
                f"'{sample_name}': extraction of {genes[marker_type]}"
                " SKIPPED (output files already exist)"
            )
            return message
        elif psl_initial_file == "0 BLAT HITS":
            message = red(f"'{sample_name}': FAILED extraction of {genes[marker_type]} (0 BLAT hits)")
            return message
        elif psl_initial_file == "0 BLAT HITS ACCEPTED":
            message = red(
                f"'{sample_name}': FAILED extraction of {genes[marker_type]} (0 BLAT hits accepted)"
            )
            return message
        else:
            final_target, final_query = filter_query_and_target(
                query_dict,
                fasta_to_dict(target_path),
                psl_initial_file,
                marker_type,
            )

        # Parse YAML to subselect only the best proteins from the 'query' and the contigs with hits
        # from the assembly ('target')
        if final_target is None or final_query is None:
            message = red(f"'{sample_name}': FAILED extraction of {genes[marker_type]}")
            return message
        else:
            final_query_parts_paths = split_refs(
                fasta_to_dict(final_query),
                Path(final_query.parent),
                marker_type,
                threads,
                final_round=True,
            )

        # Perform final Scipio's run (more exhaustive but with fewer contigs and reference proteins)
        yaml_final_file = parallel_scipio(
            scipio_params,
            final_target,
            final_query_parts_paths,
            query_info,
            ignore_depth,
            depth_tolerance,
            disable_stitching,
            overwrite,
            threads,
            ram_bytes,
            debug,
            keep_all,
            stage="final",
        )

        # Remove split final refs
        for path in final_query_parts_paths:
            path.unlink()

    # Run a single Scipio run when 'num_refs' exceeds 'max_loci_scipio_x2'
    else:
        yaml_final_file = parallel_scipio(
            scipio_params,
            target_path,
            query_parts_paths,
            query_info,
            ignore_depth,
            depth_tolerance,
            disable_stitching,
            overwrite,
            threads,
            ram_bytes,
            debug,
            keep_all,
            stage="single",
        )

    if yaml_final_file is None:
        message = dim(
            f"'{sample_name}': extraction of {genes[marker_type]} SKIPPED (output files already exist)"
        )
        return message
    elif yaml_final_file == "0 BLAT HITS":
        message = red(f"'{sample_name}': FAILED extraction of {genes[marker_type]} (0 BLAT hits)")
        return message
    else:
        yaml_final_dir = yaml_final_file.parent
        final_models = scipio_yaml_to_dict(
            yaml_final_file,
            min_score,
            min_identity,
            min_coverage,
            marker_type,
            transtable,
            max_paralogs,
            paralog_tolerance,
            predict,
        )

    # Parse final YAML to produce output FASTA files and reports
    if final_models is None:
        message = red(f"'{sample_name}': FAILED extraction of {genes[marker_type]}")
        return message
    else:
        write_gff3(
            final_models,
            marker_type,
            disable_stitching,
            tsv_comment,
            Path(yaml_final_dir, f"{marker_type}_contigs.gff"),
        )
        recovery_stats = write_fastas_and_report(
            final_models,
            sample_name,
            fasta_to_dict(target_path),
            yaml_final_dir,
            marker_type,
            max_loci_files,
            tsv_comment,
            overwrite,
        )
        message = (
            f"'{sample_name}': recovered {recovery_stats['num_loci']:,} {genes[marker_type]}"
            f" ({recovery_stats['num_loci'] / query_info['num_loci']:.1%} of {query_info['num_loci']:,}),"
            f" {recovery_stats['total_length_best_hits'] / query_info['total_length_loci']:.1%} of"
            f" total reference length, {recovery_stats['num_paralogs']:,} paralogs"
            f" found [{elapsed_time(time.time() - start)}]"
        )
        return message


def parallel_scipio(
    scipio_params: dict,
    target_path,
    query_parts_paths,
    query_info,
    ignore_depth,
    depth_tolerance,
    disable_stitching,
    overwrite,
    threads,
    ram_bytes,
    debug,
    keep_all,
    stage,
):
    # Set output directory and files according to 'sample_dir' and 'stage'
    marker_type = scipio_params["marker_type"]
    scipio_out_dir = Path(scipio_params["sample_dir"], settings.MARKER_DIRS[marker_type])
    if stage == "initial":
        scipio_out_dir = Path(scipio_out_dir, f"00_initial_scipio_{marker_type}")
        blat_prot_out_file = Path(scipio_out_dir, f"{marker_type}_scipio_initial.psl")
        scipio_out_file = Path(scipio_out_dir, f"{marker_type}_scipio_initial.yaml")
        scipio_log_file = Path(scipio_out_dir, f"{marker_type}_scipio_initial.log")

        # Adjust Scipio's 'min_score' for initial round allowing more lenient matching
        scipio_min_score = round(scipio_params["min_score"] * settings.SCIPIO_SCORE_INITIAL_FACTOR, 5)
    else:
        blat_prot_out_file = Path(scipio_out_dir, f"{marker_type}_scipio_final.psl")
        scipio_out_file = Path(scipio_out_dir, f"{marker_type}_scipio_final.yaml")
        scipio_log_file = Path(scipio_out_dir, f"{marker_type}_scipio_final.log")
        scipio_min_score = scipio_params["min_score"]

    # Set minimum BLAT score according to number of references (too many hits created when using
    # a large number of reference proteins), higher score, fewer hits
    if stage == "single":
        blat_min_score = settings.SCIPIO_ROUND1_SETTINGS[marker_type]["blat"]["min_score"]
    else:
        blat_min_score = settings.SCIPIO_ROUND2_SETTINGS[marker_type]["blat"]["min_score"]

    # Additional Scipio settings according to first (or single) or second Scipio round
    if stage == "final":
        extra_settings = settings.SCIPIO_ROUND2_SETTINGS[marker_type]["scipio"]
        tile_size = settings.SCIPIO_ROUND2_SETTINGS[marker_type]["blat"]["tile_size"]
        one_off = settings.SCIPIO_ROUND2_SETTINGS[marker_type]["blat"]["one_off"]
        max_intron = settings.SCIPIO_ROUND2_SETTINGS[marker_type]["blat"]["max_intron"]
    else:
        extra_settings = settings.SCIPIO_ROUND1_SETTINGS[marker_type]["scipio"]
        tile_size = settings.SCIPIO_ROUND1_SETTINGS[marker_type]["blat"]["tile_size"]
        one_off = settings.SCIPIO_ROUND1_SETTINGS[marker_type]["blat"]["one_off"]
        max_intron = settings.SCIPIO_ROUND1_SETTINGS[marker_type]["blat"]["max_intron"]

    scipio_out_dir, _ = make_output_dir(scipio_out_dir)

    if dir_is_empty(scipio_out_dir) is True or overwrite is True:
        # 1. Run BLAT in parallel on the assembly file split in as many parts as threads
        blat_psl = parallel_blat(
            scipio_params["blat_path"],
            "dnax",
            "prot",
            scipio_params["min_identity"],
            target_path,
            query_parts_paths,
            blat_prot_out_file,
            scipio_log_file,
            threads,
            ram_bytes,
            debug,
            tile_size,
            one_off,
            blat_min_score,
            max_intron,
        )
        if blat_psl is None:
            return "0 BLAT HITS"

        # 2. Filter hits according the depth of coverage of the contigs and by cross-loci overlaps
        blat_psl_prefiltered = prefilter_blat_psl(
            marker_type,
            blat_prot_out_file,
            stage,
            bool(query_info["separators_found"]),
            depth_tolerance,
            ignore_depth,
            threads,
            debug,
            keep_all,
        )
        if blat_psl_prefiltered is None:
            return "0 BLAT HITS ACCEPTED"
        if stage == "initial":
            return blat_prot_out_file

        # 3. Split the previous BLAT PSL according to the reference splits
        split_blat_psl_params = []
        for query_part_path in query_parts_paths:
            part = query_part_path.stem
            split_blat_psl_params.append(
                (
                    blat_prot_out_file,
                    query_parts_paths[query_part_path],
                    Path(scipio_out_dir, f"{part}_scipio.psl"),
                )
            )
        psls_to_del = []
        if debug:
            for params in split_blat_psl_params:
                blat_part_psl = split_psl_by_targets(*params)
                psls_to_del.append(blat_part_psl)
        else:
            subexecutor = ProcessPoolExecutor(max_workers=threads)
            futures = [
                subexecutor.submit(split_psl_by_targets, *params) for params in split_blat_psl_params
            ]
            for future in as_completed(futures):
                blat_part_psl = future.result()
                psls_to_del.append(blat_part_psl)
            subexecutor.shutdown()

        # Fourth, run a Scipio instance on each partial file
        blat_only = False
        run_scipio_command_params = []
        for query_part_path in query_parts_paths:
            part = query_part_path.stem
            run_scipio_command_params.append(
                (
                    scipio_params["scipio_path"],
                    Path(scipio_out_dir, f"{part}_scipio.psl"),
                    blat_only,
                    scipio_min_score,
                    scipio_params["min_identity"],
                    scipio_params["min_coverage"],
                    disable_stitching,
                    scipio_params["blat_path"],
                    blat_min_score,
                    scipio_params["transtable"],
                    extra_settings,
                    target_path,
                    query_part_path,
                    Path(scipio_out_dir, f"{part}_scipio.yaml"),
                    Path(scipio_out_dir, f"{part}_scipio.log"),
                )
            )
        yamls_to_cat = []
        if debug:
            for params in run_scipio_command_params:
                scipio_yaml = run_scipio_command(*params)
                if scipio_yaml is not None:
                    yamls_to_cat.append(scipio_yaml)
        else:
            subexecutor = ProcessPoolExecutor(max_workers=threads)
            futures = [
                subexecutor.submit(run_scipio_command, *params) for params in run_scipio_command_params
            ]
            for future in as_completed(futures):
                scipio_yaml = future.result()
                if scipio_yaml is not None:
                    yamls_to_cat.append(scipio_yaml)
            subexecutor.shutdown()
        with open(scipio_out_file, "wt") as yaml_out:
            with open(scipio_log_file, "at") as yaml_log:
                for scipio_yaml in yamls_to_cat:
                    part = scipio_yaml.stem.split("_")[1].replace("part", "")
                    with open(scipio_yaml, "rt") as yaml_in:
                        for line in yaml_in:
                            if line.startswith("  - ID: "):
                                line = line.replace("  - ID: ", f"  - ID: {part}_")
                            yaml_out.write(line)
                    scipio_log = Path(scipio_yaml.parent, f"{scipio_yaml.stem}.log")
                    with open(scipio_log, "rt") as log_in:
                        yaml_log.write(f"\n\n{scipio_log}:\n")
                        yaml_log.writelines(log_in.readlines())

        # Merge partial output files, then delete all partial files
        if not file_is_empty(scipio_out_file):
            for scipio_yaml in yamls_to_cat:
                scipio_log = Path(scipio_yaml.parent, f"{scipio_yaml.stem}.log")
                scipio_yaml.unlink()
                scipio_log.unlink()
            for psl in psls_to_del:
                psl.unlink()

        return scipio_out_file

    else:
        return None


def split_psl_by_targets(psl_path_in, seq_names_set, psl_out_path):
    with open(psl_out_path, "wt") as part_psl_out:
        with open(psl_path_in, "rt") as full_psl_in:
            for line in full_psl_in:
                if line.split()[9] in seq_names_set:
                    part_psl_out.write(line)
    return psl_out_path


def run_scipio_command(
    scipio_path,
    blat_out_file,
    blat_only,
    scipio_min_score,
    min_identity,
    min_coverage,
    disable_stitching,
    blat_path,
    blat_min_score,
    transtable,
    extra_settings,
    target_path,
    query_path,
    scipio_out_file,
    scipio_log_file,
):
    # Build the basic part of the command
    min_identity_scipio = math.floor(min_identity * settings.SEARCH_IDENT_PROP)
    min_identity_blat = math.floor(min_identity * settings.PROT_BLAT_IDENT_PROP)
    basic = [
        "perl",
        f"{scipio_path}",
        f"--blat_output={blat_out_file}",
        "--verbose",
        "--keep_blat_output",
        "--show_blatline",
        f"--min_score={scipio_min_score}",
        f"--min_identity={min_identity_scipio}",
        f"--min_coverage={min_coverage}",
        "--max_mismatch=0",  # 0 means infinite
        "--multiple_results",
        f"--blat_bin={blat_path}",
        f"--blat_score={blat_min_score}",
        f"--blat_identity={min_identity_blat}",
        f"--transtable={transtable}",
        f"--accepted_intron_penalty={settings.SCIPIO_ACCEPTED_INTRON_PENALTY}",
    ]
    if disable_stitching:
        basic += ["--single_target_hits"]  # Scipio finds matches in a single contig only

    # Finish building the command and run Scipio
    if blat_only:
        scipio_cmd = basic + extra_settings + ["--blat_only"] + [f"{target_path}", f"{query_path}"]
        with open(scipio_log_file, "wt") as yaml_log:
            try:
                subprocess.run(scipio_cmd, stderr=yaml_log)
            except Exception:
                return None
        if file_is_empty(blat_out_file):
            return None
        else:
            return blat_out_file
    else:
        scipio_cmd = basic + extra_settings + [f"{target_path}", f"{query_path}"]
        with open(scipio_out_file, "w") as yaml_out:
            with open(scipio_log_file, "w") as yaml_log:
                try:
                    subprocess.run(scipio_cmd, stdout=yaml_out, stderr=yaml_log)
                except Exception:
                    return None
        if file_is_empty(scipio_out_file):
            return None
        else:
            return scipio_out_file


def filter_query_and_target(query_dict, target_dict, psl_initial_file, marker_type):
    """
    Retain only best model proteins and contigs that were hit after initial Scipio run
    """
    psl_initial_dir = psl_initial_file.parent
    best_proteins, hit_contigs = {}, {}
    with open(psl_initial_file, "rt") as psl_ini:
        for line in psl_ini:
            p = parse_psl_record(line)
            best_proteins[p["q_name"]] = query_dict[p["q_name"]]
            hit_contigs[p["t_name"]] = target_dict[p["t_name"]]
    final_query = Path(psl_initial_dir, f"{marker_type}_best_proteins.faa")
    final_target = Path(psl_initial_dir, f"{marker_type}_hit_contigs.fna")
    dict_to_fasta(best_proteins, final_query, wrap=80)
    if marker_type in ["MIT", "PTD"]:
        filtered_contigs = {}
        for contig in hit_contigs:
            if len(hit_contigs[contig]["sequence"]) <= settings.MAX_ORGANELLE_SIZE[marker_type]:
                filtered_contigs[contig] = hit_contigs[contig]
        dict_to_fasta(filtered_contigs, final_target, wrap=80)
    else:
        dict_to_fasta(hit_contigs, final_target, wrap=80)
    if not best_proteins:
        final_query = None
    if not hit_contigs:
        final_target = None
    return final_target, final_query


def write_fastas_and_report(
    hits, sample_name, target_dict, out_dir, marker_type, max_loci_files, tsv_comment, overwrite
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

    def calc_nx_lx(total_size: int, x: float, sorted_block_sizes: list):
        x_size = total_size * x
        c_size, nx, lx = 0, 0, 0
        for i in range(len(sorted_block_sizes)):
            c_size += sorted_block_sizes[i]
            if c_size >= x_size:
                nx, lx = sorted_block_sizes[i], i + 1
                break
        return nx, lx

    def fragmentation(ref_coords, ref_size):
        block_sizes = []
        for segment in ref_coords.split("\n"):
            block_size = 0
            for coord in segment.split(","):
                start, end = int(coord.split("-")[0]), int(coord.split("-")[1])
                if start == end:
                    block_size += 1
                else:
                    block_size += end - start
            block_sizes.append(block_size)
        block_sizes = sorted(block_sizes, reverse=True)
        rec_size = sum(block_sizes)
        n50, l50 = calc_nx_lx(rec_size, 0.50, block_sizes)
        n90, l90 = calc_nx_lx(rec_size, 0.90, block_sizes)
        ng50, lg50 = calc_nx_lx(ref_size, 0.50, block_sizes)
        ng90, lg90 = calc_nx_lx(ref_size, 0.90, block_sizes)
        stats = {
            "num_contigs": len(block_sizes),
            "N50": n50,
            "N90": n90,
            "NG50": ng50,
            "NG90": ng90,
            "L50": l50,
            "L90": l90,
            "LG50": lg50,
            "LG90": lg90,
        }
        return stats

    num_loci, num_paralogs = 0, 0
    lengths_best_hits, coverages_best_hits = [], []
    flanked_seqs, gene_seqs, cds_aa_seqs, cds_nt_seqs, hit_contigs = {}, {}, {}, {}, {}
    stats_header = "\t".join(settings.EXT_STATS_HEADER)
    stats = []
    for ref in sorted(hits):
        num_loci += 1
        for h in range(len(hits[ref])):
            length = hits[ref][h]["match_len"]

            if h == 0:
                lengths_best_hits.append(length)
                coverages_best_hits.append(hits[ref][h]["coverage"])
            else:
                num_paralogs += 1

            description = (
                f"[hit={h:02}] [wscore={hits[ref][h]['wscore']:.3f}] "
                f"[cover={hits[ref][h]['coverage']:.2f}] [ident={hits[ref][h]['identity']:.2f}] "
                f"[score={hits[ref][h]['score']:.3f}] "
            )

            seq_flanked = hits[ref][h]["seq_flanked"]
            len_flanked = f"[length={len(seq_flanked)}] "
            seq_gene = hits[ref][h]["seq_gene"]
            len_gene = f"[length={len(seq_gene)}] "
            seq_nt = hits[ref][h]["seq_nt"]
            len_nt = f"[length={len(seq_nt)}] "
            seq_aa = hits[ref][h]["seq_aa"]
            len_aa = f"[length={len(seq_aa)}] "

            ref_coords = format_coords(hits[ref][h]["ref_coords"])
            # This query description sometimes becomes too long for alignment software
            # query = (
            #     f'[query={hits[ref][h]["ref_name"]}:{ref_coords}] '
            #     f'[contigs={hits[ref][h]["hit_contigs"]}]'.replace("\n", ";")
            # )
            query = f"[query={hits[ref][h]['ref_name']}] "

            shifts_flanked, shifts_gene, shifts_nt, shifts_aa = "", "", "", ""
            if marker_type in ["NUC", "PTD", "MIT"]:
                shifts_flanked = [str(p + 1) for p in range(len(seq_flanked)) if seq_flanked[p] == "N"]
                if shifts_flanked:
                    shifts_flanked = f"[frameshifts={','.join(shifts_flanked)}] "
                else:
                    shifts_flanked = ""
                shifts_gene = [str(p + 1) for p in range(len(seq_gene)) if seq_gene[p] == "N"]
                if shifts_gene:
                    shifts_gene = f"[frameshifts={','.join(shifts_gene)}] "
                else:
                    shifts_gene = ""
                shifts_nt = [str(p + 1) for p in range(len(seq_nt)) if seq_nt[p] == "N"]
                shifts_aa = [str(math.ceil(int(p) / 3)) for p in shifts_nt]
                if shifts_nt:
                    shifts_aa = f"[frameshifts={','.join(shifts_aa)}]"
                    shifts_nt = f"[frameshifts={','.join(shifts_nt)}]"
                else:
                    shifts_aa, shifts_nt = "", ""

            if len(hits[ref]) == 1:
                seq_name = settings.SEQ_NAME_SEP.join([sample_name, ref])
            else:
                seq_name = settings.SEQ_NAME_SEP.join([sample_name, ref, f"{h:02}"])

            flanked_seqs[seq_name] = {
                "description": f"{query}{description}{len_flanked}{shifts_flanked}",
                "sequence": seq_flanked,
                "ref_name": ref,
            }
            gene_seqs[seq_name] = {
                "description": f"{query}{description}{len_gene}{shifts_gene}",
                "sequence": seq_gene,
                "ref_name": ref,
            }

            if marker_type in ["NUC", "PTD", "MIT"]:
                cds_nt_seqs[seq_name] = {
                    "description": f"{query}{description}{len_nt}{shifts_nt}",
                    "sequence": seq_nt,
                    "ref_name": ref,
                }
                cds_aa_seqs[seq_name] = {
                    "description": f"{query}{description}{len_aa}{shifts_aa}",
                    "sequence": seq_aa,
                    "ref_name": ref,
                }

            # Format data for summary table
            hit_len = len(seq_flanked.replace("n", ""))
            flanks_len = len(seq_flanked.replace("n", "")) - len(seq_gene.replace("n", ""))
            if marker_type in ["NUC", "PTD", "MIT"]:
                intron_len = max(len(seq_gene.replace("n", "")) - len(seq_nt), 0)
                stats_row = {
                    "ref_type": "prot",
                    "cds_len": f"{len(seq_nt)}",
                    "intron_len": f"{intron_len}",
                    "frameshifts": shifts_nt.strip().replace("[frameshifts=", "").strip("]"),
                }
            else:
                stats_row = {
                    "ref_type": "nucl",
                    "cds_len": "NA",
                    "intron_len": "NA",
                    "frameshifts": "NA",
                }
            frag_stats = fragmentation(hits[ref][h]["ref_coords"], hits[ref][h]["ref_size"])
            stats.append(
                "\t".join(
                    [
                        sample_name,
                        marker_type,
                        ref,
                        hits[ref][h]["ref_name"],
                        ref_coords,
                        stats_row["ref_type"],
                        f"{length}",
                        f"{h:02}",
                        f"{hits[ref][h]['coverage']:.2f}",
                        f"{hits[ref][h]['identity']:.2f}",
                        f"{hits[ref][h]['score']:.3f}",
                        f"{hits[ref][h]['wscore']:.3f}",
                        f"{hit_len}",
                        stats_row["cds_len"],
                        stats_row["intron_len"],
                        f"{flanks_len}",
                        stats_row["frameshifts"],
                        f"{frag_stats['num_contigs']}",
                        f"{frag_stats['L50']}",
                        f"{frag_stats['L90']}",
                        f"{frag_stats['LG50']}",
                        f"{frag_stats['LG90']}",
                        f"{hits[ref][h]['hit_contigs']}".replace("\n", ";"),
                        f"{hits[ref][h]['strand']}".replace("\n", ";"),
                        format_coords(hits[ref][h]["hit_coords"]),
                    ]
                )
            )

            for contig in hits[ref][h]["hit_contigs"].split("\n"):
                if contig not in hit_contigs:
                    hit_contigs[contig] = dict(target_dict[contig])

    # Write statistics table
    with open(Path(out_dir, f"{marker_type}_recovery_stats.tsv"), "w") as stats_out:
        stats_out.write(tsv_comment)
        stats_out.write(stats_header + "\n" + "\n".join(stats) + "\n")

    # Write multi-sequence FASTAs and setup directories for locus-wise files, only prepare a
    # separate file per marker if the number of references is not greater than args.max_loci_files
    with open(Path(out_dir, f"{marker_type}_contigs.list.txt"), "w") as hit_list_out:
        for contig in hit_contigs:
            hit_list_out.write(f"{contig}\n")

    if marker_type in ["NUC", "PTD", "MIT"]:
        dict_to_fasta(cds_aa_seqs, Path(out_dir, f"{marker_type}{settings.FORMAT_SUFFIXES['AA']}"))
        dict_to_fasta(cds_nt_seqs, Path(out_dir, f"{marker_type}{settings.FORMAT_SUFFIXES['NT']}"))
        dict_to_fasta(gene_seqs, Path(out_dir, f"{marker_type}{settings.FORMAT_SUFFIXES['GE']}"))
        dict_to_fasta(flanked_seqs, Path(out_dir, f"{marker_type}{settings.FORMAT_SUFFIXES['GF']}"))
    elif marker_type in ["DNA", "CLR"]:
        dict_to_fasta(gene_seqs, Path(out_dir, f"{marker_type}{settings.FORMAT_SUFFIXES['MA']}"))
        dict_to_fasta(flanked_seqs, Path(out_dir, f"{marker_type}{settings.FORMAT_SUFFIXES['MF']}"))

    if len(gene_seqs) <= max_loci_files:
        if marker_type in ["NUC", "PTD", "MIT"]:
            cds_aa_seqs_dir, _ = make_output_dir(Path(out_dir, settings.FORMAT_DIRS["AA"]))
            cds_nt_seqs_dir, _ = make_output_dir(Path(out_dir, settings.FORMAT_DIRS["NT"]))
            gene_seqs_dir, _ = make_output_dir(Path(out_dir, settings.FORMAT_DIRS["GE"]))
            flanked_seqs_dir, _ = make_output_dir(Path(out_dir, settings.FORMAT_DIRS["GF"]))
            if dir_is_empty(cds_aa_seqs_dir) is True or overwrite is True:
                for faa_file in cds_aa_seqs_dir.glob("*.faa"):
                    faa_file.unlink()
            if dir_is_empty(cds_nt_seqs_dir) is True or overwrite is True:
                for fna_file in cds_nt_seqs_dir.glob("*.fna"):
                    fna_file.unlink()
        elif marker_type in ["DNA", "CLR"]:
            gene_seqs_dir, _ = make_output_dir(Path(out_dir, settings.FORMAT_DIRS["MA"]))
            flanked_seqs_dir, _ = make_output_dir(Path(out_dir, settings.FORMAT_DIRS["MF"]))

        if dir_is_empty(gene_seqs_dir) is True or overwrite is True:
            for fna_file in gene_seqs_dir.glob("*.fna"):
                fna_file.unlink()
        if dir_is_empty(flanked_seqs_dir) is True or overwrite is True:
            for fna_file in flanked_seqs_dir.glob("*.fna"):
                fna_file.unlink()

        for seq_name_full in sorted(flanked_seqs):
            if len(seq_name_full.split(settings.SEQ_NAME_SEP)) == 3:
                seq_name_short = settings.SEQ_NAME_SEP.join(
                    seq_name_full.split(settings.SEQ_NAME_SEP)[::2]
                )
            else:
                seq_name_short = seq_name_full.split(settings.SEQ_NAME_SEP)[0]
            dict_to_fasta(
                {seq_name_short: dict(flanked_seqs[seq_name_full])},
                Path(flanked_seqs_dir, f"{flanked_seqs[seq_name_full]['ref_name']}.fna"),
                append=True,
            )
            dict_to_fasta(
                {seq_name_short: dict(gene_seqs[seq_name_full])},
                Path(gene_seqs_dir, f"{gene_seqs[seq_name_full]['ref_name']}.fna"),
                append=True,
            )
            if marker_type in ["NUC", "PTD", "MIT"]:
                dict_to_fasta(
                    {seq_name_short: dict(cds_aa_seqs[seq_name_full])},
                    Path(cds_aa_seqs_dir, f"{cds_aa_seqs[seq_name_full]['ref_name']}.faa"),
                    append=True,
                )
                dict_to_fasta(
                    {seq_name_short: dict(cds_nt_seqs[seq_name_full])},
                    Path(cds_nt_seqs_dir, f"{cds_nt_seqs[seq_name_full]['ref_name']}.fna"),
                    append=True,
                )

    recovery_stats = {
        "num_loci": num_loci,
        "num_paralogs": num_paralogs,
        "mean_length_best_hits": round(statistics.mean(lengths_best_hits), 2),
        "total_length_best_hits": sum(lengths_best_hits),
        "mean_coverage_best_hits": round(statistics.mean(coverages_best_hits), 2),
    }

    return recovery_stats


def blat_misc_dna(
    blat_path,
    min_identity,
    min_coverage,
    target_path,
    sample_dir,
    sample_name,
    query_path,
    query_dict,
    query_parts_paths: dict,
    query_info,
    marker_type,
    ignore_depth,
    depth_tolerance,
    disable_stitching,
    max_loci_files,
    max_paralogs,
    paralog_tolerance,
    tsv_comment,
    threads,
    ram_bytes,
    debug,
    overwrite,
    keep_all,
):
    """
    Extract matches of miscellaneous DNA sequences by comparing the assemblies to a set of
    references, these can be formatted as the proteins references '>sample-locus_name'
    """

    start = time.time()

    # Set BLAT path in case of using the bundled version
    if blat_path == "bundled":
        blat_path = settings.BUNDLED_BLAT

    # Create output directory
    blat_dna_out_dir = Path(sample_dir, settings.MARKER_DIRS[marker_type])
    blat_dna_out_file = Path(blat_dna_out_dir, f"{marker_type}_blat_search.psl")
    blat_dna_log_file = Path(blat_dna_out_dir, f"{marker_type}_blat_search.log")
    dna_gff_file = Path(blat_dna_out_dir, f"{marker_type}_contigs.gff")
    blat_dna_out_dir, _ = make_output_dir(blat_dna_out_dir)

    if dir_is_empty(blat_dna_out_dir) is True or overwrite is True:
        dna_target = fasta_to_dict(target_path)
        # 1. Run BLAT in parallel on the assembly file split in as many parts as threads
        blat_psl = parallel_blat(
            blat_path,
            "dna",
            "dna",
            min_identity,
            target_path,
            query_parts_paths,
            blat_dna_out_file,
            blat_dna_log_file,
            threads,
            ram_bytes,
            debug,
            None,
            None,
            settings.DNA_BLAT_MIN_SCORE,
            None,
        )
        if blat_psl is None:
            message = red(
                f"'{sample_name}': FAILED extraction of miscellaneous DNA markers (0 BLAT hits)"
            )
            return message

        # 2. Filter hits according the depth of coverage of the contigs and by cross-loci overlaps
        blat_psl_prefiltered = prefilter_blat_psl(
            marker_type,
            blat_dna_out_file,
            "single",
            bool(query_info["separators_found"]),
            depth_tolerance,
            ignore_depth,
            threads,
            debug,
            keep_all,
        )
        if blat_psl_prefiltered is None:
            message = red(
                f"'{sample_name}': FAILED extraction of miscellaneous DNA markers (0 BLAT hits accepted)"
            )
            return message

        dna_hits = blat_misc_dna_psl_to_dict(
            blat_dna_out_file,
            dna_target,
            min_identity,
            min_coverage,
            marker_type,
            disable_stitching,
            max_paralogs,
            paralog_tolerance,
        )
        if not dna_hits:
            message = red(
                f"'{sample_name}': FAILED extraction of miscellaneous DNA markers (0 BLAT hits assembled)"
            )
            return message
        else:
            if not keep_all:
                Path(blat_dna_out_file).unlink()
            write_gff3(dna_hits, marker_type, disable_stitching, tsv_comment, dna_gff_file)
            recovery_stats = write_fastas_and_report(
                dna_hits,
                sample_name,
                dna_target,
                blat_dna_out_dir,
                marker_type,
                max_loci_files,
                tsv_comment,
                overwrite,
            )
            message = (
                f"'{sample_name}': recovered {recovery_stats['num_loci']:,} DNA markers"
                f" ({recovery_stats['num_loci'] / query_info['num_loci']:.1%} of"
                f" {query_info['num_loci']:,}),"
                f" {recovery_stats['total_length_best_hits'] / query_info['total_length_loci']:.1%}"
                f" of total reference length, {recovery_stats['num_paralogs']:,} paralogs found"
                f" [{elapsed_time(time.time() - start)}]"
            )
            return message

    else:
        message = dim(
            f"'{sample_name}': extraction of miscellaneous DNA markers"
            " SKIPPED (output files already exist)"
        )
        return message


def parallel_blat(
    blat_path: Path,
    target_type: str,  # options: dna dnax, type of search in assembly
    query_type: str,  # options: dna prot, if target loci are DNA or PROTEIN
    min_identity: float,
    target_path: Path,
    query_parts_paths: dict,  # dictionary where keys are Paths to target file parts (query for BLAT)
    psl_out_path: Path,
    psl_log_path: Path,
    threads: int,
    ram_bytes: int,
    debug: bool,
    tile_size: int,
    one_off: int,
    min_score: int,
    max_intron: int,
    no_head=True,
):
    # 1. Reduce number of threads if there is less than 1GB RAM per thread
    while ram_bytes / threads < (settings.RAM_FRACTION * 1024**3):
        threads -= 1
    threads = max(1, threads)

    # 2. Prepare parameters to run as many instances of BLAT as threads
    blat_params = []
    for query_part_path in query_parts_paths:
        part = query_part_path.stem.split("_part")[-1]
        psl_part_path = Path(psl_out_path.parent, f"{psl_out_path.stem}_part{part}.psl")
        blat_params.append(
            (
                blat_path,
                target_type,
                query_type,
                min_identity,
                target_path,
                query_part_path,
                psl_part_path,
                tile_size,
                one_off,
                min_score,
                no_head,
                max_intron,
            )
        )

    # 3. BLAT assembly parts in parallel
    psl_parts_paths_to_cat = []
    psl_log_parts_paths_to_cat = []
    if debug:
        for params in blat_params:
            psl_part_path, psl_log_part_path = run_blat_command(*params)
            if psl_part_path is not None:
                psl_parts_paths_to_cat.append(psl_part_path)
            if psl_log_part_path is not None:
                psl_log_parts_paths_to_cat.append(psl_log_part_path)
    else:
        subexecutor = ProcessPoolExecutor(max_workers=threads)
        futures = [subexecutor.submit(run_blat_command, *params) for params in blat_params]
        for future in as_completed(futures):
            psl_part_path, psl_log_part_path = future.result()
            if psl_part_path is not None:
                psl_parts_paths_to_cat.append(psl_part_path)
            if psl_log_part_path is not None:
                psl_log_parts_paths_to_cat.append(psl_log_part_path)
        subexecutor.shutdown()

    # 4. Concatenate partial PSLs and the logs into a single file
    with open(psl_out_path, "wb") as psl_log:
        for psl_part_path in psl_parts_paths_to_cat:
            with open(psl_part_path, "rb") as psl_part:
                while True:
                    chunk = psl_part.read(settings.CHUNK_SIZE)
                    if not chunk:
                        break
                    psl_log.write(chunk)
    with open(psl_log_path, "wb") as psl_log:
        for psl_log_part_path in psl_log_parts_paths_to_cat:
            with open(psl_log_part_path, "rb") as psl_part:
                while True:
                    chunk = psl_part.read(settings.CHUNK_SIZE)
                    if not chunk:
                        break
                    psl_log.write(chunk)
    for path in psl_parts_paths_to_cat:
        path.unlink()
    for path in psl_log_parts_paths_to_cat:
        path.unlink()

    # 5. Check if concatenated PSL is not empty
    if psl_out_path.exists():
        if psl_out_path.stat().st_size > 0:
            return psl_out_path
        else:
            psl_out_path.unlink()
            return None
    else:
        return None


def split_asm_by_contigs(
    asm_full: dict,
    contigs_parts: dict,
    part: int,
    asm_part_path: Path,
):
    asm_part = {}
    for contig in asm_full:
        if contigs_parts[contig] == part:
            asm_part[contig] = asm_full[contig]
    dict_to_fasta(asm_part, asm_part_path, write_if_empty=False)

    if asm_part_path.exists():
        return asm_part_path
    else:
        return None


def run_blat_command(
    blat_path: Path,
    target_type: str,
    query_type: str,
    min_identity: float,
    target_path: Path,
    query_part_path: Path,
    psl_part_path: Path,
    tile_size,
    one_off,
    min_score,
    no_head,
    max_intron,
):
    if query_type == "prot":
        min_identity_blat = math.floor(min_identity * settings.PROT_BLAT_IDENT_PROP)
    else:
        min_identity_blat = math.floor(min_identity * settings.SEARCH_IDENT_PROP)
    blat_cmd = [
        f"{blat_path}",
        f"-t={target_type}",
        f"-q={query_type}",
    ]
    if tile_size:
        blat_cmd += [f"-tileSize={tile_size}"]
    if one_off:
        blat_cmd += [f"-oneOff={one_off}"]
    if min_score:
        blat_cmd += [f"-minScore={min_score}"]
    blat_cmd += [f"-minIdentity={min_identity_blat:.0f}"]
    if no_head:
        blat_cmd += ["-noHead"]
    if query_type == "prot":
        if max_intron:
            blat_cmd += [f"-maxIntron={max_intron}"]
    else:
        blat_cmd += [f"-maxIntron={settings.DNA_BLAT_MAX_INSERTION}"]
    blat_cmd += [
        f"{target_path}",
        f"{query_part_path}",
        f"{psl_part_path}",
    ]
    psl_log_part_path = Path(psl_part_path.parent, psl_part_path.name.replace(".psl", ".log"))
    with open(psl_log_part_path, "w") as blat_log:
        blat_log.write(f"Captus' BLAT command:\n  {' '.join(blat_cmd)}\n\n")
    with open(psl_log_part_path, "a") as blat_log:
        subprocess.run(blat_cmd, stdout=blat_log, stderr=blat_log)
        blat_log.write("\n\n")

    if psl_part_path.exists():
        if psl_part_path.stat().st_size > 0:
            return psl_part_path, psl_log_part_path
        else:
            psl_part_path.unlink()
            return None, psl_log_part_path
    else:
        return None, psl_log_part_path


def prefilter_blat_psl(
    marker_type: str,
    blat_out_file: Path,
    stage: str,
    separators_found: bool,  # Reference targets have potentially multiple seqs per locus
    depth_tolerance: float,
    ignore_depth: bool,
    threads: int,
    debug: bool,
    keep_all: bool,
):
    size_mul = settings.SIZE_MUL[marker_type]

    # 1. Determine if contigs contain depth of coverage info in their names
    contigs_have_depth = None
    with open(blat_out_file, "rt") as psl_in:
        for line in psl_in:
            p = parse_psl_record(line)
            if "_cov_" in p["t_name"]:
                contigs_have_depth = True
            else:
                contigs_have_depth = False
            break

    # 2. Filter psl hits by min wscore, split accepted hits in as many files as threads
    psls_to_filter_paths, psl_rejected_path = split_psl_by_contigs_add_wscore(
        blat_out_file, stage, separators_found, threads, debug
    )
    if psls_to_filter_paths is None:
        return None

    # 3. Filter the parts in parallel
    filter_psl_overlaps_params = []
    for path in psls_to_filter_paths:
        filter_psl_overlaps_params.append((path, separators_found, size_mul))
    accepted_psls_to_cat = []
    rejected_psls_to_cat = []
    if debug:
        for params in filter_psl_overlaps_params:
            psl_accepted_part_path, psl_rejected_part_path = filter_psl_overlaps(*params)
            if psl_accepted_part_path is not None:
                accepted_psls_to_cat.append(psl_accepted_part_path)
            if psl_rejected_part_path is not None:
                rejected_psls_to_cat.append(psl_rejected_part_path)
    else:
        subexecutor = ProcessPoolExecutor(max_workers=threads)
        futures = [
            subexecutor.submit(filter_psl_overlaps, *params) for params in filter_psl_overlaps_params
        ]
        for future in as_completed(futures):
            psl_accepted_part_path, psl_rejected_part_path = future.result()
            if psl_accepted_part_path is not None:
                accepted_psls_to_cat.append(psl_accepted_part_path)
            if psl_rejected_part_path is not None:
                rejected_psls_to_cat.append(psl_rejected_part_path)
        subexecutor.shutdown()

    # 4. Append rejected files to initial rejected file
    with open(psl_rejected_path, "ab") as psl_rej:
        for path in rejected_psls_to_cat:
            with open(path, "rb") as psl_part_rej:
                while True:
                    chunk = psl_part_rej.read(settings.CHUNK_SIZE)
                    if not chunk:
                        break
                    psl_rej.write(chunk)
    for path in rejected_psls_to_cat:
        path.unlink()

    # 5. Filter by depth or finish writing the accepted hits and finish
    if ignore_depth is True or contigs_have_depth is False:
        with open(blat_out_file, "wt") as psl_acc:
            for path in accepted_psls_to_cat:
                with open(path, "rt") as psl_part_acc:
                    for line in psl_part_acc:
                        psl_acc.write("\t".join(line.strip().split()[0:21]) + "\n")
        for path in accepted_psls_to_cat:
            path.unlink()
        if keep_all is False:
            psl_rejected_path.unlink()
        return blat_out_file
    else:
        with open(blat_out_file, "wb") as psl_acc:
            for path in accepted_psls_to_cat:
                with open(path, "rb") as psl_part_acc:
                    while True:
                        chunk = psl_part_acc.read(settings.CHUNK_SIZE)
                        if not chunk:
                            break
                        psl_acc.write(chunk)
        for path in accepted_psls_to_cat:
            path.unlink()

        # 6. Load accepted hits and organize them by locus, add wscore and depth
        loci = {}
        with open(blat_out_file, "rt") as psl_in:
            for line in psl_in:
                p = parse_psl_record(line)
                locus = p["q_name"]
                if separators_found is True:
                    locus = p["q_name"].split(settings.REF_CLUSTER_SEP)[-1]
                hit_info = {
                    "psl_record": line.strip(),
                    "locus": locus,
                    "target": p["q_name"],
                    "contig": p["t_name"],
                    "wscore": p["wscore"],
                    "depth": float(p["t_name"].split("_cov_")[1].split("_")[0]),
                }
                if locus not in loci:
                    loci[locus] = [hit_info]
                else:
                    loci[locus].append(hit_info)

        # 7. Get highest wscore per locus across contigs and store the depth of the contig
        loci_depths = {}
        for locus in sorted(loci):
            for i in range(len(loci[locus])):
                contig = loci[locus][i]["contig"]
                wscore = loci[locus][i]["wscore"]
                depth = loci[locus][i]["depth"]
                wscoreld = wscore * math.log10(depth)
                contig_info = {"wscore": wscore, "depth": depth, "max_wscoreld": wscoreld}
                if locus not in loci_depths:
                    loci_depths[locus] = {contig: contig_info}
                else:
                    if contig not in loci_depths[locus]:
                        loci_depths[locus][contig] = contig_info
                    else:
                        if wscoreld > loci_depths[locus][contig]["max_wscoreld"]:
                            loci_depths[locus][contig] = contig_info

        # 8. Determine min_depth per locus based on depth of contig with the best wscore in the locus
        loci_min_depths = {}
        for locus in sorted(loci_depths):
            best_hit_wscoreld = 0
            depth = 0
            for contig in loci_depths[locus]:
                if loci_depths[locus][contig]["max_wscoreld"] > best_hit_wscoreld:
                    best_hit_wscoreld = loci_depths[locus][contig]["max_wscoreld"]
                    depth = loci_depths[locus][contig]["depth"]
                elif loci_depths[locus][contig]["max_wscoreld"] == best_hit_wscoreld:
                    if loci_depths[locus][contig]["depth"] > depth:
                        depth = loci_depths[locus][contig]["depth"]
            try:
                loci_min_depths[locus] = 10 ** (round(math.log10(depth) / depth_tolerance, 2))
            except ValueError:
                loci_min_depths[locus] = 0

        # 9. Write rejected and accepted hits after filtering by depth
        with open(blat_out_file, "wt") as psl_acc:
            with open(psl_rejected_path, "at") as psl_rej:
                for locus in sorted(loci):
                    for i in range(len(loci[locus])):
                        if loci[locus][i]["depth"] >= loci_min_depths[locus]:
                            psl_record = "\t".join(loci[locus][i]["psl_record"].split()[0:21])
                            psl_acc.write(f"{psl_record}\n")
                        else:
                            reason = f"min_locus_depth={loci_min_depths[locus]:.4f}"
                            psl_rej.write(f"{loci[locus][i]['psl_record']};{reason}\n")

        # 10. Write depth contig statistics for accepted hits
        accepted_contigs_depths = Path(
            blat_out_file.parent, blat_out_file.name.replace(".psl", "_depths.tsv")
        )
        depths_header = "\t".join(
            [
                "locus",
                "contig",
                "wscore",
                "depth",
                "max_wscoreld",
            ]
        )
        with open(accepted_contigs_depths, "wt") as dep_out:
            dep_out.write(f"{depths_header}\n")
            for locus in sorted(loci_depths):
                for contig in sorted(
                    loci_depths[locus],
                    key=lambda x: loci_depths[locus][x]["max_wscoreld"],
                    reverse=True,
                ):
                    if loci_depths[locus][contig]["depth"] >= loci_min_depths[locus]:
                        dep_out.write(
                            "\t".join(
                                [
                                    locus,
                                    contig,
                                    f"{loci_depths[locus][contig]['wscore']:.5f}",
                                    f"{loci_depths[locus][contig]['depth']:.4f}",
                                    f"{loci_depths[locus][contig]['max_wscoreld']:.5f}",
                                ]
                            )
                            + "\n"
                        )

    if keep_all is False:
        psl_rejected_path.unlink()
        accepted_contigs_depths.unlink()

    return blat_out_file


def split_psl_by_contigs_add_wscore(
    psl_raw_path: Path, stage: str, separators_found: bool, threads: int, debug: bool
):
    psl_rejected_path = Path(psl_raw_path.parent, psl_raw_path.name.replace(".psl", "_rejected.psl"))
    psl_wscore_path = Path(psl_raw_path.parent, psl_raw_path.name.replace(".psl", "_wscore.psl"))

    # If the number of targets per locus is high as in Mega353 for example, the potential number
    # of cross-loci combinations to check for overlaps grows enormously. To reduce the total universe
    # of combinations we can retain the hits from the targets that achieved the n-best wscores, in
    # this way we still allow Captus to choose among n targets in a locus the best full assembly
    # possible. The problem is not so bad with contigs assembled in Captus, but it gets super slow
    # when trying to test all cross-loci overlaps in a full-size chromosome sequence

    def calc_wscore(psl_record: dict):
        coverage = (
            psl_record["matches"] + psl_record["rep_matches"] + psl_record["mismatches"]
        ) / psl_record["q_size"]
        score = (
            psl_record["matches"] + psl_record["rep_matches"] - psl_record["mismatches"]
        ) / psl_record["q_size"]
        return score * coverage

    # 1. Store 'max_wscore' per target sequence
    target_max_wscores = {}
    with open(psl_raw_path, "rt") as psl_in:
        for line in psl_in:
            p = parse_psl_record(line)
            wscore = calc_wscore(p)
            target = p["q_name"]
            if target not in target_max_wscores:
                target_max_wscores[target] = wscore
            else:
                if wscore > target_max_wscores[target]:
                    target_max_wscores[target] = wscore

    # 1. Filter hits according to 'MIN_PROP_WSCORE', store contig sizes for accepted hits
    contig_sizes = {}
    target_all_wscores = {}
    with open(psl_raw_path, "rt") as psl_in:
        with open(psl_rejected_path, "wt") as psl_rej:
            with open(psl_wscore_path, "wt") as psl_wsc:
                for line in psl_in:
                    p = parse_psl_record(line)
                    wscore = calc_wscore(p)
                    locus = p["q_name"]
                    if separators_found is True:
                        locus = p["q_name"].split(settings.REF_CLUSTER_SEP)[-1]
                    target = p["q_name"]
                    line_out = f"{line.strip()}\t{wscore:.5f}"
                    max_target_wscore = target_max_wscores[target]
                    if wscore >= max(
                        settings.MIN_ABS_WSCORE, max_target_wscore * settings.MIN_PROP_WSCORE
                    ):
                        if locus not in target_all_wscores:
                            target_all_wscores[locus] = {target: wscore}
                        else:
                            if target not in target_all_wscores[locus]:
                                target_all_wscores[locus][target] = wscore
                            else:
                                if wscore > target_all_wscores[locus][target]:
                                    target_all_wscores[locus][target] = wscore
                        contig_sizes[p["t_name"]] = p["t_size"]
                        psl_wsc.write(f"{line_out}\n")
                    else:
                        psl_rej.write(f"{line_out};max_target_wscore={max_target_wscore:.5f}\n")

    # Exit if there were no accepted hits
    if len(contig_sizes) == 0:
        return None, psl_rejected_path

    if separators_found is True:
        num_targets = sum([len(target_all_wscores[locus]) for locus in target_all_wscores])
        potential_combos = (num_targets) * (num_targets - 1)
        if potential_combos > settings.MAX_CROSS_LOCI_COMP or stage == "initial":
            # 2. Sort targets by wscores in reverse order inside each locus in target_wscores
            for locus in target_all_wscores:
                target_all_wscores[locus] = dict(
                    sorted(target_all_wscores[locus].items(), key=lambda item: item[1], reverse=True)
                )
            # 3. Determine best_n_targets
            if stage == "initial" and potential_combos < settings.MAX_CROSS_LOCI_COMP:
                best_n_targets = settings.BEST_N_TARGETS_INITIAL
                if max(contig_sizes.values()) >= settings.MIN_CHROM_SIZE:
                    best_n_targets = max(1, best_n_targets // 2)
            else:
                best_n_targets = settings.BEST_N_TARGETS_GLOBAL
                if max(contig_sizes.values()) >= settings.MIN_CHROM_SIZE:
                    best_n_targets //= 10
                num_loci = len(target_all_wscores)
                for n in range(best_n_targets, 1, -1):
                    loci_target_combos = (num_loci * n) * ((num_loci - 1) * n)
                    if loci_target_combos > settings.MAX_CROSS_LOCI_COMP:
                        best_n_targets -= 1
                    else:
                        break
                best_n_targets = max(1, best_n_targets)
            # 4. Get n_best targets from every locus in target_wscores
            allowed_targets = []
            for locus in target_all_wscores:
                counter = 0
                for target in target_all_wscores[locus]:
                    if counter < best_n_targets:
                        allowed_targets.append(target)
                        counter += 1
                    else:
                        break
            del target_all_wscores

            # 5. Remove hits to targets that are not present in allowed_targets
            unfiltered_psl_wscore_path = Path(
                psl_wscore_path.parent, psl_wscore_path.name.replace(".psl", "_unfiltered.psl")
            )
            psl_wscore_path.rename(unfiltered_psl_wscore_path)
            with open(unfiltered_psl_wscore_path, "rt") as psl_in:
                with open(psl_rejected_path, "at") as psl_rej:
                    with open(psl_wscore_path, "wt") as psl_wsc:
                        for line in psl_in:
                            if line.split()[9] in allowed_targets:
                                psl_wsc.write(line)
                            else:
                                line_out = f"{line.strip()};LWT\n"  # LWT = Low Wscore Target
                                psl_rej.write(line_out)
            unfiltered_psl_wscore_path.unlink()

    # 6. Assign contigs to as many groups as threads, filling the groups with smaller contigs first
    # and increasing size after that
    contigs_lists = {i: [] for i in range(threads)}
    part = 0
    for contig in sorted(contig_sizes, key=contig_sizes.get, reverse=True):
        contigs_lists[part].append(contig)
        part += 1
        if part == threads:
            part = 0
    contigs_parts = {}
    for part in contigs_lists:
        for contig in contigs_lists[part]:
            contigs_parts[contig] = part

    # 7. Prepare parameters to run the split in parallel
    split_psl_wscore_params = []
    for part in contigs_lists:
        psl_wscore_part_path = Path(psl_raw_path.parent, f"{psl_raw_path.stem}_part{part}.psl")
        split_psl_wscore_params.append((psl_wscore_path, contigs_parts, part, psl_wscore_part_path))
    del contigs_lists

    # 8. Split PSL file with wscores
    psls_to_filter_paths = []
    if debug:
        for params in split_psl_wscore_params:
            psl_wscore_part_path = split_psl_by_contigs(*params)
            if psl_wscore_part_path is not None:
                psls_to_filter_paths.append(psl_wscore_part_path)
    else:
        subexecutor = ProcessPoolExecutor(max_workers=threads)
        futures = [
            subexecutor.submit(split_psl_by_contigs, *params) for params in split_psl_wscore_params
        ]
        for future in as_completed(futures):
            psl_wscore_part_path = future.result()
            if psl_wscore_part_path is not None:
                psls_to_filter_paths.append(psl_wscore_part_path)
        subexecutor.shutdown()

    psl_wscore_path.unlink()

    return psls_to_filter_paths, psl_rejected_path


def split_psl_by_contigs(
    psl_path_wscore: Path, contigs_parts: dict, part: int, psl_wscore_part_path: Path
):
    with open(psl_path_wscore, "rt") as full_psl_in:
        with open(psl_wscore_part_path, "wt") as part_psl_out:
            for line in full_psl_in:
                if contigs_parts[line.split()[13]] == part:
                    part_psl_out.write(line)
    if psl_wscore_part_path.exists():
        if psl_wscore_part_path.stat().st_size > 0:
            return psl_wscore_part_path
        else:
            psl_wscore_part_path.unlink()
            return None
    else:
        return None


def filter_psl_overlaps(psl_wscore_part_path: Path, separators_found: bool, size_mul: int):
    def flip_coords(t_starts_minus, t_ends_minus, t_size):
        t_starts_plus = [t_size - x for x in t_ends_minus[::-1]]
        t_ends_plus = [t_size - x for x in t_starts_minus[::-1]]
        return t_starts_plus, t_ends_plus

    psl_wscore_part_rejected_path = Path(
        psl_wscore_part_path.parent, psl_wscore_part_path.name.replace(".psl", "_rejected.psl")
    )
    contigs = {}
    with open(psl_wscore_part_path, "rt") as psl_in:
        for line in psl_in:
            p = parse_psl_record(line)
            locus = p["q_name"]
            if separators_found is True:
                locus = p["q_name"].split(settings.REF_CLUSTER_SEP)[-1]
            contig = p["t_name"]
            t_starts = p["t_starts"]
            t_ends = [
                p["t_starts"][i] + size_mul * p["block_sizes"][i] for i in range(len(p["block_sizes"]))
            ]
            if size_mul == 3 and p["strand"][-1] == "-":
                t_starts, t_ends = flip_coords(t_starts, t_ends, p["t_size"])
            hit_info = {
                "psl_record": line.strip(),
                "locus": locus,
                "target": p["q_name"],
                "contig": contig,
                "strand": p["strand"][-1],
                "t_size": p["t_size"],
                "t_starts": t_starts,
                "t_ends": t_ends,
                "wscore": p["wscore"],
                "hit_length": size_mul * sum(p["block_sizes"]),
                "accepted": None,
                "reason": None,
            }
            if contig not in contigs:
                contigs[contig] = [hit_info]
            else:
                contigs[contig].append(hit_info)

    # Mark lower wscore hits for deletion if cross-loci overlaps are detected
    with open(psl_wscore_part_path, "wt") as psl_acc:
        with open(psl_wscore_part_rejected_path, "wt") as psl_rej:
            for contig in sorted(contigs):
                contigs[contig] = sorted(contigs[contig], key=lambda i: i["wscore"], reverse=True)
                contigs[contig][0]["accepted"] = True
                if len(contigs[contig]) == 1:
                    psl_acc.write(f"{contigs[contig][0]['psl_record']}\n")
                else:
                    for i in range(len(contigs[contig])):
                        for j in range(i + 1, len(contigs[contig])):
                            if (
                                contigs[contig][i]["target"] == contigs[contig][j]["target"]
                                or contigs[contig][i]["locus"] != contigs[contig][j]["locus"]
                            ):
                                if contigs[contig][j]["accepted"] is not False:
                                    pct_overlap, overlap_types = calculate_hit_overlap(
                                        contigs[contig][i], contigs[contig][j]
                                    )
                                    if (
                                        pct_overlap > settings.HIT_MAX_PCT_OVERLAP
                                        and "allowed" not in overlap_types
                                    ):
                                        contigs[contig][j]["accepted"] = False
                                        contigs[contig][j]["reason"] = (
                                            f"pct_overlap={pct_overlap:.2f};overlap_types={overlap_types}"
                                        )
                    for i in range(len(contigs[contig])):
                        if (
                            contigs[contig][i]["accepted"] is None
                            or contigs[contig][i]["accepted"] is True
                        ):
                            psl_acc.write(f"{contigs[contig][i]['psl_record']}\n")
                        else:
                            psl_rej.write(
                                f"{contigs[contig][i]['psl_record']};{contigs[contig][j]['reason']}\n"
                            )
                            # t_ends = ",".join([str(j) for j in contigs[contig][i]["t_ends"]])
                            # printout = (f"{contigs[contig][i]['psl']}\t{t_ends}\t{contigs[contig][i]['locus']}\t"
                            #             f"{contigs[contig][i]['wscore']:.4f}\t{contigs[contig][i]['accepted']}\t")
                            # print(printout)
    psl_accepted = None
    psl_rejected = None
    if psl_wscore_part_path.exists():
        if psl_wscore_part_path.stat().st_size > 0:
            psl_accepted = psl_wscore_part_path
        else:
            psl_accepted.unlink()
    if psl_wscore_part_rejected_path.exists():
        if psl_wscore_part_rejected_path.stat().st_size > 0:
            psl_rejected = psl_wscore_part_rejected_path
        else:
            psl_wscore_part_rejected_path.unlink()
    return psl_accepted, psl_rejected


def calculate_hit_overlap(hit1: dict, hit2: dict):
    """Determine the percentage of overlap between two BLAT hits, dividing the
        number of overlapped bases by the length of the longest hit

    Args:
        hit1 (dict): BLAT hit1
        hit2 (dict): BLAT hit2
        size_mul (int): 1 for DNA, 3 for PROTEIN

    Returns:
        overlapped bp / longest hit * 100, list of overlap types found
    """

    h1s, h1e = hit1["t_starts"], hit1["t_ends"]
    h2s, h2e = hit2["t_starts"], hit2["t_ends"]

    if h1e[-1] <= h2s[0] or h2e[-1] <= h1s[0]:
        return 0.0, "N"

    if sorted([hit1["locus"], hit2["locus"]]) in settings.VALID_OVERLAPS:
        return 0.0, "allowed"

    overlapped_bp = 0
    overlap_types = []
    for i in range(len(h1s)):
        for j in range(len(h2s)):
            s1, e1 = h1s[i], h1e[i]
            s2, e2 = h2s[j], h2e[j]
            hits_coords = f"hit1:{s1}-{e1}, hit2:{s2}-{e2}"
            if s1 >= e2 or s2 >= e1:
                # print(f"{hits_coords}, NO OVERLAP")
                overlap_types.append("N")  # None
                continue
            elif s2 <= s1 and e2 >= e1:
                # print(f"{hits_coords}, OVERLAP hit1 contained in hit2")
                overlapped_bp += e1 - s1
                overlap_types.append("F")  # Full
            elif s1 <= s2 and e1 >= e2:
                # print(f"{hits_coords}, OVERLAP hit2 contained in hit1")
                overlapped_bp += e2 - s2
                overlap_types.append("F")  # Full
            elif s1 < s2 and e1 < e2:
                # print(f"{hits_coords}, OVERLAP hit1 extends to the left")
                overlapped_bp += e1 - s2
                overlap_types.append("P")  # Partial
            elif s2 < s1 and e2 < e1:
                # print(f"{hits_coords}, OVERLAP hit1 extends to the right")
                overlapped_bp += e2 - s1
                overlap_types.append("P")  # Partial
            else:
                print(f"{hits_coords}, UNDEFINED!")
    overlap_types = ",".join(overlap_types)
    pct_overlap = overlapped_bp / max(hit1["hit_length"], hit2["hit_length"]) * 100.0
    return pct_overlap, overlap_types


def cleanup_post_extraction(
    sample_name,
    sample_dir,
    assembly_path,
    tsv_comment,
    keep_all,
    overwrite,
    skip_clustering,
    cluster=False,
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

    if dir_is_empty(annotated_assembly_dir) is True or overwrite is True:
        # Concatenate all GFF from different genomes
        sample_gffs = sorted(list(sample_dir.resolve().rglob("[A-Z]*_contigs.gff")))
        sample_gffs = [gff for gff in sample_gffs if gff.parts[-2] != "06_assembly_annotated"]
        if not cluster and not skip_clustering:
            sample_gffs = [gff for gff in sample_gffs if gff.parts[-2] != settings.MARKER_DIRS["CLR"]]
        gff_lines = 0
        with open(gff_file, "wt") as gff_out:
            gff_out.write("##gff-version 3\n")
            tsv_comment = tsv_comment.split("\n")
            gff_out.write(f"#{tsv_comment[0]}\n#{tsv_comment[1]}\n")
            if sample_gffs:
                for gff in sorted(sample_gffs):
                    with open(gff, "rt") as gff_in:
                        for line in gff_in:
                            if not line.startswith("##"):
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
        stats_header = "\t".join(settings.EXT_STATS_HEADER) + "\n"
        sample_stats = sorted(list(sample_dir.resolve().rglob("[A-Z]*_recovery_stats.tsv")))
        sample_stats = [tsv for tsv in sample_stats if tsv.parts[-2] != "06_assembly_annotated"]
        if not cluster and not skip_clustering:
            sample_stats = [tsv for tsv in sample_stats if tsv.parts[-2] != settings.MARKER_DIRS["CLR"]]
        tsv_lines = 0
        with open(stats_file, "wt") as tsv_out:
            tsv_out.write(f"{tsv_comment[0]}\n{tsv_comment[1]}\n")
            tsv_out.write(stats_header)
            if sample_stats:
                for table in sorted(sample_stats):
                    with open(table, "rt") as stats_in:
                        for line in stats_in:
                            if not line.startswith("#") and line != stats_header:
                                tsv_out.write(line)
                                tsv_lines += 1
        if tsv_lines == 0:
            tsv_out.unlink()

        # Write FASTAs of contigs with and without hits
        hit_contigs = {}
        contig_lists = list(sample_dir.resolve().rglob("[A-Z]*_contigs.list.txt"))
        if not cluster and not skip_clustering:
            contig_lists = [cl for cl in contig_lists if cl.parts[-2] != settings.MARKER_DIRS["CLR"]]
        if contig_lists:
            for cl in contig_lists:
                with open(cl, "rt") as clin:
                    for line in clin:
                        hit_contigs[line.strip("\n")] = {}
        all_contigs = fasta_to_dict(assembly_path)
        leftover_contigs = {}
        for contig_name in all_contigs:
            if contig_name in hit_contigs:
                hit_contigs[contig_name] = all_contigs[contig_name]
            else:
                leftover_contigs[contig_name] = all_contigs[contig_name]
        dict_to_fasta(hit_contigs, hit_contigs_file, wrap=80)
        if cluster:
            dict_to_fasta(leftover_contigs, leftovers_clust_file, wrap=80, write_if_empty=True)
        else:
            dict_to_fasta(leftover_contigs, leftovers_file, wrap=80, write_if_empty=True)

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
            f"'{sample_name}': cleanup and merging of GFFs and stats SKIPPED"
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
        leftovers_file = Path(
            fastas_to_extract[sample]["sample_dir"],
            "06_assembly_annotated",
            "leftover_contigs.fasta.gz",
        )
        if leftovers_file.is_file():
            fastas_to_cluster[sample] = leftovers_file.resolve()
            num_leftovers += 1
        else:
            fastas_to_cluster[sample] = fastas_to_extract[sample]["assembly_path"]
    return fastas_to_cluster, num_leftovers


def rehead_and_concatenate_fastas(
    fastas_to_cluster,
    clustering_dir,
    clustering_input_file,
    clust_max_seq_len,
    threads_max,
    show_less,
):
    """
    Since all the FASTAs coming from all the samples will be joined in a single file for clustering,
    we have to include the sample names in the headers to be able to distinguish them later.
    The sample's name + '__' + original contig name will be used as new headers, so try to avoid
    the use of '__' to name the samples. The descriptions will be lost (MMseqs ignores them afaik).
    """
    start = time.time()
    rehead_params = []
    for sample in fastas_to_cluster:
        rehead_params.append((sample, fastas_to_cluster[sample], clustering_dir))
    tqdm_parallel_async_run(
        rehead_fasta_with_sample_name,
        rehead_params,
        "Preparing FASTA assemblies for clustering",
        "FASTA pre-processing completed",
        "file",
        threads_max,
        show_less,
    )
    log.log("")
    fastas_to_concatenate = list(clustering_dir.glob("*_leftover_contigs.fasta"))
    log.log(bold("Concatenating input for clustering:"))
    tqdm_cols = min(shutil.get_terminal_size().columns, 120)
    with tqdm(total=len(fastas_to_concatenate), ncols=tqdm_cols, unit="file") as pbar:
        for cat_fasta in fastas_to_concatenate:
            full_fasta = fasta_to_dict(cat_fasta)
            if clust_max_seq_len > 0:
                filtered_fasta = {
                    seq_name: full_fasta[seq_name]
                    for seq_name in full_fasta
                    if len(full_fasta[seq_name]["sequence"]) <= clust_max_seq_len
                }
                dict_to_fasta(filtered_fasta, clustering_input_file, append=True)
            else:
                dict_to_fasta(full_fasta, clustering_input_file, append=True)
            cat_fasta.unlink()
            pbar.update()
    log.log(
        bold(
            f" \u2514\u2500\u2192 File '{clustering_input_file.name}'"
            f" prepared in {elapsed_time(time.time() - start)}(s)"
        )
    )


def rehead_fasta_with_sample_name(sample_name, sample_fasta_path, clustering_dir):
    """
    Prepend sample name to sequence header followed by a '__'
    """
    start = time.time()
    sample_fasta = fasta_to_dict(sample_fasta_path)
    reheaded_fasta = {}
    reheaded_fasta_file = Path(clustering_dir, f"{sample_name}_leftover_contigs.fasta")
    for header in sample_fasta:
        reheaded_fasta[f"{sample_name}{settings.SEQ_NAME_SEP}{header}"] = dict(sample_fasta[header])
    dict_to_fasta(reheaded_fasta, reheaded_fasta_file)
    message = (
        f"'{sample_name}': file '{sample_fasta_path.name}'"
        f" reheaded [{elapsed_time(time.time() - start)}]"
    )
    return message


def cluster_and_select_refs(
    num_samples,
    clust_min_samples,
    clust_max_copies,
    clust_rep_min_len,
    mmseqs_path,
    mmseqs_method,
    cluster_mode,
    cluster_sensitivity,
    clustering_input_file,
    clustering_dir,
    min_identity,
    seq_id_mode,
    min_coverage,
    cov_mode,
    clust_tmp_dir,
    mafft_path,
    max_threads,
    ram_mb,
    debug,
):
    log.log("")
    log.log(bold(f"Initial clustering of contigs at {min_identity}% identity:"))
    clust1_prefix = f"cl{min_identity:.2f}_cov{min_coverage:.2f}"
    clust1_message = mmseqs_cluster(
        mmseqs_path,
        mmseqs_method,
        clustering_dir,
        clustering_input_file,
        clust1_prefix,
        clust_tmp_dir,
        cluster_sensitivity,
        min_identity,
        seq_id_mode,
        min_coverage,
        cov_mode,
        cluster_mode,
        max_threads,
        ram_mb,
    )
    log.log(clust1_message)
    log.log("")
    msg_p1 = bold(f"Filtering clusters with fewer than {clust_min_samples} samples,")
    msg_p2 = bold(f" more than {clust_max_copies} copies in average,")
    msg_p3 = bold(f" and with centroids shorter than {clust_rep_min_len} bp:")
    log.log(f"{msg_p1}{msg_p2}{msg_p3}")
    start = time.time()
    clust1_all_seqs_file = Path(clustering_dir, f"{clust1_prefix}_all_seqs.fasta")
    clust1_clusters = split_mmseqs_clusters_file(clust1_all_seqs_file)
    passed = []
    failed = 0
    singletons = 0
    tqdm_cols = min(shutil.get_terminal_size().columns, 120)
    with tqdm(total=len(clust1_clusters), ncols=tqdm_cols, unit="cluster") as pbar:
        for cluster in clust1_clusters:
            samples_in_cluster = len(
                set([cluster[i][1:].split(settings.SEQ_NAME_SEP)[0] for i in range(0, len(cluster), 2)])
            )
            avg_copies_in_cluster = len(cluster) / 2 / samples_in_cluster
            if samples_in_cluster == 1:
                singletons += 1
                pbar.update()
                continue
            elif (
                samples_in_cluster >= clust_min_samples
                and avg_copies_in_cluster <= clust_max_copies
                and len(cluster[1]) >= clust_rep_min_len
            ):
                passed.append(cluster)
                pbar.update()
                continue
            else:
                failed += 1
                pbar.update()
                continue
    num_clusters = len(passed)
    num_digits = len(str(num_clusters))
    num_cluster = 1
    fasta_to_recluster = {}
    cluster_lenghts = {}
    for cluster in passed:
        clr = f"captus{num_cluster:0{num_digits}}"
        for i in range(0, len(cluster), 2):
            if i == 0:
                cluster_lenghts[clr] = len(cluster[1])
            h = cluster[i][1:].split(settings.SEQ_NAME_SEP)
            smp, ctg = h[0], h[1]
            ref_sep = settings.REF_CLUSTER_SEP
            seq_sep = settings.SEQ_NAME_SEP
            seq_name = f"{smp}_C{(i // 2) + 1}{ref_sep}{clr}{seq_sep}{ctg}"
            fasta_to_recluster[seq_name] = {
                "sequence": cluster[i + 1],
                "description": "",
            }
        num_cluster += 1
    clust2_input_fasta = Path(clustering_dir, f"{clust1_prefix}_passed.fasta")
    dict_to_fasta(fasta_to_recluster, clust2_input_fasta)
    Path(clustering_dir, f"{clust1_prefix}_all_seqs.fasta").unlink()
    Path(clustering_dir, f"{clust1_prefix}_rep_seq.fasta").unlink()
    Path(clustering_dir, f"{clust1_prefix}_cluster.tsv").unlink()
    log.log(
        bold(
            f" \u2514\u2500\u2192 Clusters retained: {len(passed)} [{elapsed_time(time.time() - start)}]"
        )
    )
    log.log(
        dim(f"     Filtered out {singletons + failed} clusters, of which {singletons} were singletons")
    )
    log.log(dim(f"     Retained sequences written to {clust2_input_fasta}"))
    log.log("")
    min_id2 = min(min_identity + ((100 - min_identity) / 2), min_identity + 1)
    log.log(bold(f"Reducing passing clusters by re-clustering at {min_id2}% identity:"))
    clust2_prefix = f"cl{min_id2:.2f}_cov{min_coverage:.2f}"
    clust2_message = mmseqs_cluster(
        mmseqs_path,
        mmseqs_method,
        clustering_dir,
        clust2_input_fasta,
        clust2_prefix,
        clust_tmp_dir,
        cluster_sensitivity,
        min_id2,
        seq_id_mode,
        min_coverage,
        cov_mode,
        cluster_mode,
        max_threads,
        ram_mb,
    )
    log.log(clust2_message)
    log.log("")
    log.log(bold("Selecting final cluster representatives:"))
    start = time.time()
    clust2_all_seqs_file = Path(clustering_dir, f"{clust2_prefix}_all_seqs.fasta")
    clust2_clusters = split_mmseqs_clusters_file(clust2_all_seqs_file)
    cluster_refs = {}
    min_samples_in_cluster = max(1, math.ceil(num_samples * settings.CLR_REP_MIN_SAMPLE_PROP))
    tqdm_cols = min(shutil.get_terminal_size().columns, 120)
    with tqdm(total=len(clust2_clusters), ncols=tqdm_cols, unit="cluster") as pbar:
        for cluster in clust2_clusters:
            samples_in_cluster = len(
                set(
                    [
                        "_C".join(cluster[i][1:].split(settings.SEQ_NAME_SEP)[0].split("_C")[:-1])
                        for i in range(0, len(cluster), 2)
                    ]
                )
            )
            h = cluster[0][1:].split(settings.SEQ_NAME_SEP)
            seq_name = h[0]
            clr = seq_name.split(settings.REF_CLUSTER_SEP)[-1]
            size = len(cluster) // 2
            seq = cluster[1]
            desc = f"[cluster_size={size}] [contig={h[1]}]"
            if (
                len(seq) / cluster_lenghts[clr] >= settings.CLR_REP_MIN_LEN_PROP
                and samples_in_cluster >= min_samples_in_cluster
            ):
                if clr not in cluster_refs:
                    cluster_refs[clr] = {
                        seq_name: {
                            "sequence": seq,
                            "description": desc,
                        }
                    }
                else:
                    cluster_refs[clr][seq_name] = {
                        "sequence": seq,
                        "description": desc,
                    }
            pbar.update()
    log.log("")
    cluster_refs_fasta = Path(clustering_dir, f"{clust1_prefix}_captus_cluster_refs.fasta")
    check_strand_and_save_refs(cluster_refs, cluster_refs_fasta, mafft_path, max_threads, debug)
    Path(clustering_dir, f"{clust2_prefix}_all_seqs.fasta").unlink()
    Path(clustering_dir, f"{clust2_prefix}_rep_seq.fasta").unlink()
    Path(clustering_dir, f"{clust2_prefix}_cluster.tsv").unlink()
    shutil.rmtree(clust_tmp_dir, ignore_errors=True)
    msg_p1 = bold(f"Reference saved to {cluster_refs_fasta}")
    msg_p2 = bold(f" [{elapsed_time(time.time() - start)}]")
    log.log(f"{msg_p1}{msg_p2}")
    log.log("")
    compress_list_files([clustering_input_file, clust2_input_fasta], max_threads)
    return cluster_refs_fasta


def check_strand_and_save_refs(
    cluster_refs_unstranded: dict,
    cluster_refs_fasta: Path,
    mafft_path: str,
    max_threads: int,
    debug: bool,
):
    """
    Takes a dictionary where the main keys are cluster names with fasta dictionaries as values and
    aligns each one with MAFFT to verify each sequence in the cluster is in the same strand

    Parameters
    ----------
    cluster_refs_unstranded : dict
        Main dictionary keys are cluster names, values are Captus' fasta dictionaries
    cluster_refs_fasta : Path
        Destination path to save reference FASTA file
    mafft_path : str
        Path to MAFFT executable
    max_threads : int
        Maximum number of threads, used to run as many MAFFT alignments in parallel using 2 threads
        each
    debug : bool
        When True the MAFFT alignments will not be parallelized
    """

    # Create temporary directories for FASTAs to align
    tmp_unaligned_dir = Path(cluster_refs_fasta.parent, "tmp_unaligned")
    tmp_aligned_dir = Path(cluster_refs_fasta.parent, "tmp_aligned")
    for tmp_dir in [tmp_unaligned_dir, tmp_aligned_dir]:
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir, ignore_errors=True)
        tmp_dir.mkdir(parents=True)

    # Output FASTA dictionary with correctly stranded sequences
    cluster_refs = {}

    # Align non-singleton clusters with MAFFT to fix strands in cluster representatives
    mafft_params = []
    for clr in sorted(cluster_refs_unstranded):
        if len(cluster_refs_unstranded[clr]) == 1:
            for seq_name in cluster_refs_unstranded[clr]:
                cluster_refs[seq_name] = cluster_refs_unstranded[clr][seq_name]
        else:
            fasta_in = Path(tmp_unaligned_dir, f"{clr}.fasta")
            fasta_out = Path(tmp_aligned_dir, f"{clr}.fasta")
            dict_to_fasta(cluster_refs_unstranded[clr], fasta_in)
            mafft_params.append([mafft_path, fasta_in, fasta_out])
    show_less = True
    concurrent = max_threads // 2
    if debug:
        tqdm_serial_run(
            mafft_auto_strand,
            mafft_params,
            "Verifying cluster strands with MAFFT",
            "Strand verification completed",
            "cluster",
            show_less,
        )
    else:
        tqdm_parallel_async_run(
            mafft_auto_strand,
            mafft_params,
            "Verifying cluster strands with MAFFT",
            "Strand verification completed",
            "cluster",
            concurrent,
            show_less,
        )
    log.log("")

    # Load aligned and correctly stranded FASTAs to add to the reference, disalign first
    stranded_fastas = sorted(list(tmp_aligned_dir.resolve().glob("*.fasta")))
    for stranded_fasta_path in stranded_fastas:
        stranded_fasta = fasta_to_dict(stranded_fasta_path)
        for seq_name in stranded_fasta:
            cluster_refs[seq_name] = {
                "sequence": stranded_fasta[seq_name]["sequence"].replace("-", ""),
                "description": stranded_fasta[seq_name]["description"],
            }
    dict_to_fasta(cluster_refs, cluster_refs_fasta)

    # Remove temporary directories
    for tmp_dir in [tmp_unaligned_dir, tmp_aligned_dir]:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    return


def mafft_auto_strand(mafft_path: str, fasta_in: Path, fasta_out: Path):
    start = time.time()
    mafft_cmd = [
        mafft_path,
        "--nuc",
        "--maxiterate",
        "1000",
        "--thread",
        "2",
        "--adjustdirection",
        f"{fasta_in}",
    ]
    fasta_out_short = Path(*fasta_out.parts[-3:])
    mafft_log_file = Path(fasta_out.parent, f"{fasta_out.stem}.mafft.log")
    with open(mafft_log_file, "w") as mafft_log:
        mafft_log.write(f"Captus' MAFFT Command:\n  {' '.join(mafft_cmd)} > {fasta_out}\n\n\n")
    with open(fasta_out, "w") as mafft_out:
        with open(mafft_log_file, "a") as mafft_log:
            subprocess.run(mafft_cmd, stdout=mafft_out, stderr=mafft_log)
            if file_is_empty(fasta_out):
                message = red(f"'{fasta_out_short}': FAILED alignment, empty output file")
                fasta_out.unlink()
            else:
                rehead_root_msa(fasta_in, fasta_out, None, remove_R_=True)
                message = f"'{fasta_out_short}': aligned [{elapsed_time(time.time() - start)}]"
    return message


def collect_ext_stats(out_dir, tsv_comment):
    # Resolve each sample extraction subdirectory to follow symlinks correctly
    ext_dirs = sorted(p.resolve() for p in list(Path(out_dir).rglob("*__captus-ext")))
    samples_stats = []
    for ext_dir in ext_dirs:
        samples_stats += sorted(list(ext_dir.rglob("*_recovery_stats.tsv")))
    samples_stats = [tsv for tsv in samples_stats if tsv.parts[-2] == "06_assembly_annotated"]
    if not samples_stats:
        return None
    else:
        stats_file_out = Path(out_dir, "captus-extract_stats.tsv")
        header = "\t".join(settings.EXT_STATS_HEADER) + "\n"
        with open(stats_file_out, "wt") as tsv_out:
            tsv_out.write(tsv_comment)
            tsv_out.write(header)
            for tsv in samples_stats:
                with open(tsv, "rt") as tsv_in:
                    for line in tsv_in:
                        if not line.startswith("#") and line != header:
                            tsv_out.write(line)
        return stats_file_out
