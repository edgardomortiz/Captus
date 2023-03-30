#!/usr/bin/env python3
"""
Copyright 2020-2023 Edgardo M. Ortiz (e.ortiz.v@gmail.com)
https://github.com/edgardomortiz/Captus

This file is part of Captus. Captus is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Captus is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Captus. If
not, see <http://www.gnu.org/licenses/>.
"""

import time
from pathlib import Path

from . import log, settings
from .bioformats import cds_from_gff, dict_to_fasta, fasta_to_dict
from .misc import (bold, dim, elapsed_time, find_and_match_fastas_gffs, format_dep_msg,
                   mafft_path_version, make_output_dir, mmseqs_path_version, python_library_check,
                   quit_with_error, set_ram, set_threads, successful_exit, tqdm_parallel_async_run,
                   tqdm_serial_run, vsearch_path_version)
from .version import __version__


def cluster(full_command, args):

    captus_start = time.time()
    out_dir, out_dir_msg = make_output_dir(args.out)
    log.logger = log.Log(Path(out_dir, "captus-design_cluster.log"), stdout_verbosity_level=1)

    mar = 21  # Margin for aligning parameters and values

    ################################################################################################
    ############################################################################### STARTING SECTION
    log.log_section_header("Starting Captus-design: CLUSTER", single_newline=False)
    log.log_explanation(
        "Welcome to the marker clustering step of Captus-design. In this step, Captus will process"
        " the files provided with the '--markers_to_cluster' option. If a FASTA file is found"
        " together with a matching GFF annotation file, Captus will extract the CDS for clustering,"
        " otherwise every sequence in the FASTA file that is shorter than '--max_seq_len' will be"
        " clustered", extra_empty_lines_after=0
    )
    if len(args.markers_to_cluster) == 1 and Path(args.markers_to_cluster[0]).is_dir():
        intro_msg = (
            "Since you provided a directory name, Captus will look in that location for all the"
            " FASTA files and their matching GFF annotation files whenever posssible. To match"
            " a FASTA file with its GFF annotation file, the names should be identical except for"
            " their extension. Sample names are derived from the text before the extension. Use '_'"
            " to separate parts of the sample name, first part will be taken as Genus and second as"
            " Species, if multiple samples of the same species are used, add a identifier as third"
            " part to distinguish them."
        )
    else:
        intro_msg = (
            "Since you provided a file name or a list of file names, Captus will only process those"
            " FASTA files and their matching GFF annotation files whenever posssible. To match"
            " a FASTA file with its GFF annotation file, the names should be identical except for"
            " their extension. Sample names are derived from the text before the extension. Use '_'"
            " to separate parts of the sample name, first part will be taken as Genus and second as"
            " Species, if multiple samples of the same species are used, add a identifier as third"
            " part to distinguish them."
        )
    log.log_explanation(intro_msg, extra_empty_lines_after=0)

    intro_msg = ""
    no_clust_program = False
    _, mmseqs_version, mmseqs_status = mmseqs_path_version(args.mmseqs_path)
    _, vsearch_version, vsearch_status = vsearch_path_version(args.vsearch_path)
    if args.clust_program.lower() == "mmseqs" and mmseqs_status == "not found":
        if vsearch_status == "OK":
            args.clust_program = "vsearch"
            intro_msg = ("MMseqs2 was not found, but VSEARCH will be used for dereplicating and"
                         " clustering instead.")
        else:
            no_clust_program = True
    elif args.clust_program.lower() == "vsearch" and vsearch_status == "not found":
        if mmseqs_status == "OK":
            args.clust_program = "mmseqs"
            intro_msg = ("VSEARCH was not found, but MMseqs2 will be used for dereplicating and"
                         " clustering instead.")
        else:
            no_clust_program = True

    log.log_explanation(intro_msg, extra_empty_lines_after=0)
    log.log_explanation("For more information, please see https://github.com/edgardomortiz/Captus")

    log.log(f'{"Captus version":>{mar}}: {bold(f"v{__version__}")}')
    log.log(f'{"Command":>{mar}}: {bold(full_command)}')
    _, ram_MB, ram_GB, ram_GB_total = set_ram(args.ram)
    log.log(f'{"Max. RAM":>{mar}}: {bold(f"{ram_GB:.1f}GB")} {dim(f"(out of {ram_GB_total:.1f}GB)")}')
    threads_max, threads_total = set_threads(args.threads)
    log.log(f'{"Max. Threads":>{mar}}: {bold(threads_max)} {dim(f"(out of {threads_total})")}')
    log.log("")

    log.log(f'{"Dependencies":>{mar}}:')
    if args.clust_program.lower() == "mmseqs":
        log.log(format_dep_msg(f'{"MMseqs2":>{mar}}: ', mmseqs_version, mmseqs_status))
        log.log(format_dep_msg(f'{"VSEARCH":>{mar}}: ', "", "not used"))
    if args.clust_program.lower() == "vsearch":
        log.log(format_dep_msg(f'{"MMseqs2":>{mar}}: ', "", "not used"))
        log.log(format_dep_msg(f'{"VSEARCH":>{mar}}: ', vsearch_version, vsearch_status))
    _, mafft_version, mafft_status = mafft_path_version(args.mafft_path)
    log.log(format_dep_msg(f'{"MAFFT":>{mar}}: ', mafft_version, mafft_status))
    log.log("")

    log.log(f'{"Python libraries":>{mar}}:')
    numpy_found, numpy_version, numpy_status = python_library_check("numpy")
    pandas_found, pandas_version, pandas_status = python_library_check("pandas")
    plotly_found, plotly_version, plotly_status = python_library_check("plotly")
    log.log(format_dep_msg(f'{"numpy":>{mar}}: ', numpy_version, numpy_status))
    log.log(format_dep_msg(f'{"pandas":>{mar}}: ', pandas_version, pandas_status))
    log.log(format_dep_msg(f'{"plotly":>{mar}}: ', plotly_version, plotly_status))
    log.log("")

    log.log(f'{"Output directory":>{mar}}: {bold(out_dir)}')
    log.log(f'{"":>{mar}}  {dim(out_dir_msg)}')
    log.log("")

    if no_clust_program:
        quit_with_error(
            "Captus could not find a program for clustering, please provide a valid path with either"
            " '--mmseqs_path' or '--vsearch_path'"
        )
    if mafft_status == "not found":
        quit_with_error(
            "Captus could not find MAFFT, please provide a valid path with '--mafft_path'"
        )


    ################################################################################################
    ########################################################################## MARKER IMPORT SECTION
    log.log_section_header("Processing and import of markers for clustering")
    import_msg = (
        "Now Captus will import the markers to be clustered. It will look for valid FASTA files in"
        " the location provided with '--markers_to_cluster'. Additionally it will try to find"
        " matching GFF annotation files, names must be identical except for the extension. For"
        " example 'Poa_annua.fasta' and 'Poa_annua.gff.gz' will be taken as a valid FASTA file with"
        " a matching GFF file."
    )
    log.log_explanation(import_msg)

    log.log(f'{"Concurrent imports":>{mar}}: {bold(threads_max)}')
    log.log("")
    log.log(f'{"Bait length":>{mar}}: {args.bait_length}')
    log.log("")
    log.log(f'{"Overwrite files":>{mar}}: {bold(args.overwrite)}')
    log.log(f'{"Keep all files":>{mar}}: {bold(args.keep_all)}')
    fastas_to_import = find_and_match_fastas_gffs(args.markers_to_cluster)
    log.log(f'{"FASTAs to import":>{mar}}: {bold(len(fastas_to_import))}')
    log.log("")

    import_params = []
    for fasta in fastas_to_import:
        import_params.append((
            fasta,
            fastas_to_import[fasta]["fasta_dir"],
            fastas_to_import[fasta]["gff_path"],
            args.bait_length,
            out_dir,
            args.overwrite,
        ))

    if args.debug:
        tqdm_serial_run(import_fasta, import_params,
                        f"Importing and preprocessing FASTAs and GFFs",
                        f"Data import completed", "sample", args.show_less)
    else:
        tqdm_parallel_async_run(import_fasta, import_params,
                                f"Importing and preprocessing FASTAs and GFFs",
                                f"Data import completed", "sample", threads_max, args.show_less)
    log.log("")


    ################################################################################################
    ################################################################################# ENDING SECTION
    successful_exit(
        "Captus-design: CLUSTER -> successfully completed"
        f" [{elapsed_time(time.time() - captus_start)}]"
    )


def import_fasta(fasta_name, fasta_parent, gff_path, bait_length, out_dir, overwrite):

    start = time.time()

    fasta_path = Path(fasta_parent, fasta_name)
    sample_name = fasta_path.name.rstrip("".join(fasta_path.suffixes))
    sample_out_dir = Path(out_dir, f"{sample_name}__captus-clr")

    if overwrite is True or not sample_out_dir.exists():
        sample_out_dir, _ = make_output_dir(Path(out_dir, f"{sample_name}__captus-clr"))
        if gff_path:
            cds, long, short, data = cds_from_gff(gff_path, fasta_path, bait_length)
            cds_path =   Path(sample_out_dir, f'{sample_name}{settings.IMPORT_SUFFIXES["CDS"]}')
            long_path =  Path(sample_out_dir, f'{sample_name}{settings.IMPORT_SUFFIXES["LONG"]}')
            short_path = Path(sample_out_dir, f'{sample_name}{settings.IMPORT_SUFFIXES["SHORT"]}')
            data_path =  Path(sample_out_dir, f'{sample_name}{settings.IMPORT_SUFFIXES["DATA"]}')
            dict_to_fasta(cds, cds_path, write_if_empty=True)
            dict_to_fasta(long, long_path, write_if_empty=True)
            dict_to_fasta(short, short_path, write_if_empty=True)
            with open(data_path, "wt") as data_out:
                for line in data:
                    data_out.write("\t".join(line)+"\n")
            message = f"'{sample_name}': imported CDS and exon data [{elapsed_time(time.time() - start)}]"
        else:
            markers = fasta_to_dict(fasta_path)
            markers_path = Path(sample_out_dir, f'{sample_name}{settings.IMPORT_SUFFIXES["MARKERS"]}')
            dict_to_fasta(markers, markers_path, write_if_empty=True)
            message = f"'{sample_name}': imported FASTA [{elapsed_time(time.time() - start)}]"
    else:
        message = dim(f"'{sample_name}': skipped import (output files already exist)")

    return message
