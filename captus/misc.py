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


import argparse
import importlib
import multiprocessing
import os
import platform
import re
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from importlib import util
from pathlib import Path

from tqdm import tqdm

from . import log, settings


class ElapsedTimeThread(threading.Thread):
    """
    Stoppable thread that prints the time elapsed, from https://stackoverflow.com/a/44381654
    """
    def __init__(self):
        super(ElapsedTimeThread, self).__init__()
        self._stop_event = threading.Event()

    def stop(self):
        self._stop_event.set()

    def stopped(self):
        return self._stop_event.is_set()

    def run(self):
        thread_start = time.time()
        update_delay = 1 / 24
        while not self.stopped():
            print(f"\rElapsed time: [{time.time() - thread_start:.3f}s]", end="")
            time.sleep(update_delay)


def get_ram():
    """
    Use 'sysctl' in Mac or 'os.sysconf' in Linux to return RAM size in bytes
    """
    os_type = platform.system()
    if os_type == "Darwin":  # a.k.a. Mac
        return int(os.popen("sysctl hw.memsize").readlines()[0].split()[-1])
    elif os_type == "Linux":
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    else:
        return 0


def set_ram(ram):
    """
    Determines RAM size to be used according to Captus' '--ram' argument
    """
    memsize = get_ram()
    if ram == "auto":
        ram_B = int(memsize * settings.RAM_FRACTION)
    else:
        ram_B = min(int(float(ram) * 1024 ** 3), memsize)
    ram_MB = ram_B // 1024 ** 2
    ram_GB = round((ram_B / 1024 ** 3), 1)
    ram_GB_total = round((memsize / 1024 ** 3), 1)
    return (ram_B, ram_MB, ram_GB, ram_GB_total)


def set_threads(threads):
    """
    Parse the string given by '--threads' to return maximum threads to use
    """
    threads_total = os.cpu_count()
    threads_max = threads_total if threads == "auto" else min(int(threads), threads_total)
    return threads_max, threads_total


def tqdm_serial_run(function, params_list, description_msg, finished_msg, unit, show_less=False):
    """
    Run a function in serial mode using a `tqdm` progress bar
    """
    start = time.time()
    log.log(bold(f"{description_msg}:"))
    tqdm_cols = min(shutil.get_terminal_size().columns, 120)
    with tqdm(total=len(params_list), ncols=tqdm_cols, unit=unit) as pbar:
        for params in params_list:
            function_message = function(*params)
            log.log(function_message, print_to_screen=False)
            if not show_less:
                tqdm.write(function_message)
            pbar.update()
    log.log(bold(
        f" \u2514\u2500\u2192 {finished_msg} for {len(params_list)} {unit}(s)"
        f" [{elapsed_time(time.time() - start)}]"
    ))


def tqdm_parallel_async_run(
    function, params_list, description_msg, finished_msg, unit, threads, show_less=False
    ):
    """
    Run a function in parallel asynchronous mode updating a tqdm progress bar
    Keep in mind that the function referred as 'function_name' cannot be nested within another
    """
    def update(function_message):
        log.log(function_message, print_to_screen=False)
        if not show_less:
            tqdm.write(function_message)
        pbar.update()

    start = time.time()
    log.log(bold(f"{description_msg}:"))
    process = multiprocessing.Pool(threads)
    tqdm_cols = min(shutil.get_terminal_size().columns, 120)
    pbar = tqdm(total=len(params_list), ncols=tqdm_cols, unit=unit)
    for i in range(pbar.total):
        process.apply_async(function, params_list[i], callback=update)
    process.close()
    process.join()
    pbar.close()
    log.log(bold(
        f" \u2514\u2500\u2192 {finished_msg} for {len(params_list)} {unit}(s)"
        f" [{elapsed_time(time.time() - start)}]"
    ))


def tqdm_parallel_nested_run(
    function, params_list, description_msg, finished_msg, unit, threads, show_less=False
    ):
    """
    Run a function in parallel allowing children processes to span their own
    parallel processes
    """
    start = time.time()
    log.log(bold(f"{description_msg}:"))
    tqdm_cols = min(shutil.get_terminal_size().columns, 120)
    pbar = tqdm(total=len(params_list), ncols=tqdm_cols, unit=unit)
    executor = ProcessPoolExecutor(max_workers=threads)
    futures = [executor.submit(function, *params) for params in params_list]
    for future in as_completed(futures):
        result = future.result()
        log.log(result, print_to_screen=False)
        if not show_less:
            tqdm.write(result)
        pbar.update()
    executor.shutdown()
    pbar.close()
    log.log(bold(
        f" \u2514\u2500\u2192 {finished_msg} for {len(params_list)} {unit}(s)"
        f" [{elapsed_time(time.time() - start)}]"
    ))


def tqdm_parallel_async_write(
    function, params_list, description_msg, finished_msg, unit, threads, file_out
    ):
    start = time.time()
    log.log(bold(f"{description_msg}:"))
    tqdm_cols = min(shutil.get_terminal_size().columns, 120)
    pbar = tqdm(total=len(params_list), ncols=tqdm_cols, unit=unit)
    executor = ProcessPoolExecutor(max_workers=threads)
    futures = [executor.submit(function, *params) for params in params_list]
    for future in as_completed(futures):
        result = future.result()
        with open(file_out, "at") as fout:
            fout.write("".join(result))
        pbar.update()
    executor.shutdown()
    pbar.close()
    log.log(bold(
        f" \u2514\u2500\u2192 {finished_msg} for {len(params_list)} {unit}(s)"
        f" [{elapsed_time(time.time() - start)}]"
    ))


def elapsed_time(total_seconds):
    """
    Return minutes, hours, or days if task took more than 60 seconds
    """
    if total_seconds <= 60:
        return f"{total_seconds:.3f}s"
    days, seconds = divmod(total_seconds, 86400)
    hours, seconds = divmod(seconds, 3600)
    minutes, seconds = divmod(seconds, 60)
    values = [days, hours, minutes, seconds]
    symbols = ["d", "h", "m", "s"]
    time_str = ""
    for value, symbol in zip(values, symbols):
        if value == 0 and time_str == "":
            continue
        else:
            if symbol == "s":
                time_str += f"{value:.1f}{symbol} ({total_seconds:.3f}s)"
            else:
                time_str += f"{value:.0f}{symbol} "
    return time_str


def make_output_dir(out_dir):
    """
    Creates the output directory, if it doesn't already exist. The directory can be provided as a
    str or as a Path. Returns the created directory as a Path and a status message.
    """
    out_dir = Path(out_dir)
    if not out_dir.exists():
        try:
            out_dir.mkdir(parents=True)
        except OSError:
            quit_with_error(f"Captus was unable to make the output directory {out_dir}")
        message = "Output directory successfully created"
    elif list(Path(out_dir).glob("*")):
        message = "Output directory already exists and files may be overwritten"
    else:  # directory exists but is empty
        message = "Output directory already exists"
    return out_dir.resolve(), message


def make_tmp_dir_within(tmp_dir_in, program_tmp_subdir):
    """
    (Re)creates a temporary directory inside 'tmp_dir_in' with the name in 'program_tmp_dir',
    both arguments can be str or Path. Returns the created directory as a Path.
    """
    tmp_dir_in = Path(f"{tmp_dir_in}".replace("$HOME", "~")).expanduser()
    tmp_dir_out = Path(tmp_dir_in, program_tmp_subdir)
    if not tmp_dir_out.exists():
        try:
            tmp_dir_out.mkdir(parents=True)
        except OSError:
            quit_with_error(f"Captus was unable to make the temporary directory '{tmp_dir_out}'")
    else:
        try:
            shutil.rmtree(tmp_dir_out, ignore_errors=True)
            tmp_dir_out.mkdir(parents=True)
        except OSError:
            quit_with_error(f"Captus was unable recreate the temporary directory '{tmp_dir_out}'")
    return tmp_dir_out.resolve()


def file_is_empty(file_path):
    """
    Check if file has no contents
    """
    file_path = Path(file_path)
    if file_path.stat().st_size == 0:
        return True
    else:
        return False


def dir_is_empty(dir_path, ignore_subdirectories=True):
    """
    Checks if a directory contains files and/or subdirectories, 'dir_path' can be a string or a Path
    """
    dir_path = Path(dir_path)
    if ignore_subdirectories:
        return not bool([f for f in dir_path.glob("*")
                         if not f.name.startswith(".") and not f.is_dir()])
    else:
        return not bool([f for f in dir_path.glob("*")
                         if not f.name.startswith(".")])


def has_valid_ext(file_path, valid_extensions_list):
    """
    Checks if a filename has an extension within a list of valid extension. The argument 'file_path'
    can be a str or a Path
    """
    for ext in valid_extensions_list:
        if f"{file_path}".lower().endswith(ext.lower()):
            return True
    return False


def find_and_match_fastqs(reads, recursive=False):
    """
    Receives a list of files or a drectory name. Only FASTQ files in the list or inside the folder
    are retained. Returns a dictionary with the items formatted as:
    fastqs[file_R1] = {"containing_drectory", "file_R2"} for pairs of files, and
    fastqs[file_R1] = ["containing_drectory", None} for single-end files
    """
    valid_exts = settings.FASTQ_VALID_EXTENSIONS
    if type(reads) is not list:
        reads = [reads]
    if len(reads) == 1 and Path(reads[0]).is_dir():
        if recursive:
            reads = [file for file in Path(reads[0]).resolve().rglob("*")
                    if has_valid_ext(file, valid_exts)]
        else:
            reads = [file for file in Path(reads[0]).resolve().glob("*")
                    if has_valid_ext(file, valid_exts)]
    else:
        reads = [Path(Path(file).parent.resolve(), Path(file).name) for file in reads
                 if Path(file).resolve().is_file()
                 and has_valid_ext(file, valid_exts)
                 and " " not in Path(file).name]
    # Remove hidden files from list
    reads = [Path(file) for file in reads if not f"{file.name}".startswith(".")]
    fastqs = {}
    skipped = []
    for fastq_file in sorted(reads):
        file_name = fastq_file.name
        file_dir = fastq_file.parent
        if "_R1." in file_name or "_R1_" in file_name:
            if settings.SEQ_NAME_SEP in file_name:
                skipped.append(f"'{file_name}': SKIPPED, pattern"
                               f" '{settings.SEQ_NAME_SEP}' not allowed in filenames")
            else:
                if "_R1." in file_name:
                    file_name_r2 = file_name.replace("_R1.", "_R2.")
                elif "_R1_" in file_name:
                    file_name_r2 = file_name.replace("_R1_", "_R2_")
                if Path(file_dir, file_name_r2) in reads:
                    fastqs[file_name] = {"fastq_dir": file_dir, "fastq_r2": file_name_r2}
                else:
                    fastqs[file_name] = {"fastq_dir": file_dir, "fastq_r2": None}
        elif "_R2." not in file_name or "_R2_" not in file_name:
            skipped.append(f"'{file_name}': SKIPPED, pattern '_R1'"
                           f" or '_R2' not found in filename")
    return fastqs, skipped


def find_and_match_fastas_gffs(markers, recursive=False):
    """
    Receives a list of files or a drectory name. Only FASTA and GFF files in the list or inside the
    folder are retained. Returns a dictionary with the items formatted as:
    fastas_to_cluster[fasta_name] = {"containing_drectory", "gff_path"} for FASTA+GFF, and
    fastas_to_cluster[fasta_name] = {"containing_drectory", None} for FASTA only
    """
    valid_fasta_exts = settings.FASTA_VALID_EXTENSIONS
    valid_gff_exts   = settings.GFF_VALID_EXTENSIONS
    fastas, gffs = [], []
    if type(markers) is not list:
        markers = [markers]
    if len(markers) == 1 and Path(markers[0]).is_dir():
        if recursive:
            fastas = [file for file in Path(markers[0]).resolve().rglob("*")
                      if has_valid_ext(file, valid_fasta_exts)]
            gffs   = [file for file in Path(markers[0]).resolve().rglob("*")
                      if has_valid_ext(file, valid_gff_exts)]
        else:
            fastas = [file for file in Path(markers[0]).resolve().glob("*")
                      if has_valid_ext(file, valid_fasta_exts)]
            gffs   = [file for file in Path(markers[0]).resolve().glob("*")
                      if has_valid_ext(file, valid_gff_exts)]
    else:
        fastas = [Path(Path(file).parent.resolve(), Path(file).name) for file in markers
                  if Path(file).resolve().is_file()
                  and has_valid_ext(file, valid_fasta_exts)
                  and " " not in Path(file).name]
        gffs   = [Path(Path(file).parent.resolve(), Path(file).name) for file in markers
                  if Path(file).resolve().is_file()
                  and has_valid_ext(file, valid_gff_exts)
                  and " " not in Path(file).name]
    # Remove hidden files from list
    fastas = [Path(file) for file in fastas if not file.name.startswith(".")]
    gffs   = [Path(file) for file in gffs if not file.name.startswith(".")]
    # Get stem name (sample name) from each GFF with its corresponding path in a dictionary
    gffs   = {file.name.replace("".join(file.suffixes), ""): file for file in gffs}
    # Match if possible a FASTA with a GFF with identical stem
    fastas_to_import = {}
    for fasta_file in fastas:
        fasta_name = fasta_file.name
        fasta_stem = fasta_file.name.replace("".join(fasta_file.suffixes), "")
        fasta_dir  = fasta_file.parent
        if fasta_stem in gffs:
            fastas_to_import[fasta_name] = {"fasta_dir": fasta_dir, "gff_path": gffs[fasta_stem]}
        else:
            fastas_to_import[fasta_name] = {"fasta_dir": fasta_dir, "gff_path": None}
    return fastas_to_import


def quit_with_error(message):
    """
    Displays the given message and ends the program's execution.
    """
    log.log(red(f"\nERROR: {message}\n"), 0, stderr=True)
    sys.exit(os.EX_SOFTWARE)


def successful_exit(message):
    """
    Exit the program showing a message with a successful status for UNIX
    """
    log.log_section_header(message)
    log.log("")
    sys.exit(os.EX_OK)


def gzip_compress(file_in):
    start = time.time()
    cmd = ["gzip", f"{file_in}"]
    if Path(f"{file_in}.gz").exists():
        Path(f"{file_in}.gz").unlink()
    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    file_in = "/".join(f"{Path(file_in)}".split("/")[-3:])
    message = f"'{file_in}': successfully compressed [{elapsed_time(time.time() - start)}]"
    log.log(message, print_to_screen=False)
    return message


def pigz_compress(file_in, threads):
    start = time.time()
    cmd = ["pigz", "-p", f"{threads}", "-5", f"{file_in}"]
    if Path(f"{file_in}.gz").exists():
        Path(f"{file_in}.gz").unlink()
    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    file_in = "/".join(f"{Path(file_in)}".split("/")[-3:])
    message = f"'{file_in}': successfully compressed [{elapsed_time(time.time() - start)}]"
    log.log(message, print_to_screen=False)
    return message


def compress_list_files(files_list, threads):
    compress_params = []
    t = min(threads, settings.MAX_HDD_WRITE_INSTANCES)
    if shutil.which("pigz"):
        for ungz_file in files_list:
            if ungz_file.is_file():
                compress_params.append((ungz_file, threads))
        tqdm_serial_run(pigz_compress, compress_params, "Compressing with 'pigz'",
                        "Completed 'pigz' compression", "file")
    elif shutil.which("gzip"):
        for ungz_file in files_list:
            if ungz_file.is_file():
                compress_params.append((ungz_file, ))
        tqdm_parallel_async_run(gzip_compress, compress_params, "Compressing with 'gzip'",
                                "Completed 'gzip' compression", "file", t)


def execute_jupyter_report(out_dir, jupyter_notebook, title, prefix):
    start = time.time()

    shutil.copy(jupyter_notebook, Path(out_dir, title))
    qc_html_report = Path(out_dir, f"{prefix}.report.html")
    qc_html_report_log = Path(out_dir, f"{prefix}.report.log")
    nbconvert_cmd = [
        "jupyter", "nbconvert",
        "--to", "html",
        "--execute", f'{Path(out_dir, title)}',
        "--no-input",
        "--output", f'{qc_html_report}'
    ]
    with open(qc_html_report_log, "wt") as report_log:
        subprocess.run(nbconvert_cmd, stderr=report_log)

    Path(out_dir, title).unlink()

    if qc_html_report.exists() and qc_html_report.is_file():
        qc_html_msg = dim(f"Report generated in {elapsed_time(time.time() - start)}")
    else:
        qc_html_msg = red(f"Report not generated, verify your Jupyter installation")

    return qc_html_report, qc_html_msg


####################################################################################################
################################################################ FUNCTIONS TO VERIFY SOFTWARE STATUS
def format_dep_msg(dep_text, dep_version, dep_status):
    if dep_status == "not used":
        return f"{dep_text}{dim(dep_status)}"
    elif dep_status == "OK":
        return f'{dep_text}{bold(f"v{dep_version}")} {bold_green(dep_status)}'
    else:
        return f"{dep_text}{bold_red(dep_status)}"


def mafft_path_version(mafft_path):
    # Try "mafft" and "mafft.bat" as defaults
    if mafft_path == "mafft":
        found_mafft_path = shutil.which(mafft_path)
        if found_mafft_path is None:
            found_mafft_path = shutil.which("mafft.bat")
            if found_mafft_path is None:
                return mafft_path, "", "not found"
    else:
        found_mafft_path = shutil.which(mafft_path)
        if found_mafft_path is None:
            return mafft_path, "", "not found"
    command = [found_mafft_path, "--version"]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    version = process.communicate()[0].decode().strip("\n").split()[0][1:]
    return found_mafft_path, version, "OK"

def muscle_path_version(muscle_path):
    # Try "muscle" as default
    if muscle_path == "muscle":
        found_muscle_path = shutil.which(muscle_path)
        if found_muscle_path is None:
            return muscle_path, "", "not found"
    else:
        found_muscle_path = shutil.which(muscle_path)
        if found_muscle_path is None:
            return muscle_path, "", "not found"
    command = [found_muscle_path, "--version"]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    version = process.communicate()[0].decode().strip("\n").split()[1]
    return found_muscle_path, version, "OK"

def falco_path_version(falco_path):
    found_falco_path = shutil.which(falco_path)
    if found_falco_path is None:
        return falco_path, "", "not found"
    command = [found_falco_path, "--version"]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    version = process.communicate()[0].decode().strip("\n").split()[-1]
    return found_falco_path, version, "OK"


def fastqc_path_version(fastqc_path):
    found_fastqc_path = shutil.which(fastqc_path)
    if found_fastqc_path is None:
        return fastqc_path, "", "not found"
    command = [found_fastqc_path, "--version"]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    version = process.communicate()[0].decode().strip("\n").split("v")[-1]
    return found_fastqc_path, version, "OK"


def bbtools_path_version(bbtools_path):
    found_bbtools_path = shutil.which(bbtools_path)
    if found_bbtools_path is None:
        return bbtools_path, "", "not found"
    command = [found_bbtools_path, "-Xmx1g", "version"]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    version = process.communicate()[0].decode().splitlines()[-2].split()[-1]
    return found_bbtools_path, version, "OK"


def megahit_path_version(megahit_path):
    found_megahit_path = shutil.which(megahit_path)
    if found_megahit_path is None:
        return megahit_path, "", "not found"
    command = [found_megahit_path, "--version"]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    version = process.communicate()[0].decode().strip("\n").split("v")[-1]
    return found_megahit_path, version, "OK"


def megahit_tk_path_version(megahit_toolkit_path):
    found_megahit_tk_path = shutil.which(megahit_toolkit_path)
    if found_megahit_tk_path is None:
        return found_megahit_tk_path, "not found"
    command = [found_megahit_tk_path, "dumpversion"]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    version = process.communicate()[0].decode().strip("\n").split("v")[-1]
    return found_megahit_tk_path, version, "OK"


def scipio_path_version(scipio_path):
    if scipio_path == "bundled":
        return "Scipio", "1.4.1", "OK"
    else:
        if Path(scipio_path).is_file():
            with open(scipio_path, "wt") as pl:
                for line in pl:
                    if line.startswith("#   Version:"):
                        version = line.strip("\n").split()[1]
                        break
            return "Scipio", version, "OK"
        else:
            return "Scipio", "", "not found"


def bioperl_get_version():
    command = ["perl", "-MBio::Root::Version", "-e" , "print $Bio::Root::Version::VERSION"]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        version = process.communicate()[0].decode()
        return "BioPerl", str(version), "OK"
    except ValueError:
        return "BioPerl", "", "not found"


def yaml_perl_get_version():
    command = ["perl", "-MYAML", "-e" , "print $YAML::VERSION"]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        version = process.communicate()[0].decode()
        return "YAML", str(version), "OK"
    except ValueError:
        return "YAML", "", "not found"


def blat_path_version(blat_path):
    if blat_path == "bundled":
        blat_path = settings.BUNDLED_BLAT
    found_blat_path = shutil.which(blat_path)
    if found_blat_path is None:
        return blat_path, "", "not found"
    command = [found_blat_path]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    version = process.communicate()[0].decode().splitlines()[0].split()[5]
    return found_blat_path, version, "OK"


def mmseqs_path_version(mmseqs_path):
    found_mmseqs_path = shutil.which(mmseqs_path)
    if found_mmseqs_path is None:
        return mmseqs_path, "", "not found"
    command = [found_mmseqs_path, "version"]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    version = process.communicate()[0].decode().strip("\n")
    return found_mmseqs_path, version, "OK"


def vsearch_path_version(vsearch_path):
    found_vsearch_path = shutil.which(vsearch_path)
    if found_vsearch_path is None:
        return vsearch_path, "", "not found"
    command = [found_vsearch_path, "--version"]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    version = process.communicate()[0].decode().splitlines()[0].split()[1].lstrip("v").rstrip(",")
    return found_vsearch_path, version, "OK"


def clipkit_path_version(clipkit_path):
    found_clipkit_path = shutil.which(clipkit_path)
    if found_clipkit_path is None:
        return clipkit_path, "", "not found"
    command = [found_clipkit_path, "-v"]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    version = process.communicate()[0].decode().strip("\n").split()[-1]
    return found_clipkit_path, version, "OK"


def python_library_check(library_name):
    library_found = bool(util.find_spec(library_name))
    library_version = ""
    library_status = "not found"
    if library_found:
        library = importlib.import_module(library_name)
        library_version = library.__version__
        library_status = "OK"
    return library_found, library_version, library_status


####################################################################################################
## FUNCTIONS TAKEN FROM UNICYCLER FOR HELP AND TEXT FORMATTING (https://github.com/rrwick/Unicycler)

END_FORMATTING = "\033[0m"
BOLD = "\033[1m"
UNDERLINE = "\033[4m"
RED = "\033[31m"
GREEN = "\033[32m"
MAGENTA = "\033[35m"
YELLOW = "\033[93m"
DIM = "\033[2m"

class MyHelpFormatter(argparse.HelpFormatter):
    """
    This is a custom formatter class for argparse. It allows for some custom formatting,
    in particular for the help texts with multiple options (like bridging mode and verbosity level).
    http://stackoverflow.com/questions/3853722
    """
    def __init__(self, prog):
        terminal_width = shutil.get_terminal_size().columns
        os.environ["COLUMNS"] = str(terminal_width)
        max_help_position = min(max(24, terminal_width // 3), 40)
        try:
            self.colours = int(subprocess.check_output(["tput", "colors"]).decode().strip())
        except (ValueError, subprocess.CalledProcessError, FileNotFoundError, AttributeError):
            self.colours = 1
        super().__init__(prog, max_help_position=max_help_position)

    def _get_help_string(self, action):
        """
        Override this function to add default values, but only when 'default' is not already in the
        help text.
        """
        help_text = action.help
        if (action.default != argparse.SUPPRESS and "default" not in help_text.lower()
                and action.default is not None):
            help_text += f" (default: {action.default})"
        return help_text

    def start_section(self, heading):
        """
        Override this method to add bold underlining to section headers.
        """
        if self.colours > 1:
            heading = f"{BOLD}{heading}{END_FORMATTING}"
        super().start_section(heading)

    def _split_lines(self, text, width):
        """
        Override this method to add special behaviour for help texts that start with:
          'B|' - loop text to the column of the equals sign if found, options are indented 2 spaces
        """
        if text.startswith("B|"):
            text_lines = text[2:].splitlines()
            wrapped_text_lines = []
            first_line = True  # use different rules for the first line of help
            for line in text_lines:
                if len(line) <= width:
                    if first_line:
                        wrapped_text_lines.append(line)
                    else:
                        wrapped_text_lines.append(f"  {line}")
                else:
                    line_parts = line.split()
                    if first_line:
                        wrap_column = 0
                    else:
                        wrap_column = 2
                    current_line = f'{" " * wrap_column}{line_parts[0]}'
                    if "=" in line:
                        wrap_column += line.find("=") + 2
                    for part in line_parts[1:]:
                        if len(current_line) + 1 + len(part) <= width:
                            current_line += f" {part}"
                        else:
                            wrapped_text_lines.append(current_line)
                            current_line = f'{" " * wrap_column}{part}'
                    wrapped_text_lines.append(current_line)
                first_line = False
            return wrapped_text_lines
        else:
            return argparse.HelpFormatter._split_lines(self, text, width)

    def _fill_text(self, text, width, indent):
        if text.startswith("R|"):
            return "".join(indent + line for line in text[2:].splitlines(keepends=True))
        else:
            return argparse.HelpFormatter._fill_text(self, text, width, indent)

    def _format_action(self, action):
        """
        Override this method to make help descriptions dim.
        """
        # determine the required width and the entry label
        help_position = min(self._action_max_length + 2, self._max_help_position)
        help_width = self._width - help_position
        action_width = help_position - self._current_indent - 2
        action_header = self._format_action_invocation(action)

        # no help; start on same line and add a final newline
        if not action.help:
            tup = self._current_indent, "", action_header
            action_header = "%*s%s\n" % tup
            indent_first = 0

        # short action name; start on the same line and pad two spaces
        elif len(action_header) <= action_width:
            tup = self._current_indent, "", action_width, action_header
            action_header = "%*s%-*s  " % tup
            indent_first = 0

        # long action name; start on the next line
        else:
            tup = self._current_indent, "", action_header
            action_header = "%*s%s\n" % tup
            indent_first = help_position

        # collect the pieces of the action help
        parts = [action_header]

        # if there was help for the action, add lines of help text
        if action.help:
            help_text = self._expand_help(action)
            help_lines = self._split_lines(help_text, help_width)
            first_line = help_lines[0]
            if self.colours > 8:
                first_line = f"{DIM}{first_line}{END_FORMATTING}"
            parts.append("%*s%s\n" % (indent_first, "", first_line))
            for line in help_lines[1:]:
                if self.colours > 8:
                    line = f"{DIM}{line}{END_FORMATTING}"
                parts.append("%*s%s\n" % (help_position, "", line))

        # or add a newline if the description doesn't end with one
        elif not action_header.endswith("\n"):
            parts.append("\n")

        # if there are any sub-actions, add their help as well
        for subaction in self._iter_indented_subactions(action):
            parts.append(self._format_action(subaction))

        # return a single string
        return self._join_parts(parts)


def green(text):
    return f"{GREEN}{text}{END_FORMATTING}"


def bold_green(text):
    return f"{GREEN}{BOLD}{text}{END_FORMATTING}"


def red(text):
    return f"{RED}{text}{END_FORMATTING}"


def magenta(text):
    return f"{MAGENTA}{text}{END_FORMATTING}"


def bold_red(text):
    return f"{RED}{BOLD}{text}{END_FORMATTING}"


def bold(text):
    return f"{BOLD}{text}{END_FORMATTING}"


def bold_underline(text):
    return f"{BOLD}{UNDERLINE}{text}{END_FORMATTING}"


def underline(text):
    return f"{UNDERLINE}{text}{END_FORMATTING}"


def dim(text):
    return f"{DIM}{text}{END_FORMATTING}"


def dim_underline(text):
    return f"{DIM}{UNDERLINE}{text}{END_FORMATTING}"


def bold_yellow(text):
    return f"{YELLOW}{BOLD}{text}{END_FORMATTING}"


def bold_yellow_underline(text):
    return f"{YELLOW}{BOLD}{UNDERLINE}{text}{END_FORMATTING}"


def bold_red_underline(text):
    return f"{RED}{BOLD}{UNDERLINE}{text}{END_FORMATTING}"


def remove_formatting(text):
    return re.sub('\033.*?m', '', text)
