# Captus
Tools for hybridization capture-based target enrichment: Probe Design, Data Assembly, Marker Extraction and Alignment


## Installation

The simplest way to install Captus is to first create an isolated software environment with the dependencies with `miniconda`:

```bash
conda create -n captus -c conda-forge -c bioconda \
"python>=3.6" jupyter matplotlib numpy pandas tqdm \
"perl-bioperl-core>=1.007002" \
bbmap fastqc mafft mmseqs2 megahit pigz vsearch
```

Then clone the Captus repository:

```bash
git clone https://github.com/edgardomortiz/Captus
```

Now install Captus within the `conda` environment you created:

```bash
$ conda activate captus
$ cd Captus
$ pip install .
```
Run the following command:
```bash
$ captus_assembly
```
And if the program was correctly installed you will see the main help page of Captus:
```
usage: captus_assembly command [options]

Captus' pipeline for targeted sequence capture assembly

Captus-assembly commands:
  command     Pipeline commands (in typical order of execution)
                clean = Trim adaptors and quality filter reads with BBTools
                assemble = Perform de novo assembly with MEGAHIT: Assembling reads that were
                           cleaned with the 'clean' command is recommended, but any other reads
                           are also allowed
                extract = Recover targeted markers with BLAT and Scipio: Extracting markers from
                          the assembly obtained with the 'assemble' command is recommended, but
                          any other assemblies in FASTA format are also allowed.
                align = Align extracted markers across samples with MAFFT: Marker alignment
                        depends on the directory structure created by the 'extract' command

Help:
  -h, --help  Show this help message and exit
  --version   Show Captus' version number

For help on a particular command: captus_assembly command -h

ERROR: Missing command
```
Alternatively, you can skip the installation with `pip` and simply run the following command in your `conda` environment (you will need to add Captus' directory to your `$PATH` to make it available from anywhere in the computer):
```bash
$ captus_assembly-runner.py
```
