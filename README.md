# **Captus**
### *Assembly of Phylogenomic Datasets from High-Throughput Sequencing data*
[https://edgardomortiz.github.io/captus.docs/](https://edgardomortiz.github.io/captus.docs/)  

[![Install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/captus/README.html) [![Bioconda downloads](https://anaconda.org/bioconda/captus/badges/downloads.svg)](https://anaconda.org/bioconda/captus) [![Version in Bioconda](https://anaconda.org/bioconda/captus/badges/version.svg)](https://anaconda.org/bioconda/captus) [![Last updated](https://anaconda.org/bioconda/captus/badges/latest_release_date.svg)](https://github.com/edgardomortiz/Captus/releases)
___
## Installation

The simplest way to install `Captus` is to create an isolated software environment using `conda`, if you don't have `conda` we recommend to install `miniconda` from [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html). Once you have `conda` installed in your system you need to configure your channels:
```bash
conda config --prepend channels bioconda
conda config --prepend channels conda-forge
conda config --show channels
```

The last command should show your current channels, the order matters:
```bash
channels:
  - conda-forge
  - bioconda
  - defaults
```

Now we are ready to create a separate environment for Captus:
```bash
conda create -n captus -c bioconda captus
```

`conda` sometimes takes too long to find and configure dependencies, if that happens we recommend installing `mamba` first, and installing `Captus` with it:
```bash
conda install mamba
mamba create -n captus -c bioconda captus
```

Finally, test that `Captus` was correctly installed:
```bash
conda activate captus
captus_assembly
```

And if the program was correctly installed you will see the main help page of Captus:
```text
usage: captus_assembly command [options]

Captus 0.9.90: Assembly of Phylogenomic Datasets from High-Throughput Sequencing data

Captus-assembly commands:
  command     Program commands (in typical order of execution)
                clean = Trim adaptors and quality filter reads with BBTools, run
                        FastQC on the raw and cleaned reads
                assemble = Perform de novo assembly with MEGAHIT: Assembling reads
                           that were cleaned with the 'clean' command is
                           recommended, but reads cleaned elsewhere are also allowed
                extract = Recover targeted markers with BLAT and Scipio: Extracting
                          markers from the assembly obtained with the 'assemble'
                          command is recommended, but any other assemblies in FASTA
                          format are also allowed.
                align = Align extracted markers across samples with MAFFT or MUSCLE:
                        Marker alignment depends on the directory structure created
                        by the 'extract' command. This step also performs paralog
                        filtering and alignment trimming using ClipKIT

Help:
  -h, --help  Show this help message and exit
  --version   Show Captus' version number

For help on a particular command: captus_assembly command -h

ERROR: Missing command
```

## Usage

Documentation and tutorials available at https://edgardomortiz.github.io/captus.docs/
