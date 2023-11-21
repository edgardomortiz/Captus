# **Captus**
### *Assembly of Phylogenomic Datasets from High-Throughput Sequencing data*
[https://edgardomortiz.github.io/captus.docs/](https://edgardomortiz.github.io/captus.docs/)  

[![Install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/captus/README.html) [![Bioconda downloads](https://anaconda.org/bioconda/captus/badges/downloads.svg)](https://anaconda.org/bioconda/captus) [![Version in Bioconda](https://anaconda.org/bioconda/captus/badges/version.svg)](https://anaconda.org/bioconda/captus) [![Last updated](https://anaconda.org/bioconda/captus/badges/latest_release_date.svg)](https://github.com/edgardomortiz/Captus/releases)
___
>[!WARNING]
>### VERY IMPORTANT FOR MAC USERS!
>`Bioconda` has broken the latest build of `MEGAHIT` for macOS. Please follow the instructions below, either **Using micromamba** or **Using conda**, to install `Captus` and **AFTER** `Captus` is installed please run the following two commands to get a functional `MEGAHIT` (if you used `conda` just replace `micromamba` by `conda`):
>```text
>micromamba activate captus
>micromamba install -c bioconda megahit=1.2.9=hfbae3c0_0
>```
---
>[!TIP]
>### USE BUSCO DATABASES AS REFERENCE TARGETS
>Now Captus can parse any BUSCO lineage database and use it as reference targets. Just download one of the `.tar.gz` files from [https://busco-data.ezlab.org/v5/data/lineages/](https://busco-data.ezlab.org/v5/data/lineages/) and provide the path to `captus extract`, for example:
>```text
>captus extract --captus_assemblies_dir my_assemblies --nuc_refs ~/Downloads/aves_odb10.2021-02-19.tar.gz
>```
---
## INSTALLATION
### Using micromamba
The fastest way to install `Captus` is to create an isolated software environment using `micromamba` (https://mamba.readthedocs.io/en/latest/installation.html), if you don't have `micromamba` it can very easily be installed:
For linux with `bash` shell:
```
curl micro.mamba.pm/install.sh | bash
```
For macOS with `zsh` shell:
```
curl micro.mamba.pm/install.sh | zsh
```

. Once you have `micromamba` installed in your system you need to configure your channels:
```bash
micromamba config prepend channels bioconda
micromamba config prepend channels conda-forge
micromamba config list
```

The last command should show your current channels, the order matters:
```bash
channels:
  - conda-forge
  - bioconda
show_banner: false
```

Now we are ready to create a separate environment for Captus:
```bash
micromamba create -n captus captus
```

Finally, test that `Captus` was correctly installed:
```bash
micromamba activate captus
captus_assembly
```

### Using conda
A more "mainstream" but slower way to install `Captus` is to create an isolated software environment using `conda`, if you don't have `conda` we recommend to install `miniconda` from [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html). Once you have `conda` installed in your system you need to configure your channels:
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

Captus 0.9.98: Assembly of Phylogenomic Datasets from High-Throughput Sequencing data

Captus-assembly commands:
  command     Program commands (in typical order of execution)
                clean = Trim adaptors and quality filter reads with BBTools,
                        run FastQC on the raw and cleaned reads
                assemble = Perform de novo assembly with MEGAHIT: Assembling
                           reads that were cleaned with the 'clean' command is
                           recommended, but reads cleaned elsewhere are also
                           allowed
                extract = Recover targeted markers with BLAT and Scipio:
                          Extracting markers from the assembly obtained with
                          the 'assemble' command is recommended, but any other
                          assemblies in FASTA format are also allowed
                align = Align extracted markers across samples with MAFFT or
                        MUSCLE: Marker alignment depends on the directory
                        structure created by the 'extract' command. This step
                        also performs paralog filtering and alignment trimming
                        using ClipKIT

Help:
  -h, --help  Show this help message and exit
  --version   Show Captus' version number

For help on a particular command: captus_assembly command -h

ERROR: Missing command
```
## USAGE

Documentation and tutorials available at https://edgardomortiz.github.io/captus.docs/

## CITATION

Ortiz, E.M., A. Hoewener, G. Shigita, M. Raza, O. Maurin, A. Zuntini, F. Forest, W.J. Baker, H. Schaefer. (2023). _A novel phylogenomics pipeline revels complex pattern of reticulate evolution in Cucurbitales_. bioRxiv [https://doi.org/10.1101/2023.10.27.564367](https://doi.org/10.1101/2023.10.27.564367)
