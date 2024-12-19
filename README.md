<p align="center">
  <img src=docs/static/images/logo.svg alt=captus_logo width=500>
</p>

### *Assembly of Phylogenomic Datasets from High-Throughput Sequencing data*

Documentation: [https://edgardomortiz.github.io/captus.docs](https://edgardomortiz.github.io/captus.docs)  

[![Install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/captus/README.html)
[![Bioconda downloads](https://anaconda.org/bioconda/captus/badges/downloads.svg)](https://anaconda.org/bioconda/captus)
[![Version in Bioconda](https://anaconda.org/bioconda/captus/badges/version.svg)](https://anaconda.org/bioconda/captus)
[![Last updated](https://anaconda.org/bioconda/captus/badges/latest_release_date.svg)](https://github.com/edgardomortiz/Captus/releases)

## INSTALLATION

### Using conda/mamba (recommended)

`Captus` is available as a [conda package](https://anaconda.org/bioconda/captus). If you have [`conda`](https://docs.conda.io/projects/conda/en/latest/index.html) or [`mamba`](https://mamba.readthedocs.io/en/latest/index.html) installed, you can easily create a new environment and install `Captus` with all dependencies using the following command:

```shell
conda create -n captus -c bioconda -c conda-forge captus
```

> [!WARNING]
>
> ### Important for macOS users
>
> One of the builds of the latest version of `MEGAHIT` (v.1.2.9) on `Bioconda` is broken. To ensure you install a functional one, please specify the build number as follows:
>
> ```shell
> conda create -n captus -c bioconda -c conda-forge captus megahit=1.2.9=hfbae3c0_0
> ```

Check that `Captus` was correctly installed:

```shell
conda activate captus
captus -h
```

If the program was correctly installed, you will see the following help message:

```text
usage: captus command [options]

Captus 1.1.0: Assembly of Phylogenomic Datasets from High-Throughput Sequencing data

Captus-assembly commands:
  command     Program commands (in typical order of execution)
                clean = Trim adaptors and quality filter reads with BBTools,
                        run FastQC on the raw and cleaned reads
                assemble = Perform de novo assembly with MEGAHIT and estimate
                           contig depth of coverage with Salmon: Assembling
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
```

### Manual installation

If you are unable to use conda/mamba for any reason, you will need to install all dependencies manually. Please refer to the [documentation](https://edgardomortiz.github.io/captus.docs/basics/installation/#manual-installation) for details.

## USAGE

Documentation and tutorials are available at [https://edgardomortiz.github.io/captus.docs](https://edgardomortiz.github.io/captus.docs)

> [!TIP]
>
> ### Extract BUSCO
>
> Since v1.0.1, `Captus` supports the use of any BUSCO lineage dataset as reference targets for extraction. Just download one of the `.tar.gz` files from [https://busco-data.ezlab.org/v5/data/lineages](https://busco-data.ezlab.org/v5/data/lineages) and provide the path to `captus extract`.  
> For example:
>
> ```text
> captus extract -a 02_assemblies -n ~/Downloads/aves_odb10.2021-02-19.tar.gz
> ```

## CITATION

Ortiz, E.M., A. Höwener, G. Shigita, M. Raza, O. Maurin, A. Zuntini, F. Forest, W.J. Baker, H. Schaefer. (2023). *A novel phylogenomics pipeline reveals complex pattern of reticulate evolution in Cucurbitales*. bioRxiv [https://doi.org/10.1101/2023.10.27.564367](https://doi.org/10.1101/2023.10.27.564367)
