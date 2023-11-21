+++
title = "Installation"
weight = 15
+++

### *The easy way (recommended)*

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
conda create -n captus captus
```

`conda` sometimes takes too long to find and configure dependencies, if that happens we recommend installing `mamba` first, and installing `Captus` with it:
```bash
conda install mamba
mamba create -n captus captus
```

Once `Captus` is installed in its own environment, let's activate it:
```console
conda activate captus
```
And that is all! Notice that the beginning of your prompt should have changed from`(base)$` to `(captus)$` as we activate the environment.  

Just to verify it is correctly installed try typing `captus_assembly --help`, if everything went OK you should see the following output in the terminal:
```console
usage: captus_assembly command [options]

Captus 1.0.0: Assembly of Phylogenomic Datasets from High-Throughput Sequencing data

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
```

___
### *The manual way*

You will have to install all the the dependencies separately yourself:

`Captus` was written for `python >= v3.7`, the only required library is `tqdm` but if you want to produce the HTML reports you will also need `pandas` and `plotly`

- `BBTools` (https://jgi.doe.gov/data-and-tools/bbtools/)

- `Falco` (https://github.com/smithlabcode/falco) or `FastQC` (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

- `MEGAHIT` (https://github.com/voutcn/megahit)

- **\*** `Scipio` (https://www.webscipio.org/)

- `BioPerl` (https://bioperl.org/)

- `YAML` (https://metacpan.org/pod/YAML)

- **\*** `BLAT >= 36x7` (http://hgdownload.soe.ucsc.edu/admin/exe/)

- `MMseqs2` (https://github.com/soedinglab/MMseqs2)

- `MAFFT` (https://mafft.cbrc.jp/alignment/software/)

- `MUSCLE` (https://www.drive5.com/muscle/)

- `ClipKIT` (https://github.com/JLSteenwyk/ClipKIT)

- `pigz` (https://zlib.net/pigz/)

**\*** Bundled with `Captus`

Once you have all the dependencies installed you can proceed to clone the repository and install `Captus` as described before:

```console
git clone https://github.com/edgardomortiz/Captus.git
cd Captus
pip install .
```

{{% notice tip %}}
If you don't want to install `Captus` you can simply add the directory where you cloned the repository to your system `$PATH` and use `captus_assembly-runner.py` instead of `captus_assembly`
{{% /notice %}}
___
Created by [Edgardo M. Ortiz]({{< ref "../more/credits/#edgardo-m-ortiz">}}) (06.08.2021)  
Last modified by [Edgardo M. Ortiz]({{< ref "../more/credits/#edgardo-m-ortiz">}}) (21.11.2023)