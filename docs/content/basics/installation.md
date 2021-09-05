---
title: "Installation"
weight: 15
---

### *The easy way (recommended)*

The easiest way (until we eventually get `Captus` installing directly from Bioconda) is to create a separate environment using [miniconda](https://docs.conda.io/en/latest/miniconda.html) for all the dependencies and *THEN* clone and install `Captus` within that environment.

First, let's create the `conda` environment:
```console
conda create -n captus -c bioconda -c conda-forge -c jlsteenwyk clipkit "python>=3.6" pandas plotly tqdm "perl-bioperl-core>=1.007002" bbmap falco fastqc mafft mmseqs2 megahit pigz vsearch
```

Once the environment is ready, let's activate it:
```console
conda activate captus
```

Notice that the beginning of your prompt should have changed from`(base)$` to `(captus)$` as we activate the environment. Now let's clone the `Captus` repository and install the package:

```console
git clone https://github.com/edgardomortiz/Captus.git
cd Captus
pip install .
```

And that is all!, just to verify it is installed try typing `captus_assembly --help`, if everything went OK you should see the following output in the terminal:

```console
(captus)$ captus_assembly --help
usage: captus_assembly command [options]

Captus 0.0.17: Assembly of Phylogenomic Datasets from High-Throughput Sequencing data

Captus-assembly commands:
  command     Program commands (in typical order of execution)
                clean = Trim adaptors and quality filter reads with BBTools, run FastQC on the raw
                        and cleaned reads
                assemble = Perform de novo assembly with MEGAHIT: Assembling reads that were
                           cleaned with the 'clean' command is recommended, but reads cleaned
                           elsewhere are also allowed
                extract = Recover targeted markers with BLAT and Scipio: Extracting markers from
                          the assembly obtained with the 'assemble' command is recommended, but
                          any other assemblies in FASTA format are also allowed.
                align = Align extracted markers across samples with MAFFT: Marker alignment
                        depends on the directory structure created by the 'extract' command. This
                        step also performs paralog filtering and alignment trimming using ClipKIT

Help:
  -h, --help  Show this help message and exit
  --version   Show Captus' version number

For help on a particular command: captus_assembly command -h
```

{{% notice tip %}}
If you don't want to install `Captus` you can simply add the directory where you cloned the repository to your system `$PATH` and use `captus_assembly-runner.py` instead of `captus_assembly`
{{% /notice %}}
___
### *The manual way*

You will have to install all the the dependencies separately yourself:

`Captus` was written for `python >= v3.6`, the only required library is `tqdm` but if you want to produce the HTML reports you will also need `pandas` and `plotly`

- `BBTools` (https://jgi.doe.gov/data-and-tools/bbtools/)

- `Falco` (https://github.com/smithlabcode/falco) or `FastQC` (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

- `MEGAHIT` (https://github.com/voutcn/megahit)

- `Scipio` (https://www.webscipio.org/) *

- `BioPerl` (https://bioperl.org/)

- `YAML` (https://metacpan.org/pod/YAML)

- `BLAT >= 36x7` (http://hgdownload.soe.ucsc.edu/admin/exe/) *

- `MMseqs2` (https://github.com/soedinglab/MMseqs2)

- `MAFFT` (https://mafft.cbrc.jp/alignment/software/)

- `ClipKIT` (https://github.com/JLSteenwyk/ClipKIT)

- `pigz` (https://zlib.net/pigz/)

\* Bundled with `Captus`

Once you have all the dependencies installed you can proceed to clone the repository and install `Captus` as described before:

```console
git clone https://github.com/edgardomortiz/Captus.git
cd Captus
pip install .
```

___
Created by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (2021-08-06)  
Last modified by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (2021-09-05)