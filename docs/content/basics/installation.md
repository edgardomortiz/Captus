+++
title = "Installation"
weight = 15
+++

### Using conda/mamba (recommended)

`Captus` is available as a [conda package](https://anaconda.org/bioconda/captus). If you have `conda` [<i class="fas fa-question-circle fa-sm"></i>](https://docs.conda.io/projects/conda/en/latest/index.html) or `mamba` [<i class="fas fa-question-circle fa-sm"></i>](https://mamba.readthedocs.io/en/latest/index.html) installed, you can easily create a new environment and install `Captus` with all dependencies using the following command:

##### Linux:
```shell
conda create -n captus -c bioconda -c conda-forge captus
```
or
```shell
mamba create -n captus -c bioconda -c conda-forge captus
```

##### Mac computers with Intel processors:
```shell
conda create -n captus -c bioconda -c conda-forge captus megahit=1.2.9=hfbae3c0_0
```
or
```shell
mamba create -n captus -c bioconda -c conda-forge captus megahit=1.2.9=hfbae3c0_0
```

##### Mac computers with Apple silicon (M processors):
```shell
conda create --platform osx-64 -n captus -c bioconda -c conda-forge captus megahit=1.2.9=hfbae3c0_0 mmseqs2=16.747c6
```
or
```shell
mamba create --platform osx-64 -n captus -c bioconda -c conda-forge captus megahit=1.2.9=hfbae3c0_0 mmseqs2=16.747c6
```

Then check that `Captus` was correctly installed:
```shell
conda activate captus
captus -h
```

If the program was correctly installed, you will see the following help message:

```console
usage: captus command [options]

Captus 1.6.2: Assembly of Phylogenomic Datasets from High-Throughput Sequencing data

Captus-assembly commands:
  command     Program commands (in typical order of execution)
                clean = Trim adaptors and quality filter reads with BBTools,
                        run FastQC on the raw and cleaned reads
                assemble = Perform de novo assembly with MEGAHIT and estimate
                           contig depth of coverage with Salmon. Assembling
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
                        using TAPER and ClipKIT

Help:
  -h, --help  Show this help message and exit
  --version   Show Captus' version number

For help on a particular command: captus_assembly command -h
```

___

### Manual installation

If you are unable to use `conda`/`mamba` for any reason, you will need to manually install all the dependencies listed below:

|Dependency|Version|URL|
|-|-|-|
|`Python` |>=3.6|<https://www.python.org>|
|`BBTools`||<https://jgi.doe.gov/data-and-tools/bbtools>|
|`BioPerl`||<https://bioperl.org>|
|`ClipKIT`|>=1.3.0|<https://github.com/JLSteenwyk/ClipKIT>|
|`Falco`|>=0.3.0|<https://github.com/smithlabcode/falco>|
|`FastQC`||<https://www.bioinformatics.babraham.ac.uk/projects/fastqc>|
|`MAFFT`||<https://mafft.cbrc.jp/alignment/software>|
|`MEGAHIT`|1.2.9|<https://github.com/voutcn/megahit>|
|`MMseqs2`||<https://github.com/soedinglab/MMseqs2>|
|`MUSCLE`|>=5.1|<https://www.drive5.com/muscle>|
|`pandas`|>=2.1.0|<https://pandas.pydata.org>|
|`Plotly`||<https://github.com/plotly/plotly.py>|
|`pigz`||<https://zlib.net/pigz>|
|`Salmon`|>=1.10.0|<https://github.com/COMBINE-lab/salmon>|
|`tqdm`||<https://github.com/tqdm/tqdm>|
|`VSEARCH`||<https://github.com/torognes/vsearch>|
|`YAML`||<https://metacpan.org/pod/YAML>|

**\*** The following two dependencies are bundled with `Captus`, so no additional installation is required.

|Dependency|Version|URL|
|-|-|-|
|`BLAT`|37x1|<http://hgdownload.soe.ucsc.edu/admin/exe>|
|`Scipio`|1.4.1|<https://www.webscipio.org>|

Once you have all the dependencies installed, you can proceed to clone the repository and install `Captus` as follows:

```shell
git clone https://github.com/edgardomortiz/Captus.git
cd Captus
pip install .
captus -h
```

Alternatively, you can run `Captus` using the wrapper script `captus_assembly-runner.py` as follows:

```shell
git clone https://github.com/edgardomortiz/Captus.git
./Captus/captus_assembly-runner.py -h
```

___
Created by [Edgardo M. Ortiz]({{< ref "../more/credits/#edgardo-m-ortiz">}}) (06.08.2021)<br>
Last modified by [Edgardo M. Ortiz]({{< ref "../../more/credits/#gentaro-shigita">}}) (02.21.2026)
