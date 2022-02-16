---
title: "Basic"
weight: 1
pre: '<i class="fas fa-school"></i> '
plotly: true
summary: through analysing target-capture sequence data.
---
This basic tutorial takes you through
In this basic tutorial, you will master the most common and minimal usage of the `Captus` assembly pipeline through analyzing four sample data.

### Preparation

---

#### Installation

To run this tutorial, you need to install the following programs on your system:

- `Captus` and its dependencies (see [Installation]({{< ref "basics/installation">}}))
- `IQ-TREE` (see [IQ-TREE documentation](http://www.iqtree.org/doc/Quickstart#installation))

#### Getting data

Download the file below and place it in a directory where you want to run this tutorial.

After moving to the directory using `cd` command, run the following command to decompress the archive and then delete it.

```shell
tar -zxvf 00_raw_reads.tar.gz && rm 00_raw_reads.tar.gz
```

Now you should have a directory named `00_raw_reads` containing eight gzip-compressed FASTQ files:

```console
00_raw_reads
├── GenusA_speciesA_CAP_R1.fq.gz
├── GenusA_speciesA_CAP_R2.fq.gz
├── GenusB_speciesB_CAP_R1.fq.gz
├── GenusB_speciesB_CAP_R2.fq.gz
├── GenusC_speciesC_CAP_R1.fq.gz
├── GenusC_speciesC_CAP_R2.fq.gz
├── GenusD_speciesD_CAP_R1.fq.gz
└── GenusD_speciesD_CAP_R2.fq.gz
```

representing paired-end reads (R1 and R2) from four samples (`GenusA_speciesA_CAP`, `GenusB_speciesB_CAP`, `GenusC_speciesC_CAP`, `GenusD_speciesD_CAP`)
These data are obtained by following a method described in [Johnson *et al*. (2018)](https://academic.oup.com/sysbio/article/68/4/594/5237557).  
Briefly, genomic DNA were extracted from harbarium specimens and fragmented by sonication into ca. 350-bp, Illumina library prepared, enriched universal (conserved) coding sequences among angiosperms using [Angiosperms353](https://github.com/mossmatters/Angiosperms353) probe set, PCR amplified, pool sequenced with paired-end mode.
{{% notice note %}}
When you analyze your own data, please be sure that all files are named according to the [naming convention]({{< ref "assembly/clean/preparation.md">}}).
All files are named according to the [naming convention]({{< ref "assembly/clean/preparation.md">}}).
{{% /notice %}}

### 1. Clean

---

Now we would like to trim adapters and low-quality bases out from the raw reads.  
`Captus` will automatically recognize sample names and library layout(s) (single-end or paired-end) based on the file names, and process them appropriately.
Run the following command to clean all FASTQ files in the `00_raw_reads` directory.
You can customize cleaning criteria in detail, see [options]({{< ref "assembly/clean/options">}}).

```shell
captus_assembly clean -r 00_raw_reads
```

- **`-r`** Path to directory containing FASTQ files

This command will end up creating a directory, `01_clean_reads`, containing the following directories and files:

```console
01_clean_reads
├── 00_adaptors_trimmed
├── 01_qc_stats_before
├── 02_qc_stats_after
├── 03_qc_extras
├── GenusA_speciesA_CAP.cleaning.log
├── GenusA_speciesA_CAP.cleaning.stats.txt
├── GenusA_speciesA_CAP_R1.fq.gz
├── GenusA_speciesA_CAP_R2.fq.gz
├── GenusB_speciesB_CAP.cleaning.log
├── GenusB_speciesB_CAP.cleaning.stats.txt
├── GenusB_speciesB_CAP_R1.fq.gz
├── GenusB_speciesB_CAP_R2.fq.gz
├── GenusC_speciesC_CAP.cleaning.log
├── GenusC_speciesC_CAP.cleaning.stats.txt
├── GenusC_speciesC_CAP_R1.fq.gz
├── GenusC_speciesC_CAP_R2.fq.gz
├── GenusD_speciesD_CAP.cleaning.log
├── GenusD_speciesD_CAP.cleaning.stats.txt
├── GenusD_speciesD_CAP_R1.fq.gz
├── GenusD_speciesD_CAP_R2.fq.gz
├── captus-assembly_clean.log
└── captus-assembly_clean.report.html
```

Of these, the files with extension `*.fq.gz` are the files containing the cleaned reads, and going to be used for next step.  
For descriptions of the other output files, see [here]({{< ref "assembly/clean/output">}}).

### 2. Assemble

---
Next, we would like to assemble contigs using the clean reads from the previous step.
Run the following command to perform *de novo* assembly using all the cleaned reads for each sample.

```shell
captus_assembly assemble -r 01_clean_reads
```

- **`-r`** Path to directory containing clean reads

You can customize assembler settings, see [options]({{< ref "assembly/assemble/options">}}).

This command will end up creating a directory, `02_assemblies`, containing the following directories and files:

```console
02_assemblies
├── GenusA_speciesA_CAP__captus-asm
├── GenusB_speciesB_CAP__captus-asm
├── GenusC_speciesC_CAP__captus-asm
├── GenusD_speciesD_CAP__captus-asm
├── captus-assembly_assemble.log
├── captus-assembly_assemble.report.html
└── captus-assembly_assemble.stats.tsv
```

`assembly.fasta` in each directory are the contigs assembled.  
For descriptions of the other output files, see [here]({{< ref "assembly/assemble/output">}}).

### 3. Extract

---
Now let's find out the target from the contigs
we would like to find target from the assembled contigs.
the following command to extract

```shell
captus_assembly extract -a 02_assemblies -n Angiosperms353
```

- **`-a`** Path to 
- **`-n`** nuclear protein reference. Here we use `Angiosperms353` built-in reference dataset.
Set of nuclear protein references,
options are:
Angiosperms353 = The original set of target proteins from Angiosperms353

Alternatively, provide a path to a FASTA file containing your reference protein sequences in either nucleotide or aminoacid. When the FASTA file is in nucleotides, '--nuc_transtable' will be used to translate it to aminoacids.

For descriptions of the other available options, see [options]({{< ref "assembly/extract/options">}}).

The command will end up with creating a directory, `03_extractions`, containing the following directories and files:

```console
03_extractions
├── GenusA_speciesA_CAP__captus-ext
├── GenusB_speciesB_CAP__captus-ext
├── GenusC_speciesC_CAP__captus-ext
├── GenusD_speciesD_CAP__captus-ext
├── captus-assembly_extract.log
├── captus-assembly_extract.refs.json
├── captus-assembly_extract.report.html
└── captus-assembly_extract.stats.tsv
```

For descriptions of the other output files, see [here]({{< ref "assembly/extract/output">}}).

### 4. Align

---
align across all samples

```shell
captus_assembly align -e 03_extractions
```

- **`-e`** Path to the output directory that contains the assemblies and extractions from previous steps of Captus-assembly. This directory is called '02_assemblies' if you did not specify a different name during the 'assemble' or 'extract' steps (default: ./03_extractions)

For descriptions of other available options, see [here]({{< ref "assembly/align/options">}}).  

The command will end up creating a directory, `04_alignments`, containing the following directories and files:

```console
04_alignments
├── 01_unaligned
├── 02_aligned_untrimmed
├── 03_aligned_trimmed
├── captus-assembly_align.log
├── captus-assembly_align.paralog_stats.tsv
├── captus-assembly_align.report.html
└── captus-assembly_align.stats.tsv
```

`02_aligned_untrimmed`, `03_aligned_trimmed` directories contain alignments in multi-FASTA format.
For descriptions of the other output files, see [here]({{< ref "assembly/align/output">}}).

### 5. Phylogenetic inference

---


Infer a phylogenetic tree using the alignments using [IQ-TREE](http://www.iqtree.org).
Captus creates several versions of alignment, users have to choose 
Here we show an example using the trimmed amino-acid alignments

Run the following commands to infer a concatenation-based species tree with an edge-linked proportional partition model:

```shell
# Create a directory and enter
mkdir 05_phylogeny && cd 05_phylogeny
# Create symbolic link
ln -s ../04_alignments/03_aligned_trimmed/06_careful_no_refs/01_coding_NUC/01_AA ALN_DIR
# Run IQ-TREE
iqtree -p ALN_DIR --prefix concat -B 1000 -T AUTO
```

- `-p` NEXUS/RAxML partition file or directory with 1000 ultrafast bootstrap and an Edge-linked proportional partition model
- `--prefix` Prefix for all output files (default: aln/partition)
- `-B` No. of ultrafast bootstrap replicates
- `-T` No. cores/threads or AUTO-detect (default: 1)

```console
05_phylogeny
├── ALN_DIR -> ../04_alignments/03_aligned_trimmed/06_careful_no_refs/01_coding_NUC/01_AA
├── concat.best_model.nex
├── concat.best_scheme
├── concat.best_scheme.nex
├── concat.bionj
├── concat.ckp.gz
├── concat.contree
├── concat.iqtree
├── concat.log
├── concat.mldist
├── concat.model.gz
├── concat.splits.nex
└── concat.treefile
```

Of these, `concat.treefile` is the inferred maximum likelihood phylogenetic tree.
For descriptions of the other output files, see [here]().
You can open this file with a tree viewer software such as [FigTree](http://tree.bio.ed.ac.uk/software/figtree).

tree rooted using `GenusC_speciesC_CAP` should look like:
![figtree](/images/tutorial_basic_figtree.png?height=800))

This is a **minimal** tutorial, [Advanced tutorial]({{< ref "tutorials/assembly/advanced">}}) will show you more deeper world.

---
Created by [Gentaro Shigita]({{< ref "../../credits/#gentaro-shigita">}}) (01.10.2021)  
Last modified by [Gentaro Shigita]({{< ref "../../credits/#gentaro-shigita">}}) (16.11.2021)
