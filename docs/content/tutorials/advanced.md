+++
title = "Advanced Tutorial"
menuTitle = "Advanced"
weight = 2
pre = '<i class="fas fa-graduation-cap"></i> '
+++

Under construction...

1. Handling different data types
2. Assembling reads cleaned outside Captus
3. Combining all assembled data types
4. Importing pre-assembled sample (e.g., GenBank ref genome)
5. Extracting all marker types simultaneously
6. Discovering new markers by clustering in capture data
7. Adding new markers by clustering unused contigs (capture data)
8. Phylogenetic tree reconstruction using palalogs - ASTRAL-Pro


<!--### Preparation

---

#### Installation

Install `Captus` and dependencies on your system by following the [instruction]({{< ref "basics/installation">}}).

#### Getting data

Download the example data below, and place it in a directory where you want to run this tutorial.

Run the following command to decompress the archive and then delete it.

```shell
cd {path/to/directory}
tar -xf 00_raw_reads.tar.gz && rm 00_raw_reads.tar.gz
```

Now you have a directory named `00_raw_reads` containing 16 gzip-compressed FASTQ files:

```console
00_raw_reads
├── GenusA_speciesA_CAP_R1.fq.gz
├── GenusA_speciesA_CAP_R2.fq.gz
├── GenusA_speciesA_RNA_R1.fq.gz
├── GenusA_speciesA_RNA_R2.fq.gz
├── GenusB_speciesB_CAP_R1.fq.gz
├── GenusB_speciesB_CAP_R2.fq.gz
├── GenusC_speciesC_CAP_R1.fq.gz
├── GenusC_speciesC_CAP_R2.fq.gz
├── GenusC_speciesC_RNA_R1.fq.gz
├── GenusC_speciesC_RNA_R2.fq.gz
├── GenusD_speciesD_CAP_R1.fq.gz
├── GenusD_speciesD_CAP_R2.fq.gz
├── GenusD_speciesD_RNA_R1.fq.gz
├── GenusD_speciesD_RNA_R2.fq.gz
├── GenusH_speciesH_RNA_R1.fq.gz
└── GenusH_speciesH_RNA_R2.fq.gz
```

CAP: Illumina universal coding sequences for angiosperms are enriched using [Angiosperms353](https://github.com/mossmatters/Angiosperms353) probe set.
RNA: mRNA library enriched with poly-T probe
paired-end
They are named following the [naming convention]({{< ref "assembly/clean/preparation.md">}}).

### Clean

---

Run the following commands to clean.
`--rna` flag enable to trim poly-A tails as well as Illumina adapters.

```shell
captus_assembly clean -r 00_raw_reads/*CAP_R?.fq.gz -o 01_clean_reads_CAP
captus_assembly clean -r 00_raw_reads/*RNA_R?.fq.gz -o 01_clean_reads_RNA --rna
```

For descriptions of the other output files, see [here]({{< ref "assembly/clean/output">}}).

### Assemble

---
you can easily integrate samples assembled with another tool.

```shell
captus_assembly assemble -r 01_clean_reads_CAP -o 02_assemblies_CAP
captus_assembly assemble -r 01_clean_reads_RNA -o 02_assemblies_RNA --preset RNA
captus_assembly assemble -r 01_clean_reads_WGS -o 02_assemblies_WGS --preset WGS
```

For descriptions of the other output files, see [here]({{< ref "assembly/assemble/output">}}).

### Extract

---
Download the example data below, and place it in a directory where you want to run this tutorial.

Run the following command to extract targets from the contigs.
built-in reference dataset

```shell
captus_assembly extract -a 02_assemblies_CAP -f 02_assemblies_WGS -n Angiosperms353 -p SeedPlantsPTD -m SeedPlantsMIT -c
```

For descriptions of the other output files, see [here]({{< ref "assembly/extract/output">}}).

### Align

---

```shell
captus_assembly align
```

For descriptions of the other output files, see [here]({{< ref "assembly/align/output">}}).

### Phylogenetic inference

---

```shell
iqtree
``` -->

---
Created by [Gentaro Shigita]({{< ref "../../more/credits/#gentaro-shigita">}}) (01.10.2021)  
Last modified by [Gentaro Shigita]({{< ref "../../more/credits/#gentaro-shigita">}}) (31.10.2022)
