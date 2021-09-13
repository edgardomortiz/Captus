---
title: "Output Files"
weight: 14
pre: '<i class="fas fa-dna"></i> '
---

Imagine we start with a directory called `00_raw_reads` with the following content:

![Raw reads](/images/raw_reads.png?width=604)

We have a samples with different data types, to distinguish them we added `_CAP` to the samples where hybridization-capture was used, `_WGS` for high-coverage whole genome sequencing, `_RNA` for RNA-Seq reads, and `_GSK` for genome skimming data (notice also the difference in file sizes). For this example, we only want to clean the samples in red rectangles correspondig to capture data. We run the following `Captus` command:

```console
captus_assembly clean --reads ./00_raw_reads/*_CAP_R?.fq.gz
```

Notice we are using default settings, the only required argument is the location of the raw reads. The output was written to a new directory called `01_clean_reads`. Let's take a look at the contents:

![Clean reads](/images/clean_reads.png?width=604)

### 1. **`[sample]_R1.fq.gz`, `[sample]_R2.fq.gz`**
In case of paired-end input we will have a pair of file like in the image, the forward reads are indicated by _R1 and the reverse reads by _R2. Single-end input will only return forward reads.
___
### 2. **`[sample].cleaning.log`**
This file contains the cleaning command used for `bbduk.sh` as well the data shown as screen output, this and other information is compiled in the [Cleaning report]({{< ref "assembly/clean/report">}}).
___
### 3. **`[sample].cleaning.stats.log`**
List of contaminants found by `bbduk.sh` in the input reads, sorted by abundance.
___
### 4. **`captus-assembly_clean.report.html`**
This is the final [Cleaning report]({{< ref "assembly/clean/report">}}), summarizing statistics across all samples analyzed.
___
### 5. **`captus-assembly_clean.log`**
This is the log from `Captus`, it contains the command used and all the information shown during the run. If the option `--show_less` was enabled the log will also contain all the extra detailed information that was hidden during the run.
___
### 6. **`00_adaptors_trimmed`**
This is an intermediate directory that contains the FASTQ files without adaptors, prior to quality-trimming and filtering. The directory also stores `bbduk.sh` commands and logs for the adaptor trimming stage. If the option `--keep_all` was enabled the FASTQs from this intermediate are kept after the run, otherwise they are deleted.
___
### 7. **`01_qc_stats_before`, `02_qc_stats_after`**
These directories contain the results from either `Falco` or `FastQC`, organized in a subdirectory per FASTQ file analyzed.
___
### 8. **`03_qc_extras`**
This directory contains all the tab-separated-values tables needed to build the [Cleaning report]({{< ref "assembly/clean/report">}}). We provide them separately to allow the user more detailed analyses.

___
Created by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (2021-08-06)  
Last modified by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (2021-09-12)