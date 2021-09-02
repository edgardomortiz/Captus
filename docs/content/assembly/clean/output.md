---
title: "Output Files"
weight: 14
pre: '<i class="fas fa-dna"></i> '
---

Imagine we start with a directory called `00_raw_reads` with the following content:

![Raw reads](/images/raw_reads.png?width=500)

We have a samples with different data types, to recognize them we added `_CAP` to the samples where hybridization-capture was used, `_WGS` for high-coverage whole genome sequencing, `_RNA` for RNA-Seq reads, and `_GSK` for genome skimming data (notice also the difference in file sizes). For this example, we only want to clean the samples in red rectanlges correspondig to capture data. We run the following `Captus` command:

```console
(captus) $ captus_assembly clean --reads ./00_raw_reads/*_CAP_R?.fq.gz
```

Notice we are using default settings, the only required argument is the localization of the raw reads. The time it took to process these four samples was around 24 seconds and the output was written to anew directory called `01_clean_reads`. Let's take a look at the contents:

![Clean reads](/images/clean_reads.png?width=500)


