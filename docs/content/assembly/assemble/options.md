---
title: "Options"
weight: 13
pre: '<i class="fas fa-cog"></i> '
---
# assemble
___
To show all available options and their default values you can type in your terminal:
```console
captus_assembly assemble --help
```

___
### *Input*
___
#### **`-r, --reads`**
With this option you provide the location of your clean FASTQ files, there are several ways to list them:

- _**Directory:**_ providing the path to the directory containing your cleaned FASTQ files is the typical way to tell `Captus` which files to analyze. Subdirectories will not be searched in this case. If you cleaned your reads with `Captus` just use `-r 01_clean_reads` or the name you gave to the output directory.

- _**List of files:**_ you can also provide the individual path to each of your FASTQ files separated by a single space. This is useful if you only want to analyze only a couple of samples within a directory with many other samples. Another use for lists is when your FASTQ files are located in different directories.

- _**UNIX pattern:**_ another easy way to provide lists of files is using the wildcards `*` and `?` to match many or just one character respectively.

This argument is **required** <i class="fas fa-exclamation-triangle"></i>, the default is **./01_clean_reads/**

{{% expand "Read this if you want to import FASTQ files cleaned elsewhere" %}}
Imagine that you have a directory `clean_reads` with the following structure:
```console
clean_reads
├── A_R1.fastq.gz
├── A_R2.fastq.gz
├── B_R1.fastq.gz
├── B_R2.fastq.gz
├── C_R1.fq
├── C_R2.fq
├── D_R1.fq
└── D_R2.fq
```
- If you want to analyze all the FASTQ files in `clean_reads` you can simply use `-r clean_reads`
- To analyze only samples `A` and `D`: `-r cleaned_reads/A_R1.fastq.gz cleaned_reads/A_R2.fastq.gz cleaned_reads/D_R1.fq cleaned_reads/D_R2.fq` or more easily `-r cleaned_reads/A_R?.fastq.gz cleaned_reads/D_R?.fq`
{{% /expand %}}

___
#### **`--sample_reads_target`**
In case that you want to subsample a fixed amount of reads (e.g. if your FASTQ files have hundreds of millions of reads) you can indicate the number with this option. For example, `--sample_reads_target 10_000_000` will randomly sample 10 million reads (if single-end) or 10 million pairs (if paired-end).

This argument is optional, the default is **0** (= no subsampling).
___
### *Output*
___
#### **`-o, --out`**
With this option you can redirect the output directory to a path of your choice, that path will be created if it doesn't already exist.

Inside this directory, the assembly of each sample will be stored in a subdirectory with the ending `__captus-asm`.

This argument is optional, the default is **./02_assemblies/**
___
#### **`--keep_all`**
Many intermediate files are created by MEGAHIT during assembly, some are large (like the individual FASTA files from each kmer size), `Captus` deletes all the unnecesary intermediate files unless you enable this flag.
___
#### **`--overwrite`**
Use this flag with caution, this will replace any previous result within the output directory (for the sample names that match).
___
### *MEGAHIT*
___
#### **`--k_list`**
Comma-separated list (no spaces) of kmer sizes, all must be odd values in the range 15-255, in increments of at most 28. The final kmer size will be adjusted automatically so it never exceeds the mean read length of the sample by more than 31 (so don't worry if the largest kmer size in the default list seems too big for your read length).

This argument is optional, the default is **31,39,47,63,79,95,111,127,143,159,175** (= optimized for hybridization capture and/or genome skimming data).
___
#### **`--min_count`**
Minimum contig depth (called multiplicity in MEGAHIT). Acceptable values are integers >= 1. If your FASTQ files have many tens of million reads, it is recommended to use `--min_count 3` or higher to avoid low-coverage contigs resulting from erroneous reads (the higher the sequencing depth, the higher the number of reads with errors that will get assembled into spurious contigs).

This argument is optional, the default is **2**.
___
#### **`--prune_level`**
Prunning strength for low-coverage edges during graph cleaning. Increasing the value beyond 2 can speed up the assembly at the cost of losing low-depth contigs. Acceptable values are integers between 0 and 3.

This argument is optional, the dedfault is **2**.
___
#### **`--merge_level`**
Thresholds for merging complex bubbles, the first number multiplied by the kmer size represents the maximum bubble length to merge, the second number represents the minimum similarity required to merge bubbles. For example, `--merge_level 20,0.97` at kmer 175 will merge any two paths of at least 175 bp x 20 (3500 bp) with a similarity >= 97%.

This argument is optional, the default is **20,0.95**.
___
#### **`--preset`**
We have optimized other MEGAHIT parameter combinations specific to high-coverage (>= 30x depth) `WGS` data or `RNA`-Seq data. These modes will only work well with at least 8 GB of RAM. If, in addition to a `preset`, you provide your own `k_list`, `min_count`, or `prune_level`, your settings take priority over the preset's.
- `RNA` = `--k-list 27,47,67,87,107,127,147,167 --min-count 2 --prune-level 2`
- `WGS` = `--k-list 31,39,51,71,91,111,131,151,171 --min-count 3 --prune-level 2`

This argument is optional and has no default.
___
#### **`--min_contig_len`**
Minimum contig length in bp in output assembly.

This argument is optional, the default is **auto** (= mean read length + smallest kmer in `k_list`)
___
#### **`-tmp_dir`**
MEGAHIT needs a temporary directory in an internal hard drive, otherwise it refuses to run.

This argument is optional, the default is **$HOME**
___
### *Other*
___
#### **`--reformat_path`**, **`--megahit_path`**, **`--megahit_toolkit_path`**
If you have installed your own copies of `reformat.sh` (from `BBTools`) or `MEGAHIT` (and its `megahit_toolkit`) you can provide the full path to those copies.

These arguments are optional, the defaults are **reformat.sh**, **megahit**, and **megahit_toolkit** respectively.
___
#### **`--ram`**, **`--threads`**, **`--concurrent`**, **`--debug`**, **`--show_less`**
See [Parallelization (and other common options)]({{< ref "parallelization">}})

___
Created by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (2021-08-06)  
Last modified by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (2021-08-15)