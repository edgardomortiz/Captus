---
title: "Options"
weight: 13
pre: '<i class="fas fa-cog"></i> '
---
# clean
___
To show all available options and their default values you can type in your terminal:
```console
captus_assembly clean --help
```

___
### *Input*
___
#### **`-r, --reads`**
With this option you provide the location of your raw FASTQ files, there are several ways to list them:

- _**Directory:**_ the path to the directory containing your FASTQ files is usually the easiest way to tell `Captus` which files to analyze. When you provide a directory, `Captus` searches within all its subdirectories for files with [valid FASTQ extensions]({{< ref "assembly/clean/preparation">}}).

- _**List of files:**_ you can also provide the individual path to each of your FASTQ files separated by spaces. This is useful if you only want to analyze only a couple of samples within a directory with many other samples for example. Another use for lists is when your FASTQ files are located in different directories.

- _**UNIX pattern:**_ another easy way to provide lists of files is using the wildcards `*` and `?` to match many or just one character respectively.

This argument is **required** <i class="fas fa-exclamation-triangle"></i>.

{{% expand "Examples" %}}
Imagine that your directory `raw_reads` has the following structure:
```console
raw_reads
├── batch_1
│   ├── C_R1.fastq
│   ├── C_R2.fastq
│   ├── D_R1.fastq
│   └── D_R2.fastq
├── A_R1.fq.gz
├── A_R2.fq.gz
├── B_R1.fq.gz
└── B_R2.fq.gz
```
- If you want to analyze all the FASTQ files in `raw_reads` (including the subdirectory `batch_1`) you can simply use `-r raw_reads`
- To analyze only samples `A` and `D`: `-r raw_reads/A_R1.fq.gz raw_reads/A_R2.fq.gz raw_reads/batch_1/D_R1.fastq raw_reads/batch_1/D_R2.fastq` or more easily `-r raw_reads/A_R?.fq.gz raw_reads/batch_1/D_R?.fastq`
- To analyze only the samples inside `raw_reads` but not in `batch_1`: `-r raw_reads/*.*`
{{% /expand %}}
___
### *Output*
___
#### **`-o, --out`**
With this option you can redirect the output directory to a path of your choice, that path will be created if it doesn't already exist.

This argument is optional, the default is **./01_clean_reads/**
___
#### **`--keep_all`**
Many intermediate files are created during the read cleanup, some are large (like FASTQ files) while others small (like temporary logs). `Captus` deletes all the unnecesary intermediate files unless you enable this flag.
___
#### **`--overwrite`**
Use this flag with caution, this will replace any previous result within the output directory (for the sample names that match).
___
### *Adaptor trimming*
___
#### **`--adaptor_set`**
We have bundled with `Captus` adaptor sequences, these options are available:

- `Illumina` = Adaptor set copied from `BBTools`.
- `BGI` = Including BGISEQ, DNBSEQ, and MGISEQ.
- `ALL` = If you are unsure of the technology used for your sequences this combines both sets of adaptors.

This argument is optional, the default is **Illumina**.
___
#### **`--rna`**
Enable this flag to trim poly-A tails from RNA-Seq reads.
___
### *Quality trimming and filtering*
Here you can control [PHRED](https://drive5.com/usearch/manual/quality_score.html) quality score thresholds. `BBTools` uses the [PHRED algorithm](http://seqanswers.com/forums/showpost.php?p=144154&postcount=17) to trim low-quality bases or to discard low-quality reads.
___
#### **`--trimq`**
Leading and trailing read regions with average PHRED quality score below this value will be trimmed.

Many people raise this value to 20 or even higher but that usually [discards lots of useful data for *de novo* assembly](https://www.biostars.org/p/124207/). In general, unless you have really high sequencing depth, don't increase this threshold beyond ~16.

This argument is optional, the default is **13**.
___
#### **`--maq`**
Once the trimming of low-quality bases from both ends of the reads has been completed, the average PHRED score of the entire read is recalculated and reads that do not have at least this **m**inimum **a**verage **q**uality are discarded.

Again, very high thresholds will throw away useful data. In general, set it to at least `trimq` or just a couple numbers higher.

This argument is optional, the default is **16**.
___
#### **`--ftl`**
Trim any base to the left of this position. For example, if you want to remove 4 bases from the left of the reads set this number to 5.

This argument is optional, the default is **0** (no `ftl` applied).
___
#### **`--ftr`**
Trim any base to the right of this position. For example, if you want to truncate your reads length to 100 bp set this number to 100

This argument is optional, the default is **0** (no `ftr` applied).
___
### *QC Statistics*
___
#### **`--qc_program`**
Select the program for obtaining the statistics from your FASTQ files. Both programs should return identical results, but `Falco` is much faster. Valid options are:
- `Falco`
- `FastQC`

This argument is optional, the default is **Falco**.
___
#### **`--skip_qc_stats`**
This flag disables the `Falco` or `FastQC` analysis, keep in mind that the final HTML report can't be created without the results from this analysis.
___
### *Other*
___
#### **`--bbduk_path`**, **`--falco_path`**, **`--fastqc_path`**
If you have installed your own copies of `bbduk.sh`, `Falco`, or `FastQC` you can provide the full path to those copies.

These arguments are optional, the defaults are **bbduk.sh**, **falco**, and **fastqc** respectively.

___
#### **`--ram`**, **`--threads`**, **`--concurrent`**, **`--debug`**, **`--show_less`**
See [Parallelization (and other common options)]({{< ref "parallelization">}})

___
Created by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (2021-08-06)  
Last modified by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (2021-09-15)