---
title: "Options"
weight: 13
pre: '<i class="fas fa-cog"></i> '
---

To show all available options and their default values you can type in your terminal:
```console
captus_assembly clean --help
```

___
### *Input*
___
#### **`-r, --reads`**
With this option you provide the location of your FASTQ files, there are several ways to list them:

- _**Directory:**_ the path to the directory containing your FASTQ files is usually the easiest way to tell `Captus` which files to analyze. When you provide a directory, `Captus` searches within all its subdirectories for files with [valid FASTQ extensions]({{< ref "assembly/clean/preparation">}}).

- _**List of files:**_ you can also provide the individual path to each of your FASTQ files separated by a single space. This is useful if you only want to analyze only a couple of samples within a directory with many other samples. Another use for lists is when your FASTQ files are located in different directories.

- _**UNIX pattern:**_ another easy way to provide lists of files is using the wildcards `*` and `?` to match many or just one character respectively.

{{% expand "Examples:" %}}
Imagine that your directory `raw_reads` has the following structure:
```console
raw_reads
├── batch_1
│   ├── C_R1.fq
│   ├── C_R2.fq
│   ├── D_R1.fq
│   └── D_R2.fq
├── A_R1.fq.gz
├── A_R2.fq.gz
├── B_R1.fq.gz
└── B_R2.fq.gz
```
- If you want to analyze all the FASTQ files in `raw_reads` (including the subdirectory `batch_1`) you can simply use `-r raw_reads`
- To analyze only samples `A` and `D`: `-r raw_reads/A_R1.fq.gz raw_reads/A_R2.fq.gz raw_reads/batch_1/D_R1.fq raw_reads/batch_1/D_R2.fq` or more easily `-r raw_reads/A_R?.fq.gz raw_reads/batch_1/D_R?.fq`
- To analyze only the samples inside `raw_reads` but not in `batch_1`: `-r raw_reads/*.*`
{{% /expand %}}

___
### *Output*
___
#### **`-o, --out`**
With this option you can redirect the output directory to a path of your choice, that path will be created if it doesn't already exist. When you don't provide an output directory name, `Captus` creates a directory called `01_clean_reads` to contain the processed reads.
___
#### **`--keep_all`**
Many intermediate files are created during the read cleanup, some are large (like FASTQ files) while others small (like temporary logs), `Captus` deletes all the unnecesary intermediate files unless you enable this flag.
___
#### **`--overwrite`**
Use this flag with caution, this will replace any previous result within the output directory only when the sample names match.
___
### *Adaptor trimming*
___
#### **`--adaptor_set`**
We have bundled with `Captus` the adaptor sequences from `Illumina` (copied from `BBTools`) and `BGI` (which include BGISEQ, DNBSEQ, and MGISEQ), if you are unsure of the technology used for your sequences just use `--adator_set ALL`. Alternatively, if you have some custom adaptors you can provide a path to the FASTA file with your adaptor sequences (e.g. `--adaptor_set custom.fa`).
___
#### **`--rna`**
Enable this flag to trim poly-A tails from RNA-Seq reads.
___
### *Quality trimming and filtering*
Here you can control [PHRED](https://drive5.com/usearch/manual/quality_score.html) quality score thresholds. `BBTools` uses the [PHRED algorithm](http://seqanswers.com/forums/showpost.php?p=144154&postcount=17) to trim low-quality bases or to discard low-quality reads.
___
#### **`--trimq`**
Leading and trailing read regions with average PHRED quality score below this value will be trimmed.
___
#### **`--maq`**
Once the the trimming of low-quality bases from both ends of the reads has been completed, the average PHRED score of the entire read is recalculated and reads that do not have at least this **m**inimum **a**verage **q**uality are discarded.
___
### *Other*
___
#### **`--bbduk_path`**, **`--fastqc_path`**
If you have installed your own copies of `bbduk.sh` or `FastQC` you can provide the full path to those copies.
___
#### **`--skip_fastqc`**
This flag disables the `FastQC` analysis, keep in mind that the final HTML report can't be created without the results from `FastQC`.
___
#### **`--ram`**, **`--threads`**, **`--concurrent`**, **`--debug`**, **`--show_less`**
See [Parallelization (and other common options)]({{< ref "parallelization">}})
___