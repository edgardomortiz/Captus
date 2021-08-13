---
title: "Options"
weight: 13
pre: '<i class="fas fa-cog"></i> '
---

To show all available options and their default values you can type in your terminal:
```console
captus_assembly align --help
```

___
### *Input*
___
#### **`-e, --captus_extractions_dir`**
Path to the output directory from the `extract` command, (e.g. `03_extractions` iy you used the default name). The `align` command depends entirely on the output from the `extract` step, in other words, you can't provide your unaligned or aligned FASTA files for processing.
___
#### **`-k, --markers`**
Which type(s) of markers to align, you can provide a comma-separated list (no spaces). These are the available marker types:
- _**NUC**_ = **Nuc**lear proteins inside directories '01_coding_NUC'
- _**PTD**_ = **P**las**t**i**d**ial proteins inside directories '02_coding_PTD'
- _**MIT**_ = **Mit**ochondrial proteins inside directories '03_coding_MIT'
- _**DNA**_ = Miscellaneous **DNA** markers inside directories '04_misc_DNA'
- _**CLR**_ = **Cl**uste**r**-derived DNA markers inside directories '05_clusters'
- _**ALL**_ = Shortcut for NUC,PTD,MIT,DNA,CLR
___
#### **`-f, --formats`**
For each marker type, `Captus` creates several different formats. You can provide a comma-separated list (no spaces) of the formats you wish to align. These are the available formats:
- _**AA**_ = Coding sequences in **a**mino**a**cids
- _**NT**_ = Coding sequences in **n**ucleo**t**ides
- _**GE**_ = Complete **ge**ne sequences (exons + introns) without extra flanking sequence
- _**GF**_ = Complete **g**ene sequences with **f**lanking upstream and downstream basepairs
- _**MA**_ = **Ma**tched sequences without extra flanking sequence
- _**MF**_ = **M**atched sequences with **f**lanking upstream and downstream basepairs
- _**ALL**_ = Shortcut for AA,NT,GE,GF,MA,MF

\* AA, NT, GE, and GF are valid only for NUC, PTD, and MIT markers, while MA and MF are valid only for DNA and CLR
___
### *Output*
___
#### **`-o, --out`**
With this option you can redirect the output directory to a path of your choice, that path will be created if it doesn't already exist. If you don't provide an output directory name, `Captus` creates a directory called `04_alignments` to contain the output.
___
#### **`--keep_all`**
Many intermediate log files are created by `MAFFT` and `ClipKIT` during assembly, `Captus` deletes all the unnecesary intermediate files unless you enable this flag.
___
#### **`--overwrite`**
Use this flag with caution, this will replace any previous results.
___
### *MAFFT*
___
#### **`--mafft_algorithm`**
Select [MAFFT's alignment algorithm](https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html). Valid algorithm names are:
- _**auto**_ 
- _**genafpair**_ = E-INS-i
- _**localpair**_ = L-INS-i
- _**globalpair**_ = G-INS-i
- _**retree1**_ = FFT-NS-1
- _**retree2**_ = FFT-NS-2
___
#### **`--mafft_timeout`**
Modify the waiting time in seconds for an individual alignment to complete. When using more exhaustive MAFFT algorithm (e.g. `genafpair`), alignment can take very long (up to hours depending on sample number an length of the sequences). `Captus` terminates an alignment if it doesn't finish in 21600 seconds (6 hours).
___
### *Paralog filtering*
___
#### **`--filter_method`**
We provide two filtering methods for paralog removal, you can select either or both:
- _**fast**_ = Only the best hit for each sample (marked as hit=00) is retained, when the reference only contains a single sequence per locus it is equivalent to the `carefu` method.
- _**careful**_ = Only keep the copy (regardless of hit ranking) that is most similar to the reference sequence that was
 chosen most frequently among all other samples in the alignment. This method was designed to take advantage of references that contain several sequences per locus (like `Angiosperms353`), if the reference only contains a single reference per locus the result will be identical to the `fast` method.
- _**both**_ = Two separate folders will be created, each containing the results from each filtering method.
___
### *ClipKIT*
___
#### **`--clipkit_algorithm`**
Select [ClipKIT's trimming mode](https://jlsteenwyk.com/ClipKIT/advanced/index.html#modes). Valid trimming modes are:
- _**smart-gap**_
- _**gappy**_
- _**kpic**_
- _**kpic-smart-gap**_
- _**kpic-gappy**_
- _**kpi**_
- _**kpi-smart-gap**_
- _**kpi-gappy**_
___
#### **`--clipkit_gaps`**
Gappynes threshold per position. Accepted values between 0 and 1. This argument is ignored when using the `kpi` and `kpic` algorithms or intermediate steps that use `smart-gap`.
___
### *Other*
___
#### **`--redo_from`**
You can repeat the analysis without undoing all the steps. These are the points from which you ca restart the `align` command:
- _**alignment**_ = Delete all subdirectories with alignments and restart.
- _**filtering**_ = Delete all subdirectories with paralog-filtered alignments and restart.
- _**removal**_ = Delete all subdirectories with alignments whose references have been removed and restart.
- _**trimming**_ = Delete all subdirectories with trimmed alignments and restart.
___
#### **`--mafft_path`**, **`--clipkit_path`**
If you have installed your own copies of `MAFFT` or `ClipKIT` you can provide the full path to those copies.
___
#### **`--ram`**, **`--threads`**, **`--concurrent`**, **`--debug`**, **`--show_less`**
See [Parallelization (and other common options)]({{< ref "parallelization">}})

___
Created by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (11.08.2021)  
Last modified by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (11.08.2021)