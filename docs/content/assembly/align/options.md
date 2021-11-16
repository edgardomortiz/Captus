---
title: "Options"
weight: 13
pre: '<i class="fas fa-cog"></i> '
---
# align
___
To show all available options and their default values you can type in your terminal:
```console
captus_assembly align --help
```

___
### *Input*
___
#### **`-e, --captus_extractions_dir`**
Path to the output directory from the `extract` command, (e.g. `03_extractions` iy you used the default name). The `align` command depends entirely on the output from the `extract` step, in other words, you can't provide your unaligned or aligned FASTA files for processing.

This argument is **required** <i class="fas fa-exclamation-triangle"></i>, the default is **./03_extractions/**
___
#### **`-m, --markers`**
Which type(s) of markers to align, you can provide a comma-separated list (no spaces). These are the available marker types:
- `NUC` = **Nuc**lear proteins inside directories '01_coding_NUC'
- `PTD` = **P**las**t**i**d**ial proteins inside directories '02_coding_PTD'
- `MIT` = **Mit**ochondrial proteins inside directories '03_coding_MIT'
- `DNA` = Miscellaneous **DNA** markers inside directories '04_misc_DNA'
- `CLR` = **Cl**uste**r**-derived DNA markers inside directories '05_clusters'
- `ALL` = Shortcut for NUC,PTD,MIT,DNA,CLR

This argument is optional, the default is **ALL**.
___
#### **`-f, --formats`**
For each marker type, `Captus` creates several different formats. You can provide a comma-separated list (no spaces) of the formats you wish to align. These are the available formats:
- `AA` = Coding sequences in **a**mino**a**cids
- `NT` = Coding sequences in **n**ucleo**t**ides
- `GE` = Complete **ge**ne sequences (exons + introns) without extra flanking sequence
- `GF` = Complete **g**ene sequences with **f**lanking upstream and downstream basepairs
- `MA` = **Ma**tched sequences without extra flanking sequence
- `MF` = **M**atched sequences with **f**lanking upstream and downstream basepairs
- `ALL` = Shortcut for AA,NT,GE,GF,MA,MF

\* AA, NT, GE, and GF are valid only for NUC, PTD, and MIT markers, while MA and MF are valid only for DNA and CLR

This argument is optional, the default is **AA,NT,GE,MA**
___
#### **`-p, --max_paralogs`**
Maximum number of marker copies (paralogs) allowed per sample in an alignment. Large numbers of marker copies per sample can increase alignment times. Copies are ranked from best to worst during the extraction step, this number selects the top _n_ copies to align.

This argument is optional, the default is **5**
___
### *Output*
___
#### **`-o, --out`**
With this option you can redirect the output directory to a path of your choice, that path will be created if it doesn't already exist.

This argument is optional, the default is **./04_alignments/**
___
#### **`--keep_all`**
Many intermediate log files are created by `MAFFT` and `ClipKIT` during assembly, `Captus` deletes all the unnecesary intermediate files unless you enable this flag.
___
#### **`--overwrite`**
Use this flag with caution, this will replace any previous result within the output directory (for the sample names that match).
___
### *MAFFT*
___
#### **`--mafft_algorithm`**
Select [MAFFT's alignment algorithm](https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html). Valid algorithm names are:
- `auto` = Automatic selection by `MAFFT` based on amount of data
- `genafpair` = E-INS-i
- `localpair` = L-INS-i
- `globalpair` = G-INS-i
- `retree1` = FFT-NS-1
- `retree2` = FFT-NS-2

This argument is optional, the default is **auto**.
___
#### **`--mafft_timeout`**
Modify the waiting time in seconds for an individual alignment to complete. When using more exhaustive MAFFT algorithm (e.g. `genafpair`), alignment can take very long (up to hours depending on sample number an length of the sequences).

This argument is optional, the default is **21600** (= 6 hours).
___
### *Paralog filtering*
___
#### **`--filter_method`**
We provide two filtering methods for paralog removal, you can select either or both:
- `fast` = Only the best hit for each sample (marked as hit=00) is retained, when the reference only contains a single sequence per locus it is equivalent to the `carefu` method.
- `careful` = Only keep the copy (regardless of hit ranking) that is most similar to the reference sequence that was
 chosen most frequently among all other samples in the alignment. This method was designed to take advantage of references that contain several sequences per locus (like `Angiosperms353`), if the reference only contains a single reference per locus the result will be identical to the `fast` method.
- `both` = Two separate folders will be created, each containing the results from each filtering method.

This argument is optional, the default is **both**.
___
### *ClipKIT*
___
#### **`--clipkit_algorithm`**
Select [ClipKIT's trimming mode](https://jlsteenwyk.com/ClipKIT/advanced/index.html#modes). Valid trimming modes are:
- `smart-gap`
- `gappy`
- `kpic`
- `kpic-smart-gap`
- `kpic-gappy`
- `kpi`
- `kpi-smart-gap`
- `kpi-gappy`

This argument is optional, the default is **smart-gap**.
___
#### **`--clipkit_gaps`**
Gappynes threshold per position. Accepted values between 0 and 1. This argument is ignored when using the `kpi` and `kpic` algorithms or intermediate steps that use `smart-gap`.

This argument is optional, the default is **0.9**.
___
### *Other*
___
#### **`--redo_from`**
You can repeat the analysis without undoing all the steps. These are the points from which you ca restart the `align` command:
- `alignment` = Delete all subdirectories with alignments and restart.
- `filtering` = Delete all subdirectories with paralog-filtered alignments and restart.
- `removal` = Delete all subdirectories with alignments whose references have been removed and restart.
- `trimming` = Delete all subdirectories with trimmed alignments and restart.

This argument is optional and has no default.
___
#### **`--mafft_path`**, **`--clipkit_path`**
If you have installed your own copies of `MAFFT` or `ClipKIT` you can provide the full path to those copies.

These arguments are optional, the defaults are **mafft** and **clipkit** respectively.
___
#### **`--ram`**, **`--threads`**, **`--concurrent`**, **`--debug`**, **`--show_less`**
See [Parallelization (and other common options)]({{< ref "parallelization">}})

___
Created by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (2021-08-06)  
Last modified by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (2021-08-15)