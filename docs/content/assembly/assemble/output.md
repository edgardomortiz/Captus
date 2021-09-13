---
title: "Output Files"
weight: 14
pre: '<i class="fas fa-dna"></i> '
---

For this example we will use the directory `01_clean_reads` previously created with the [`clean` module]({{< ref "assembly/clean/output">}}). We run the following `Captus` command to assemble our cleaned reads:

```console
captus_assembly assemble --reads 01_clean_reads --sample_reads_target 1000000
```

We are including the option `--sample_reads_target 1000000` to show the expected output even though [this option]({{< ref "assembly/assemble/options#--sample_reads_target" >}}) will not be very commonly used.

After the run is finished we should see a new directory called `02_assemblies` with the following structure and files:

![Assemblies](/images/assemblies.png?width=604)

### 1. **`[sample]__captus-asm`**
A subdirectory ending in `__captus-asm` is created to contain the assembly of each sample separately (S1, S2, S3, and S4 in the image).
___
### 2. **`00_subsampled_reads`**
This directory is **only** created when the option `--sample_reads_target` is used. It contains the subsampled reads that were used for the assembly.
___
### 3. **`01_assembly`**
This directory contains the FASTA and FASTG assembly files as well as assembly statistics and logs.
___
### 4. **`assembly.fasta`**
The main assembly file in FASTA format, this file contains the contigs assembled by `MEGAHIT`. The sequence header have been modified by `Captus` and will have the following format:
```text
>NODE_230_length_748_cov_16.0000_k_175_flag_1
AAGCAGCCTTTAGAATTTGACTTTTTATTTGTCTTTGTTTTTTATTTATTTATTTATAATTTAAAAAAACAAAAAACAAA
AAACAAAACATGTGTTTGCTAACTATTTTTTTTTATTAAAAAAAAAATGTGTCTAATATATATTTACCGAATTAATAAAA
CAATGCTCCCAATCCATGCAAAAAAAAAGAAAAAAAAAAGAAAAAAAAAGGGGTGAGGCTGCGAAAGAACTCAGCTGTGC
CATATTCCGCATGTCATCATTAGCTATCACAATGAATAGAAGTAGAATACACATTCCATTGGTTGAAATATTATATTTCA
AAACCAAAACCAAAACCAACACCCCATCATCATCATCATCATCACCATCATCTGTTTTTCCAGTTCTAAATCTTTCTCTT
CTCTTTCCTTCACAGATCATCAGGTTTCGACAAACACCAACTTGCCTTCCACCATTTCTTCATTACTCAAATCTATCACT
CTCGTCTCCACAGACTCACACAAAGAGAGATTCATCTCATCTCATTACTACTGAAGATGGAAACTTTTACACATACACAT
GAATGAAGGAAACAAAAGTAAGCAAGCTTCTTCTTGCCTAAATATTGACAACTCAACTACAACAACAATATGAGCATTAT
ACTACTATATATGGTGATCACAACTGACCCAACCGAACGTCTCTATCAAGGCGGTCAGCGTTCAGGGTCCCCGGAAGCTC
TCTGAAACCTTCATTGAAGATGGAATAT
>NODE_1150_length_795_cov_29.0000_k_175_flag_1
GATAACCGAGGTGTCCGGAGGGGTAAGCCTCCAACTAGTAGGGGGAACCCATCTTCTTCGGGGTCCTGGCGACTCGACCG
AGAGACTGCCCCAGTGGCCCCAAAATTTTGGCACTAGAGTGGATCGAACTCAGGACGTGCGCCTAACCCACGCGTCCCAG
AATTCACCCTTACCACTAGGCCAAACCCTTGGGGTTACTTTCCGCTACTTTTTTAGAATGATTATCCTAAATCAAGAAAG
GAATAGCATTGAGAAAAACATATCATGAAAATAAAAGTTTCTGTTTGCCTAAGTATGCATCACCTGTCGGATTATTACAC
....
```
___
### 5. **`assembly_graph.fastg`**
The assembly graph in [FASTG format](http://fastg.sourceforge.net/FASTG_Spec_v1.00.pdf). This file can be explored in [Bandage](https://rrwick.github.io/Bandage/) or similar software which are able to plot the connections between contigs, loops, circular segments, etc.
![FASTG in Bandage](/images/fastg_in_bandage.png?classes=shadow)
___
### 6. **`assembly.stats.tsv`, `assembly.stats.t.tsv`**
Assembly statistics, the `assembly.stats.t.tsv` is just a transposed version of `assembly.stats.tsv`:
```text
              sample : GenusA_speciesA_CAP
           n_contigs : 1830
 pct_contigs_>=_1kbp : 33.661
 pct_contigs_>=_2kbp : 11.694
 pct_contigs_>=_5kbp : 1.749
pct_contigs_>=_10kbp : 0.984
      longest_contig : 26676
     shortest_contig : 183
        total_length : 2041738
  pct_length_>=_1kbp : 70.061
  pct_length_>=_2kbp : 42.875
  pct_length_>=_5kbp : 18.15
 pct_length_>=_10kbp : 13.836
          avg_length : 1116
       median_length : 659
                 N50 : 1648
          GC_content : 39.028
           avg_depth : 24.1
   pct_contigs_>=_1x : 100.0
   pct_contigs_>=_2x : 70.109
   pct_contigs_>=_5x : 49.563
  pct_contigs_>=_10x : 38.47
```
___
### 7. **`megahit.brief.log`, `megahit.full.log`**
`MEGAHIT` program logs, the _brief_ version contains just the screen output from each `MEGAHIT` run.
___
### 8. **`captus-assembly_assemble.stats.tsv`**
Statistics tab-separated-values table compiled across all assembled samples.
___
### 9. **`captus-assembly_assemble.report.html`**
This is the final [Assembly report]({{< ref "assembly/assemble/report">}}), summarizing statistics across all samples assembled.
___
### 10. **`captus-assembly_assemble.log`**
This is the log from `Captus`, it contains the command used and all the information shown during the run. If the option `--show_less` was enabled the log will also contain all the extra detailed information that was hidden during the run.

___
Created by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (2021-08-06)  
Last modified by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (2021-09-12)