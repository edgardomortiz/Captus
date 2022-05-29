+++
title = "Output Files"
weight = 14
pre = '<i class="fas fa-dna"></i> '
+++

For this example we will use the directory `01_clean_reads` previously created with the [`clean` module]({{< ref "assembly/clean/output">}}). We run the following `Captus` command to assemble our cleaned reads:

```console
captus_assembly assemble --reads 01_clean_reads --sample_reads_target 1000000 --max_contig_gc 60.0
```

We are including the option `--sample_reads_target 1000000` to show the expected output even though [this option]({{< ref "assembly/assemble/options#--sample_reads_target" >}}) will not be very commonly used. Additionally, the option `--max_contig_gc 60.0` is used to filter contigs with GC content over 60% after assembly, only when this option is used the file `filtered_contigs.fasta` is produced.

After the run is finished we should see a new directory called `02_assemblies` with the following structure and files:

![Assemblies](/images/assemblies.png?width=640&classes=shadow)

### 1. **`[sample]__captus-asm`**
A subdirectory ending in `__captus-asm` is created to contain the assembly of each sample separately (**S1**, **S2**, **S3**, and **S4** in the image).
___
### 2. **`00_subsampled_reads`**
This directory is **only** created when the option `--sample_reads_target` is used. It contains the subsampled reads in FASTQ format that were used for the assembly.

{{% expand "Example" %}}
![FASTQ format](/images/fastq_format.png?width=1000&classes=shadow)
{{% /expand %}}
___
### 3. **`01_assembly`**
This directory contains the FASTA and FASTG assembly files as well as assembly statistics and logs.
___
### 4. **`assembly.fasta`**
The main assembly file in FASTA format, this file contains the contigs assembled by `MEGAHIT`. The sequence headers are modified by `Captus` to resemble the headers produced by the assembler `Spades`.

{{% expand "Example" %}}
![FASTA format](/images/fasta_format.png?width=1000&classes=shadow)
{{% /expand %}}
___
### 5. **`assembly_graph.fastg`**
The assembly graph in [FASTG format](http://fastg.sourceforge.net/FASTG_Spec_v1.00.pdf). This file can be explored in [Bandage](https://rrwick.github.io/Bandage/) or similar software which are able to plot the connections between contigs, loops, circular segments, etc.

{{% expand "Example" %}}
![FASTG in Bandage](/images/fastg_in_bandage.png?width=1000&classes=shadow)
{{% /expand %}}
___
### 6. **`filtered_contigs.fasta`**
This file in **only** created when the option `--max_contig_gc` is used, `Captus` will place any contig exceeding the specified GC content in this FASTA-formatted file (same format as in **4**).
___
### 7. **`assembly.stats.tsv`**, **`assembly.stats.t.tsv`**
Assembly statistics, the `assembly.stats.t.tsv` is just a transposed version of `assembly.stats.tsv`:

{{% expand "Example" %}}
**`assembly.stats.tsv`**
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
{{% /expand %}}
___
### 8. **`megahit.brief.log`**, **`megahit.full.log`**
`MEGAHIT` program logs, the _brief_ version contains just the screen output from each `MEGAHIT` run.

{{% expand "Example" %}}
**`megahit.brief.log`**
```text
Captus' MEGAHIT Command:
  megahit -1 /Volumes/Shuttle500G/for_docs_output/02_assemblies/GenusA_speciesA_CAP__captus-asm/00_subsampled_reads/GenusA_speciesA_CAP_R1.fq.gz -2 /Volumes/Shuttle500G/for_docs_output/02_assemblies/GenusA_speciesA_CAP__captus-asm/00_subsampled_reads/GenusA_speciesA_CAP_R2.fq.gz --min-count 2 --k-list 31,39,47,63,79,95,111,127,143,159,175 --merge-level 20,0.95 --prune-level 2 --memory 8504035246 --num-cpu-threads 4 --out-dir /Volumes/Shuttle500G/for_docs_output/02_assemblies/GenusA_speciesA_CAP__captus-asm/01_assembly --min-contig-len 182 --tmp-dir /Users/emortiz/captus_megahit_tmp


2022-04-08 01:09:30 - MEGAHIT v1.2.9
2022-04-08 01:09:30 - Using megahit_core with POPCNT and BMI2 support
2022-04-08 01:09:30 - Convert reads to binary library
2022-04-08 01:09:31 - b'INFO  sequence/io/sequence_lib.cpp  :   77 - Lib 0 (/Volumes/Shuttle500G/for_docs_output/02_assemblies/GenusA_speciesA_CAP__captus-asm/00_subsampled_reads/GenusA_speciesA_CAP_R1.fq.gz,/Volumes/Shuttle500G/for_docs_output/02_assemblies/GenusA_speciesA_CAP__captus-asm/00_subsampled_reads/GenusA_speciesA_CAP_R2.fq.gz): pe, 720090 reads, 151 max length'
2022-04-08 01:09:31 - b'INFO  utils/utils.h                 :  152 - Real: 1.3603\tuser: 0.6603\tsys: 0.2130\tmaxrss: 64000000'
2022-04-08 01:09:31 - Start assembly. Number of CPU threads 4 
2022-04-08 01:09:31 - k list: 31,39,47,63,79,95,111,127,143,159,175 
2022-04-08 01:09:31 - Memory used: 8504035246
2022-04-08 01:09:31 - Extract solid (k+1)-mers for k = 31 
2022-04-08 01:09:41 - Build graph for k = 31 
2022-04-08 01:09:43 - Assemble contigs from SdBG for k = 31
2022-04-08 01:09:46 - Local assembly for k = 31
2022-04-08 01:09:55 - Extract iterative edges from k = 31 to 39 
2022-04-08 01:09:57 - Build graph for k = 39 
2022-04-08 01:09:58 - Assemble contigs from SdBG for k = 39
2022-04-08 01:10:00 - Local assembly for k = 39
2022-04-08 01:10:15 - Extract iterative edges from k = 39 to 47 
2022-04-08 01:10:17 - Build graph for k = 47 
2022-04-08 01:10:18 - Assemble contigs from SdBG for k = 47
2022-04-08 01:10:20 - Local assembly for k = 47
2022-04-08 01:10:39 - Extract iterative edges from k = 47 to 63 
2022-04-08 01:10:41 - Build graph for k = 63 
2022-04-08 01:10:42 - Assemble contigs from SdBG for k = 63
2022-04-08 01:10:44 - Local assembly for k = 63
2022-04-08 01:11:03 - Extract iterative edges from k = 63 to 79 
2022-04-08 01:11:05 - Build graph for k = 79 
2022-04-08 01:11:06 - Assemble contigs from SdBG for k = 79
2022-04-08 01:11:08 - Local assembly for k = 79
2022-04-08 01:11:28 - Extract iterative edges from k = 79 to 95 
2022-04-08 01:11:30 - Build graph for k = 95 
2022-04-08 01:11:31 - Assemble contigs from SdBG for k = 95
2022-04-08 01:11:33 - Local assembly for k = 95
2022-04-08 01:11:49 - Extract iterative edges from k = 95 to 111 
2022-04-08 01:11:50 - Build graph for k = 111 
2022-04-08 01:11:51 - Assemble contigs from SdBG for k = 111
2022-04-08 01:11:53 - Local assembly for k = 111
2022-04-08 01:12:08 - Extract iterative edges from k = 111 to 127 
2022-04-08 01:12:09 - Build graph for k = 127 
2022-04-08 01:12:10 - Assemble contigs from SdBG for k = 127
2022-04-08 01:12:12 - Local assembly for k = 127
2022-04-08 01:12:26 - Extract iterative edges from k = 127 to 143 
2022-04-08 01:12:26 - Build graph for k = 143 
2022-04-08 01:12:27 - Assemble contigs from SdBG for k = 143
2022-04-08 01:12:29 - Local assembly for k = 143
2022-04-08 01:12:40 - Extract iterative edges from k = 143 to 159 
2022-04-08 01:12:40 - Build graph for k = 159 
2022-04-08 01:12:41 - Assemble contigs from SdBG for k = 159
2022-04-08 01:12:43 - Local assembly for k = 159
2022-04-08 01:12:53 - Extract iterative edges from k = 159 to 175 
2022-04-08 01:12:54 - Build graph for k = 175 
2022-04-08 01:12:54 - Assemble contigs from SdBG for k = 175
2022-04-08 01:12:56 - Merging to output final contigs 
2022-04-08 01:12:56 - 1850 contigs, total 2031858 bp, min 182 bp, max 26295 bp, avg 1098 bp, N50 1616 bp
2022-04-08 01:12:56 - ALL DONE. Time elapsed: 206.394586 seconds 
```
{{% /expand %}}
___
### 9. **`captus-assembly_assemble.stats.tsv`**
Statistics tab-separated-values table compiled across all assembled samples.
___
### 10. **`captus-assembly_assemble.report.html`**
This is the final [Assembly report]({{< ref "assembly/assemble/report">}}), summarizing statistics across all samples assembled.
___
### 11. **`captus-assembly_assemble.log`**
This is the log from `Captus`, it contains the command used and all the information shown during the run. If the option `--show_less` was enabled, the log will also contain all the extra detailed information that was hidden during the run.

___
Created by [Edgardo M. Ortiz]({{< ref "../../more/credits/#edgardo-m-ortiz">}}) (06.08.2021)  
Last modified by [Edgardo M. Ortiz]({{< ref "../../more/credits/#edgardo-m-ortiz">}}) (29.05.2022)