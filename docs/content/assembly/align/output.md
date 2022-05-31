+++
title = "Output Files"
weight = 14
pre = '<i class="fas fa-dna"></i> '
+++

For this example we will use the directory `03_extractions` previously created with the [`extract` module]({{< ref "assembly/extract/output">}}). We run the following `Captus` command to collect markers across samples and align them:

```console
captus_assembly align align -e 03_extractions_CAP/ -o 04_alignments_CAP -k ALL -f ALL
```

After the run is finished we should see a new directory called `04_alignments` with the following structure and files:

![Alignments](/images/alignments.png?width=640&classes=shadow)
___
{{% expand "Complete structure of the alignment output directory" %}}
![Alignment directory](/images/alignment_stages.png?width=1200&classes=shadow)
{{% /expand %}}
___
### 1. **`01_unaligned`**
This directory contains the unaligned FASTA files corresponding to each marker that were gathered from the extractions directory. The files are organized in subdirectories, first by [marker type]({{< relref "assembly/align/options#-m---markers" >}}) and then by [format]({{< relref "assembly/align/options#-f---formats" >}}).
___
### 2. **`02_untrimmed`**
This directory contains the aligned FASTA files corresponding to each file in the `01_unaligned` directory. The files are organized in subdirectories, first by filtering strategy, then by [marker type]({{< relref "assembly/align/options#-m---markers" >}}), and finally by [format]({{< relref "assembly/align/options#-f---formats" >}}). The subdirectory structure is identical to the one inside the `03_trimmed` directory (see **4** to **15** below).
![Untrimmed alignments](/images/alignment_untrimmed_stages.png?width=1200&classes=shadow)
___
### 3. **`03_trimmed`**
All the files present in the `02_untrimmed` directory are trimmed using `ClipKIT` which removes columns that are mostly empty (see options [`--clipkit_algorithm`]({{< relref "assembly/align/options#--clipkit_algorithm" >}}), [`--clipkit_gaps`]({{< relref "assembly/align/options#--clipkit_gaps" >}})), then `Captus` removes sequences that are too short after trimming ([`--min_coverage`]({{< relref "assembly/align/options#--min_coverage" >}})). The files are organized in subdirectories, first by filtering strategy, then by [marker type]({{< relref "assembly/align/options#-m---markers" >}}), and finally by [format]({{< relref "assembly/align/options#-f---formats" >}}). The subdirectory structure is identical to the one inside the `02_untrimmed` directory (see **4** to **15** below).
![Trimmed alignments](/images/alignment_trimmed_stages.png?width=1200&classes=shadow)
___
### 4. **`01_unfiltered_w_refs`**
This directory contains the alignments before performing any filtering. All the reference sequences selected by at least a sample will be present as well as all the paralogs per sample. The files are organized in subdirectories, first by [marker type]({{< relref "assembly/align/options#-m---markers" >}}) and then by [format]({{< relref "assembly/align/options#-f---formats" >}}).
___
### 5. **`02_naive_w_refs`**
This directory contains the alignments where paralogs have been filtered by the `naive` method, which consists in simply keeping the best hit per sample (hit ranked as `00`). All the reference sequences selected by at least a sample will still be present. The files are organized in subdirectories, first by [marker type]({{< relref "assembly/align/options#-m---markers" >}}) and then by [format]({{< relref "assembly/align/options#-f---formats" >}}).
![Naive paralog filter](/images/paralog_filter_naive.png?width=1200&classes=shadow)
___
### 6. **`03_informed_w_refs`**
This directory contains the alignments where paralogs have been filtered by the `informed` method. Under this strategy, `Captus` compares every copy to the most commonly used reference sequence (sequence `ABCD-3400` in the figure) and retains the copy with the highest similarity to that reference, regardless of its paralog ranking (in the figure, `Sample1` and `Sample4` whose selected copies had paralog rankings of `01` and `02` respectively). All the reference sequences selected by at least a sample will still be present. The files are organized in subdirectories, first by [marker type]({{< relref "assembly/align/options#-m---markers" >}}) and then by [format]({{< relref "assembly/align/options#-f---formats" >}}).
![Informed paralog filter](/images/paralog_filter_informed.png?width=1200&classes=shadow)
___
### 7. **`04_unfiltered`**, **`05_naive`**, **`06_informed`**
These contain equivalent alignments to directories `01_unfiltered_w_refs`, `02_naive_w_refs`, and `03_informed_w_refs` respectively, but excluding the reference sequences. *In most cases you will estimate phylogenies from the trimmed versions of these alignments.*
___
### 8. **`01_coding_NUC`**, **`02_coding_PTD`**, **`03_coding_MIT`**
These directories contain the aligned **coding** markers from the **NUC**lear, **P**las**T**i**D**ial, and **MIT**ochondrial genomes respectively.  
The alignments are presented in four formats: protein sequence (**coding_AA**), coding sequence in nucleotide (**coding_NT**), exons and introns concatenated (**genes**), and the concatenation of exons and introns flanked by a fixed length of sequence (**genes_flanked**):

![Protein extraction formats](/images/protein_extraction.png?width=600&classes=shadow)
___
### 9. **`01_AA`**
This directory contains the protein alignments (`AA` in the figure above) of the extracted markers gathered across samples. One FASTA file per marker, with extension `.faa`.
___
### 10. **`02_NT`**
This directory contains the alignments of coding sequence in nucleotides (`NT` in the figure above) of the extracted markers gathered across samples. One FASTA file per marker, with extension `.fna`.
___
### 11. **`03_genes`**
This directory contains the alignments of gene sequence (exons + introns) in nucleotides (`GE` in the figure above) of the extracted markers gathered across samples. One FASTA file per marker, with extension `.fna`.
___
### 12. **`04_genes_flanked`**
This directory contains the alignments of flanked gene sequence in nucleotides (`GF` in the figure above) of the extracted markers gathered across samples. One FASTA file per marker, with extension `.fna`.
___
### 13. **`04_misc_DNA`**, **`05_clusters`**
These directories contain the aligned **miscellaneous DNA** markers, either from a **DNA** custom set of references or from the **CL**uste**R**ing resulting from using the option `--cluster_leftovers` during the extraction step.  
The alignments are presented in two formats: matching DNA segments (**matches**), and the matched segments including flanks and other intervening segments not present in the reference (**matches_flanked**).

![Miscellaneous DNA extraction formats](/images/misc_dna_extraction.png?width=600&classes=shadow)
___
### 14. **`01_matches`**
This directory contains the alignments of DNA sequence matches (`MA` in the figure above) for the extracted markers gathered across samples. One FASTA file per marker, with extension `.fna`.
___
### 15. **`02_matches_flanked`**
This directory contains the alignments of DNA sequence matches (`MF` in the figure above) including flanks and intervening segments not present in the references for the extracted markers gathered across samples. One FASTA file per marker, with extension `.fna`.
___
### 16. **`captus-assembly_align.paralogs.tsv`**
A tab-separated-values table recording which copy was selected during the `informed` filtering of paralogs.

{{% expand "Information included in the table" %}}
|Column|Description|
|-|-|
|**marker_type**|Type of marker. Possible values are `NUC`, `PTD`, `MIT`, `DNA`, or `CLR`.|
|**format_filtered**|Alignment format used for identity comparison and filtering. For protein references, if the reference is provided both in aminoacids and nucleotides, the nucleotide format is used.|
|**locus**|Name of the locus.|
|**ref**|Name of most common the reference observed in the alignment.|
|**sample**|Name of the sample.|
|**hit**|Paralog ranking.|
|**sequence**|Name of the sequence as it appears in the alignments.|
|**length**|Lenght of the sequence in aminoacids if `format_filtered` is `AA`, or in nucleotides if `format_filtered` is `NT`|
|**identity**|Identity of the sequence to the `ref` sequence, over the overlapping segment, ignoring internal gaps.|
|**paralog_score**|(`length` / length of `ref`) * (`identity` / 100), the highest score decides the selected copy.|
|**accepted**|Whether the copy is accepted (`TRUE`) or not (`FALSE`).|
{{% /expand %}}
___
### 17. **`captus-assembly_align.alignments.tsv`**
A tab-separated-values table recording alignment statistics for each of the alignments produced.

{{% expand "Information included in the table" %}}
|Column|Description|
|-|-|
|**path**|Absolute path to the alignment file.|
|**trimmed**|Whether the alignment has been trimmed (`TRUE`) or not (`FALSE`).|
|**paralog_filter**|Filtering strategy applied to the file. Possible values are `unfiltered`, `naive`, or `informed`.|
|**with_refs**|Whether the file still contains the reference sequences (`TRUE`) or they have been removed (`FALSE`).|
|**marker_type**|Type of marker. Possible values are `NUC`, `PTD`, `MIT`, `DNA`, or `CLR`.|
|**format**|Alignment format. Possible values are `AA`,`NT`,`GE`,`GF`,`MA`,`MF`.|
|**locus**|Name of the locus.|
|**seqs**|Number of sequences in the alignment.|
|**samples**|Number of samples represented in the alignment. The number can be smaller than `seqs` if the alignment has paralogs.|
|**avg_copies**|`seqs` / `samples`|
|**sites**|Alignment length.|
|**informative**|Number of parsimony-informative-sites in the alignment.|
|**informativeness**|(`informative` / `sites`) * 100|
|**uninformative**|`constant` + `singleton`|
|**constant**|Number of invariant sites in the alignment.|
|**singleton**|Number of sites with variaton in a single sequence.|
|**paterns**|Number of unique columns, for a detailed explanation see IQ-TREE's F.A.Q.: [What are the differences between alignment columns/sites and patterns?](http://www.iqtree.org/doc/Frequently-Asked-Questions).|
|**avg_pid**|Average pairwise identity between sequences in the alignment.|
|**missingness**|Proportion of missing data (`-`, `N`, `X`, `?`, `.`, `~`) in the alignment.|
{{% /expand %}}
___
### 18. **`captus-assembly_align.samples.tsv`**
A tab-separated-values table recording sample statistics across the different filtering and trimming stages, as well as marker types and formats.

{{% expand "Information included in the table" %}}
|Column|Description|
|-|-|
|**sample**|Sample name.|
|**stage_marker_format**|Trimming stage / Filtering strategy / Marker type / Format.|
|**locus**|Name of the locus.|
|**cov_gapped**|Coverage of the sequence relative to alignment length, including internal gaps.|
|**cov_ungapped**|Coverage of the sequence relative to alignment length, excluding internal gaps.|
|**pct_ambig**|Percentage of ambiguities in the sequence, excluding gaps.|
|**num_copies**|Number of copies in the alignment.|
{{% /expand %}}
___
### 19. **`captus-assembly_align.astral-pro.tsv`**
ASTRAL-Pro requires a tab-separated-values file for mapping the names of the paralog sequence names (first column) to the name of the sample (second column). `Captus` produces this file automatically.

{{% expand "Example" %}}
```text
GenusA_speciesA_CAP	GenusA_speciesA_CAP
GenusA_speciesA_CAP__00	GenusA_speciesA_CAP
GenusA_speciesA_CAP__01	GenusA_speciesA_CAP
GenusA_speciesA_CAP__02	GenusA_speciesA_CAP
GenusA_speciesA_CAP__03	GenusA_speciesA_CAP
GenusA_speciesA_CAP__04	GenusA_speciesA_CAP
GenusA_speciesA_CAP__05	GenusA_speciesA_CAP
GenusB_speciesB_CAP	GenusB_speciesB_CAP
GenusB_speciesB_CAP__00	GenusB_speciesB_CAP
GenusB_speciesB_CAP__01	GenusB_speciesB_CAP
GenusB_speciesB_CAP__02	GenusB_speciesB_CAP
GenusC_speciesC_CAP	GenusC_speciesC_CAP
GenusC_speciesC_CAP__00	GenusC_speciesC_CAP
GenusC_speciesC_CAP__01	GenusC_speciesC_CAP
GenusC_speciesC_CAP__02	GenusC_speciesC_CAP
GenusC_speciesC_CAP__03	GenusC_speciesC_CAP
GenusC_speciesC_CAP__04	GenusC_speciesC_CAP
GenusC_speciesC_CAP__05	GenusC_speciesC_CAP
GenusD_speciesD_CAP	GenusD_speciesD_CAP
GenusD_speciesD_CAP__00	GenusD_speciesD_CAP
GenusD_speciesD_CAP__01	GenusD_speciesD_CAP
GenusD_speciesD_CAP__02	GenusD_speciesD_CAP
GenusD_speciesD_CAP__03	GenusD_speciesD_CAP
GenusD_speciesD_CAP__04	GenusD_speciesD_CAP
GenusD_speciesD_CAP__05	GenusD_speciesD_CAP
```
{{% /expand %}}
___
### 20. **`captus-assembly_align.report.html`**
This is the final [Aligment report]({{< ref "assembly/align/report">}}), summarizing alignment statistics across all processing stages, marker types, and formats.
___
### 21. **`captus-assembly_align.log`**
This is the log from `Captus`, it contains the command used and all the information shown during the run. Even if the option `--show_more` was disabled, the log will contain all the extra detailed information that was hidden during the run.

___
Created by [Edgardo M. Ortiz]({{< ref "../../more/credits/#edgardo-m-ortiz">}}) (06.08.2021)  
Last modified by [Edgardo M. Ortiz]({{< ref "../../more/credits/#edgardo-m-ortiz">}}) (31.05.2022)