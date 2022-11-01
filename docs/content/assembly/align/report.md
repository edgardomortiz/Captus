+++
title = "HTML Report"
weight = 15
pre = '<i class="fas fa-chart-bar"></i> '
plotly = true
+++

## Concept

---
This `align` module generates several sets of alignments that are ready-to-use in popular phylogenetic tree inference programs (e.g., [`IQ-TREE`](http://www.iqtree.org), [`RAxML`](https://cme.h-its.org/exelixis/web/software/raxml), [`ASTRAL`](https://github.com/smirarab/ASTRAL)).
Each alignment set differs from one another in the following four respects: 1) whether they are trimmed, 2) which [paralog filter <i class="fas fa-question-circle fa-sm"></i>]({{< relref "assembly/align/options#paralog-filtering" >}}) is applied, 3) whether they contain reference sequences, and 4) in which [formats <i class="fas fa-question-circle fa-sm"></i>]({{< relref "assembly/align/options#-f---formats" >}}).
Thus, **it is important to understand the differences between each alignment set and carefully evaluate their quality in order to decide which alignment set to use for subsequent analyses**.

Open the report `captus-assembly_align.report.html` with your browser (internet connection required) to explore and compare general alignment statistics for each locus and each sample!
{{% notice style="tip" title="Tips" %}}

- The entire report is based on data stored in the following two files:
  - [`captus-assembly_align.alignments.tsv`]({{< relref "assembly/align/output#17-captus-assembly_alignalignmentstsv" >}})
  - [`captus-assembly_align.samples.tsv`]({{< relref "assembly/align/output#18-captus-assembly_alignsamplestsv" >}})
- All tables and plots in the report are interactive powered by [`Plotly`](https://plotly.com/python).  
Visit the following sites once to take full advantage of its interactivity:

  - <https://plotly.com/chart-studio-help/getting-to-know-the-plotly-modebar>
  - <https://plotly.com/chart-studio-help/zoom-pan-hover-controls>
{{% /notice %}}

## Contents

---

### 1. Stats Comparison at Each Processing Step

This plot shows distributions of general alignment statistics at each processing step.

Features:

- Switch the `Marker Type` dropdown to change the marker type<br>(appeared only when you have more than one marker type).
- Switch the dropdown on the *x*-axis to change the variable to show.
- Click on the legend to toggle hide/show of each format.

{{< plotly json="/plotly/alignment_report_ridgeline.json" height="700px" >}}
{{% expand "Description of each processing step" %}}
Depending on [`--filter_method`]({{< relref "assembly/align/options#--filter_method" >}}) argument, you will have up to 12 processing steps as follows:
|Processing step (Path to alignments)|Trimmed|Paralog filter|With references|
|-|-|-|-|
|**02_untrimmed/01_unfiltered_w_refs**|No|None|**Yes**|
|**02_untrimmed/02_naive_w_refs**|No|**Naive**|**Yes**|
|**02_untrimmed/03_informed_w_refs**|No|**Informed**|**Yes**|
|**02_untrimmed/01_unfiltered**|No|None|No|
|**02_untrimmed/02_naive**|No|**Naive**|No|
|**02_untrimmed/03_informed**|No|**Informed**|No|
|**03_trimmed/01_unfiltered_w_refs**|**Yes**|None|**Yes**|
|**03_trimmed/02_naive_w_refs**|**Yes**|**Naive**|**Yes**|
|**03_trimmed/03_informed_w_refs**|**Yes**|**Informed**|**Yes**|
|**03_trimmed/01_unfiltered**|**Yes**|None|No|
|**03_trimmed/02_naive**|**Yes**|**Naive**|No|
|**03_trimmed/03_informed**|**Yes**|**Informed**|No|

For more explanations, read [<i class="fab fa-readme"></i> Output Files]({{< relref "assembly/align/output" >}}).
{{% /expand %}}
{{% expand "Description of each variable" %}}
|Variable|Description|Unit|
|-|-|-|
|**Sequences**|Number of sequences in the alignment|-|
|**Samples**|Number of samples in the alignment|-|
|**Sequences Per Sample**|= `Sequences` / `Samples`|-|
|**Alignment Length**|Length of the alignment|aa/bp|
|**Informative Sites**|Number of parsimony-informative sites that have at least two different characters and at least two of which appear in at least two sequences|-|
|**Informativeness**|= (`Informative Sites` / `Alignment Length`) * 100|%|
|**Uninformative Sites**|= `Alignment Length` - `Informative Sites`<br>= `Constant Sites` + `Singleton Sites`|-|
|**Constant Sites**|Number of invariant sites in the alignment|-|
|**Singleton Sites**|Number of variable sites where one character appears in multiple sequences while other characters appear in only one sequence|-|
|**Patterns**|Number of unique sites that have different character configurations|-|
|**Mean Pairwise Identity**|Mean pairwise sequence identity in the alignment|%|
|**Missingness**|Proportion of `-`, `N`, `X`, `?`, `.`, and `~` in the alignment|%|
|**GC Content**|GC content of the alignment (inapplicable to `AA` format)|%|
|**GC Content at 1st Codon Position**|GC content at 1st codon position in the alignment<br>(only applicable to `NT` format)|%|
|**GC Content at 2nd Codon Position**|GC content at 2nd codon position in the alignment<br>(only applicable to `NT` format)|%|
|**GC Content at 3rd Codon Position**|GC content at 3rd codon position in the alignment<br>(only applicable to `NT` format)|%|

For more explanation about sites and patterns, see IQ-TREE's FAQ: [<i class="fas fa-question-circle"></i> What are the differences between alignment columns/sites and patterns?](http://www.iqtree.org/doc/Frequently-Asked-Questions#how-does-iq-tree-treat-identical-sequences)
{{% /expand %}}

---

### 2. Bivariate Relationships and Distributions

This plot shows general alignment statistics for each alignment (locus).  
When your result contains more than one marker type, the report will include separate plots for each marker type.

Features:

- Switch the `Processing Step` dropdown to change the processing step to show the statistics.
- Switch the dropdowns on the *x*- and *y*- axes to change variables to plot on each axis.
- Click on the legend to toggle hide/show of each format.

{{< plotly json="/plotly/alignment_report_scatter_loci.json" height="700px" >}}
{{% expand "Description of each processing step" %}}
Depending on [`--filter_method`]({{< relref "assembly/align/options#--filter_method" >}}) argument, you will have up to 12 processing steps as follows:
|Processing step (Path to alignments)|Trimmed|Paralog filter|With references|
|-|-|-|-|
|**02_untrimmed/01_unfiltered_w_refs**|No|None|**Yes**|
|**02_untrimmed/02_naive_w_refs**|No|**Naive**|**Yes**|
|**02_untrimmed/03_informed_w_refs**|No|**Informed**|**Yes**|
|**02_untrimmed/01_unfiltered**|No|None|No|
|**02_untrimmed/02_naive**|No|**Naive**|No|
|**02_untrimmed/03_informed**|No|**Informed**|No|
|**03_trimmed/01_unfiltered_w_refs**|**Yes**|None|**Yes**|
|**03_trimmed/02_naive_w_refs**|**Yes**|**Naive**|**Yes**|
|**03_trimmed/03_informed_w_refs**|**Yes**|**Informed**|**Yes**|
|**03_trimmed/01_unfiltered**|**Yes**|None|No|
|**03_trimmed/02_naive**|**Yes**|**Naive**|No|
|**03_trimmed/03_informed**|**Yes**|**Informed**|No|

For more explanations, read [<i class="fab fa-readme"></i> Output Files]({{< relref "assembly/align/output" >}}).
{{% /expand %}}
{{% expand "Description of each variable" %}}
|Variable|Description|Unit|
|-|-|-|
|**Sequences**|Number of sequences in the alignment|-|
|**Samples**|Number of samples in the alignment|-|
|**Sequences Per Sample**|= `Sequences` / `Samples`|-|
|**Alignment Length**|Length of the alignment|aa/bp|
|**Informative Sites**|Number of parsimony-informative sites that have at least two different characters and at least two of which appear in at least two sequences|-|
|**Informativeness**|= (`Informative Sites` / `Alignment Length`) * 100|%|
|**Uninformative Sites**|= `Alignment Length` - `Informative Sites`<br>= `Constant Sites` + `Singleton Sites`|-|
|**Constant Sites**|Number of invariant sites in the alignment|-|
|**Singleton Sites**|Number of variable sites where one character appears in multiple sequences while other characters appear in only one sequence|-|
|**Patterns**|Number of unique sites that have different character configurations|-|
|**Mean Pairwise Identity**|Mean pairwise sequence identity in the alignment|%|
|**Missingness**|Proportion of `-`, `N`, `X`, `?`, `.`, and `~` in the alignment|%|
|**GC Content**|GC content of the alignment (inapplicable to `AA` format)|%|
|**GC Content at 1st Codon Position**|GC content at 1st codon position in the alignment<br>(only applicable to `NT` format)|%|
|**GC Content at 2nd Codon Position**|GC content at 2nd codon position in the alignment<br>(only applicable to `NT` format)|%|
|**GC Content at 3rd Codon Position**|GC content at 3rd codon position in the alignment<br>(only applicable to `NT` format)|%|

For more explanation about sites and patterns, see IQ-TREE's FAQ: [<i class="fas fa-question-circle"></i> What are the differences between alignment columns/sites and patterns?](http://www.iqtree.org/doc/Frequently-Asked-Questions#how-does-iq-tree-treat-identical-sequences)
{{% /expand %}}

---

### 3. Stats Per Sample

This plot shows general alignment statistics for each sample.  
When your result contains more than one marker type, the report will include separate plots for each marker type.

Features:

- Switch the `Sort Samples by` dropdown to re-sort the *x*-axis by sample name or mean of the variable.
- Switch the dropdown on the *y*-axis to change the variable to show.
- Click on the legend to toggle hide/show of each data series.

{{< plotly json="/plotly/alignment_report_scatter_samples.json" height="700px" >}}
{{% expand "Description of each variable" %}}
|Variable|Description|Unit|
|-|-|-|
|**Number of Loci**|Number of alignments (loci) containing the sample|-|
|**Mean Ungapped Length**|Mean sequence length of the sample excluding gaps (`-`)|aa/bp|
|**Total Ungapped Length**|Cumulative sequence length of the sample excluding gaps (`-`)|aa/bp|
|**Mean Gaps**|Mean length of internal gaps (`-`)|-|
|**Total Gaps**|Cumulative length of internal gaps (`-`)|-|
|**Mean Ambiguities**|Mean count of ambiguous characters of the sample|-|
|**Mean GC Content**|Mean GC content of the sample (inapplicable to `AA` format)|%|
|**Mean GC Content at 1st Codon Position**|Mean GC content at 1st codon position of the sample (only applicable to `NT` format)|%|
|**Mean GC Content at 2nd Codon Position**|Mean GC content at 2nd codon position of the sample (only applicable to `NT` format)|%|
|**Mean GC Content at 3rd Codon Position**|Mean GC content at 3rd codon position of the sample (only applicable to `NT` format)|%|
|**Mean Copies**|Mean number of sequences per alignment (always `1` for alignments with paralog filter applied)|%|
{{% /expand %}}

---
Created by [Gentaro Shigita]({{< ref "../../more/credits/#gentaro-shigita">}}) (11.08.2021)  
Last modified by [Gentaro Shigita]({{< ref "../../more/credits/#gentaro-shigita">}}) (17.10.2022)
