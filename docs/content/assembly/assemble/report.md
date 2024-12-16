+++
title = "HTML Report"
weight = 15
pre = '<i class="fas fa-chart-bar"></i> '
plotly = true
+++

## Concept

---
**No successful marker extractions can be achieved without successful assemblies**. Even though this `assemble` module offers presets tuned for different data types, it is recommendable to repeat this step some times with different parameters to find optimal settings for your own data.

`Captus` assists you in this tedious process by automatically generating a useful report for assembly evaluation.  
Just open `captus-assembly_assemble.report.html` with your browser (internet connection required) to get general assembly statistics across all your samples!
{{% notice tip %}}

- The entire report is based on data stored in the following three files:
  - [`captus-assembly_assemble.assembly_stats.tsv`]({{< relref "assembly/assemble/output#9-captus-assembly_assemblestatstsv" >}})
  - [`captus-assembly_assemble.depth_stats.tsv`]({{< relref "assembly/assemble/output" >}})
  - [`captus-assembly_assemble.length_stats.tsv`]({{< relref "assembly/assemble/output" >}})
- All tables and plots in the report are interactive powered by [`Plotly`](https://plotly.com/python).  
Visit the following sites once to take full advantage of its interactivity:

  - <https://plotly.com/chart-studio-help/getting-to-know-the-plotly-modebar>
  - <https://plotly.com/chart-studio-help/zoom-pan-hover-controls>
{{% /notice %}}

## Contents

---

### 1. Summary Table

This table shows the general assembly statistics for each sample.  
All values shown in this table are calculated **after** filtering by GC content ([--max_contig_gc]({{< relref "assembly/assemble/options#--max_contig_gc" >}})) and/or depth ([--min_contig_depth]({{< relref "assembly/assemble/options#--min_contig_depth" >}})).

Features:

- Switch the `Sort by` dropdown to re-sort the table by any column value.
- Cells are color-coded according to value (high = green; low = pink).

{{< plotly json="/captus.docs/plotly/assemble_report_summary_table.json" height="230px" >}}
{{% expand "Description of each column" %}}
|Column|Description|Unit|
|-|-|-|
|**Sample**|Sample name|-|
|**#Contigs**|Number of contigs|-|
|**Total Length**|Total length of all contigs|bp|
|**Shortest Contig**|Shortest contig length|bp|
|**Longest Contig**|Longest contig length|bp|
|**N50**|Weighted average of contig lengths that 50% of total assembly length consists of contigs over this length|bp|
|**L50**|Least number of contigs that contain 50% of total assembly length|-|
|**GC Content**|Overall GC content|%|
|**Mean Depth**|Mean contig depth|x|
{{% /expand %}}

---

### 2. Visual Stats

In addition to the general statistics shown in the `Summary Table`, this plot shows more detailed statistics before and after filtering for each sample as a bar graph.  

Features:

- Switch the dropdown at the *x*-axis to change the variable to show.
- Click on the legend to toggle hide/show of each data series (only applicable to some variables).

{{< plotly json="/captus.docs/plotly/assemble_report_visual_stats.json" height="300px" >}}

{{% expand "Description of each variable" %}}
|Variable<sup>a</sup>|Description|Unit|
|-|-|-|
|**#Contigs**|Number of contigs|-|
|**Total Length**|Total length of all contigs|bp|
|**Shortest Contig**|Shortest contig length|bp|
|**Longest Contig**|Longest contig length|bp|
|**N50**|Weighted average of contig lengths that 50% of total assembly length consists of contigs over this length|bp|
|**N75**|Weighted average of contig lengths that 75% of total assembly length consists of contigs over this length|bp|
|**L50**|Least number of contigs that contain 50% of total assembly length|-|
|**L75**|Least number of contigs that contain 75% of total assembly length|-|
|**Mean Length**|Mean contig length|bp|
|**Median Length**|Median contig length|bp|
|**Contig Breakdown by Length**|Percentage of contigs over 1, 2, 5, 10, 20, and 50 kbp in total number of contigs|%|
|**Length Breakdown by Contig Length**|Percentage of contigs over 1, 2, 5, 10, 20, and 50 kbp in total length of contigs|%|
|**GC Content**|Overall GC content|%|
|**Mean Depth**|Mean contig depth|x|
|**Median Depth**|Median contig depth|x|
{{% /expand %}}

---

### 3. Length Distribution

This plot shows the distribution of contig lengths before and after filtering for each sample as a heatmap.

Feature:

- Switch the `Variable` dropdown at the colorscale to change the variable to show.

{{< plotly json="/captus.docs/plotly/assemble_report_length_distribution.json" height="300px" >}}

{{% expand "Description of each variable" %}}
|Variable|Description|Unit|
|-|-|-|
|**Length**|Total contig length in the bin|bp|
|**Fraction**|(Total contig length in the bin) / (Total assembly length) * 100|%|
|**#Contigs**|Number of contigs in the bin|-|
{{% /expand %}}

---

### 4. Depth Distribution

This plot shows the distribution of contig depths before and after filtering for each sample as a heatmap.

Feature:

- Switch the `Variable` dropdown at the colorscale to change the variable to show.

{{< plotly json="/captus.docs/plotly/assemble_report_depth_distribution.json" height="300px" >}}

{{% expand "Description of each variable" %}}
|Variable|Description|Unit|
|-|-|-|
|**Length**|Total contig length in the bin|bp|
|**Fraction**|(Total contig length in the bin) / (Total assembly length) * 100|%|
|**#Contigs**|Number of contigs in the bin|-|
{{% /expand %}}

---
Created by [Gentaro Shigita]({{< ref "../../more/credits/#gentaro-shigita">}}) (11.08.2021)  
Last modified by [Gentaro Shigita]({{< ref "../../more/credits/#gentaro-shigita">}}) (16.12.2024)
