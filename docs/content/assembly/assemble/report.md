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

- The entire report is based on data stored in [`captus-assembly_assemble.stats.tsv`]({{< relref "assembly/assemble/output#9-captus-assembly_assemblestatstsv" >}}).
- All tables and plots in the report are interactive powered by [`Plotly`](https://plotly.com/python).  
Visit the following sites once to take full advantage of its interactivity:

  - <https://plotly.com/chart-studio-help/getting-to-know-the-plotly-modebar>
  - <https://plotly.com/chart-studio-help/zoom-pan-hover-controls>
{{% /notice %}}

## Contents

---

### 1. Summary Table

This table shows the general assembly statistics for each sample.  

Features:

- Switch the `Sort by` dropdown to re-sort the table by any column value.
- Cells are color-coded according to value (high = green; low = pink).

{{< plotly json="/plotly/assemble_report_summary_table.json" height="250px" >}}
{{% expand "Description of each column" %}}
|Column<sup>a</sup>|Description|Unit|
|-|-|-|
|**Sample**|Sample name|-|
|**Total Length Passed**|Cumulative length of all contigs passed GC content filter|bp|
|**Total Length Removed**|Cumulative length of all contigs removed by GC content filter<br>(shown only when [`--max_contig_gc`]({{< relref "assembly/assemble/options#--max_contig_gc" >}}) argument is enabled)|bp|
|**Number of Contigs Passed**|Number of contigs passed GC content filter|-|
|**Number of Contigs Removed**|Number of contigs removed by GC content filter<br>(shown only when [`--max_contig_gc`]({{< relref "assembly/assemble/options#--max_contig_gc" >}}) argument is enabled)|-|
|**\*Mean Length**|Mean contig length|bp|
|**\*Contig N50**|Weighted average of contig lengths that 50% of total assembly length consists of contigs over this length|bp|
|**\*Longest Contig**|Longest contig length|bp|
|**\*Shortest Contig**|Shortest contig length|bp|
|**\*≥1 kbp Contigs**|Percentage of contigs over 1 kbp in total number of contigs|%|
|**GC Content Passed**|Mean GC content of contigs passed GC content filter|%|
|**GC Content Removed**|Mean GC content of contigs removed by GC content filter<br>(shown only when [`--max_contig_gc`]({{< relref "assembly/assemble/options#--max_contig_gc" >}}) argument is enabled)|%|
|**\*Mean Depth**|Mean contig depth|x|

<sup>a</sup> Statistics denoted by `*` are calculated only from contigs passed GC content filter.
{{% /expand %}}

---

### 2. Visual Stats

This plot shows the statistics shown in the table above in addition to more detailed statistics as a bar graph.  

Features:

- Switch the dropdown at the *y*-axis to switch the variable to show.
- Click on the legend to toggle hide/show of each data series.
{{< plotly json="/plotly/assemble_report_visual_stats.json" height="500px" >}}

{{% expand "Description of each variable" %}}
|Variable<sup>a</sup>|Description|Unit|
|-|-|-|
|**Total Length**|Cumulative length of all contigs|bp|
|**Number of Contigs**|Number of contigs|-|
|**Mean Length**|Mean contig length|bp|
|**\*Median Length**|Median contig length|bp|
|**\*Contig N50**|Weighted average of contig lengths that 50% of total assembly length consists of contigs over this length|bp|
|**\*Longest Contig**|Longest contig length|bp|
|**\*Shortest Contig**|Shortest contig length|bp|
|**\*Contig Breakdown by Length**|Percentage of contigs over 1, 2, 5, and 10 kbp in total number of contigs|%|
|**\*Length Breakdown by Contig Length**|Percentage of contigs over 1, 2, 5, and 10 kbp in total length of contigs|%|
|**GC Content**|Mean GC content of contigs|%|
|**\*Mean Depth**|Mean contig depth|x|
|**\*Contig Breakdown by Depth**|Percentage of contigs with over 1, 2, 5, and 10 x depth in total number of contigs|%|

<sup>a</sup> Statistics denoted by `*` are calculated only from contigs passed GC content filter.
{{% /expand %}}

---
Created by [Gentaro Shigita]({{< ref "../../more/credits/#gentaro-shigita">}}) (11.08.2021)  
Last modified by [Gentaro Shigita]({{< ref "../../more/credits/#gentaro-shigita">}}) (17.10.2022)
