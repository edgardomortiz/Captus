+++
title = "HTML Report"
weight = 15
pre = '<i class="fas fa-chart-bar"></i> '
plotly = true
+++

## Concept

---
**No successful marker extractions can be achieved without successful assemblies**.
Therefore, it is advisable to control your impatience and carefully evaluate your assembly results before moving on to the next step.
(Even though `Captus` provides some presets for different data types, it would be the most prudent way to repeat this step several times with different parameter combinations and find optimal assembler settings for your dataset.)
A HTML report, `captus-assembly_assemble.report.html` can help you with that.
Just open this file in your browser (Microsoft Edge, Google Chrome, Mozilla Firefox, Safari, etc., internet connection required) and you can compare general assembly statistics across all your samples!
{{% notice info %}}

- All original data for the report is stored in `captus-assembly_assemble.stats.tsv`.
- Since all tables and plots in the report are created using [`Plotly`](https://plotly.com/python), you can use some interactive functions such as zoom in/out, pan, hover, and download plot as a SVG.
For more information, please visit the following sites:

  - <https://plotly.com/chart-studio-help/zoom-pan-hover-controls>
  - <https://plotly.com/chart-studio-help/getting-to-know-the-plotly-modebar>
{{% /notice %}}

## Contents

---

### 1. Summary Table

This table shows the general assembly statistics for each sample.  

##### Features:

- The cells are colored according to their relative value within each column (green = high; pink = low).
- By switching `Sort by` dropdown, you can re-sort the table in descending order according to the values in each column (by default, the table is sorted alpha-numerically by sample name).

{{< plotly json="/plotly/assemble_report_summary_table.json" height="250px" >}}
{{% expand "Description of each column" %}}
|Column|Description|
|-|-|
|**Sample**|Sample name|
|**Total Length (bp)**|Total length of all contigs|
|**Number of Contigs**|Total number of contigs|
|**Mean Length (bp)**|Mean contig length|
|**Contig N50 (bp)**|Weighted average of contig length such that 50% of total assembly length is represented by contigs equal to or longer than this value|
|**Longest Contig (bp)**|Longest contig length|
|**Shortest Contig (bp)**|Shortest contig length|
|**≥1kbp Contig (%)**|Percentage of contigs ≥ 1 kbp in total number of contigs|
|**GC Content (%)**|Mean GC content of all contigs|
|**Mean Depth (x)**|Mean depth of all contigs|
{{% /expand %}}

---

### 2. Visual Stats

In addition to the above statistics, this plot shows some more detailed statistics as a simple or stacked bar plot.

##### Feature:

- By switching the dropdown at the X-axis, you can change the variable to show.
- Samples are sorted in descending order by each variable.
{{< plotly json="/plotly/assemble_report_visual_stats.json" height="260px" >}}

{{% expand "Description of each variable" %}}
|Variable|Description|
|-|-|
|**Total Length (bp)**|Total length of all contigs|
|**Number of Contigs**|Total number of contigs|
|**Mean Length (bp)**|Mean contig length|
|**Median Length (bp)**|Median contig length|
|**Contig N50 (bp)**|Weighted average of contig length such that 50% of total assembly length is represented by contigs equal to or longer than this value|
|**Longest Contig (bp)**|Longest contig length|
|**Shortest Contig (bp)**|Shortest contig length|
|**Contig Breakdown by Length (%)**|Percentage of contigs ≥ 1, 2, 5, and 10 kbp in the total number of contigs|
|**Length Breakdown by Contig Length (%)**|Percentage of contigs ≥ 1, 2, 5, and 10 kbp in the total length of all contigs|
|**GC Content (%)**|Mean GC content of all contigs|
|**Mean Depth (x)**|Mean depth of all contigs|
|**Contig Breakdown by Depth (%)**|Percentage of contigs with ≥ 1, 2, 5, and 10x depth in the total number of contigs|
{{% /expand %}}

---
Created by [Gentaro Shigita]({{< ref "../../more/credits/#gentaro-shigita">}}) (11.08.2021)  
Last modified by [Gentaro Shigita]({{< ref "../../more/credits/#gentaro-shigita">}}) (31.05.2022)
