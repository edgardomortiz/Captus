+++
title = "HTML Report"
weight = 15
pre = '<i class="fas fa-chart-bar"></i> '
plotly = true
+++

## Concept

---
**Proper cleaning is the first step to perform proper analyses** on high-throughput sequencing data.
To assess the quality of raw reads and how it is improved by the cleaning, the `clean` module internally runs the famous quality check program, [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc), or its faster emulator, [`Falco`](https://github.com/smithlabcode/falco), on the reads before and after cleaning.
Although both programs generate informative reports, they are in separate files for each sample, each read direction (for paired-end), and before and after cleaning.
This makes it tedious to review every report, and can lead to overlook some serious problems, such as residual low-quality bases or adaptor sequences, contamination of different samples, and improper setting of cleaning parameters.

`Captus` summarizes the information in those disparate reports into a single HTML file. All you need to do is open `captus-assembly_clean.report.html` with your browser (internet connection required) to get a quick overview on all your samples, both reads (for paired-end), and before and after cleaning!

{{% notice tip %}}

- The entire report is based on tables stored in the [`03_qc_extras`]({{< relref "assembly/clean/output#8-03_qc_extras" >}}) directory.
- All tables and plots in the report are interactive powered by [`Plotly`](https://plotly.com/python).  
Visit the following sites once to take full advantage of its interactivity:

  - <https://plotly.com/chart-studio-help/getting-to-know-the-plotly-modebar>
  - <https://plotly.com/chart-studio-help/zoom-pan-hover-controls>
{{% /notice %}}

## Contents

---
The report comprises the following nine sections:

1. [Summary Table](#1-summary-table)
2. [Stats on Reads/Bases](#2-stats-on-readsbases)
3. [Per Base Quality](#3-per-base-quality)
4. [Per Read Quality](#4-per-read-quality)
5. [Read Length Distribution](#5-read-length-distribution)
6. [Per Base Nucleotide Content](#6-per-base-nucleotide-content)
7. [Per Read GC Content](#7-per-read-gc-content)
8. [Sequence Duplication Level](#8-sequence-duplication-level)
9. [Adaptor Content](#9-adaptor-content)

A brief description and interactive example of each section is given below.  
By switching the tabs at the top of each plot, you can compare the plot produced by `Captus` with the corresponding plot from [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc).

---

### 1. Summary Table

This table shows general cleaning statistics for each sample.

Features:

- Switch the `Sort by` dropdown to re-sort the table by any column value.
- Cells are color-coded according to value (high = green; low = pink).

{{< tabs groupId="Summary Table" >}}
{{% tab name="Captus" %}}
{{< plotly json="/captus.docs/plotly/cleaning_report_summary_table.json" height="230px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
![fastqc_basic_statistics](/captus.docs/images/fastqc_basic_statistics.png?height=200px)
{{% /tab %}}
{{< /tabs >}}

{{% expand "Description of each column" %}}
|Column|Description|Unit|
|-|-|-|
|**Sample**|Sample name|-|
|**Input Reads**|Number of reads before cleaning|-|
|**Input Bases**|Number of bases before cleaning|bp|
|**Output Reads**|Number of reads passed cleaning|-|
|**Output Reads%**|= (`Output Reads` / `Input Reads`) * 100|%|
|**Output Bases**|Number of bases passed cleaning|bp|
|**Output Bases%**|= (`Output Bases` / `Input Bases`) * 100|%|
|**Mean Read Length%**|= (Mean read length after cleaning / Mean read length before cleaning) * 100|%|
|**≥Q20 Reads%**|Percentage of reads with mean Phred quality score over 20 after cleaning|%|
|**≥Q30 Reads%**|Percentage of reads with mean Phred quality score over 30 after cleaning|%|
|**GC%**|Mean GC content in the reads after cleaning|%|
|**Adapter%**|Percentage of reads containing adaptor sequences before cleaning|%|
{{% /expand %}}

---

### 2. Stats on Reads/Bases

`Captus` cleans reads through two consecutive rounds of adaptor trimming (`Round1`, `Round2`) followed by quality trimming and filtering.  
This plot shows changes in the number of reads (left panel) and bases (right panel) at each step of the cleaning process.

Features:

- Switch the buttons at the top to choose whether to show counts or percentages.
- Samples are sorted by the number or percentage of bases passed cleaning.
- Click on the legend to toggle hide/show of each data series.

{{< tabs groupId="Stats on Reads/Bases" >}}
{{% tab name="Captus" %}}
{{< plotly json="/captus.docs/plotly/cleaning_report_stats_on_reads_bases.json" height="240px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
There is no corresponding plot.
{{% /tab %}}
{{< /tabs >}}

---

### 3. Per Base Quality

This plot shows the range of Phred quality score at each position in the reads before and after cleaning.  
For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html).  

Feature:

- Switch the dropdown at the top to change the variable to show, these variables represent the elements of the boxplots in the `FastQC` report.

{{< tabs groupId="Per Base Quality" >}}
{{% tab name="Captus" %}}
{{< plotly json="/captus.docs/plotly/cleaning_report_per_base_quality.json" height="300px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
![fastqc_per_base_quality](/captus.docs/images/fastqc_per_base_quality.png?height=300px)
{{% /tab %}}
{{< /tabs >}}

---

### 4. Per Read Quality

This plot shows the distribution of mean Phred quality score for each read before and after cleaning.  
For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/3%20Per%20Sequence%20Quality%20Scores.html).
{{< tabs groupId="Per Read Quality" >}}
{{% tab name="Captus" %}}
{{< plotly json="/captus.docs/plotly/cleaning_report_per_read_quality.json" height="300px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
![fastqc_per_sequence_quality](/captus.docs/images/fastqc_per_sequence_quality.png?height=300px)
{{% /tab %}}
{{< /tabs >}}

---

### 5. Read Length Distribution

This plot shows the distribution of read lengths before and after cleaning.  
For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/7%20Sequence%20Length%20Distribution.html).
{{< tabs groupId="Read Length Distribution" >}}
{{% tab name="Captus" %}}
{{< plotly json="/captus.docs/plotly/cleaning_report_seq_len_dist.json" height="300px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
![fastqc_sequence_length_distribution](/captus.docs/images/fastqc_sequence_length_distribution.png?height=300px)
{{% /tab %}}
{{< /tabs >}}

---

### 6. Per Base Nucleotide Content

This plot shows the composition of each nucleotide (A, T, G, C) at each position in the reads before and after cleaning.  
If a particular nucleotide is overrepresented at a certain position in the reads, you will see the color corresponding to that nucleotide; otherwise, the plot will be a uniform grayish color.  
For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html).
{{< tabs groupId="Per Base Nucleotide Content" >}}
{{% tab name="Captus" %}}
{{< plotly json="/captus.docs/plotly/cleaning_report_per_base_seq_content.json" height="300px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
![fastqc_per_base_sequence_content](/captus.docs/images/fastqc_per_base_sequence_content.png?height=300px)
{{% /tab %}}
{{< /tabs >}}

---

### 7. Per Read GC Content

This plot shows the frequency of GC content in the reads before and after cleaning.  
Broader or bimodal peaks may indicate contamination with DNA from different organisms.  
For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html).
{{< tabs groupId="Per Read GC Content" >}}
{{% tab name="Captus" %}}
{{< plotly json="/captus.docs/plotly/cleaning_report_per_seq_gc_content.json" height="300px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
![fastqc_per_sequence_gc_content](/captus.docs/images/fastqc_per_sequence_gc_content.png?height=300px)
{{% /tab %}}
{{< /tabs >}}

---

### 8. Sequence Duplication Level

This plot shows the percentage of sequences with different degrees of duplication before and after cleaning.  
For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html).  

Feature:

- Click on the legend to toggle hide/show of each data series.

{{< tabs groupId="Sequence Duplication Level" >}}
{{% tab name="Captus" %}}
{{< plotly json="/captus.docs/plotly/cleaning_report_seq_dup_level.json" height="300px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
![fastqc_duplication_levels](/captus.docs/images/fastqc_duplication_levels.png?height=300px)
{{% /tab %}}
{{< /tabs >}}

---

### 9. Adaptor Content

This plot shows the cumulative adaptor content at each position in the reads before and after cleaning.  
For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html).
{{< tabs groupId="Adaptor Content" >}}
{{% tab name="Captus" %}}
{{< plotly json="/captus.docs/plotly/cleaning_report_adaptor_content.json" height="300px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
![fastqc_adaptor_content](/captus.docs/images/fastqc_adaptor_content.png?height=300px)
{{% /tab %}}
{{< /tabs >}}

---
Created by [Gentaro Shigita]({{< ref "../../more/credits/#gentaro-shigita">}}) (11.08.2021)  
Last modified by [Gentaro Shigita]({{< ref "../../more/credits/#gentaro-shigita">}}) (22.02.2023)
