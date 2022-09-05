+++
title = "HTML Report"
weight = 15
pre = '<i class="fas fa-chart-bar"></i> '
plotly = true
+++

## Concept

---
**Proper cleaning is the first step to perform proper analyses** on high-throughput sequencing data.
To evaluate the quality of the input reads and how it was improved by the cleaning, `Captus` internally runs the well-known quality check program, [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc), or its faster emulator, [`Falco`](https://github.com/smithlabcode/falco), on the reads before and after cleaning.
Although both programs generate nice visual reports, those reports are generated as separate files for each sample, for each read direction (for paired-end), and before and after cleaning.
This makes it tedious to check every report, and may lead to overlooking some serious problems, such as residual low-quality bases and adaptor sequences, contamination of different samples, and inappropriate settings of cleaning parameters.

`Captus` summarizes the information in those disparate reports into a single HTML file, `captus-assembly_clean.report.html`. All you need to do is open this single file in your browser (e.g., Microsoft Edge, Google Chrome, Mozilla Firefox, Safari, etc.; internet connection required) to get a quick overview on all your samples, both reads (for paired-end), and before and after cleaning!

{{% notice tip %}}
Since all tables and plots in the report are created using [`Plotly`](https://plotly.com/python), you can use some useful functions such as zoom in/out, pan, hover, and download plots as an image (SVG format).
For more information, please visit the following sites:

- <https://plotly.com/chart-studio-help/zoom-pan-hover-controls>
- <https://plotly.com/chart-studio-help/getting-to-know-the-plotly-modebar>

{{% /notice %}}

## Contents

---
The report comprises the following nine plots:

1. [Summary Table](#1-summary-table)
2. [Stats on Reads/Bases](#2-stats-on-readsbases)
3. [Per Base Quality](#3-per-base-quality)
4. [Per Read Quality](#4-per-read-quality)
5. [Read Length Distribution](#5-read-length-distribution)
6. [Per Base Nucleotide Content](#6-per-base-nucleotide-content)
7. [Per Read GC Content](#7-per-read-gc-content)
8. [Sequence Duplication Level](#8-sequence-duplication-level)
9. [Adaptor Content](#9-adaptor-content)

A brief description and interactive example for each plot is given below.  
By switching the tabs at the top of each plot, you can compare the plot produced by `Captus` with the corresponding plot from [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc).

---

### 1. Summary Table

This table shows the general statistics for each sample.

##### Features:

- The cells are color-coded according to their relative value within each column (green = high; pink = low).
- By switching `Sort by` dropdown, you can re-sort the table in descending order according to the values in each column (by default, the table is sorted by sample name).

Original data for this table is stored in the files, `03_qc_extras/reads_bases.tsv`, `03_qc_extras/seq_len_dist.tsv`, `03_qc_extras/per_seq_qual_scores.tsv`, `03_qc_extras/per_seq_gc_content.tsv`, `03_qc_extras/adaptor_content.tsv`.

{{< tabs groupId="Summary Table" >}}
{{% tab name="Captus" %}}
{{< plotly json="/plotly/cleaning_report_summary_table.json" height="230px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
![fastqc_basic_statistics](/images/fastqc_basic_statistics.png?height=200px)
{{% /tab %}}
{{< /tabs >}}

{{% expand "Description of each column" %}}
|Column|Description|
|-|-|
|**Sample**|Sample name|
|**Input Reads**|Number of reads before cleaning|
|**Input Bases**|Number of bases before cleaning|
|**Output Reads**|Number of reads after cleaning|
|**Output Reads%**|Percentage of reads passed cleaning|
|**Output Bases**|Number of bases after cleaning|
|**Output Bases%**|Percentage of bases passed cleaning|
|**Mean Read Length%**|Percentage of mean read length after cleaning to the mean read length before cleaning|
|**≥Q20 Reads%**|Percentage of reads with mean Phred quality score ≥ 20 after cleaning|
|**≥Q30 Reads%**|Percentage of reads with mean Phred quality score ≥ 30 after cleaning|
|**GC%**|Mean GC content of reads after cleaning|
|**Adapter%**|Percentage of reads containing adapters before cleaning|
{{% /expand %}}

---

### 2. Stats on Reads/Bases

`Captus` cleans the reads through two consecutive rounds of adaptor trimming (`Round1`, `Round2`) followed by quality trimming and filtering.  
This plot shows the change in the number of reads (left panel) and bases (right panel) at each step of the cleaning process.

##### Features:

- By switching the buttons at the top, you can choose whether to show counts or percentages.
- Samples are sorted in descending order by number or percentage of bases passed cleaning.
- By clicking on the legend, you can toggle between showing and hiding each series.

Original data for this plot is stored in `03_qc_extras/reads_bases.tsv`.
{{< tabs groupId="Stats on Reads/Bases" >}}
{{% tab name="Captus" %}}
{{< plotly json="/plotly/cleaning_report_stats_on_reads_bases.json" height="240px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
There is no corresponding plot.
{{% /tab %}}
{{< /tabs >}}

---

### 3. Per Base Quality

This plot shows the range of Phred quality score at each position in the reads before and after cleaning.  
For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html).  

##### Features:

- By switching the dropdown at the top, you can change the variable to show, these variables represent the elements of the boxplots in the `FastQC` graph.

Original data for this plot is stored in `03_qc_extras/per_base_seq_qual.tsv`.
{{< tabs groupId="Per Base Quality" >}}
{{% tab name="Captus" %}}
{{< plotly json="/plotly/cleaning_report_per_base_quality.json" height="300px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
![fastqc_per_base_quality](/images/fastqc_per_base_quality.png?height=300px)
{{% /tab %}}
{{< /tabs >}}

---

### 4. Per Read Quality

This plot shows the distribution of the mean Phred quality score of the reads before and after cleaning.  
For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/3%20Per%20Sequence%20Quality%20Scores.html).  
Original data for this plot is stored in `03_qc_extras/per_seq_qual_scores.tsv`.
{{< tabs groupId="Per Read Quality" >}}
{{% tab name="Captus" %}}
{{< plotly json="/plotly/cleaning_report_per_read_quality.json" height="300px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
![fastqc_per_sequence_quality](/images/fastqc_per_sequence_quality.png?height=300px)
{{% /tab %}}
{{< /tabs >}}

---

### 5. Read Length Distribution

This plot shows the distribution of the read length before and after cleaning.  
For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/7%20Sequence%20Length%20Distribution.html).  
Original data for this plot is stored in `03_qc_extras/seq_len_dist.tsv`.
{{< tabs groupId="Read Length Distribution" >}}
{{% tab name="Captus" %}}
{{< plotly json="/plotly/cleaning_report_seq_len_dist.json" height="300px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
![fastqc_sequence_length_distribution](/images/fastqc_sequence_length_distribution.png?height=300px)
{{% /tab %}}
{{< /tabs >}}

---

### 6. Per Base Nucleotide Content

This plot shows the composition of the four nucleotides (A, T, G, C) at each position in the reads before and after cleaning.  
For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html).  
Original data for this plot is stored in `03_qc_extras/per_base_seq_content.tsv`.
{{< tabs groupId="Per Base Nucleotide Content" >}}
{{% tab name="Captus" %}}
{{< plotly json="/plotly/cleaning_report_per_base_seq_content.json" height="750px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
![fastqc_per_base_sequence_content](/images/fastqc_per_base_sequence_content.png?height=300px)
{{% /tab %}}
{{< /tabs >}}

---

### 7. Per Read GC Content

This plot shows the distribution of GC content in the reads before and after cleaning.
Broader or multiple peaks might be a sign of contamination with a different sample.  
For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html).  
Original data for this plot is stored in `03_qc_extras/per_seq_gc_content.tsv`.  
{{< tabs groupId="Per Read GC Content" >}}
{{% tab name="Captus" %}}
{{< plotly json="/plotly/cleaning_report_per_seq_gc_content.json" height="300px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
![fastqc_per_sequence_gc_content](/images/fastqc_per_sequence_gc_content.png?height=300px)
{{% /tab %}}
{{< /tabs >}}

---

### 8. Sequence Duplication Level

This plot shows the percentage of sequences with different degrees of duplication before and after cleaning.  
For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html).  

##### Features:

- By clicking on the legend, you can toggle between showing and hiding each series.

Original data for this plot is stored in `03_qc_extras/seq_dup_levels.tsv`.  
{{< tabs groupId="Sequence Duplication Level" >}}
{{% tab name="Captus" %}}
{{< plotly json="/plotly/cleaning_report_seq_dup_level.json" height="300px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
![fastqc_duplication_levels](/images/fastqc_duplication_levels.png?height=300px)
{{% /tab %}}
{{< /tabs >}}

---

### 9. Adaptor Content

This plot shows the cumulative adaptor content at each position in the reads before and after cleaning.  
For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html).  
Original data for this plot is stored in `03_qc_extras/adaptor_content.tsv`.
{{< tabs groupId="Adaptor Content" >}}
{{% tab name="Captus" %}}
{{< plotly json="/plotly/cleaning_report_adaptor_content.json" height="300px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
![fastqc_adaptor_content](/images/fastqc_adaptor_content.png?height=300px)
{{% /tab %}}
{{< /tabs >}}

---
Created by [Gentaro Shigita]({{< ref "../../more/credits/#gentaro-shigita">}}) (11.08.2021)  
Last modified by [Gentaro Shigita]({{< ref "../../more/credits/#gentaro-shigita">}}) (02.09.2022)
