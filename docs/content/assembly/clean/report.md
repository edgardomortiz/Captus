---
title: "HTML Report"
weight: 15
pre: '<i class="fas fa-chart-bar"></i> '
plotly: true
---

`captus-assembly_clean.report.html`
you can open the file with browsers (internet connection required)
[Plotly](https://plotly.com/python/)

Plotly zoom-in/out, pan, export as image
for details, please visit following sites
<https://plotly.com/chart-studio-help/getting-to-know-the-plotly-modebar/>
<https://plotly.com/chart-studio-help/zoom-pan-hover-controls/>

Captus runs [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) internally, 
This report enables you to have an overview of the quality of your data.

the report comprise one table and eight plots below.

- [1. *Summary Table*](#1-summary-table)
- [2. *Stats on Reads/Bases*](#2-stats-on-readsbases)
- [3. *Per Base Quality*](#3-per-base-quality)
- [4. *Per Read Quality*](#4-per-read-quality)
- [5. *Read Length Distribution*](#5-read-length-distribution)
- [6. *Per Base Nucleotide Content*](#6-per-base-nucleotide-content)
- [7. *Per Read GC Content*](#7-per-read-gc-content)
- [8. *Sequence Duplication Level*](#8-sequence-duplication-level)
- [9. *Adapter Content*](#9-adapter-content)

By switching tabs, you can compare with corresponding plot from [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

### 1. *Summary Table*

---
This table shows the general statistics.
{{< tabs >}}
{{% tab name="Captus" %}}
{{< plotly json="/plotly/cleaning_report_1_summary_table.json" height="225px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
![fastqc_per_base_qual](/images/fastqc_per_base_qual.png)
{{% /tab %}}
{{< /tabs >}}

Description of columns

- **Sample** = Sample name
- **Input Reads** = Number of reads before cleaning
- **Input Bases** = Number of bases before cleaning
- **Output Reads** = Number of reads after cleaning
- **Output Reads%** = Percentage of reads passed cleaning
- **Output Bases** = Number of bases after cleaning
- **Output Bases%** = Percentage of bases passed cleaning
- **Mean Read Length%** = Mean cleaned read length / Original read length * 100
- **≥Q20 Reads%** = Percentage of reads with mean Phred quality score ≥ 20 after cleaning
- **≥Q30 Reads%** = Percentage of reads with mean Phred quality score ≥ 30 after cleaning
- **GC%** = Mean GC content after cleaning

cell colors represent normalized values for each column (Green = high, Red = low).  
By switching `Sort by` dropdown at the top of the table, you can sort descending  
The original data of this table are stored in  

- `03_qc_extras/reads_bases.tsv`
- `seq_len_dist.tsv`
- `per_seq_qual_scores.tsv`
- `per_seq_gc_content.tsv`

### 2. *Stats on Reads/Bases*

---
This plot shows the changes on read and bases due to cleaning.  
{{< tabs >}}
{{% tab name="Captus" %}}
{{< plotly json="/plotly/cleaning_report_2_stats_on_reads_bases.json" height="240px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
There is no corresponding plot.
{{% /tab %}}
{{< /tabs >}}

The original data of this plot is stored in `03_qc_extras/reads_bases.tsv`

### 3. *Per Base Quality*

---
This plot shows the range of Phred quality score at each base position in the reads as a heatmap.
color scale shows 


- Mean
- 90<sub>th</sub> Percentile
- 75<sub>th</sub> Percentile
- 50<sub>th</sub> Percentile (= Median)
- 25<sub>th</sub> Percentile
- 10<sub>th</sub> Percentile

For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html)

The original data of this plot is stored in `03_qc_extras/per_base_seq_qual.tsv`
{{< tabs >}}
{{% tab name="Captus" %}}
{{< plotly json="/plotly/cleaning_report_3_per_base_quality.json" height="300px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
GenusC_speciesC_CAP After cleaning
![per_base_quality](/images/per_base_quality.png?height=300)
![per_base_quality](/images/per_base_quality.png)
{{% /tab %}}
{{< /tabs >}}

### 4. *Per Read Quality*

---
This plot shows distribution of Phred quality score in your data as heatmap.  
The x-axis shows the Phred quality score.  
The y-axis shows the sample and stage (before or after cleaning).  
The color scale shows the proportion of reads.  
The number of reads with average quality scores. Shows if a subset of reads has poor quality.

For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/3%20Per%20Sequence%20Quality%20Scores.html)
The original data of this plot is stored in `03_qc_extras/per_seq_qual_scores.tsv`
{{< tabs >}}
{{% tab name="Captus" %}}
{{< plotly json="/plotly/cleaning_report_4_per_read_quality.json" height="300px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}

{{% /tab %}}
{{< /tabs >}}

### 5. *Read Length Distribution*

---
This plot shows the read length distribution after cleaning.  
For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/7%20Sequence%20Length%20Distribution.html)  
The original data of this plot is stored in `03_qc_extras/seq_len_dist.tsv`
{{< tabs >}}
{{% tab name="Captus" %}}
{{< plotly json="/plotly/cleaning_report_5_seq_len_dist.json" height="240px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}

{{% /tab %}}
{{< /tabs >}}

### 6. *Per Base Nucleotide Content*

---
The proportion of the four nucleotides (A, T, G, C) at each base position in the read.
If there is a distortion, that region will turn red (high frequency) or blue (low frequency).

For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html)
The original data of this plot is stored in `03_qc_extras/per_base_seq_content.tsv`
{{< tabs >}}
{{% tab name="Captus" %}}
{{< plotly json="/plotly/cleaning_report_6_per_base_seq_content.json" height="750px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}

{{% /tab %}}
{{< /tabs >}}

### 7. *Per Read GC Content*

---
The mean GC content of reads. fit normal distribution. if multiple peaks are observed, it may be a sign of contamination.
Broader or multiple peaks may represent contamination with different species.
For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html)
The original data of this plot is stored in `03_qc_extras/per_seq_gc_content.tsv`
{{< tabs >}}
{{% tab name="Captus" %}}
{{< plotly json="/plotly/cleaning_report_7_per_seq_gc_content.json" height="300px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}

{{% /tab %}}
{{< /tabs >}}

### 8. *Sequence Duplication Level*

---
For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html)
The original data of this plot is stored in `03_qc_extras/seq_dup_levels.tsv`
{{< tabs >}}
{{% tab name="Captus" %}}
{{< plotly json="/plotly/cleaning_report_8_seq_dup_level.json" height="300px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}

{{% /tab %}}
{{< /tabs >}}
{{% notice note %}}
in non-random library (RNA-Seq, Target sequencing, etc.).
{{% /notice %}}

### 9. *Adapter Content*

---
For more details, read [<i class="fab fa-readme"></i> FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html)  
The original data of this plot is stored in `03_qc_extras/adapter_content.tsv`
{{< tabs >}}
{{% tab name="Captus" %}}
{{< plotly json="/plotly/cleaning_report_9_adapter_content.json" height="300px" >}}
{{% /tab %}}
{{% tab name="FastQC" %}}
![fastqc_adapt_content](/images/adapter_content.png?height=200px)
{{% /tab %}}
{{< /tabs >}}

---
Created by [Gentaro Shigita]({{< ref "../../credits/#gentaro-shigita">}}) (11.08.2021)  
Last modified by [Gentaro Shigita]({{< ref "../../credits/#gentaro-shigita">}}) (11.08.2021)
