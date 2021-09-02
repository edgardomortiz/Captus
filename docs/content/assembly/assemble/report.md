---
title: "HTML Report"
weight: 15
pre: '<i class="fas fa-chart-bar"></i> '
plotly: true
---

`captus-assembly_assemble.report.html`
two sections below

### 1. *Summary Table*

---
This table shows general statistics.  
The original data is stored in `captus-assembly_assemble.stats.tsv`
{{< plotly json="/plotly/assemble_report_1_summary_table.json" height="250px" >}}

- **Sample** = Sample name
- **Total Length (bp)** = Total length of all contigs
- **Number of Contigs** = Number of contigs assembled
- **Mean Length (bp)** = Mean contig length
- **Contig N50 (bp)** = 50% of total assembly size is represented by contigs with length equal to or longer than this value
- **Longest Contig (bp)** = Longest contig length
- **Shortest Contig (bp)** = Shortest contig length
- **≥1kbp Contig (%)** = Percentage of contigs with length ≥ 1 kbp
- **GC Content (%)** = overall GC content of contigs
- **Mean Depth (x)** = Mean depth of all contigs

**Sample**
: Sample name

**Total Length (bp)**
: Total length of all contigs

**Number of Contigs**
: Number of contigs assembled

**Mean Length (bp)**
: Mean contig length

**Contig N50 (bp)**
: 50% of total assembly size is represented by contigs with length equal to or longer than this value

**Longest Contig (bp)**
: Longest contig length

**Shortest Contig (bp)**
: Shortest contig length

**≥1kbp Contig (%)**
: Percentage of contigs with length ≥ 1 kbp

**GC Content (%)**
: overall GC content of contigs

**Mean Depth (x)**
: Mean depth of all contigs

### 2. *Detailed Stats*

---
This plot shows with some more detailed as simple or stacked bar plot.
{{< plotly json="/plotly/assemble_report_2_detailed_stats.json" height="400px" >}}

by switching the dropdown menu on the top right
cell colors represent normalized values in each column (Green: high, Red: low).

- **Median Length (bp)** = Median contig length
- **Contig breakdown by length (%)** = Proportion of 1 kbp, 2 kbp, 5 kbp, and 10 kbp
- **Length breakdown by contig length (%)** = 
- **Contig breakdown by depth (%)** = 

---
Created by [Gentaro Shigita]({{< ref "../../credits/#gentaro-shigita">}}) (11.08.2021)  
Last modified by [Gentaro Shigita]({{< ref "../../credits/#gentaro-shigita">}}) (11.08.2021)
