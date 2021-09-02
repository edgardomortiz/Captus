---
title: "HTML Report"
weight: 15
pre: '<i class="fas fa-chart-bar"></i> '
plotly: true
---
`captus-assembly_align.report.html`
You can open with browser (e.g. Microsoft Edge <i class="fab fa-edge"></i>, Google Chrome <i class="fab fa-chrome"></i>, Mozilla Firefox <i class="fab fa-firefox"></i>, Safari <i class="fab fa-safari"></i>).
`captus-assembly_align.stats.tsv`
If your result contains multiple marker types, separated table for each marker type.

{{< plotly json="/plotly/alignment_report.json" height="400px" >}}
![fastqc_per_base_qual](/images/align_format.drawio.svg?width=800)
- Marker type:
- Format:
- Locus:

column names tell you that in which directory the alignment located in
the tables are sorted in descending order of average value of rows.
colors represent normalized values in each row.

read [`-k, --markers`](http://localhost:1313/assembly/align/options/#-k---markers)

Table has two dropdowns
Format
You can filter the table
For description of each format, read [`-f, --formats`](http://localhost:1313/assembly/align/options/#-f---formats) option.

Variable
You can switch the variable to show from the following 9 variables.

- Number of sequences: Number of sequences in the alignment.
- Alignment length: in bp or aa
- Informaive sites: Number of sites that have at least two different nucleotides or amino acids and each character should occur in at least two samples.
- Constant sites: Number of sites that containing a single nucleotide or amino acid over all sequences.
- Singleton sites: 
- Patterns: Number of sites have distinct patterns.
- Mean pairwise identity (%): 
- Missingness (%): number of gaps / (alignment length \* number of sequences) \* 100 (i.e. gappyness)

|Variable|Description|
|:-|:-|
|Number of sequences|Number of sequences in the alignment|
|Alignment length|in bp or aa|
|Informaive sites|Number of sites that have at least two different nucleotides or amino acids and each character should occur in at least two samples|
|Constant sites|Number of sites that containing a single nucleotide or amino acid over all sequences|
|Singleton sites||
|Patterns|Number of sites have distinct patterns|
|Mean pairwise identity (%)||
|Missingness (%)||number of gaps / (alignment length \* number of sequences) \* 100 (i.e. gappyness)|

<!-- Put the link to IQ-TREE documentation -->
For sites, we refers [<i class="fab fa-readme"></i> IQ-TREE documentation]()

---
Created by [Gentaro Shigita]({{< ref "../../credits/#gentaro-shigita">}}) (11.08.2021)  
Last modified by [Gentaro Shigita]({{< ref "../../credits/#gentaro-shigita">}}) (11.08.2021)
