+++
title = "HTML Report"
weight = 15
pre = '<i class="fas fa-chart-bar"></i> '
plotly = true
+++

## Concept

---
The output from this `extract` module, such as **how many loci are recovered, in how many samples, and to what extent, would be the most direct indication of whether your analysis is successful or not**, and thus would be of most interest to many users.
However, collecting, summarizing, and visualizing such important information can be backbreaking, especially in a phylo"genomic" project which typically employs hundreds or even thousands of samples and loci.  

Don't worry, `Captus` automatically generates an informative report!
Open `captus-assembly_extract.report.html` with your browser (internet connection required) to explore your extraction result at various scales, from the global level to the single sample or single locus level.
{{% notice tip %}}

- The entire report is based on data stored in [`captus-assembly_extract.stats.tsv`]({{< relref "assembly/extract/output#26-captus-assembly_extractstatstsv" >}}).
- All tables and plots in the report are interactive powered by [`Plotly`](https://plotly.com/python).  
Visit the following sites once to take full advantage of its interactivity:

  - <https://plotly.com/chart-studio-help/getting-to-know-the-plotly-modebar>
  - <https://plotly.com/chart-studio-help/zoom-pan-hover-controls>
{{% /notice %}}

## Example

---
Here is a small example of the report you can play with!  
The heatmap shows a extraction result of the `Angiosperms353` ([Johnson *et al*., 2019](https://doi.org/10.1093/sysbio/syy086)) loci from targeted-capture data of four plant species.
The blue bars along with *x*- and *y*-axes indicate how many loci are recovered in each sample and how many samples each locus is recovered in, respectively.  
{{< plotly json="/plotly/extraction_report.json" height="500px" >}}
{{% notice note %}}

- When your result contains more than one marker type, the report will include separate plots for each marker type.
- For loci with more than one copy found in a sample, information on best hit (hit with the highest `weighted score`) will be shown.
- Information on loci with no samples recovered and samples with no loci recovered will not be shown.
{{% /notice %}}

## Features

---

### 1. Hover information

Hover mouse cursor over the heatmap to see detailed information about each single data point.  
{{% expand "List of the information to be shown" %}}
|Field|Description|Unit|
|-|-|-|
|**Sample**|Sample name|-|
|**Marker type**|Marker type (`NUC` = Nuclear proteins; `PTD` = Plastidial proteins; `MIT` = Mitochondrial proteins; `DNA` = Miscellaneous DNA markers; `CLR` = Cluster-derived DNA markers)|-|
|**Locus**|Locus name|-|
|**Ref name**|Reference sequence name selected|-|
|**Ref coords**|Matched coordinates with respect to the reference sequence (Consecutive coordinates separated by `,` indicate partial hits on the same contig; coordinates separated by `;` indicate hits on different contigs)|-|
|**Ref type**|Reference sequence format (`nucl` = nucleotides; `prot` = amino acids)|-|
|**Ref len matched**|Total length of `Ref coords`|aa/bp|
|**Total hits (copies)**|Number of hits found (Values greater than 1 imply the presence of paralogs)|-|
|**Recovered length**|Percentage of reference sequence length recovered, calcurated as (`Ref len matched` / Reference sequence length) * 100|%|
|**Identity**|Sequence identity of the recovered sequence to the reference sequence|%|
|**Score**|Score inspired by [`Scipio`](https://www.webscipio.org/help/webscipio#setting), calculated as (matches - mismatches) / reference sequence length|-|
|**Weighted score**|Weighted `score` to address multiple reference sequences per locus<br>(for details, read [<i class="fab fa-readme"></i> Information included in the table]({{< relref "assembly/extract/output#26-captus-assembly_extractstatstsv" >}}))|-|
|**Hit length**|Length of sequence recovered|bp|
|**CDS length**|Total length of coding sequences (CDS) recovered (always `NA` when the `ref_type` is `nucl`)|bp|
|**Intron length**|Total length of introns recovered (always `NA` when the `ref_type` is `nucl`)|bp|
|**Flanking length**|Total length of flanking sequences recovered|bp|
|**Number of frameshifts**|Number of corrected frameshifts in the extracted sequence<br>(always `0` when `ref_type` is `nucl`)|-|
|**Position of frameshifts**|Positions of corrected frameshifts in the extracted sequence<br>(`NA` when no frameshift is detected or `ref_type` is `nucl`)|-|
|**Contigs in best hit**|Number of contigs used to assemble the best hit|-|
|**Best hit L50**|Least number of contigs in best hit that contain 50% of the best hit's recovered length|-|
|**Best hit L90**|Least number of contigs in best hit that contain 90% of the best hit's recovered length|-|
|**Best hit LG50**|Least number of contigs in best hit that contain 50% of the reference locus length|-|
|**Best hit LG90**|Least number of contigs in best hit that contain 90% of the reference locus length|-|

\* When your data is huge (number of samples * number of loci > 500k), only `Sample`, `Locus`, and the variable selected in the `Variable` dropdown will be shown.
{{% /expand %}}

---

### 2. `Variable` dropdown

Switch this dropdown to change the variable to be shown as a heatmap among the following options:  
|Variable|Description|Unit|
|-|-|-|
|**Recovered Length**|Percentage of reference sequence length recovered|%|
|**Identity**|Sequence identity of the recovered sequence to the reference sequence|%|
|**Total Hits (Copies)**|Number of hits found (Values greater than 1 imply the presence of paralogs)|-|
|**Score**|Score inspired by [`Scipio`](https://www.webscipio.org/help/webscipio#setting), calculated as (matches - mismatches) / reference sequence length|-|
|**Weighted Score**|Weighted `score` to address multiple reference sequences per locus<br>(for details, read [<i class="fab fa-readme"></i> Information included in the table]({{< relref "assembly/extract/output#26-captus-assembly_extractstatstsv" >}}))|-|
|**Number of Frameshifts**|Number of corrected frameshifts in the extracted sequence<br>(always `0` if the reference sequence is in nucleotide)|-|
|**Contigs in Best Hit**|Number of contigs used to assemble the best hit|-|
|**Best Hit L50**|Least number of contigs in best hit that contain 50% of the best hit's recovered length|-|
|**Best Hit L90**|Least number of contigs in best hit that contain 90% of the best hit's recovered length|-|
|**Best Hit LG50**|Least number of contigs in best hit that contain 50% of the reference locus length|-|
|**Best Hit LG90**|Least number of contigs in best hit that contain 90% of the reference locus length|-|

---

### 3. `Sort by Value` dropdown

Switch this dropdown to change the sorting manner of each axis as follow:
|Label|Locus (*x*-axis)|Sample (*y*-axis)|
|-|-|-|
|**None**|Sort by name|Sort by name|
|**Mean X**|Sort by **mean** value|Sort by name|
|**Mean Y**|Sort by name|Sort by **mean** value|
|**Mean Both**|Sort by **mean** value|Sort by **mean** value|
|**Total X**|Sort by **total** value|Sort by name|
|**Total Y**|Sort by name|Sort by **total** value|
|**Total Both**|Sort by **total** value|Sort by **total** value|

---
Created by [Gentaro Shigita]({{< ref "../../more/credits/#gentaro-shigita">}}) (11.08.2021)  
Last modified by [Gentaro Shigita]({{< ref "../../more/credits/#gentaro-shigita">}}) (16.09.2022)
