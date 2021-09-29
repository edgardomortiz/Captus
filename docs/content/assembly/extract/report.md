---
title: "HTML Report"
weight: 15
pre: '<i class="fas fa-chart-bar"></i> '
plotly: true
---
### Concept

---
In the `Captus` workflow, the results of this step, such as **how many loci/markers were recovered, in how many samples, and to what extent, would be the most direct indication of whether the analysis was successful or not**, and thus would be the most interesting for many users.
However, it is not easy to collect, organize, and visualize such important information, especially in a phylo"genomic" project which employs hundreds or thousands of samples and loci/markers.
Don't worry, `Captus` will automatically generate an informative report, `captus-assembly_extract.report.html`.
Please take a look just by opening it in your browser (Microsoft Edge, Google Chrome, Mozilla Firefox, Safari, etc., internet connection required).
Interactive heatmaps allow you to explore results at various scales, from the comprehensive level to the single sample or single locus/marker level!
{{% notice info %}}

- All original data for the report is stored in `captus-assembly_extract.stats.tsv`.
- Since all plots in the report are created using [`Plotly`](https://plotly.com/python), you can use some interactive functions such as zoom in/out, pan, hover, and download plot as a PNG.
For more information, please visit the following sites:

  - <https://plotly.com/chart-studio-help/zoom-pan-hover-controls>
  - <https://plotly.com/chart-studio-help/getting-to-know-the-plotly-modebar>

{{% /notice %}}

### Example

---
Here is a small example of the report you can play with.  
This heatmap shows a extraction result of the [`Angiosperms353`](https://github.com/mossmatters/Angiosperms353) loci from four hybridization-capture samples.
{{< plotly json="/plotly/extraction_report.json" height="400px" >}}
{{% notice note %}}

- If your result contains more than one marker type, the report will include a separete heatmap for each marker type as well as a global heatmap for all marker types.
- If there are multiple overlapping hits, only the information about the best hit (the hit with the highest `Length-weighted score`) will be shown, except for the `Hit count (paralogs)`.
- Information about loci/markers that were not recovered in all samples, or samples where not all loci/markers were recovered, will not appear in the report.
{{% /notice %}}

### Features

---

#### Hover information

By hovering the mouse cursor over the heat map, you can see detailed information of each single data point.  
{{% expand "List of the information to be shown" %}}
|Field|Description|
|-|-|
|**Sample**|Sample name|
|**Marker type**|Marker type (`NUC` = Nuclear proteins \| `PTD` = Plastidial proteins \| `MIT` = Mitochondrial proteins \| `DNA` = Miscellaneous DNA markers \| `CLR` = Cluster-derived DNA markers)|
|**Locus**|Locus/marker name|
|**Ref name**|Reference sequence name (If your reference dataset only contains single sequence per locus, this field will be identical to the `Locus` field)|
|**Ref coords**|Coordinates of the matched parts in the reference sequence (Consecutive coordinates separated by `,` (commas) indicate multiple partial hits on a single contig, whereas coordinates separated by `;` (semicolons) indicate hits on different contigs)|
|**Ref type**|Reference sequence format (`nucl` = nucleotides \| `prot` = amino acids)|
|**Ref len matched**|Total length of the matched parts in the reference sequence|
|**Hit count (paralogs)**|Number of hits found (If the value is greater than 1, it means that there are multiple overlapping hits, which implies the existence of paralogs)|
|**Recovered length**|Percentage of recovered length in the reference sequence length, calcurated as `Ref len matched` / Reference sequence length * 100|
|**Identity**|Sequence identity between the recovered sequence and the reference sequence|
|**Score**|`Scipio` score, calculated as (Number of matched residues - Number of mismatched residues) / Reference sequence length|
|**Length-weighted score**|Modified `Scipio` score to take account into the proportion recovered, calculated as (Number of matched residues - Number of mismatched residues) / Reference sequence length * Proportion recovered|
|**Hit length**|Total length of the contigs where the hit was found (This value should be equal to the sum of subsequent three fields)|
|**CDS length**|Total length of recovered coding sequences (CDS)|
|**Intron length**|Total length of recovered introns|
|**Flanking length**|Total length of recovered flanking sequences|
|**Frameshift**|Positions of additonal or missing bases have been found that would lead to frameshifts during translation (These are most probable due to sequencing/assembly errors, but might also hint to the existance of pseudogenes)|
{{% /expand %}}

---

#### Variable dropdown

By default, the heatmaps show the `Recovered length (%)` of each sample and each locus/marker.  
This dropdown allows you to switch the variable to be shown as a heatmap.  
There are five options:
|Variable|Description|
|-|-|
|**Recovered Length (%)** (default)|Percentage of recovered length in the reference sequence length|
|**Identity (%)**|Sequence identity between the recovered sequence and the reference sequence|
|**Hit Count (Paralog)**|Number of hits found (If the value is greater than 1, it means that there are multiple overlapping hits, which implies the existence of paralogs)|
|**Score**|`Scipio` score, calculated as (Number of matched residues - Number of mismatched residues) / Reference sequence length|
|**Length-weighted Score**|Modified `Scipio` score to take account into the proportion recovered, calculated as (Number of matched residues - Number of mismatched residues) / Reference sequence length * Proportion recovered|

---

#### Sort by Value dropdown

By default, the X and Y axes are sorted alpha-numerically by locus/marker name and sample name, respectively.  
This dropdown allows you to change the sorting manner of each axis as follows:
|Label|X axis|Y axis|
|-|-|-|
|**None** (default)|Sort by name|Sort by name|
|**Mean X**|Sort by **mean** value|Sort by name|
|**Mean Y**|Sort by name|Sort by **mean** value|
|**Mean Both**|Sort by **mean** value|Sort by **mean** value|
|**Total X**|Sort by **total** value|Sort by name|
|**Total Y**|Sort by name|Sort by **total** value|
|**Total Both**|Sort by **total** value|Sort by **total** value|

---
Created by [Gentaro Shigita]({{< ref "../../credits/#gentaro-shigita">}}) (11.08.2021)  
Last modified by [Gentaro Shigita]({{< ref "../../credits/#gentaro-shigita">}}) (28.09.2021)
