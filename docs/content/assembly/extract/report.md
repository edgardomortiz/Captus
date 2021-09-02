---
title: "HTML Report"
weight: 15
pre: '<i class="fas fa-chart-bar"></i> '
plotly: true
---
`captus-assembly_extract.report.html`
interactive heatmap
you need internet connection
All the information are also available in a table format file `captus-assembly_extract.stats.tsv`.

you can get a quick overview across all samples and all markers.
of best hit.
If you extracted multiple marker types, you will get one global heatmap with all marker types data in addition to separeted heatmaps for each marker type.
locus recovered on at least one sample will be shown.
loci no single sample recovered will not be shown.
here is a small example with 4 samples and 353 [`Angiosperms353`](https://github.com/mossmatters/Angiosperms353) loci.
note that This heatmap contains information on loci recovered in at least one sample
This heatmap does not include information on loci that were not recovered in all samples
{{< plotly json="/plotly/extraction_report.json" height="300px" >}}

#### Hover information

By hovering your mouse cursor on the heatmap, you can check detailed information listed below at each single data point.  

**Sample**
: Sample name

**Marker type**
: marker type read [`-k, --markers`](http://localhost:1313/assembly/align/options/#-k---markers)  
    - NUC = Nucler proteins
    - PTD = Plastidial proteins
    - MIT = Mitochondrial proteins
    - DNA = Miscellaneous DNA markers
    - CLR = Cluster-derived DNA markers

**Locus**
: Locus name

**Ref name**
: Reference sequence name.
If your reference dataset only contains single sequence per locus, this field will be identical with `Locus` field.

**Ref coords**
: Coordinates of matched region

**Ref type**
: Format of the reference (nucl = nucleotide | prot = amino acid).

**Ref len matched**
: Total matched length.

**Hit count**
: Number of Scipio hit. This value can be used as a paralogy

**Recovered length**
: Percentage of recovered contig length to reference length

**Identity**
: Sequence identity between recovered contig and reference sequence

**Score**
: `Scipio` score. (match - mismatch) / length. For details, read [<i class="fab fa-readme"></i> Scipio's settings](https://www.webscipio.org/help/webscipio#setting)

**Length-weighted score**
: Modified `Scipio` score taking account into 

**Hit length**
: Total 

**CDS length**
: Length of CDS recovered

**Intron length**
: Length of intron(s) recovered

**Flanking length**
: Length of flanking sequence(s) recovered

**Frameshift**
: 

#### Variable dropdown

You can switch the variable to be displayed as a heatmap from following 5 variables to show as heatmap.

- Recovered Length (%)
- Identity (%)
- Hit Number (Paralog): Number of hits
- Score
- Length-weighted score

#### Sort by Value dropdown

Axis sorting buttons
toggle buttons
sort the order of by total value of each sample/locus.

- None: Sort both axis by name
- X only: Sort loci by value
- Y only: Sort samples by value
- Both: Sort loci and samples by value

(called as modebar)
![fastqc_per_base_qual](/images/plotly_modebar.png)
You can zoom in/out, download plot as a PNG
for details, please visit following sites
<https://plotly.com/chart-studio-help/getting-to-know-the-plotly-modebar/>
<https://plotly.com/chart-studio-help/zoom-pan-hover-controls/>

---
Created by [Gentaro Shigita]({{< ref "../../credits/#gentaro-shigita">}}) (11.08.2021)  
Last modified by [Gentaro Shigita]({{< ref "../../credits/#gentaro-shigita">}}) (11.08.2021)
