+++
title = "Align"
weight = 25
pre = "<b>4. </b>"
+++

The last step in the `Captus` workflow is to align the extracted markers so you can estimate phylogenetic trees with your favorite tool (e.g., `IQ-TREE`, `RAxML`, `MrBayes`, `SVDQuartets`, etc.).
`Captus` starts this step by gathering all the markers across your extracted samples and building a FASTA file per marker. Then, it will add the references used during extraction, these are useful to improve alignment since they serve as guides. Then `Captus` aligns the files using `MAFFT` or `MUSCLE`, however, if you choose to align coding sequence in aminoacid (`AA`) and nucleotide (`NT`) in the same run, `Captus` will first align the `AA` version with `MAFFT` or `MUSCLE` and then use that alignment as template to align the `NT` version, thus producing an codon-aware alignment of the coding sequence in nucleotides. Once alignment is completed, paralogs are filtered. Finally, during the trimming stage, Captus removes erroneous stretches of aligned sequences using `TAPER` and then the gappy columns are removed with `ClipKIT` (but you can any other `ClipKIT`'s filtering strategy) as well as exceptionally short sequences which tend to decrease phylogenetic accuracy.
At the end, `Captus` provides several alignment formats from which you can choose the most appropriate for your needs as well as a comprehensive HTML report summarizing alignment statistics along the multiple stages of the `align` step.

- [<i class="fas fa-cog"></i> Options]({{< relref "assembly/align/options" >}})  
- [<i class="fas fa-dna"></i> Output Files]({{< relref "assembly/align/output" >}})  
- [<i class="fas fa-chart-bar"></i> HTML Report]({{< relref "assembly/align/report" >}})

___
Created by [Edgardo M. Ortiz]({{< ref "../../more/credits/#edgardo-m-ortiz">}}) (06.08.2021)  
Last modified by [Gentaro Shigita]({{< ref "../../more/credits/#gentaro-shigita">}}) (12.05.2025)
