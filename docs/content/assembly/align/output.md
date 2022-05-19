---
title: "Output Files"
weight: 14
pre: '<i class="fas fa-dna"></i> '
---

For this example we will use the directory `03_extractions` previously created with the [`extract` module]({{< ref "assembly/extract/output">}}). We run the following `Captus` command to collect markers across samples and align them:

```console
captus_assembly align align -e 03_extractions_CAP/ -o 04_alignments_CAP -k ALL -f ALL
```

After the run is finished we should see a new directory called `04_alignments` with the following structure and files:

![Alignments](/images/alignments.png?width=640&classes=shadow)

### 1. **`01_unaligned`**
This directory contains the unaligned FASTA files corresponding to each marker that were gathered from the extractions directory. The files are organized in subdirectories, first by [marker type]({{< relref "assembly/align/options#-m---markers" >}}) and then by [format]({{< relref "assembly/align/options#-f---formats" >}}).
___
### 2. **`02_untrimmed`**
This directory contains the aligned FASTA files corresponding to each file in the `01_unaligned` directory. The files are organized in subdirectories, first by filtering strategy, then by [marker type]({{< relref "assembly/align/options#-m---markers" >}}), and finally by [format]({{< relref "assembly/align/options#-f---formats" >}}). The subdirectory structure is identical to the one inside the `03_trimmed` directory (see **4** to **15** below).
![Untrimmed alignments](/images/alignment_untrimmed_stages.png?width=1200&classes=shadow)
___
### 3. **`03_trimmed`**
All the files present in the `02_untrimmed` directory are trimmed using `ClipKIT` which removes columns that are mostly empty (see options [**--clipkit_algorithm**]({{< relref "assembly/align/options#--clipkit_algorithm" >}}) and [**--clipkit_gaps**]({{< relref "assembly/align/options#--clipkit_gaps" >}})). The files are organized in subdirectories, first by filtering strategy, then by [marker type]({{< relref "assembly/align/options#-m---markers" >}}), and finally by [format]({{< relref "assembly/align/options#-f---formats" >}}). The subdirectory structure is identical to the one inside the `02_untrimmed` directory (see **4** to **15** below).
![Trimmed alignments](/images/alignment_trimmed_stages.png?width=1200&classes=shadow)
___
### 4. **`01_unfiltered`**
This directory contains the alignments before performing any filtering. All the reference sequences selected by at least a sample will be present as well as all the paralogs per sample. The files are organized in subdirectories, first by [marker type]({{< relref "assembly/align/options#-m---markers" >}}) and then by [format]({{< relref "assembly/align/options#-f---formats" >}}).
___
### 5. **`02_fast`**
This directory contains the alignments where paralogs have been filtered by the `fast` method, which consists in simply keeping the best hit per sample (hit ranked as `00`). All the reference sequences selected by at least a sample will still be present. The files are organized in subdirectories, first by [marker type]({{< relref "assembly/align/options#-m---markers" >}}) and then by [format]({{< relref "assembly/align/options#-f---formats" >}}).
___
### 6. **`03_careful`**
This directory contains the alignments where paralogs have been filtered by the `careful` method. All the reference sequences selected by at least a sample will still be present. The files are organized in subdirectories, first by [marker type]({{< relref "assembly/align/options#-m---markers" >}}) and then by [format]({{< relref "assembly/align/options#-f---formats" >}}).
___
### 7. **`04_unfiltered_no_refs`**, **`05_fast_no_refs`**, **`06_careful_no_refs`**
A subdirectory.
___
### 8. **`01_coding_NUC`**, **`02_coding_PTD`**, **`03_coding_MIT`**
A subdirectory.
___
### 9. **`01_AA`**
A subdirectory.
___
### 10. **`02_NT`**
A subdirectory.
___
### 11. **`03_genes`**
A subdirectory.
___
### 12. **`04_genes_flanked`**
A subdirectory.
___
### 13. **`04_misc_DNA`**, **`05_clusters`**
A subdirectory.
___
### 14. **`01_matches`**
A subdirectory.
___
### 15. **`02_matches_flanked`**
A subdirectory.
___
### 16. **`captus-assembly_align.astral-pro.tsv`**
A subdirectory.
___
### 17. **`captus-assembly_align.paralog_stats.tsv`**
A subdirectory.
___
### 18. **`captus-assembly_align.stats.tsv`**
A subdirectory.
___
### 19. **`captus-assembly_align.report.html`**
A subdirectory.
___
### 20. **`captus-assembly_align.log`**
A subdirectory.

___
Created by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (06.08.2021)  
Last modified by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (24.04.2022)