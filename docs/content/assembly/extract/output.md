---
title: "Output Files"
weight: 14
pre: '<i class="fas fa-dna"></i> '
---

For this example we will use the directory `02_assemblies` previously created with the [`assemble` module]({{< ref "assembly/assemble/output">}}). We run the following `Captus` command to search and extract our reference markers from the assemblies:

```console
captus_assembly extract -a 02_assemblies \
-n Angiosperms353 -p SeedPlantsPTD -m SeedPlantsMIT \
-d noncoding_DNA.fasta --cluster_leftovers
```

After the run is finished we should see a new directory called `03_extractions` with the following structure and files:

![Extractions](/images/extractions.png?width=604)

### 1. **`[SAMPLE_NAME]__captus-ext`**
A subdirectory ending in `__captus-ext` is created to contain the extracted markers of each sample separately (S1, S2, S3, and S4 in the image).
___
### 2. **`01_coding_NUC`**, **`02_coding_PTD`**, **`03_coding_MIT`**
These directories contain the extracted **coding** markers from the **NUC**lear, **P**las**T**i**D**ial, and **MIT**ochondrial genomes respectively.
___
### 3. **`[MARKER_TYPE]_coding_AA.faa`**, **`01_AA`**
Coding sequence in **aminoacids**. Prefixes can be `NUC`, `PTD`, or `MIT`.
___
### 4. **`[MARKER_TYPE]_coding_NT.fna`**, **`02_NT`**
Coding sequence in **nucleotides**. Prefixes can be `NUC`, `PTD`, or `MIT`.
___
### 5. **`[MARKER_TYPE]_genes.fna`**, **`03_genes`**
Gene sequence (exons in capital letters + introns in lowercase letters) in **nucleotides**. Prefixes can be `NUC`, `PTD`, or `MIT`.
___
### 6. **`[MARKER_TYPE]_genes_flanked.fna`**, **`04_genes_flanked`**
Gene sequence (exons in capital letters + introns in lowercase letters) plus additional flanking sequence in lowercase **nucleotides**. Prefixes can be `NUC`, `PTD`, or `MIT`.
___
### 7. **`[MARKER_TYPE]_contigs_list.txt`**
List of contig names that had protein hits. Prefixes can be `NUC`, `PTD`, or `MIT`.
___
### 8. **`[MARKER_TYPE]_contigs.gff`**
Annotation track in GFF format for protein hits to contigs in assembly. Prefixes can be `NUC`, `PTD`, or `MIT`.
___
### 9. **`[MARKER_TYPE]_recovery_stats.tsv`**
Tab-separated-values table with marker recovery statistics, these are summarized in the final [Marker Recovery report]({{< ref "assembly/extract/report">}}). Prefixes can be `NUC`, `PTD`, or `MIT`.
___
### 10. **`[MARKER_TYPE]_scipio_final.log`**
Log of the second Scipio's run, where best references have already been selected (when using multi-sequence per locus references) and only the contigs that had hits durin Scipio's initial run are used. Prefixes can be `NUC`, `PTD`, or `MIT`.
___
### 11. **`00_initial_scipio_[MARKER_TYPE]`**
Directory for Scipio's initial run results. The directory contains the set of filtered protein references `[MARKER_TYPE]_best_proteins.faa` (when using multi-sequence per locus references) and the log of Scipio's initial run `[MARKER_TYPE]_scipio_initial.log`. Suffixes can be `NUC`, `PTD`, or `MIT`.
___
### 12. **`04_misc_DNA`**, **`05_clusters`**
These directories contain the extracted **miscellaneous DNA** markers, either from a **DNA** custom set of references or from the **CL**uste**R**ing resulting from using the option `--cluster_leftovers`.
___
### 13. **`[MARKER_TYPE]_matches.fna`**, **`01_matches`**
Matches per miscellaneous DNA marker in **nucleotides**. Prefixes can be `DNA` or `CLR`.
___
### 14. **`[MARKER_TYPE]_matches_flanked.fna`**, **`02_matches_flanked`**
Matches plus additional flanking sequence per miscellaneous DNA marker in **nucleotides**. Prefixes can be `DNA` or `CLR`.
___
### 15. **`[MARKER_TYPE]_contigs_list.txt`**
List of contig names that had miscellaneous DNA marker hits. Prefixes can be `DNA` or `CLR`.
___
### 16. **`[MARKER_TYPE]_contigs.gff`**
Fill
___
### 17. **`[MARKER_TYPE]_recovery_stats.tsv`**
Fill
___
### 18. **`[MARKER_TYPE]_blat_search.log`**
Fill
___
### 19. **`06_assembly_annotated`**
Fill
___
### 20. **`[SAMPLE_NAME]_hit_contigs.fasta`**
Fill
___
### 21. **`[SAMPLE_NAME]_hit_contigs.gff`**
Fill
___
### 22. **`[SAMPLE_NAME]_recovery_stats.tsv`**
Fill
___
### 23. **`leftover_contigs.fasta.gz`**
Fill
___
### 24. **`leftover_contigs_after_custering.fasta.gz`**
Fill
___
### 25. **`captus-assembly_extract.refs.json`**
Fill
___
### 26. **`captus-assembly_extract.stats.tsv`**
Fill
___
### 27. **`captus-assembly_extract.report.html`**
Fill
___
### 28. **`captus-assembly_extract.log`**
Fill
___
### 29. **`clust_id##.##_cov##.##_captus_clusters_refs.fasta`**
Fill

___
Created by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (2021-08-06)  
Last modified by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (2021-09-13)