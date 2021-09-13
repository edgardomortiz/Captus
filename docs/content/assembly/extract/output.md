---
title: "Output Files"
weight: 14
pre: '<i class="fas fa-dna"></i> '
---

For this example we will use the directory `02_assemblies` previously created with the [`assemble` module]({{< ref "assembly/assemble/output">}}). We run the following `Captus` command to search and extract our reference markers from the assemblies:

```console
captus_assembly extract \
-a 02_assemblies -o 03_extractions \
-n Angiosperms353 \
-p SeedPlantsPTD \
-m SeedPlantsMIT \
-d noncoding_DNA.fasta \
--cluster_leftovers
```

After the run is finished we should see a new directory called `03_extractions` with the following structure and files:

![Extractions](/images/extractions.png?width=604)

### 1. **`[sample]__captus-ext`**
A subdirectory ending in `__captus-ext` is created to contain the extracted markers of each sample separately (S1, S2, S3, and S4 in the image).
___
### 2. **`01_coding_NUC`, `02_coding_PTD`, `03_coding_MIT`**
These directories contain the extracted **coding** markers from the **NUC**lear, **P**las**T**i**D**ial, and **MIT**ochondrial genomes respectively.
___
### 3. **`_coding_AA.faa`**
Coding sequence in **aminoacids**. Prefixes can be `NUC`, `PTD`, or `MIT`.
___
### 4. **`_coding_NT.fna`**
Coding sequence in **nucleotides**. Prefixes can be `NUC`, `PTD`, or `MIT`.
___
### 5. **`_genes.fna`**
Gene sequence (exons in capital letters + introns in lowercase letters) in **nucleotides**. Prefixes can be `NUC`, `PTD`, or `MIT`.
___
### 6. **`_genes_flanked.fna`**
Gene sequence (exons in capital letters + introns in lowercase letters) plus additional flanking sequence in lowercase **nucleotides**. Prefixes can be `NUC`, `PTD`, or `MIT`.
___
### 7. **`_contigs_list.txt`**
List of contig names that had protein hits. Prefixes can be `NUC`, `PTD`, or `MIT`.
___
### 8. **`_contigs.gff`**
Fill
___
### 9. **`_recovery_stats.tsv`**
Fill
___
### 10. **`_scipio_final.log`**
Fill
