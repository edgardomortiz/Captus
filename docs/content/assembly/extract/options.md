+++
title = "Options"
weight = 13
pre = '<i class="fas fa-cog"></i> '
+++

# extract
___
To show all available options and their default values you can type in your terminal:
```console
captus_assembly extract --help
```

___
## *Input*
___
### **`-a, --captus_assemblies_dir`**
Path to the output directory from the `assemble` command, (e.g. `02_assemblies` if you used the default name). If you didn't assemble any sample in `Captus` this directory will be created in order to import your contigs assembled elsewhere.

This argument is **required** <i class="fas fa-exclamation-triangle"></i>, the default is **./02_assemblies/**
___
### **`-f, --fastas`**
With this option you can provide the location of your own FASTA assembly files, there are several ways to list them:

- _**Directory:**_ the path to the directory containing your FASTA files assembled elsewhere. This will be the typical way to tell `Captus` which assemblies to import prior to marker extraction. All subdirectories will be searched for FASTA (`.fa`, `.fna`, `.fasta`, `.fa.gz`, `.fna.gz`, `.fasta.gz`) files in this case.

- _**List of files:**_ you can also provide the individual path to each of your FASTA file separated by a single space. This is useful if you only want to analyze only a couple of samples within a directory with many other samples. Another use for lists is when your FASTA assembly files are located in different directories.

- _**UNIX pattern:**_ another easy way to provide lists of files is using the wildcards `*` and `?` to match many or just one character respectively.

Remember, all the FASTA files provided with this option will be imported to the path set by `--captus_assemblies_dir`

This argument is optional and has no default.
___
## *Output*
___
### **`-o, --out`**
With this option you can redirect the output directory to a path of your choice, that path will be created if it doesn't already exist.

Inside this directory, the extracted markers for each sample will be stored in a subdirectory with the ending `__captus-ext`.

This argument is optinal, the default is **./03_extractions/**
___
### **`--max_paralogs`**
Maximum number of secondary hits (copies) of any particular reference marker allowed in the output. We recommend disabling the removal of paralogs (secondary hits/copies) during the `extract` step because the `align` step uses a more sophisticated filter for paralogs. This can be useful for exploratory runs, for example: if after an initial run allowing all paralogs we found out that the average number of secondary hits across samples is 5, we could use this number to get rid of outliers.

This argument is optional, the default is **-1** (include all paralogs in the output).
___
#### **`--max_loci_files`**
When the number of loci in the reference exceeds this value, `Captus` will not write a separate FASTA file per sample per marker, otherwise the hard drive fills up with tons of small files. The file that includes all the extracted markers grouped per sample is still written (this is the only file needed by the final step `align` to produce the marker alignments across all samples).

This argument is optional, the default is **0** (do not write separate loci files).
___
### **`--keep_all`**
Many intermediate files are created during the marker extraction, some are large (like `BLAT`'s `.psl` files) while others small temporary logs of intermediate steps, `Captus` deletes all the unnecesary intermediate files unless you enable this flag.
___
### **`--overwrite`**
Use this flag with caution, this will replace any previous result within the output directory (for the sample names that match).
___
## *Proteins extraction global options (Scipio)*
___
### **`--max_loci_scipio_x2`**
When the number of loci in a protein reference file exceeds this number, `Captus` will not run a second, more exhaustive round of `Scipio`. Usually the results from the first round are extremely similar and sufficient, the second round can become extremely slow as the number of reference proteins grows.

This argument is optional, the default is **2000**.
___
### **`--predict`**
Scipio flags introns as dubious when the splice signals are not found at the exon edges, this may indicate that there are additional aminoacids in the recovered protein that are not present in the reference protein. Enable this flag to attempt translation of these dubious introns, if the translation does not introduce premature stop codons they will be added to the recovered protein. Recommended if you are extracting from RNA assemblies (introns would have been removed in this type of data).
___
## *Nuclear proteins (Scipio)*
___
### **`-n, --nuc_refs`**
The reference set of nuclear proteins to search and extract from the assemblies. `Captus` includes two sets:
- `Angiosperms353` = extract the original Angiosperms353.
- `Mega353` = extract the later expanded version of Angiosperms353; because `Mega353` includes many more sequences, processing times become longer. 
- Alternatively, you can provide the path to a FASTA file (e.g. `-n my_refs/my_nuclear_prots.fa`) that includes your own coding references either in nucleotide or aminoacid.

If you provide a nucleotide file, please also specify the translation table to be used, otherwise `Captus` will translate it using `--nuc_transtable 1`, the "Standard code". See the [complete list of translation tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).

This argument is optional and has no default.
___
### **`--nuc_transtable`**
If you provide a nucleotide file for extraction you can specify the genetic code to translate it.

This argument is optional, the default is **1** (= the "Standard code". See the [complete list of translation tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)).
___
### **`--nuc_min_score`**
Keep hits to the reference proteins that have at least this `Scipio` score. The default has been optimized to perform cross-species extraction in fragmented assemblies like the ones obtained from hybridization capture data. Accepted values are decimals from 0 to 1. For more details, read [<i class="fab fa-readme"></i> Scipio's settings](https://www.webscipio.org/help/webscipio#setting).

This argument is optional, the default is **0.13**.
___
### **`--nuc_min_identity`**
Minimum percentage of identity to the reference protein for a hit to be retained. Accepted values are any number between 0 and 100. For more details, read [<i class="fab fa-readme"></i> Scipio's settings](https://www.webscipio.org/help/webscipio#setting).

This argument is optional, the default is **65**.
___
### **`--nuc_min_coverage`**
Minimum percentage of coverage of the reference protein for a hit to be retained. Accepted values are any number between 0 and 100. For more details, read [<i class="fab fa-readme"></i> Scipio's settings](https://www.webscipio.org/help/webscipio#setting).

This argument is optional, the default is **20**.
___
## *Plastidial proteins (Scipio)*
___
### **`-p, --ptd_refs`**
The reference set of plastidial proteins to search and extract from the assemblies. The options available are:
- `SeedPlantsPTD` = Curated set of chloroplast proteins that spans all Seed Plants.
- Alternatively, you can provide the path to a FASTA file (e.g. `-p my_refs/my_plastome_prots.fa`) that includes your own coding references either in nucleotide or aminoacid.

If you provide a nucleotide file, please also specify the translation table to be used, otherwise `Captus` will translate it using `--ptd_transtable 11`, the "Bacterial, Archaeal and Plant Plastid code". See the [complete list of translation tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).

This argument is optional and has no default.
___
### **`--ptd_transtable`**
If you provide a nucleotide file for extraction you can specify the genetic code to translate it.

This argument is optional, the default is **11** (= the "Bacterial, Archaeal and Plant Plastid code". See the [complete list of translation tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)).
___
### **`--ptd_min_score`**
Keep hits to the reference proteins that have at least this `Scipio` score. The default has been optimized to perform cross-species extraction in fragmented assemblies like the ones obtained from hybridization capture data. Accepted values are decimals from 0 to 1. For more details, read [<i class="fab fa-readme"></i> Scipio's settings](https://www.webscipio.org/help/webscipio#setting).

This argument is optional, the default is **0.2**.
___
### **`--ptd_min_identity`**
Minimum percentage of identity to the reference protein for a hit to be retained. Accepted values are any number between 0 and 100. For more details, read [<i class="fab fa-readme"></i> Scipio's settings](https://www.webscipio.org/help/webscipio#setting).

This argument is optional, the default is **65**.
___
### **`--ptd_min_coverage`**
Minimum percentage of coverage of the reference protein for a hit to be retained. Accepted values are any number between 0 and 100. For more details, read [<i class="fab fa-readme"></i> Scipio's settings](https://www.webscipio.org/help/webscipio#setting).

This argument is optional, the default is **20**.
___
## *Mitochondrial proteins (Scipio)*
___
### **`-m, --mit_refs`**
The reference set of mitochondrial proteins to search and extract from the assemblies. The options available are:
- `SeedPlantsMIT` = Curated set of mitochondrial proteins that spans all Seed Plants.
- `InvertebratesMIT` = Curated set of mitochondrial proteins that spans all Invertebrates (Arthropods, Molluscs, and Nematodes), uses translation table 5.
- Alternatively, you can provide the path to a FASTA file (e.g. `-p my_refs/my_mitome_prots.fa`) that includes your own coding references either in nucleotide or aminoacid.

If you provide a nucleotide file, please also specify the translation table to be used, otherwise `Captus` will translate it using `--mit_transtable 1`, the "Standard code" which is the genetic code used by plant mitochondria. See the [complete list of translation tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).

This argument is optional and has no default.
___
### **`--mit_transtable`**
If you provide a nucleotide file for extraction you can specify the genetic code to translate it.

This argument is optional, the default is **1** (= the "Standard code". See the [complete list of translation tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)).
___
### **`--mit_min_score`**
Keep hits to the reference proteins that have at least this `Scipio` score. The default has been optimized to perform cross-species extraction in fragmented assemblies like the ones obtained from hybridization capture data. Accepted values are decimals from 0 to 1. For more details, read [<i class="fab fa-readme"></i> Scipio's settings](https://www.webscipio.org/help/webscipio#setting).

This argument is optional, the default is **0.2**.
___
### **`--mit_min_identity`**
Minimum percentage of identity to the reference protein for a hit to be retained. Accepted values are any number between 0 and 100. For more details, read [<i class="fab fa-readme"></i> Scipio's settings](https://www.webscipio.org/help/webscipio#setting).

This argument is optional, the default is **65**.
___
### **`--mit_min_coverage`**
Minimum percentage of coverage of the reference protein for a hit to be retained. Accepted values are any number between 0 and 100. For more details, read [<i class="fab fa-readme"></i> Scipio's settings](https://www.webscipio.org/help/webscipio#setting).

This argument is optional, the default is **20**.
___
## *Miscellaneous DNA markers (BLAT)*
___
### **`--dna_refs`**
Path to a FASTA nucleotide file of miscellaneous DNA references. You can include non-coding regions, complete genes (exons + introns), complete mRNAS (including UTRs), etc.

This argument is optional and has no default.
___
### **`--dna_min_identity`**
Minimum percentage of identity to the reference sequence for ahit to be retained. Accepted values are any number between 0 and 100.

This argument is optional, the default is **80**.
___
### **`--dna_min_coverage`**
Minimum percetange of coverage of the reference sequence for a hit to be retained. Accepted values are any number between 0 and 100.

This argument is optional, the default is **20**.
___
## *Assemblies clustering (MMSeqs2)*
Most contigs in your assemblies will not contain hits to your references, that means that many potentially useful markers usually remain unused and unexplored. With this feature you can cluster your unused contigs across samples to find new markers that can complement your phylogenomic dataset.
___

### **`-c, --cluster_leftovers`**
Enable [`MMseqs2`](https://github.com/soedinglab/MMseqs2) clustering across samples for the contigs that had no hits to the reference markers. A new miscellaneous DNA reference is built from the best representative of each cluster in order to perform a miscellaneous DNA marker extraction.
___
### **`--mmseqs2_method`**
Select [MMseqs2's clustering algorithm](https://github.com/soedinglab/mmseqs2/wiki#easy-workflows). Valid options are:
- `easy-linclust` = Fast linear time (for huge datasets), less sensitive clustering
- `easy-cluster` = Sensitive homology search (slower)

This argument is optional, the default is **easy-linclust**.
___
### **`--cluster_mode`**
Select [MMseqs2's clustering mode](https://github.com/soedinglab/mmseqs2/wiki#clustering-modes). Valid options are:
- `0` = Greedy set cover
- `1` = Connected component
- `2` = Greedy incremental (analogous to CD-HIT)

This argument is optional, the default is **2**.
___
### **`--cl_min_identity`**
Minimum percentage of similarity between sequences within a cluster. Accepted values are any number between 75 and 100. Since `Captus` will perform a [Miscellaneous DNA  marker extraction]({{< ref "#miscellaneous-dna-markers-blat">}}) using as reference the best representative of each cluster, it is convenient to set `--cl_min_identity` a little lower than `--dna_min_identity`.

This identity threshold only affects the clustering used for the creation of the new reference of DNA markers, the actual marker extraction still depends of `dna_min_identity`.

This argument is optional, the default is **auto** (= 99% of `dna_min_identity`).
___
### **`--cl_seq_id_mode`**
Select MMseqs2's sequence identity mode. Valid options are:
- `0` = Alignment length
- `1` = Shorter sequence
- `2` = Longer sequence

This argument is optional, the default is **1**.
___
### **`--cl_min_coverage`**
For a sequence to be included in a cluster, this percentage of its length has to be matched by the longest sequence in the cluster. Accepted values are number between 0 and 100.

This only affects the clustering used for the creation of the new reference of DNA markers, the actual marker extraction still depends of `dna_min_coverage`.

This argument is optional, the default is **80**. 
___
### **`--cl_cov_mode`**
Select [MMseqs2's sequence coverage mode](https://github.com/soedinglab/mmseqs2/wiki#how-to-set-the-right-alignment-coverage-to-cluster). Valid options are:
- `0` = Bidirectional (query and target)
- `1` = Target
- `2` = Query
- `3` = Target-in-query

This argument is optional, the default is **1**.
___
### **`--cl_max_seq_len`**
Do not cluster sequences longer than this length in bp, the maximum allowed by MMseqs2 is 65535. Use the value `0` to disable this filter.

This argument is optional, the default is **20000**.
___
### **`--cl_rep_min_len`**
After clustering is finished, only accept cluster representatives of at least this length to be part of the new miscellaneous DNA reference. This avoids the creation of very short locus alignments. Use the value `0` to disable this filter.

This argument is optional, the default is **500**.
___
### **`--cl_min_samples`**
Minimum number of samples per cluster.

This argument is optional, the default is **auto** (= 30% of the total number of samples or at least 4).
___
### **`--cl_tmp_dir`**
Path where to create the temporary directory for `MMseqs2`. Clustering can become slow when done on external drives, set this location to an internal drive.

This argument is optional, the default is **$HOME**.
___
## *Other*
___
### **`--scipio_path`**, **`--blat_path`**
If you have installed your own copies of `Scipio`, `BLAT` you can provide the full path to those copies.

These arguments are optional, the default is to use the **bundled** copies of these programs.

### **`--mmseqs2_path`**
If you have installed your own copy of `MMSeqs2` you can provide the full path to that copy.

This argument is optional, the default is **mmseqs**
___
### **`--ram`**, **`--threads`**, **`--concurrent`**, **`--debug`**, **`--show_less`**
See [Parallelization (and other common options)]({{< ref "parallelization">}})

___
Created by [Edgardo M. Ortiz]({{< ref "../../more/credits/#edgardo-m-ortiz">}}) (06.08.2021)  
Last modified by [Edgardo M. Ortiz]({{< ref "../../more/credits/#edgardo-m-ortiz">}}) (23.02.2023)