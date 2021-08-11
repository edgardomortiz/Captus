---
title: "Input Preparation"
weight: 12
pre: '<i class="fas fa-check-square"></i> '
---

At this step `Captus` is also flexible allowing you to provide your own assemblies in FASTA format (e.g. NCBI genomes or transcriptomes). This can be done in addition to other samples assembled with `Captus` so you can complement your dataset with published assemblies. If you want to provide sequences assembled elsewhere put them together in a directory and please read the following note about renaming them for `Captus`:

{{% expand "How to rename your FASTA assemblies for Captus" %}}
{{% notice note %}}
In order to add your own FASTA assembly files, a <i class="fas fa-exclamation-triangle"></i> **VERY IMPORTANT** step is to rename them so they are clearly identified throughout the rest of the analysis: 
- Valid extension for assemblies are: `.fa`, `.fna`, `.fasta`, `.fa.gz`, `.fna.gz`, `.fasta.gz`.
- The only special characters that are safe to use in the sample name are `-`, and `_` (`_` is commonly used to replace spaces in many phylogenetic programs). Otherwise, do not use spaces other special characters (``! " # $ % & ( ) * + , . / : ; < = > ? @ [ \ ] ^ ` { | } ~``) or accented letters (like `á`, `è`, `ü`, `ç`, `ñ`), they are just guaranteed to give you headaches at some point.
- In general, a good tip for renaming your samples is to think on how you want the names in your final phylogenetic tree.

These are examples of **valid** FASTA filenames for marker extraction with `Captus`:

- `Arabidopsis_arenosa.fna.gz`
- `Mus_cervicolor_GY747683.fasta`
- `ERI_Vaccinium_macrocarpon.fa.gz`
- `Z.fa`

And here, some examples or **invalid** FASTA filenames:

- `ERR246535.fast` will be ignored because of the invalid extension `.fast`
- `Octomeles_sp.8.fasta`, it is better to replace the `.` in the sample name by a `-` to get `Octomeles_sp-8.fasta`
- `Malus_spontánea.fna.gz`, the sample name contains and accent `á`, it is better to change it to `Malus_spontanea.fna.gz`
{{% /notice %}}
{{% /expand %}}
