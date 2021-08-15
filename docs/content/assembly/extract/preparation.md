---
title: "Input Preparation"
weight: 12
pre: '<i class="fas fa-clipboard-check"></i> '
---

At this point you should have *de novo* assemblies from your samples ready. However, `Captus` also gives you the flexibility of starting the analysis from this point by providing your own assemblies as FASTA files or complementing you newly aasembled samples with other assemblies (e.g. genomes or transcriptomes from GenBank). If you want to do so, please read the following note:

{{% expand "How to rename your FASTA assemblies for Captus" %}}
{{% notice note %}}
In order to add your own FASTA assembly files, a <i class="fas fa-exclamation-triangle"></i> **VERY IMPORTANT** step is to rename them so they are clearly identified throughout the rest of the analysis.
{{% /notice %}}

The same tips apply for renaming you FASTA assemblies. In general, a good tip for renaming your samples is to think on how you want the names in your final phylogenetic tree.

The only special characters that are safe to use in the sample name are `-`, and `_` (`_` is commonly used to replace spaces in many phylogenetic programs). Otherwise, do not use spaces other special characters (``! " # $ % & ( ) * + , . / : ; < = > ? @ [ \ ] ^ ` { | } ~``) or accented letters (like `á`, `è`, `ü`, `ç`, `ñ`), they are just guaranteed to give you headaches at some point.

Also, please use this naming convention for your FASTA files:

![Naming convention for FASTA files](/images/fasta.png?width=600)

- Any text before the **extension** will become your **sample name**.
- Valid extension for assemblies are: `.fa`, `.fna`, `.fasta`, `.fa.gz`, `.fna.gz`, `.fasta.gz`.

These are examples of **valid** FASTA filenames for marker extraction with `Captus`:

- `Arabidopsis_arenosa.fna.gz`
- `Mus_cervicolor_GY747683.fasta`
- `ERI_Vaccinium_macrocarpon.fa.gz`
- `Z.fa`

And here, some examples or **invalid** FASTA filenames:

- `ERR246535.fast` will be ignored because of the invalid extension `.fast`
- `Octomeles_sp.8.fasta`, it is better to replace the `.` in the sample name by a `-` to get `Octomeles_sp-8.fasta`
- `Malus_spontánea.fna.gz`, the sample name contains and accent `á`, it is better to change it to `Malus_spontanea.fna.gz`
{{% /expand %}}
___
### *Reference dataset formatting*

Most importantly, in order to extract markers, the sequences in your reference datasets have to follow some naming conventions if you want to take advantage of using multiple reference sequences per locus and a more careful method for paralog filtering. When multiple reference sequences per locus are found in the reference dataset, `Captus` will decide during the extraction which of those references matches your assembly best based on similarity and total recovered length percentage.

Here is an example of a reference protein dataset that has **two loci** (called *accD* and *cemA*) with **five** reference sequences each (probably coming from different taxa to expand phylogenetic coverage). Coding sequences can be provided in either aminoacid or nucleotide. [Miscellaneous DNA markers]({{< ref "/assembly/extract/options#--dna_refs">}}) can only contain nucleotide sequences.

```text
>AA-S46062.1-accD [cluster_size=80]
MALQSLRGSMRSVVGKRICPLIEYAIFPPLPRIIVYASRRARMQRGNYSLIKKPKKVSTLRQYQSTKSPMYQSLQRICGVREWLNKYCMWKEVDEKDFG*
>AAZ94660.1-accD [cluster_size=17]
MEKRWLNSMLSKGELEYRCRLSKSINSLGPIESEGSIINNMNKNIPSHSDSYNSSYSTVDDLVGIRNFVSYDTFLVRDSNSSSYSIYLDIENQIFEIDN*
>ABH88096.1-accD [cluster_size=3]
MQKWRFNSMLLNRELEYGCEFKESLGPIENTSLNEEPKILSDIHKKINRWDDSDNSSYNSLDYLVGADNIQDFLSDKTFLVRDNKRNSYSIYLDIEKKT*
>ABW20568.1-accD [cluster_size=7]
MQNWINNSFQAEFEQESYFGSLGENSMNPRSGGDRYPEALIIRDITGETSAIYFDITDQILENDTHQTILASPIENDLWAEKDISIDIYRYINELIFYD*
>ACU46588.1-accD [cluster_size=1]
MAKYWFNLMLSYKMLSYNKLEHRCGLSKSMDNLNDLGHIGGNEELILNENDAKKNILGLENYNTHSINYLFDSRNIYNLIYNETFLVRNSNGYHYFVYF*
>QNK04966.1-cemA [cluster_size=1]
MKNKKAFIPLLHITFIVFLPWWIAFLFNKGLESWVINWWNTSKSEIFLNDIQEKNILEKFIELEELLLLDEMIKEYPET*
>QNP0849-5.1-cemA [cluster_size=4]
MTKKKAFTPIFYLSFLLFLPWWIDLLFNKCLRSWPTHWWNTRQSEMFLTTLQEKSLLEKFLELEEFLFLDKIIKKEFET*
>QNP08626.1-cemA [cluster_size=2]
MIKNKVFTPLFYLAFIVFLPWGIYFLLNKCMGSWTTNWWTTRESEILSTNINENSLLEKFIQFEEFLLLDEIIKKDTET*
>QNQ64689.1-cemA [cluster_size=9]
MAKNKICIPFISIVFLPWWISFLFKKDFESWVTNWWNTSKSEILLNDIQEKSILKTFIELEELFLLDEMLKEYPETRLQ*
>QPZ48083.1-cemA [cluster_size=7]
MAKKKAFISLIYLASIVFLPWWLSFTFNKSMESWVKNCWNTGPSENFLNDIEEKIIIKKFIELEELSLFDEILKDYTQD*
```

So, if you want to format your reference dataset to include multiple sequences per locus you have to use this naming convention:

![Naming convention for reference sequences](/images/multi_seq_per_locus.png?width=600)

- The **sequence name** (any text found before the first space) can contain multiple **`-`** characters, but only the *last one* will become the **separator**.
- Any text found before the **separator** will be considered as **sequence ID**.
- Any text found after the the **separator** will become the **locus name**.
- Any text found after the first space is considered the **description** and this text is optional.

[`Angiosperms353`](https://github.com/mossmatters/Angiosperms353) and [`HybPiper`](https://github.com/mossmatters/HybPiper) use this format, therefore, in order to mantain compatibility, we also used it to build our reference datasets for plastome proteins `SeedPlantsPTD` and mitochondrial proteins `SeedPlantsMIT`. At the same time it is also sompatible with the references needed by other pipelines (e.g. [`SeCaPr`](content.com/AntonelliLab/seqcap_processor/master/docs/documentation/subdocs/extract_contigs.html#Extracting-target-contigs)) which only consider a single sequence per locus.

{{% notice warning %}}
When the **separator** is not present (even in a single sequence in the whole reference dataset), the entire **sequence name** will be used as the **locus name**. This is fine if your reference dataset only contains a single sequence per locus for example, but if you intend to use the "muti-sequence per locus" format, please ensure that every single sequence contains a **`-`** as **separator**.
{{% /notice %}}