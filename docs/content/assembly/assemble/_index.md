+++
title = "Assemble"
weight = 15
pre = "<b>2. </b>"
+++

The next step in the workflow is to perform *de novo* assembly with [`MEGAHIT`](https://github.com/voutcn/megahit) on the cleaned reads. Once the reads have been assembled into contigs, depth of coverage is calculated using [`Salmon`](https://github.com/COMBINE-lab/salmon). `Captus` makes it easy to process many samples in a consistent manner, automatically, and providing a comprehensive [Assembly Control HTML report]({{< ref "assembly/assemble/report" >}}).

`MEGAHIT` is an extremely fast and robust *de novo* assembler that is very memory-efficient thanks to the use of a new data structure called succint De Bruijn Graph (sDBG). Most importantly, accuracy is not sacrificed because it was originally designed for assembling metagenomes. Since then, it has been shown to be an excellent generalist assembler as well.

If you are analyzing data with extremely high sequencing depth, or you are trying out parameters and want a quick result, `Captus` can subsample a fixed number of reads prior to assembly.

`MEGAHIT` only provides a quick estimation of contig depth of coverage, therefore we use `Salmon` to rapidly and accurately calculate the depth of coverage of each contig. Afterwards, contigs with little evidence (by default <1.5x) are automatically filtered.

If you assemble your reads in `Captus` you can also filter the contigs by GC content. For example, remove contigs with >60% GC will remove mostly bacterial contigs. The filtering by GC content and/or by depth of coverage can be repeated multiple times without needing to repeat the assembly and depth of coverage estimation steps.

`Captus` allows you the flexibility to also provide pre-assembled samples. However, we recommend that, whenever you have read data, to assemble it using `Captus`. For transcriptome assemblies for example, other assemblers will produce lots of redundant contigs (due to isoforms, alleles, etc.) which `MEGAHIT` tends to collapse into a single contig. This is ideal for phylogenomics (and perhaps also to build non-redundant reference transcriptomes).

- [<i class="fas fa-clipboard-check"></i> Input Preparation]({{< relref "assembly/assemble/preparation" >}})  
- [<i class="fas fa-cog"></i> Options]({{< relref "assembly/assemble/options" >}})  
- [<i class="fas fa-dna"></i> Output Files]({{< relref "assembly/assemble/output" >}})  
- [<i class="fas fa-chart-bar"></i> HTML Report]({{< relref "assembly/assemble/report" >}})

{{% notice note %}}
In case you assembled your reads elsewhere or you want to use only pre-assembled genomes (e.g., downloaded from GenBank), you can jump ahead to the [`extract` command]({{< ref "extract">}}) page.
{{% /notice %}}

___
Created by [Edgardo M. Ortiz]({{< ref "../../more/credits/#edgardo-m-ortiz">}}) (06.08.2021)  
Last modified by [Gentaro Shigita]({{< ref "../../more/credits/#gentaro-shigita">}}) (12.05.2025)
