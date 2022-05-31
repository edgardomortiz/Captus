+++
title = "Assemble"
weight = 15
pre = "<b>2. </b>"
+++

The next step in the workflow is to perform *de novo* assembly with [`MEGAHIT`](https://github.com/voutcn/megahit) on the cleaned reads. `Captus` makes it easy to process many samples in a consistent manner, automatically, and providing a comprehensive [Assembly Control HTML report]({{< ref "assembly/assemble/report" >}}).

`MEGAHIT` is an extremely fast and robust *de novo* assembler that is very memory-efficient thanks to the use of a new data structure called succint De Bruijn Graph (sDBG). Most importantly, accuracy is not sacrificed because it was originally designed for assembling metagenomes. Since then, it has been shown to be an excellent generalist assembler as well.

`Captus` allows you the flexibility to also provide pre-assembled samples. However, we recommend that, whenever you have read data, to assemble it using `Captus`. For transcriptome assemblies for example, other assemblers will produce lots of redundant contigs (due to isoforms, alleles, etc.) which `MEGAHIT` tends to collapse into a single contig. This is ideal for phylogenomics (and perhaps also to build non-redundant reference transcriptomes).

If you are analyzing data with extremely high sequencing depth, or you are trying out parameters and want a quick result, `Captus` can subsample a fixed number of reads prior to assembly.

{{% notice note %}}
In case you assembled your reads elsewhere or you want to use only pre-assembled genomes (e.g., downloaded from GenBank), you can jump ahead to the [`extract` command]({{< ref "extract">}}) page.
{{% /notice %}}

___
Created by [Edgardo M. Ortiz]({{< ref "../../more/credits/#edgardo-m-ortiz">}}) (06.08.2021)  
Last modified by [Edgardo M. Ortiz]({{< ref "../../more/credits/#edgardo-m-ortiz">}}) (30.05.2022)