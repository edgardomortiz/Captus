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

___
Created by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (06.08.2021)  
Last modified by [Edgardo M. Ortiz]({{< ref "../../credits/#edgardo-m-ortiz">}}) (27.03.2022)