+++
title = "Overview"
weight = 10
+++

Welcome to `Captus`, a toolkit for the assembly of phylogenomic datasets (many samples and loci) from High-Throughput Sequencing (a.k.a. Next Generation Sequencing) data. `Captus` was initially designed specifically for the assembly of target-enriched (e.g. via hybridization of RNA or DNA probes) sequencing data, but has since been expanded to accomodate other common types of HTS data such as Genome Skimming, Hyb-Seq (Target Enrichment + Genome Skimming), RNA-seq, and Whole Genome Sequencing. The toolkit will also include a module for the design of probes for target enrichment experiments.

We wanted to provide as much as possible the same advantages that [`ipyrad`](https://ipyrad.readthedocs.io/en/master/index.html) provides for RAD-seq data, following the same ["Ethos"](https://ipyrad.readthedocs.io/en/master/1-ethos.html) we wanted `Captus` to be:
- _**Simple:**_ Easy to install and easy to use
- _**Fast:**_ We optimized the speed of every step without sacrificing sensitivity
- _**Reproducible:**_ Commands are clear, output directories are well organized, and extensive logs are always kept.
- _**Flexible:**_ You can start the analysis with raw data or with your own cleaned or even assembled reads. You can add samples to your existing datasets without having to reanalyze everything or you can redo analysis on only a subset of samples.
- _**Transparent:**_ We provide informative, well organized, and easy to read output as well as neat HTML reports so you can quickly assess the status of any sample at any stage of the workflow.
___