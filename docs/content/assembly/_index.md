+++
title = "Assembly"
weight = 2
+++

Captus' first module is called `captus_assembly` but can also be executed as simply `captus` and it aims to create marker alignments across samples starting from the raw sequencing reads of each sample (or alternatively from previously cleaned reads or even previously assembled reads).

To accomplish this the module has four commands: `clean`, `assemble`, `extract`, and `align` which are tipically run in that order:
___
## 1. *clean*
This command will perform adaptor trimming (plus poly-A trimming if you are cleaning RNAseq reads) followed by quality trimming using `bbduk.sh` from the [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/) suite. Once the cleaning is completed, [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) or [`Falco`](https://github.com/smithlabcode/falco) is run on the raw and cleaned reads and a HTML report is generated summarizing the results from all the samples.
___
## 2. *assemble*
Using the cleaned reads produced by the previous step, `Captus` will perform *de novo* assembly using [`MEGAHIT`](https://github.com/voutcn/megahit). The default assembly parameters are tuned for hybridization capture or genome skimming data or a combination of both (`CAPSKIM` preset), additionally we provide two more presets for RNA-Seq (`RNA`) or high-coverage Whole Genome Shotgun (`WGS`) data. `MEGAHIT` only outputs a rough estimation of depth of coverage for the assembled contigs, therefore we run `Salmon` right after assembly to rapidly and accurately estimate contig depth of coverage and automatically remove contigs with depth <1.5x (but this default threshold can be changed). `Captus` can also remove contigs exceeding a given percentage of GC content. This is particularly useful if, for example, you are working with Eukaryotes and want to remove bacterial contamination, whose contigs typically have GC contents above 60 %. An HTML report summarizing several assembly statistics is also produced after this step. After examining the report, the filtering by GC content and/or minimum depth of coverage can be repeated multiple times without the need to reassemble the reads. Even though we recommend using the cleaned reads produced by the `clean` command you can also provide your own previously cleaned reads.
___
## 3. *extract*
During this step `Captus` will search the assemblies produced by the previous step for the loci contained in the provided reference target sequence datasetsets (aminoacids or nucleotides) and then extract them. Proteins can be provided in either aminoacid or nucleotide, these are searched and extracted using [`Scipio`](https://www.webscipio.org/). Additionally, you can provide as reference targets any other DNA sequences (e.g., ribosomal genes, individual exons, entire genes with introns, non-coding regions, RAD loci, etc.), in this case `Captus` uses [`BLAT`](http://hgdownload.soe.ucsc.edu/admin/exe/) for searching and our own code for extracting and stitching partial hits if needed. Finally, since most of the contigs in the assembly are tipically not used because the targets will not be found in most contigs, we provide the option of clustering those unused contigs across samples using [`MMseqs2`](https://github.com/soedinglab/MMseqs2) in order to discover new putative homologous regions that can be used for phylogenomics. If you have your own assemblies in FASTA format you can use them instead of the assemblies produced by the `assemble` command (e.g., downloaded genomes from NCBI). Like in the previous steps, `Captus` will produce an HTML report summarizing the marker recovery statistics across all samples and extracted markers.
___
## 4. *align*
In this step `Captus` will process the results from the `extract` command. First, it will collect all the markers across samples and create a separate FASTA file per marker. Then, the reference sequences used for extraction will be added to their corresponding FASTA marker file to aid as an alignment guide. This is followed by alignment using [`MAFFT`](https://mafft.cbrc.jp/alignment/software/) or [`MUSCLE`](https://www.drive5.com/muscle/). If you are aligning coding sequences, `Captus` will codon-align the nucleotide version using as template the aminoacid alignment of the locus. `Captus` extracts all the copies (hits) of a marker that are found in the assembly and ranks them by their similarity to the reference sequence, once the sequences are aligned, the program filters the paralogs using either the `naive` method which retains the best hit as the ortholog or the `informed` method which takes into account the references and the frequency with which they were selected across all samples to decide which of the copies most likely represents the ortholog (which is not necessarily the best hit). After paralogs have been filtered, the references used for guiding the alignment are removed. Finally, the alignments are trimmed first using [`TAPER`](https://github.com/chaoszhang/TAPER) followed by [`ClipKIT`](https://github.com/JLSteenwyk/ClipKIT). As in previous steps, `Captus` will summarize the alignment statistics of all the markers (e.g. length, mean pairwise identity, missingness, number of informative sites, etc.) and produce an HTML report.
___
To show the main help of the `captus_assembly` module just type `captus --help`:
```console
(captus)$ captus --help
usage: captus command [options]

Captus 1.3.3: Assembly of Phylogenomic Datasets from High-Throughput Sequencing data

Captus-assembly commands:
  command     Program commands (in typical order of execution)
                clean = Trim adaptors and quality filter reads with BBTools,
                        run FastQC on the raw and cleaned reads
                assemble = Perform de novo assembly with MEGAHIT and estimate
                           contig depth of coverage with Salmon: Assembling
                           reads that were cleaned with the 'clean' command is
                           recommended, but reads cleaned elsewhere are also
                           allowed
                extract = Recover targeted markers with BLAT and Scipio:
                          Extracting markers from the assembly obtained with
                          the 'assemble' command is recommended, but any other
                          assemblies in FASTA format are also allowed
                align = Align extracted markers across samples with MAFFT or
                        MUSCLE: Marker alignment depends on the directory
                        structure created by the 'extract' command. This step
                        also performs paralog filtering and alignment trimming
                        using TAPER and ClipKIT

Help:
  -h, --help  Show this help message and exit
  --version   Show Captus' version number

For help on a particular command: captus_assembly command -h
```

So, for example, if you want to show the help of the `extract` command you can type:
```console
captus extract --help
```

___
Created by [Edgardo M. Ortiz]({{< ref "../../more/credits/#edgardo-m-ortiz">}}) (06.08.2021)<br>
Last modified by [Edgardo M. Ortiz]({{< ref "../../more/credits/#edgardo-m-ortiz">}}) (09.04.2025)