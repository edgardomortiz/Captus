---
title: "Options"
weight: 13
pre: '<i class="fas fa-cog"></i> '
---

To show all available options and their default values you can type in your terminal:
```console
captus_assembly align --help
```

___
### *Input*
  -e CAPTUS_EXTRACTIONS_DIR, --captus_extractions_dir CAPTUS_EXTRACTIONS_DIR
                                 Path to the output directory that contains the assemblies and
                                 extractions from previous steps of Captus-assembly. This
                                 directory is called '02_assemblies' if you did not specify a
                                 different name during the 'assemble' or 'extract' steps (default:
                                 ./03_extractions)
  -k MARKERS, --markers MARKERS  Which markers to align, you can provide a comma-separated list,
                                 no spaces (default: all)
                                   NUC = Nuclear proteins inside directories '01_coding_NUC'
                                   PTD = Plastidial proteins inside directories '02_coding_PTD'
                                   MIT = Mitochondrial proteins inside directories '03_coding_MIT'
                                   DNA = Miscellaneous DNA markers inside directories '04_misc_DNA'
                                   CLR = Cluster-derived DNA markers inside directories
                                         '05_clusters'
                                   ALL = Shortcut for NUC,PTD,MIT,DNA,CLR
  -f FORMATS, --formats FORMATS  Which alignment format(s) to prepare for each marker category,
                                 you can provide a comma-separated list, no spaces (default:
                                 AA,NT,GE,MA)
                                   Valid types for NUC, PTD, and MIT markers:
                                   AA = Coding sequences in aminoacids
                                   NT = Coding sequences in nucleotides
                                   GE = Complete gene sequences (exons + introns) without flanking
                                        upstream or downstream basepairs
                                   GF = Complete gene sequences with flanking upstream and
                                        downstream basepairs
                                   Valid types for miscellaneous DNA and CLusteR-derived markers:
                                   MA = Matched sequences without flanking upstream or downstream
                                        basepairs
                                   MF = Matched sequences with flanking upstream and downstream
                                        basepairs
                                   ALL = Shortcut for AA,NT,GE,GF,MA,MF
___
### *Output*
___
#### **`-o, --out`**
With this option you can redirect the output directory to a path of your choice, that path will be created if it doesn't already exist. If you don't provide an output directory name, `Captus` creates a directory called `04_alignments` to contain the output.
___
#### **`--keep_all`**
Many intermediate log files are created by `MAFFT` and `ClipKIT` during assembly, `Captus` deletes all the unnecesary intermediate files unless you enable this flag.
___
#### **`--overwrite`**
Use this flag with caution, this will replace any previous results.
___
### *Extraction references*
  -n NUC_REFS, --nuc_refs NUC_REFS
                                 Set of nuclear protein references. These will be used as guides
                                 for alignment and removed from the final alignment files. The
                                 references are also used when the 'careful' method for paralog
                                 removal is chosen. Options are:
                                   Angiosperms353 = The original set of target proteins from
                                                    Angiosperms353
                                   Mega353 = The improved set of target proteins from Angiosperms353
                                   PathAA,PathNT = Paths to protein reference set both in
                                                   aminoacids and nucleotides separated by a
                                                   comma, no spaces (e.g.: ref.faa,ref.fna)
                                   PathNT,transtable = Path to protein reference set in
                                                       nucleotides, followed by the Genetic Code
                                                       number to translate it, separated by a
                                                       comma, no spaces. If omitted, GC defaults
                                                       to 1 (e.g.: ref.fna,1)
                                   PathAA = Path to protein reference set in aminoacids, in this
                                            case it is not possible to add references to the
                                            nucleotide CDS alignments (e.g.: ref.faa)
  -p PTD_REFS, --ptd_refs PTD_REFS
                                 Set of plastidial protein references. These will be used as
                                 guides for alignment and removed from the final alignment files.
                                 The references are also used when the 'careful' method for
                                 paralog removal is chosen. Options are:
                                   SeedPlantsPTD = A set of plastidial proteins for Seed Plants,
                                                   curated by us
                                   PathAA,PathNT = paths to protein reference set both in
                                                   aminoacids and nucleotides separated by a
                                                   comma, no spaces (e.g.: ref.faa,ref.fna)
                                   PathNT,transtable = Path to protein reference set in
                                                       nucleotides, followed by the Genetic Code
                                                       number to translate it, separated by a
                                                       comma, no spaces. If omitted, GC defaults
                                                       to 11 (e.g.: ref.fna,11)
                                   PathAA = path to protein reference set in aminoacids, in this
                                            case it is not possible to add references to the
                                            nucleotide CDS alignments (e.g.: ref.faa)
  -m MIT_REFS, --mit_refs MIT_REFS
                                 Set of mitochondrial protein references. These will be used as
                                 guides for alignment and removed from the final alignment files.
                                 The references are also used when the 'careful' method for
                                 paralog removal is chosen. Options are:
                                   SeedPlantsMIT = A set of mitochondrial proteins for Seed
                                                   Plants, curated by us
                                   PathAA,PathNT = paths to protein reference set both in
                                                   aminoacids and nucleotides separated by a
                                                   comma, no spaces (e.g.: ref.faa,ref.fna)
                                   PathNT,transtable = Path to protein reference set in
                                                       nucleotides, followed by the Genetic Code
                                                       number to translate it, separated by a
                                                       comma, no spaces. If omitted, GC defaults
                                                       to 1 (e.g.: ref.fna,1)
                                   PathAA = path to protein reference set in aminoacids, in this
                                            case it is not possible to add references to the
                                            nucleotide CDS alignments (e.g.: ref.faa)
  -d DNA_REFS, --dna_refs DNA_REFS
                                 Path to a FASTA nucleotide file of miscellaneous DNA references.
                                 These will be used as guides for alignment and removed from the
                                 final alignment files.
___
### *MAFFT*
___
#### **`--mafft_algorithm`**
Select [MAFFT's alignment algorithm](https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html)
___
#### **`--mafft_timeout`**
Modify the waiting time in seconds for an individual alignment to complete. When using more exhaustive MAFFT algorithm (e.g. `genafpair`), alignment can take very long (up to hours depending on sample number an length of the sequences). `Captus` terminates an alignment if it doesn't finish in 21600 seconds (6 hours).
___
### *Paralog filtering*
___
#### **`--filter_method`**
We provide two filtering methods for paralog removal, you can select either or both:
- _**fast**_ = Only the best hit for each sample (marked as hit=00) is retained, when the reference only contains a single sequence per locus it is equivalent to the `carefu` method.
- _**careful**_ = Only keep the copy (regardless of hit ranking) that is most similar to the reference sequence that was
 chosen most frequently among all other samples in the alignment. This method was designed to take advantage of references that contain several sequences per locus (like `Angiosperms353`), if the reference only contains a single reference per locus the result will be identical to the `fast` method.
- _**both**_ = Two separate folders will be created, each containing the results from each filtering method.
___
### *ClipKIT*
___
#### **`--clipkit_algorithm`**
Select [ClipKIT's trimming mode](https://jlsteenwyk.com/ClipKIT/advanced/index.html#modes)
___
#### **`--clipkit_gaps`**
Gappynes threshold per position. Accepted values between 0 and 1. This argument is ignored when using the `kpi` and `kpic` algorithms or intermediate steps that use `smart-gap`.
___
### *Other*
___
#### **`--redo_from`**
You can repeat the analysis without undoing all the steps. These are the points from which you ca restart the `align` command:
- _**alignment**_ = Delete all subdirectories with alignments and restart.
- _**filtering**_ = Delete all subdirectories with paralog-filtered alignments and restart.
- _**removal**_ = Delete all subdirectories with alignments whose references have been removed and restart.
- _**trimming**_ = Delete all subdirectories with trimmed alignments and restart.
___
#### **`--mafft_path`**, **`--clipkit_path`**
If you have installed your own copies of `MAFFT` or `ClipKIT` you can provide the full path to those copies.
___
#### **`--ram`**, **`--threads`**, **`--concurrent`**, **`--debug`**, **`--show_less`**
See [Parallelization (and other common options)]({{< ref "parallelization">}})
___