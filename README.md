# Captus
Tools for hybridization capture-based target enrichment: Probe Design, Data Assembly, Marker Extraction and Alignment

Provisional installation of dependencies, now you need at least Pyton v3.6:

```bash
conda create -n captus -c conda-forge -c bioconda "python>=3.6" tqdm "perl-bioperl-core>=1.007002" mmseqs2 vsearch mafft bbmap fastqc megahit pigz
```