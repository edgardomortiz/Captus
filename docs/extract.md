---
sort: 3
---
# Extract

[Scipio](https://www.webscipio.org)

[BLAT](http://www.kentinformatics.com/products.html)

[MMseqs2](https://github.com/soedinglab/MMseqs2)

## Usage

```shell
captus_assembly extract -a CAPTUS_ASSEMBLIES_DIR [options]
```

## Options

### Input

`-a` , `--captus_assemblies_dir`
: Path to the output directory from the `assemble` step of Captus-assembly. 

`-f` , `--fastas`
: 

### Output

`-o` , `--out`
: 

`--max_loci_files`
: 

`--max_loci_scipip2x`
: 

`--keep_all`
: 

`--overwrite`
: 

### Nuclear proteins extraction (Scipio)

`-n` , `--nuc_refs`
:

`--nuc_min_score`
:

`--nuc_min_identity`
:

`--nuc_min_coverage`
:

`--nuc_transtable`
:

### Plastidial proteins extraction (Scipio)

`-p` , `--ptd_refs`
:

`--ptd_min_score`
:

`--ptd_min_identity`
:

`--ptd_min_coverage`
:

`--ptd_transtable`
:

### Mitochondrial proteins extraction (Scipio)

`-m` , `--mit_refs`
:

`--mit_min_score`
:

`--mit_min_identity`
:

`--mit_min_coverage`
:

`--mit_transtable`
:

### Miscellaneous DNA extraction (BLAT)

`-d` , `--dna_refs`
:

`--dna_min_identity`
:

`--dna_min_coverage`
:

### Assemblies clustering (MMseqs2)

`-c` , `--cluster_leftovers`
:

`--cl_min_identity`
:

`--cl_min_coverage`
:

`--cl_gap_open`
:

`--cl_gap_extend`
:

`--cl_max_seq_len`
:

`--cl_temp_dir`
:

`--cl_min_samples`
:

### Other

`--scipio_path`
:

`blat_path`
:

`--mmseqs2_path`
:

`--ram`
:

`--threads`
:

`--concurrent`
:

`--show_less`
:

### Help

`-h` , `--help`
:

`--version`
:
