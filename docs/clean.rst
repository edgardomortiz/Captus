---
sort: 1
---
# Clean

[BBTools](https://jgi.doe.gov/data-and-tools/bbtools/)

## Usage

```shell
captus_assembly clean -r READS [options]
```

## Options

### Input

`-r` , `--reads`
: Input FASTQ files.
  Valid file name extensions are: .fq, .fastq, .fq.gz, and .fastq.gz.
  The file names must include the string '_R1' (and '_R2' if your data is paired-end).
  Everything before the string '_R1' (or '_R2') will used as sample name.
  
  There is three ways to provide the FASTQ files:
  - path to directory containing FASTQ files (e.g. `-r ./raw_reads`)
  - list
  - pattern


### Output

`-o` , `--out` *`path`*
: Output directory name (dafault: ./01_clean_reads)

`--overwrite`
: Enable to overwrite previous results (default: False)

### Adaptor trimming

`--adaptor_set <Illumina|BGI>`
: aaaa  

`--rna`
: Enable if cleaning RNA-seq reads to additionally trim poly-A tails (default: False)

### Quality trimming and filtering

`--trimq`
: Leading

`--maq`
: After quality trimming, reads with average PHRED quality score 

### Other

`--bbduk_path`
: 

`--fastqc_path`
: 

`--skip_fastqc`
: 

`--ram`
: 

`--threads`
: 

`--show_less`
: 

### Help

`-h` , `--help`
: Show help message and exit

`--version`
: 
