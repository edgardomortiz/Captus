---
title: "Why BBTools?"
weight: 11
hidden: true
---

We needed something *FAST*, but most importantly *ACCURATE* and *RELIABLE*. Years ago, Brian Bushnell, the author of [`BBTools`](https://jgi.doe.gov/data-and-tools/bbtools), posted a comparison between his own `bbduk.sh`, and the popular [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) and [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). You can read the experiment setup and his results [here](http://seqanswers.com/forums/showpost.php?p=138702&postcount=2).  

Since then, there have been new versions of these programs, and a new one, [`fastp`](https://github.com/OpenGene/fastp), was published. Therefore, we updated the experiment including `fastp`, increasing the minimum length to 21 bp, and processing paired-end reads which is closer to the real use of these tools.  

First, we created a `conda` environment for these three programs:
```console
conda create -n rem_adaptors -c bioconda bbmap cutadapt fastp trimmomatic
conda activate rem_adaptors
```

Then, we grabbed a million reads (PE 2x150 bp) from a plant genome project, `reformat.sh` is also part of `BBTools`:
```console
reformat.sh in=SB12_R#.fq.gz out=raw_R#.fq.gz reads=1000000
```

These are the same "rotated" adaptors that were included in the file `gruseq.fa` that Bushnell made for his test:
```text
>Gruseq_Adapter_Index_1_6
CTGACCTTCTCATATACGAGCTTAGAATCGATATGATACTGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_2
CTGACCTTCTCATATACGAGCTTAGAATCGATAACTGCGTGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_3
CTGACCTTCTCATATACGAGCTTAGAATCGATAGGTCCATGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_4
CTGACCTTCTCATATACGAGCTTAGAATCGATAGCTAATTGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_5
CTGACCTTCTCATATACGAGCTTAGAATCGATATATCGCTGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_6
CTGACCTTCTCATATACGAGCTTAGAATCGATACAATTGTGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_7
CTGACCTTCTCATATACGAGCTTAGAATCGATAATCTGATGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_8
CTGACCTTCTCATATACGAGCTTAGAATCGATATAGGCTTGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_9
CTGACCTTCTCATATACGAGCTTAGAATCGATACTGATCTGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_10
CTGACCTTCTCATATACGAGCTTAGAATCGATAGTCAGGTGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_11
CTGACCTTCTCATATACGAGCTTAGAATCGATACCAGTATGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_12
CTGACCTTCTCATATACGAGCTTAGAATCGATAAGGCGTTGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_13
CTGACCTTCTCATATACGAGCTTAGAATCGATATCGATTATTGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_14
CTGACCTTCTCATATACGAGCTTAGAATCGATATCGGAACGTGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_15
CTGACCTTCTCATATACGAGCTTAGAATCGATATGCGATCTTGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_16
CTGACCTTCTCATATACGAGCTTAGAATCGATAAACGAAACTGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_18_7
CTGACCTTCTCATATACGAGCTTAGAATCGATACGAACATATGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_19
CTGACCTTCTCATATACGAGCTTAGAATCGATACGCTTTACTGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_20
CTGACCTTCTCATATACGAGCTTAGAATCGATACGCCAAGGTGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_21
CTGACCTTCTCATATACGAGCTTAGAATCGATACGGGACCTTGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_22
CTGACCTTCTCATATACGAGCTTAGAATCGATAACGTACGTTGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_23
CTGACCTTCTCATATACGAGCTTAGAATCGATACTCGCCTGTGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_25
CTGACCTTCTCATATACGAGCTTAGAATCGATATAGCTGTGTGAGACGTGCAACGAGGAGCAGGC
>Gruseq_Adapter_Index_27
CTGACCTTCTCATATACGAGCTTAGAATCGATATGGAAGGGTGAGACGTGCAACGAGGAGCAGGC
```

Now, we add those adaptors to the raw reads with another program from `BBTools`:
```console
addadapters.sh in=raw_R#.fq.gz out=dirty_R#.fq.gz qout=33 ref=gruseq.fa right
```

`addadapters.sh` modifies the read headers adding the exact information of where was the adaptor inserted. For example, the adaptors were added from position **124** in this read:
```console
@4_150_124_150_124 /1
ANTAAAAAGAGTAGTGTCAGATAGCTTATATGGAGAAAGCCATAGCAATTTTATCAGTGCTGTAGAGGAATTAAAAATAGAATATGCAGTGGGAATCTGGAGCAATCATGGGGTCTGGCTTCCACTGACCTTCTCATATACGAGCTTAGA
+
F!FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FF:F:FFFFFFFFFFF
```

Now we proceed to remove these synthetic adaptors, these are the commands used for each software:
```console
# Cutadapt 1.18:
time cutadapt -m 21 -j 0 -b CTGACCTTCTCATATACGAGCTTAGAATCGATATGATACTGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATAACTGCGTGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATAGGTCCATGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATAGCTAATTGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATATATCGCTGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATACAATTGTGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATAATCTGATGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATATAGGCTTGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATACTGATCTGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATAGTCAGGTGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATACCAGTATGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATAAGGCGTTGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATATCGATTATTGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATATCGGAACGTGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATATGCGATCTTGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATAAACGAAACTGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATACGAACATATGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATACGCTTTACTGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATACGCCAAGGTGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATACGGGACCTTGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATAACGTACGTTGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATACTCGCCTGTGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATATAGCTGTGTGAGACGTGCAACGAGGAGCAGGC -b CTGACCTTCTCATATACGAGCTTAGAATCGATATGGAAGGGTGAGACGTGCAACGAGGAGCAGGC -o cutadapt_R1.fq.gz -p cutadapt_R2.fq.gz dirty_R1.fq.gz dirty_R2.fq.gz

# Trimmomatic 0.39:
time trimmomatic PE -phred33 dirty_R1.fq.gz dirty_R2.fq.gz trimmomatic_R1.fq.gz trimmomatic_U1.fq.gz trimmomatic_R2.fq.gz trimmomatic_U2.fq.gz ILLUMINACLIP:gruseq.fa:2:28:10:2:keepBothReads MINLEN:21

# fastp 0.20.1:
time fastp -w 8 -Q -l 21 --adapter_fasta gruseq.fa --detect_adapter_for_pe --in1 dirty_R1.fq.gz --in2 dirty_R2.fq.gz --out1 fastp_R1.fq.gz --out2 fastp_R2.fq.gz

# bbduk.sh 38.18:
time bbduk.sh in=dirty_R#.fq.gz out=bbduk_R#.fq.gz ref=gruseq.fa ktrim=r mink=12 hdist=1 minlen=21 tpe tbo

# bbduk.sh 38.18 in Captus 0.0.11:
time bbduk.sh ktrim=r minlength=21 interleaved=f tpe tbo ref=gruseq.fa in=dirty_R#.fq.gz out=stdout.fq k=21 mink=11 hdist=2 | bbduk.sh ktrim=r minlength=21 interleaved=f tpe tbo ref=gruseq.fa in=stdin.fq out=captus_R#.fq.gz k=19 mink=9 hdist=1
```

Finally, we use `addapters.sh` to grade the presence of adaptor of each set of reads:
```console
addadapters.sh in=dirty_R#.fq.gz grade
addadapters.sh in=cutadapt_R#.fq.gz grade
addadapters.sh in=trimmomatic_R#.fq.gz grade
addadapters.sh in=fastp_R#.fq.gz grade
addadapters.sh in=bbduk_R#.fq.gz grade
addadapters.sh in=captus_R#.fq.gz grade
```

The tests were made in a MacbookPro 2015 with 4 CPU cores (8 threads) and 16 GB of RAM, we can see in this table the summary of the times and the results from the grading of `addadapters.sh`:

|Metric                       |dirty  |Cutadapt |Trimmomatic|fastp    |bbduk   |Captus   |
|-----------------------------|-------|---------|-----------|---------|--------|---------|
|*Time to clean*              |*NA*   |2m16.692s|0m58.693s  |1m28.551s|0m7.782s|0m13.810s|
|*Reads retained*             |100.000|93.346   |92.514     |92.997   |93.002  |92.994   |
|*Bases retained*             |100.000|83.939   |74.436     |73.970   |74.268  |74.186   |
|*Perfectly correct (Reads)*  |49.970 |75.428   |80.922     |96.256   |94.849  |95.784   |
|*Perfectly correct (Bases)*  |49.970 |72.780   |86.426     |96.035   |93.900  |95.099   |
|*Incorrect (Reads)*          |50.030 |24.572   |19.078     |3.744    |5.151   |4.216    |
|*Incorrect (Bases)*          |50.030 |27.220   |13.574     |3.965    |6.100   |4.901    |
|*Adaptors remaining (Reads)* |50.030 |23.786   |5.846      |1.830    |3.866   |2.798    |
|*Adaptors remaining (Bases)* |25.182 |11.439   |0.422      |0.049    |0.193   |0.105    |
|*Non-adaptor removed (Reads)*|0.000  |0.786    |13.231     |1.914    |1.285   |1.418    |
|*Non-adaptor removed (Bases)*|0.000  |0.019    |0.218      |0.566    |0.308   |0.325    |

All numbers (except times) represent percentages, Reads and Bases retained are as percentage of the input "dirty" reads, while all the rest are as percentage of each method's number of reads and bases retained. 

Even though `fastp` is slightly more accurate than `bbduk.sh`, it takes >6x longer to finish. More importantly, when cleaning real unaltered data (without extra adapters added) we noticed that `fastp` sometimes mistakes real genomic sequence as adaptor, and we are not the only that noticed this, you can read about it in this [GitHub Issue](https://github.com/OpenGene/fastp/issues/160).

Otherwise, `bbduk.sh` still beats `Cutadapt` and `Trimmomatic` on every metric. The settings we use in `Captus` increase the sensitivity of `bbduk.sh` to a level very close to `fastp` while not increasing too much runtime. In the future, if the developers of `fastp` fix the false detection problem we could incorporate it in the pipeline.
___