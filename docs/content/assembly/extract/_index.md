+++
title = "Extract"
weight = 20
pre = "<b>3. </b>"
+++

The next step in the workflow is to search markers within the assemblies and extract them. For that you need to provide `Captus` reference sequence(s) of the marker(s) you want to extract. You can provide coding sequences in either nucleotide or aminoacid, in this case, protein extraction will be performed with [`Scipio`](https://www.webscipio.org/). If you want to extract any other DNA marker (entire genes with introns, full mRNAs, non-coding regions, etc.) you provide your reference(s) in nucleotide format, in this case the extraction is performed using 
[`BLAT`](http://hgdownload.soe.ucsc.edu/admin/exe/) and our own code to stitch and extract partial hits if necessary. 
