# denovoTF
Calculating change in binding affinity for predicted TF binding sites as a result of de novo mutations in non-coding regions. 

Uses JASPAR database for TF binding prediction.

## Getting Set UP
Running build.R will install hg19 annotation (for pulling sequence context) and install JASPAR2014 and TFBSTools Bioconductor packages.

Dependencies: BSGenome, Biostrings, TFBSTools, JASPAR2014, optparse

## Analyzing a set of De Novos
A tab-delimited file of de novo (non-coding) mutations is the only required input. The only required columns are: chr, pos, ref, alt and for now only SNPs are supported (no indels). Unless specified, hg19 is assumed to be the human genome annotation used.

The genomic coordinates will be used to pull the sequence context which will be run against ALL of the human genome transcription factor binding site position weight matrices (PWMs) in the JASPAR database. Binding affinity for transcription factors (expressed in -log10 scale) is calculated in the neighborhood of each of the de novos and any predicted binding event is returned in the output file with the ref_score (binding affinity with the ref nucleotide) and alt_score (binding affinity with the alt nucleotide).

The output file will have one row for each de novo + TFBS predicted binding 'event'. As some de novos may be predicted to disrupt multiple binding sites (and some to disrupt none) the output file will not look identical in row/column space to the input file. the --min_score parameter can be adjusted (defaults to 95%) to accept more or less stringent binding predictions.

Output file example:

unique_id chr       pos ref alt tfbs_name jaspar_internal ref_score alt_score
3:180462583T>C   3 180462583   T   C     GATA2        MA0036.1  5.931088 -2.124195
13:95600362T>A  13  95600362   T   A       YY1        MA0095.1  7.934169  1.507904

