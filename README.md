# denovoTF
Calculating change in binding affinity for predicted TF binding sites as a result of de novo mutations in non-coding regions. 

Uses JASPAR database for TF binding prediction.

## Setting up
Running build.R will install hg19 annotation (for pulling sequence context) and install JASPAR2014 and TFBSTools Bioconductor packages.

## Analyzing a set of De Novos
A tab-delimited file of de novo (non-coding) mutations is the only required input. The only required columns are: chr, pos, ref, alt and for now only SNPs are supported (no indels). Unless specified, hg19 is assumed to be the human genome annotation used.

The genomic coordinates will be used to pull the sequence context which will be run against ALL of the human genome transcription factor binding site position weight matrices (PWMs) in the JASPAR database. Binding affinity for transcription factors (expressed in -log10 scale) is calculated in the neighborhood of each of the de novos and any predicted binding event is returned in the output file with the ref_score (binding affinity with the ref nucleotide) and alt_score (binding affinity with the alt nucleotide).

The output file will have one row for each de novo + TFBS predicted binding 'event'. As some de novos may be predicted to disrupt multiple binding sites (and some to disrupt none) the output file will not look identical in row/column space to the input file. the --min_score parameter can be adjusted (defaults to 95%) to accept more or less stringent binding predictions.

Output file columns:
unique_id, chr, pos, ref, alt, tfbs_name, jaspar_internal, ref_score, alt_score

