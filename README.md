# denovoLOBGOB
denovoTF has been replaced by denovoLOBGOB (short for de novo Loss of Binding/Gain of Binding).

The reference and alternate sequence for a list of variants (chr, pos, ref, alt) is generated and scanned against all JASPAR human transcription factor binding site position weight matrices. A set of TFs with predicted binding affinity will be generated for both the ref and alt sequence. As a result, loss of binding (LOB), gain of binding (GOB), and silent mutations can be determined.

Uses JASPAR2014 database (and associated R package) is used for TF binding affinity prediction.

## Getting Set UP
Running build.R will install hg19 annotation (for retrieving sequence context) and install JASPAR2014 and TFBSTools Bioconductor packages.

Dependencies: BSGenome, Biostrings, TFBSTools, JASPAR2014, optparse

## Analyzing a set of De Novos
A tab-delimited file of de novo (non-coding) mutations is the only required input. The only required columns are: chr, pos, ref, alt and for now only SNPs are supported (no indels). Unless specified, hg19 is assumed to be the human genome annotation used.

The genomic coordinates will be used to retrieve the sequence context which will be run against ALL of the human genome transcription factor binding site position weight matrices (PWMs) in the JASPAR database. Binding affinity for transcription factors (expressed in -log10 scale) is calculated for both the reference and alterate (mutated) sequence around the de novos and any predicted binding event is returned in the output file with the ref_score (binding affinity with the ref nucleotide) and alt_score (binding affinity for the alt nucleotide).

The output file will have one row for each de novo + TFBS predicted binding 'event'. As some de novos may be predicted to disrupt multiple binding sites (and some to disrupt none) the output file will not look identical in row/column space to the input file. the --min_score parameter can be adjusted (defaults to 95%) to accept more or less stringent binding predictions.

Output file example:

|   unique_id    | chr |    pos    | ref | alt | tfbs_name | jaspar_internal | ref_score | alt_score |
| -------------- | --- | --------- | --- | --- | --------- | --------------- | --------- | --------- |
| 3:180462583T>C |  3  | 180462583 |  T  |  C  |   GATA2   |     MA0036.1    |  5.931088 | -2.124195 |
| 13:95600362T>A |  13 |  95600362 |  T  |  A  |   YY1     |     MA0095.1    |  7.934169 | 1.507904  |

Any additional columns that are passed i.e. patient_is_diagnosed will be preserved.


