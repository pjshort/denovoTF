# denovoLOBGOB
denovoTF has been replaced by denovoLOBGOB (short for de novo Loss of Binding/Gain of Binding).

The reference and alternate sequence for a list of variants (chr, pos, ref, alt) is generated and scanned against all JASPAR human transcription factor binding site position weight matrices. A set of TFs with predicted binding affinity will be generated for both the ref and alt sequence. As a result, loss of binding (LOB), gain of binding (GOB), and silent mutations can be determined.

Uses JASPAR2014 database (and associated R package) is used for TF binding affinity prediction.

```
Rscript denovoLOBGOB.R --de_novos=/path/to/denovos --min_score=95% --verbose --out=/path/to/annotated_denovos
```

## Getting Set Up
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

## Generating Simulation Data

simulateDN.R (in scripts folder) can be use to generate simulated de novos. The random selection of location, alt allele is non-uniform - instead, it is based on a background trinucleotide mutation rate first described by Samocha et. al (http://www.nature.com/ng/journal/v46/n9/abs/ng.3050.html). For any trinucleotide, the poisson lambda parameter of a mutation at the middle base is determined by the trinucleotide mutation table (i.e. Lambda(ATG is mutated) = Lambda(ATG -> ACG) + Lambda(ATG -> AGG) + Lambda(ATG -> AAG).

simulateDN.R requires:
--n_snps -> number of snps to simulate
--n_probands -> number of probands to assign snps to
--regions -> regions in which the SNPs will be simulated (given by columns chr, start, stop)
--iterations -> how many synthetic data sets to create
--n_chunks -> number of different files to save them to
--base_name -> base name that the chunks will have (or simply file name if only 1 chunk) e.g. --base_name=~/results/sim_data will save as ~/results/sim_data1.txt, ~/results/sim_data2.txt ...

```
Rscript simulateDN.R --n_snps=200 --n_probands=150 --regions=/path/to/regions --iterations=10 --n_chunks=2 --base_name=/path/to/chunk --verbose
```

--n_chunks will only be useful if you plan to do annotation in parallel (which is advised if a cluster is available). 1000 iterations in a single file with 453 snps in 425 probands is approximately 16MB in size. The output file will have an 'iterations' column that can be used to split the data for comparing with a observations from real data, for instance.

