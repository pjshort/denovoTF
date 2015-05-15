# the main script to run in order to annotate a list of de novos with predicted TF binding sites

# INPUT:
# de novo file with mandatory columns: "unique_id", "chr", "pos", "ref", "alt" - if no unique_id is provided, then one will be made in the form chr:posref>alt
# optional: --tf_list -> list of transcription factor motifs (by jaspar_internal id e.g. MA0098.1) in single column with NO HEADER
# optional: --hg_version -> build for human genome (this is essential to get right as it determines the sequence context used to scan)
# optional: --min_score -> minimum score cutoff to call binding event (default to 95% i.e. adjusted p-val <0.05)
# optional: --summary -> creates an additional summary file summarizing hits per TFBS, global TFBS hit rate, etc. (NOT YET IMPLEMENTED)

# OUTPUT:
# de novo output file with columns "unique_id", "chr", "pos", "ref", "alt", "tf_name", "jaspar_internal", "ref_score", "alt_score" 
# with ONE ROW PER TF binding event. the output file will likely have more rows than the input file (many more if score threshold is low)

# documentation notes:
# a triple hash (### description xyz) denotes a 'section header' in the code while (# comments..) denotes a more simple comment

### dependencies
library(optparse)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TFBSTools)
library(JASPAR2014)

source("../R/core.R")


### command line options
option_list <- list(
  make_option("--de_novos", default="../data/de_novo_filtered.txt",
              help="Pass the genomic regions that should be annotated with predicted TF binding sites."),
  make_option("--tf_list", default=FALSE,
              help="Pass a list of TFs to be run against regions."),
  make_option("--out", default="../results/",
              help="Set location to save the output of JASPAR-annotated de novos."),
  make_option("--min_score", default="95%",
              help="Set min score for the JASPAR annotations - less than 95% may yield unrealistic number of predicted binding sites."),
  make_option("--verbose", action="store_true", default=FALSE,
              help="Print extra output advising the user of progression through the code.")
)

args <- parse_args(OptionParser(option_list=option_list))

### check that the input file has columns "unique_id", "chr", "pos", "ref", "alt". if no "unique_id", create one

dn <- read.table(args$de_novos, sep = "\t", header = TRUE)

reqd_columns <- c("chr", "pos", "ref", "alt")

# throw error if any required column is missing
if (!all(reqd_columns %in% colnames(dn))){
  stop("One or more of the required column names missing from input de novo file. Requires: \"chr\", \"pos\", \"ref\", \"alt\"")
}

# unique_id not present, add column in the form chr:posREF>ALT
if (!("unique_id" %in% colnames(dn))){
  dn$unique_id <- paste0(dn$chr, ":", dn$pos, dn$ref, ">", dn$alt)
}

### load the full JASPAR position-weight matrix list using all sources (SELEX, ChIP-seq, protein-binding microarray (PBM))
# load full PWM list
if ( args$verbose ) {
  write("Loading JASPAR position weight matrices from database...", stderr())
}

# NOTE: db is initialized to ../data/myMatrixDb.sqlite after build.R is run
pwm_options = list("species" = 9606, "all_versions" = TRUE, "matrixtype" = "PWM") # 9606 = "homo sapiens"
pwm_list = getMatrixSet(JASPAR2014, pwm_options)

# TODO: update this to take names of transcription factors and slice PWM accordingly
if (args$tf_list != FALSE){  # switch to reduced set of TFs if requested
  TFBS_to_scan = read.table(args$tf_list, header = FALSE) # get single column of requested IDs
  pwm_list = pwm_list[unique(TFBS_to_scan[,1])] # assumed to be a single column
}

# get longest motif in dataset (to set bound for sequence query)
max_motif_length = max(sapply(pwm_list, function(t) ncol(t@profileMatrix)))

### get the sequence context for each de novo +/- 20 base pairs upstream and downstream. There is no need to specify strandedness as the JASPAR PWM will
### account for both + and - strand predicted binding events.

# use BSgenomes getSeq method to get +/- max_motif_length (around 20 bp) from each de novo
if ( args$verbose ) { write("Getting sequence context for each de novo...", stderr()) }

m = max_motif_length
seqs = mapply(function(chr, start, stop) getSeq(Hsapiens, chr, start, stop), paste0("chr", dn$chr), dn$pos - m, dn$pos + m)
names(seqs) = dn$unique_id
#TODO: put in test that makes sure ref/alt match the sequence returned by BSgenome getSeq

### scan every sequence against the full JASPAR list and keep track of any binding events that intersect with the de novo position.

if ( args$verbose ) { write("Scanning sequence surrounding each de novo for predicted transcription factor binding event...", stderr()) }

rel_positions <- rep(m+1, length(seqs)) # relative position should be the middle of each seq object (max motif length +1)
scanned_regions = scan_regions(seqs, rel_positions, pwm_list, min.score = args$min_score)
hits_per_de_novo = sapply(scanned_regions, length)

### calculate the increase/decrease in binding affinity
if ( args$verbose ) { write(sprintf("Analyzing change in information content for ref vs. alt on predicting TF binding events..."), stderr()) }

# flatten the list of hits - used rep with hits_per_de_denovo to repeat alt and ref allele as much as necessary for de novos
# which perturb multiple motifs
r = unlist(scanned_regions[hits_per_de_novo > 0])
rel_positions = rep(m+1, sum(hits_per_de_novo))
ref = as.character(rep(dn$ref, hits_per_de_novo))
alt = as.character(rep(dn$alt, hits_per_de_novo))

scores <- mapply(binding_change, r, rel_positions, ref, alt, MoreArgs = list("min.score" = "95%"))
scores <- t(scores)  # flip to columns (ref_score, alt_score)

### reformat the results into annotated de novo output file and exit
# results will have similar form as input with additional columns and additional rows where a de novo hits more than one TF binding site
# unique_id, chr, pos, ref, alt, tfbs_name, jaspar_internal, ref_score, alt_score

# process remaining columns for data frame
unique_id <- rep(dn$unique_id, hits_per_de_novo)
chr <- rep(dn$chr, hits_per_de_novo)
pos <- rep(dn$pos, hits_per_de_novo)
tfbs_name <- unlist(sapply(r, function(s) s@pattern@name))
jaspar_internal <- unlist(sapply(r, function(s) s@pattern@ID))
ref_score <- scores[,1]
alt_score <- scores[,2]

# create annotated de novo data frame
annotated_dn <- data.frame(unique_id, chr, pos, ref, alt, tfbs_name, jaspar_internal, ref_score, alt_score)
rownames(annotated_dn) <- NULL

if ( args$verbose ) { write(sprintf("Number of de novos passed to input: %i", nrow(dn)), stderr()) }
if ( args$verbose ) { write(sprintf("Number of de novos intersecting at least one TFBS: %i", sum(hits_per_de_novo > 0), stderr())) }
if ( args$verbose ) { write(sprintf("Total number of predicted TFBS perturbation: %i", nrow(annotated_dn), stderr())) }

write.table(annotated_dn, file = paste0(args$out, "/JASPAR_tfbs_annotated_de_novos.txt"), row.names = FALSE, sep = "\t", col.names = TRUE)

if ( args$verbose ) { write(sprintf("Finished!"), stderr()) }