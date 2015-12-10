# the main script to run in order to annotate a list of de novos with predicted TF binding sites
# author: Patrick Short <pjs90@cam.ac.uk>

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
  make_option("--de_novos", default="../data/DDD_noncoding_for_denovoTF.txt",
              help="Pass the genomic regions that should be annotated with predicted TF binding sites."),
  make_option("--pwms", default=FALSE,
              help="Pass an RData file with list of manually curated PWMs."),
  make_option("--out", default="../results/JASPAR_tfbs_annotated_de_novos.txt",
              help="Set location to save the output of JASPAR-annotated de novos."),
  make_option("--min_score", default="95%",
              help="Set min score for the JASPAR annotations - less than 95% may yield unrealistic number of predicted binding sites."),
  make_option("--start_row", default=0,
              help="If de novos is very large, chunking can be used (pass row to start analysis on)."),
  make_option("--chunk_size", default=0,
                  help="If de novos is very large, chunking can be used (pass row to stop analysis on)."),
  make_option("--verbose", action="store_true", default=FALSE,
              help="Print extra output advising the user of progression through the code.")
)

args <- parse_args(OptionParser(option_list=option_list))
### check that the input file has columns "unique_id", "chr", "pos", "ref", "alt". if no "unique_id", create one

print(args)

de_novos <- read.table(args$de_novos, sep = "\t", header = TRUE)

if (args$start_row != 0 & args$chunk_size != 0){
  write(sprintf("Analyzing chunk of de novos/rare variant file that has been passed (row %i to row %i inclusive).", args$start_row, args$start_row + args$chunk_size - 1), stderr())
  last_row = min((args$start_row + args$chunk_size - 1), nrow(de_novos))
  de_novos = de_novos[args$start_row:last_row, ]
} else {
  write("Analyzing entire input file of de novos/rare variants. If you would prefer to chunk, pass --start_row and --chunk_size arguments.", stderr())
}
  
# remove indels from de novo file - TODO: add support to analyze indels
write("Removing any indels from the de novos/rare variant file that has been passed...", stderr())
de_novos = de_novos[nchar(as.character(de_novos$ref)) == 1 & nchar(as.character(de_novos$alt)) == 1,]

reqd_columns <- c("chr", "pos", "ref", "alt")

# throw error if any required column is missing
if (!all(reqd_columns %in% colnames(de_novos))){
  stop("One or more of the required column names missing from input de novo file. Requires: \"chr\", \"pos\", \"ref\", \"alt\"")
}

# unique_id not present, add column in the form chr:posREF>ALT
if (!("unique_id" %in% colnames(de_novos))){
  de_novos$unique_id <- paste0(de_novos$chr, ":", de_novos$pos, de_novos$ref, ">", de_novos$alt)
}

### load the full JASPAR position-weight matrix list using all sources (SELEX, ChIP-seq, protein-binding microarray (PBM))
# load full PWM list
if ( args$verbose ) {
  write("Loading JASPAR position weight matrices...", stderr())
}


# updated 8th of December 2015 to take pre-curated PWM list
if (args$pwms != FALSE){  # switch to reduced set of TFs if requested
  load(args$pwms)  # expects to load as pwm_list
} else {
  # NOTE: db is initialized to ../data/myMatrixDb.sqlite after build.R is run
  pwm_options = list("species" = 9606, "all_versions" = TRUE, "matrixtype" = "PWM") # 9606 = "homo sapiens"
  pwm_list = getMatrixSet(JASPAR2016, pwm_options)
}

# get longest motif in dataset (to set bound for sequence query)
max_motif_length = max(sapply(pwm_list, function(t) ncol(t@profileMatrix)))

### get the sequence context for each de novo +/- 20 base pairs upstream and downstream. There is no need to specify strandedness as the JASPAR PWM will
### account for both + and - strand predicted binding events.

# use BSgenomes getSeq method to get +/- max_motif_length (around 20 bp) from each de novo
if ( args$verbose ) { write("Getting sequence context for each de novo...", stderr()) }

m = max_motif_length
ref_seqs = get_sequence(paste0("chr", de_novos$chr), de_novos$pos - m, de_novos$pos + m)
names(ref_seqs) = de_novos$unique_id

### scan every sequence against the full JASPAR list and keep track of any binding events that intersect with the de novo position.

if ( args$verbose ) { write("Scanning sequence surrounding each de novo for predicted transcription factor binding event...", stderr()) }

scans = mapply(LOBGOB_scan, ref_seqs, m+1, as.character(de_novos$ref), as.character(de_novos$alt), MoreArgs = list("pwm_list" = pwm_list, "min.score" = args$min_score))
# count number of hits per de novo

if ( args$verbose ) { write("LOBGOB scan complete.", stderr()) }

hits_per_de_novo = sapply(scans, function(s) ifelse(is.null(s), 0, nrow(s)))
hits_TFBS = hits_per_de_novo > 0

if ( args$verbose ) { write("Finished scanning, binding dataframes from scan results.", stderr()) }

# rbind the dataframes from scans
scores = do.call(rbind, scans)
tfbs_name = sapply(scores$jaspar_internal, function(j) pwm_list[j][[1]]@name)
scores = cbind(tfbs_name, scores)

### calculate the increase/decrease in binding affinity
if ( args$verbose ) { write(sprintf("Analyzing change in information content for ref vs. alt on predicting TF binding events..."), stderr()) }

# combine de novo IDs and scores
annotated_dn = cbind(de_novos[rep(seq(nrow(de_novos)), hits_per_de_novo),], scores)

# annotated_dn motif start is determined in LOBGOB_scan as relative to de novo position
annotated_dn$motif_start = annotated_dn$motif_start + annotated_dn$pos
annotated_dn$motif_end = annotated_dn$motif_end + annotated_dn$pos

if ( args$verbose ) { write(sprintf("Number of de novos passed to input: %i", nrow(de_novos)), stderr()) }
if ( args$verbose ) { write(sprintf("Number of de novos intersecting at least one TFBS: %i", sum(hits_TFBS), stderr())) }
if ( args$verbose ) { write(sprintf("Total number of predicted TFBS perturbations: %i", nrow(annotated_dn), stderr())) }

write.table(annotated_dn, file = args$out, row.names = FALSE, sep = "\t", col.names = TRUE, quote = FALSE)

if ( args$verbose ) { write(sprintf("Finished!"), stderr()) }
