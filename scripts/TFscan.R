# the main script to run in order to annotate a set of regions with predicted TF binding sites
# annotate a set of regions with all instances of a single motif over specified --min_score

# this will primarily be used as training data for msCentipede

# author: Patrick Short <ps14@sanger.ac.uk>

### dependencies
library(optparse)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TFBSTools)
library(JASPAR2014)

source("../R/core.R")

### command line options
option_list <- list(
  make_option("--regions", default="../data/noncoding_regions.txt",
              help="Pass the genomic regions that should be annotated with predicted TF binding sites."),
  make_option("--jaspar_motif_id", default="MA0885.1",
              help="Specify which jaspar motif to scan against."),
  make_option("--pwms", default="../data/brain_tf_pwms.RData",
              help="Specify which jaspar motif to scan against."),
  make_option("--out", default=NULL,
              help="Set directory to save the output of scan (looks like MA####.#)"),
  make_option("--min_score", default="95%",
              help="Set min score for the JASPAR annotations - less than 95% may yield unrealistic number of predicted binding sites."),
  make_option("--verbose", action="store_true", default=FALSE,
              help="Print extra output advising the user of progression through the code.")
)

args <- parse_args(OptionParser(option_list=option_list))
### check that the input file has columns "unique_id", "chr", "pos", "ref", "alt". if no "unique_id", create one

regions <- read.table(args$regions, sep = "\t", header = TRUE)

reqd_columns <- c("chr", "start", "end")

# throw error if any required column is missing
if (!all(reqd_columns %in% colnames(regions))){
  stop("One or more of the required column names missing from input de novo file. Requires: \"chr\", \"start\", \"end\"")
}


### load the full JASPAR position-weight matrix list using all sources (SELEX, ChIP-seq, protein-binding microarray (PBM))
# load full PWM list
if ( args$verbose ) {
  write("Loading JASPAR position weight matrices from database...", stderr())
}

if (args$pwms != FALSE){  # switch to reduced set of TFs if requested
  load(args$pwms)  # expects to load as pwm_list
} else {
  # NOTE: db is initialized to ../data/myMatrixDb.sqlite after build.R is run
  pwm_options = list("species" = 9606, "all_versions" = TRUE, "matrixtype" = "PWM") # 9606 = "homo sapiens"
  pwm_list = getMatrixSet(JASPAR2016, pwm_options)
}

# restrict to only jaspar motif specified
pwm_list = pwm_list[[args$jaspar_motif_id]]

# use BSgenomes getSeq method to get sequence in each region
if ( args$verbose ) { write("Getting sequence context for each region...", stderr()) }

if ("seq" %in% colnames(regions)){
  seqs = sapply(as.character(regions$seq), function(s) DNAString(s))
  names(seqs) = paste(regions$chr, regions$start, regions$end, sep = ".")
} else {
  seqs = as.character(get_sequence(paste0("chr", regions$chr), regions$start, regions$end))
  names(seqs) = paste(regions$chr, regions$start, regions$end, sep = ".")
}

### scan every sequence against the full JASPAR list and keep track of any predicted binding events
if ( args$verbose ) { write("Scanning each sequence for predicted transcription factor binding events...", stderr()) }

scanned_regions = mapply(searchSeq, seqs, MoreArgs = list("seqname"="ref_sequence", "x" = pwm_list, "min.score"=args$min_score, "strand"="*"))
hit_count = as.numeric(sapply(scanned_regions, function(r) length(r@score)))
contains_hit = hit_count > 0

motif_start = as.numeric(unlist(sapply(scanned_regions, function(r) r@views@ranges@start)))
motif_start = rep(regions$start, hit_count) + motif_start - 1
motif_width = as.numeric(unlist(sapply(scanned_regions, function(r) r@views@ranges@width)))
motif_end = motif_start + motif_width - 1

motif_strand = as.character(unlist(sapply(scanned_regions, function(r) r@strand)))
motif_score = as.numeric(unlist(sapply(scanned_regions, function(r) r@score)))

# output will contain chromosome, motif start, motif end, strand, and score
motif_df = data.frame(Chr = rep(regions$chr, hit_count), Start = motif_start, End = motif_end, Strand = motif_strand, PwmScore = motif_score)

if (all(grepl("^chr", motif_df$Chr)) == FALSE){
  motif_df$Chr = paste0("chr", motif_df$Chr)
}

if ( args$verbose ) { write(sprintf("Finished!"), stderr()) }

if (!is.null(args$out)){
  write.table(motif_df, file = args$out, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
} else {
  write.table(motif_df, file = sprintf("./JASPAR_scan_%s.txt", args$jaspar_motif_id), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
}
