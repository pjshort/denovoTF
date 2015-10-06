# combine set of annotated dataframes with the same headers (simulation data)

library(optparse)

### command line options
option_list <- list(
  make_option("--base_name", default="/tmp/denovoLOBGOB_sim_chunks/sim_data",
              help="Pass the base name of the chunks."),
  make_option("--n_chunks", default=1000,
              help="Number of chunks that should be combined."),
  make_option("--out", default="../results/simulated_de_novos_JASPAR_tfbs_annotation.txt",
              help="Set location to save the output of combined chunks.")
)

args <- parse_args(OptionParser(option_list=option_list))

chunks = vector("list", args$n_chunks)
for (i in seq(args$n_chunks)){
  fname = sprintf("%s.%i.txt", args$base_name, i) # will look for fname e.e.g /tmp/denovoLOBGOB_sim_chunks/sim_data.n.out for n 1-n_chunks
  if (file.exists(fname)){
    chunks[[i]] = read.table(fname, sep = "\t", header = TRUE)
  } else {
    chunks[[i]] = NULL
  }
}

chunks = chunks[!sapply(chunks, is.null)]

annotated_sim <- do.call(rbind, chunks)
annotated_sim$diff <- annotated_sim$ref_score - annotated_sim$alt_score

write.table(annotated_sim, file = args$out, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

