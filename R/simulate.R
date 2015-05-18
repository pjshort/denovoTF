# create de novo files by simulation

relative_seq_probabilities <- function(regions){
  
  # get probability of mutation at each point in sequence (based on sequence context) and normalize
  seq_probabilities = sapply(as.character(regions$seq), function(s) p_position(s, normalize = TRUE), USE.NAMES = FALSE)
  names(seq_probabilities) = regions$region_id
  return(seq_probabilities)
}

simulate_de_novos <- function(regions, seq_probabilities, snp_total, iteration) {
  
  # takes vector of region_ids to consider, list matching region_id to sequence relative probability, and total snps to simulate
  # samples snp_total region ids (with replacement). for each region, specific relative position will be randomly sampled
  
  region_ids = sample(regions$region_id, snp_total, replace = TRUE, prob = regions$p_relative)
  rel_pos = sapply(as.character(region_ids), function(id) sample(seq(1,length(seq_probabilities[id][[1]])), 1, prob = seq_probabilities[id][[1]]))
  coords = do.call(rbind, strsplit(region_ids, "\\.")) # chr, start, stop
  pos = as.integer(coords[,2]) + as.integer(rel_pos)
  
  # once the position within the sequence has been chosen, want to choose the mutated (alt) base based on null model
  ref_tri = as.character(mapply(function(s, p) substr(s, p-1, p+1), regions[match(region_ids, regions$region_id), "seq"], rel_pos))
  
  # replace any di nucleotides if simulation fell on edge TODO write test to check this works
  di_idx = as.integer(which(sapply(ref_tri, nchar) < 3))
  if (length(di_idx) > 0) {
    ref_tri[di_idx] = as.character(get_sequence(chr = paste0("chr", as.character(coords[di_idx, 1])), start = pos[di_idx] - 1, stop = pos[di_idx] + 1))
  }
  
  ref = as.character(sapply(ref_tri, function(s) substr(s,2,2)))
  alt = as.character(sapply(ref_tri, sample_alt))
  
  sim_dn = data.frame("chr" = coords[,1], "pos" = pos, "ref" = ref, "alt" = alt, "iteration" = iteration)
  
  return(sim_dn)  # data frame
}
