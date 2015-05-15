# core routines for TF binding analysis

get_sequence <- function(chr, start, stop, version = "hg19") {
  
  # input: (multiple) chr, start, stop, hg version (defaults to hg19)
  # output: list of sequences as DNAStrings object for each input
  
  if (version == "hg19"){
    library(BSgenome.Hsapiens.UCSC.hg19)
  } else if (version == "hg18"){
    library(BSgenome.Hsapiens.UCSC.hg18) # TODO need to add download of hg18 to build.R
  }
  
  if (!all(grepl(pattern = "^chr", chr))){  # assert that chromosome column have chr in front
    warning("Not all entries in the chromosome column start with \"chr\" - try reformatting this column e.g. \"chrX\" instead of \"X\" with paste0(\"chr\",chr_number")
    chr = paste0("chr", chr)
  }
  
  seqs = getSeq(Hsapiens, chr, start, stop)
  return(seqs)
}

### Annotated sequences with predicted TF binding

single_sequence_coverage <- function(seq, rel_pos, pwm_list, min.score = "95%"){
  
  # input: DNAString sequence, list of PWMs to query, min.score (optional)
  # returns: site
  # returns all regions predicted to have TFB affinity >= min.score
  
  # TODO: write test for this section
  
  # scan full list of PWMs against the sequence provided
  site_seq_list = searchSeq(pwm_list, seq, seqname="ref_sequence", min.score=min.score, strand="*")
  
  # keep only the TFs that have a hit greater than min score
  interval_hits = site_seq_list[which(sapply(site_seq_list, length) > 0)]
  
  # keep only the TFs that have a hit in region that overlaps with the de novo
  overlaps_dn = sapply(interval_hits, function(t) t[(rel_pos >= start(t@views@ranges)) & (rel_pos <= end(t@views@ranges))])
  
  # filter list to remove the empty TFs
  pos_hits = overlaps_dn[which(sapply(overlaps_dn, length) > 0)]

  return(pos_hits) # returns a (possibly empty) list of SiteSet objects
}

scan_regions <- function(sequences, rel_positions, pwm_list, min.score = "95%"){
  
  # input: vector of sequences, vector of relative positions of de novo within sequence, list of PWMs to query, min.score (optional)
  # output: list with one element for each pair of seq, rel_pos that contains predicted de novo binding events (if any)
  
  scan_results = mapply(single_sequence_coverage, sequences, rel_positions, MoreArgs = list("pwm_list" = pwm_list))
  
  return(scan_results)
}


### calculate ref vs. alt change in binding

binding_change <- function(site_set, rel_pos, ref, alt, min.score = "95%"){
  
  # input: SiteSet object (TFBSTools), relative position of de novo, ref, alt
  # returns: ref_score, alt_score in 2x1 vector
  
  # TODO: reformulate to allow indels
    
  # get the position of the de novo within the motif
  motif_pos = rel_pos - start(site_set@views@ranges) + 1
  
  # get position weight matrix
  pwm = site_set@pattern@profileMatrix
  
  # get score with reference allele
  ref_score = site_set@score
  
  # compute change with alt allele and add to ref to calculate new alt allele binding score
  change = as.numeric(pwm[alt, motif_pos] - pwm[ref, motif_pos])
  alt_score = ref_score + change
  
  return(cbind(ref_score, alt_score))
}