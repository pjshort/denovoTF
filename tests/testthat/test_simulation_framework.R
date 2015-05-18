# tests for core routines underlying simulateDN

source("../R/core.R")
source("../R/simulate.R")
source("../R/mutation_null_model.R")
library(testthat)

context("Testing simulation framework")

test_that("probabilities generated from regions match expectations", {
  
  regions = read.table(header=TRUE, colClasses=c("character", "numeric", "numeric"), text="
      chr start stop
      chr16 85662461 85662511
      chr19 32662337 32662387
      chr9 26205098 26205198
      chr3 180462583 180462633")
  
  seqs = get_sequence(regions$chr, regions$start, regions$stop, version = "hg19")
  regions$p_snp_null <- sapply(seqs, p_sequence)
  
  expect_more_than(object = regions$p_snp_null[3], expected = regions$p_snp_null[4]) # twice as many bp region should have greater chance of mut
  expect_less_than(regions$p_snp_null[3], 1e-5)
})

test_that("simulation works", {
  
  regions = read.table(header=TRUE, colClasses=c("character", "numeric", "numeric"), text="
      chr start stop
      chr16 85662461 85662511
      chr19 32662337 32662387
      chr9 26205098 26205198
      chr3 180462583 180462633")
  
  # note that these extras are necessary, otherwise simulate_de_novos fails
  regions$region_id <- paste(regions$chr, regions$start, regions$stop, sep = ".")
  regions$seq <- as.character(get_sequence(regions$chr, regions$start, regions$stop, version = "hg19"))
  regions$p_snp_null <- sapply(regions$seq, p_sequence)
  regions$p_relative <- regions$p_snp_null/sum(regions$p_snp_null)
  
  seq_probabilities = relative_seq_probabilities(regions) # run these two lines to regenerate
  
  set.seed(42) 
  s <- simulate_de_novos(regions, seq_probabilities, 2, 1)
  
  expect_equal(s[1,"pos"], 180462622) # first simulated de novo is that position
  expect_equal(as.character(s[2,"chr"]), "chr3") # second simulated de novo is on chr3
  
  # full result should be
  #> s
  #chr       pos ref alt iteration
  #1 chr3 180462622   C   A         1
  #2 chr3 180462591   G   A         1
  
})