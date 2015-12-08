# create the brain transcription factor position-weight-matrix

# steps to prepare data for script
# prepare list of genes with GO terms for transcription factor and brain development
# manually download PFM from JASPAR for genes present in homo sapiens data base
# combine these PFMs into single file
# for each entry first line is header, next four lines describe the PFM

# these PWMs will be saved as an .RData file for later use

library(TFBSTools)
library(JASPAR2014)

pfms = "../data/brain_transcription_factors.pfm"
text = readLines(pfms, n=1000)

i = 1 # count lines - first is header, 2-5 are matrix
pwm_list = vector(mode = "list", length = length(text)/5)
for (l in text) {
  if (i %% 5 == 1){ # header line - look for i mod 5 == 1
    info = strsplit(l, "\t")
    jaspar_id = substr(info[[1]][1], 2, 9)  # drop the > and look for 8 letter motif in form MA####.#
    gene_name = info[[1]][2]
  } else if (i %% 5 == 2) { # first line of matrix, order is A C G T
    entries = substr(l, 5, nchar(l))  # take off the "A "
    nums = as.integer(strsplit(entries, "[[:space:]]")[[1]])
    nums = nums[!is.na(nums)]
    mat = nums
  } else {
    entries = substr(l, 5, nchar(l))  # take off the "A "
    nums = as.integer(strsplit(entries, "[[:space:]]")[[1]])
    nums = nums[!is.na(nums)]
    mat = c(mat, nums)
  }
  
  if (i %% 5 == 0){  
    # build position-frequency-matrix
    pfm = PFMatrix(ID=jaspar_id, name=gene_name, bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                profileMatrix=matrix(mat, byrow=TRUE, nrow=4, dimnames=list(c("A", "C", "G", "T"))))
    
    # convert to position-weight-matrix
    pwm = toPWM(pfm, type="log2probratio", pseudocounts=0.8, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
    
    pwm_list[[i/5]] = pwm
    
  }
  
  i = i + 1
}

save(pwm_list, file = "../data/brain_tf_pwms.Rdata")
  