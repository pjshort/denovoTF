# build script for the denovo_TF pipeline that will install dependencies if missing and build database for the JASPAR position weight matrices

# source bioconductor and load BSgenomes
source("http://bioconductor.org/biocLite.R")
biocLite("BSgenomes")
biocLite("BSgenome.Hsapiens.UCSC.hg19") # load the hg19 genome build - takes a little while as it is approx. 800 MB
biocLite("TFBSTools")
biocLite("JASPAR2014")
library(TFBSTools)
library(JASPAR2014)

# initialize JASPAR database and save in data
db = "../data/myMatrixDb.sqlite"
initializeJASPARDB(db)


