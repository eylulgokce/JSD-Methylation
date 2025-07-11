#library(shannon)
source("R/core.R")
source("R/gpf_utils.R")
source("R/io.R")
source("R/constants.R")
source("R/estimators.R")
source("R/preprocessing.R")

meta <- read.csv("/Users/eylul/Desktop/JSD-Methylation/shannon/test_data/CHG_chr1.csv")

divergence(
  sample = meta,
  chrom = "1",
  data_columns = list(c(5, 6)),
  outfile = "output_chr1_r.txt",
  chunksize = 1000
)
