library(EnrichShuf)
library(dplyr)
library(data.table)

answer_dat <- readRDS(
  system.file("extdata/ATAC_within.rds.gz", package = "EnrichShuf")
)

target_dat <- ObsExpObj(
  factor   = system.file("extdata/N6mA_positive.bed.gz", package = "EnrichShuf"), 
  element  = system.file("extdata/ATAC.bed.gz", package = "EnrichShuf"), 
  genome   = system.file("extdata/hg38_no_chrYM.genomesizes", package = "EnrichShuf"), 
  excl     = system.file("extdata/shuffle_exclude.bed.gz", package = "EnrichShuf"), 
  parallel = 12
)

if (!identical(answer_dat$observe, target_dat$observe)) {
  stop("something wrong in observed data")
}

if (!length(answer_dat$expect)==length(target_dat$expect)) {
  stop("the length of expected data is different")
}

for (i in 1:length(answer_dat$expect)) {
  if (!identical(answer_dat$expect[[i]], target_dat$expect[[i]])) {
    stop("something wrong in expected data")
  }
}

message("complete the test procedure")