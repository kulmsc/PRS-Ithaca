library(bigsnpr)

args = commandArgs(trailingOnly=TRUE)

snp_readBed(paste0("new.",args[1],".filter.bed"))
