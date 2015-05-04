#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# R script to be run by HTCondor in parallel 
# that wraps add_meth_cpg_cov.R

# Usage:
# Rscript ./wrap_add_meth_cpg_cov.R
#   <feat_gr_and_all_bs.dat>  
#   <bin_dir>
#   <out_dir>
#   <number of cores>

.libPaths('/home/users/lazar/R/x86_64-redhat-linux-gnu-library/3.1/')

library(foreach)
library(doMC)
library(bsseq)

args <- commandArgs(TRUE)
#Load  feat.gr, all.bs
if(length(args) > 0) {
  load(args[1])
  source(paste0(args[2], '/R_meth_functions.R'))
  cores <- as.numeric(args[3])
  registerDoMC(cores)
  outdir <- args[4]

  feat.gr <- add_meth_cpg_cov(feat.gr, all.bs)
  # Save data to file
  save(feat.gr, file=paste0(outdir, 'feat.dat'))

  print(1)
}
