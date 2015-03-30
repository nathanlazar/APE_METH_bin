#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# R script to be run by HTCondor in parallel 
# that wraps par_rand.R

.libPaths('/home/users/lazar/R/x86_64-redhat-linux-gnu-library/3.1/')

library(foreach)
library(doMC)
library(bsseq)
source('/mnt/lustre1/users/lazar/GIBBONS/gibbon_meth/R_meth_functions.R')

# Usage:
# Rscript ./wrap_par_rand.R
#   <file_of_R_data.dat>
#   <type = all, gene, etc.>
#   <size of ends of chroms to be excluded>
#   <number of cores>

# Example:
# Rscript ./wrap_par_rand.R
#   par_permute.dat
#   all
#   1000
#   16

args <- commandArgs(TRUE)
#Load  feat.gr, bp.lr.gr, all.bs, breaks, sizes, lengths
load(args[1])

type <- args[2]
end.exclude <- as.numeric(args[3])
reps <- as.numeric(args[4])
registerDoMC(reps)

sizes <- width(bp.lr.gr)

rand <- par_rand(all.bs, feat.gr, breaks, sizes, lengths, end.exclude, type, reps)
print(rand)


