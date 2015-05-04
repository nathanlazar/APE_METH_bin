#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# R script to be run by HTCondor in parallel 
# that wraps par_rand_sides.R

.libPaths('/home/users/lazar/R/x86_64-redhat-linux-gnu-library/3.1/')

library(foreach)
library(doMC)
library(bsseq)
source('/mnt/lustre1/users/lazar/APE_METH/POST_CRASH/APE_METH_bin/R_meth_functions.R')

# Usage:
# Rscript ./wrap_par_rand.R
#   <file_of_R_data.dat>
#   <size of ends of chroms to be excluded>
#   <adjacent=1 or 0>
#   <number of cores>
#   <out file>

# Example:
# Rscript ./wrap_par_rand.R
#   par_permute_sides.dat
#   1000
#   1
#   16
#   par_sides.Rdata

args <- commandArgs(TRUE)
print(args)

load(args[1])             #Load all.bs, bp.lr.mad, breaks, lengths
end.exclude <- as.numeric(args[2])
adjacent <- as.logical(as.numeric(args[3]))
reps <- as.numeric(args[4])
registerDoMC(reps)

rand <- par_rand_sides(all.bs, bp.lr.mad, breaks, lengths, end.exclude, adjacent, reps)
save(rand, file=args[5])

print(paste('completed run data written to', args[5]))