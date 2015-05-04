#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# Make tracks showing the absolute difference in mean methylation,
# coverage and CpG count between first and second half of 20kb regions.
# The 20kb window slides by <slide> base pairs

# Usage: sweep_cpg_meth.R \
#   /mnt/lustre1/users/lazar/APE_METH/POST_CRASH/ \
#   APE_METH_bin/ \
#   Slide_analysis/Gibbon \
#   NomLeu1_0/seq_len.txt \
#   MAPPED_GOOD/Gibbon/DNA111101LC_62_HSA_normal_NoIndex_L006_1_val_1_trim.fq.gz/cpg10/ \
#   500

library(bsseq)

args <- commandArgs(TRUE)
dir <- args[1]
bindir <- args[2]
outdir <- args[3]
len_file <- args[4]
cpg_drive <- args[5]
slide <- args[6]

#bp_file <- args[6]
#gene_file <- args[7]
#rep_file <- args[8]
#cpg_isl_file <- args[9]

source(paste0(bindir, 'make_seqinfo.R'))
source(paste0(bindir, 'mcread.bsmooth.R'))
source(paste0(bindir, 'make_all_bs.R'))

# Read in lengths of chromsomes and make seqinfo object
#######################################################
seqinfo <- make_seqinfo(len_file, strsplit(len_file, "/")[[1]][1])

# Read in CpG data and make BSeq object for all CpGs
####################################################
all.bs <- make_all_bs(cpg_drive, strsplit(len_file, "/")[[1]][1], seqinfo, 0)

# Subset CpGs with coverage of at least 4
#########################################
cov4.bs <- all.bs[getCoverage(all.bs) > 3]


