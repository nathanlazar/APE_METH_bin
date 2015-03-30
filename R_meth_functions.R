# Functions used in analysis of bisulfite reads using BSmooth
# nathan dot lazar at gmail dot com

library(bsseq)
library(parallel)
library(dplyr)
library(reshape2)

bindir <- '/mnt/lustre1/users/lazar/APE_METH/POST_CRASH/APE_METH_bin/'
#bindir <- '~/gibbon_meth/'}

source(paste0(bindir, 'make_seqinfo.R'))
source(paste0(bindir, 'mcread.bsmooth.R'))
source(paste0(bindir, 'make_tracks.R'))
source(paste0(bindir, 'mcgetMeth.R'))
source(paste0(bindir, 'mcgetCoverage.R'))
source(paste0(bindir, 'make_all_bs.R'))
source(paste0(bindir, 'read_bp.R'))
source(paste0(bindir, 'read_cpg_island.R'))
source(paste0(bindir, 'make_cpg_shore.R'))
source(paste0(bindir, 'make_lr_bp.R'))
source(paste0(bindir, 'make_rep_gr.R'))
source(paste0(bindir, 'gtf2GRanges.R'))
source(paste0(bindir, 'gff2GRanges.R'))
source(paste0(bindir, 'add_meth_cpg_cov.R'))
source(paste0(bindir, 'get_region.R'))
source(paste0(bindir, 'permute.R'))
source(paste0(bindir, 'par_permute.R'))
source(paste0(bindir, 'par_rand.R'))
source(paste0(bindir, 'make_per_submit.R'))
source(paste0(bindir, 'condor_add_meth_cpg_cov.R'))
source(paste0(bindir, 'combine_per.R'))
source(paste0(bindir, 'par_rand_sides.R'))
source(paste0(bindir, 'par_permute_sides.R'))