#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# Walk through a genome looking at the absolute difference (AD) of mean values
# in methylatin, CpG counts and coverage between adjacent regions of the
# given size

# Input is a directory of CpG level measurements, output is a 
# bedgraph file with an AD value for each window that can be 
# visualized as a track in a genome viewer.
# Size is the size of each region to be compared
# Step is the step size for walking through the genome

# Usage: bp_scan.R \
#   bin_directory/ \
#   genome_chr_length.txt \
#   methylation_cpg_drive/ \
#   size=10000 step=1000 min_cov=4 \
#   out_directory/

# Example: bp_scan.R \
#  APE_METH_bin/ \
#  NomLeu1_0/seq_len.txt \
#  MAPPED_GOOD/Gibbon/Nleu3.0/DNA111101LC_62_HSA_normal_NoIndex_L006_1_val_1_trim.cpg10/ \
#  size=10000 step=1000 min_cov=4 \
#  BP_SCAN/Gibbon/Nleu3.0/

.libPaths("/home/users/lazar/R/x86_64-redhat-linux-gnu-library/3.1")
library(bsseq)
condor <- T

args <- commandArgs(TRUE)
bin_dir <- args[1]
chr_len_file <- args[2]
cpg_dir <- args[3]
size <- as.numeric(sub('size=', '', args[4]))
step <- as.numeric(sub('step=', '', args[5]))
min_cov <- as.numeric(sub('min_cov=', '', args[6]))
out_dir <- args[7]
dir <- paste0(getwd(), "/")

source(paste0(bin_dir, 'make_seqinfo.R'))
source(paste0(bin_dir, 'make_all_bs.R'))
source(paste0(bin_dir, 'mcread.bsmooth.R'))
if(condor) source(paste0(bin_dir, 'condor_add_meth_cpg_cov.R'))

# Check that size is a multiple of step
#######################################
if((size %% step) != 0) {
  print("size must be a multiple of step")
  break
}

# Make a seq info object storing the chromosome lengths
#######################################################
seqinfo <- make_seqinfo(chr_len_file, strsplit(chr_len_file, "/")[[1]][1])

# Read in CpG data and make BSeq object for all CpGs
####################################################
# If the 'all_meth.Rdata' file is already there, just load it.
if(file.exists(paste0(out_dir, 'all_meth.Rdata'))) {
  load(paste0(out_dir, 'all_meth.Rdata'))
} else {
  all.bs <- make_all_bs(cpg_dir, strsplit(chr_len_file, "/")[[1]][1], seqinfo, 0)
  save(all.bs, file=paste0(out_dir, 'all_meth.Rdata'))
}

# Make BSeq object of CpGs with coverage above <min_cov>
########################################################
cov4.bs <- all.bs[getCoverage(all.bs) >=4]

# Make a GRanges object to store the AD values for each location
#################################################################
ad.gr <- tileGenome(seqinfo, tilewidth=step, cut.last.tile.in.chrom=T) 
ad.gr$meth <- 0
ad.gr$cov <- 0
ad.gr$cpgs <- 0
ad.gr$cpgs_w_cov <- 0

# Make a GRanges object storing the left & right regions of size <size>
#######################################################################
left.gr <- flank(ad.gr, size)
right.gr <- resize(ad.gr, size)

# Calculate methylation and coverage for left and right side of 
# each region in ad.gr
###################################################################
if(condor){
  # split up adding info 10,000 ranges at a time
  splits <- rep(1:(length(left.gr)/10000+1), each=10000)[1:length(left.gr)]
  left.list <- split(left.gr, splits)
  for(i in 1:length(left.list)) {
    left.list[[i]] <- condor_add_meth_cpg_cov(left.list[[i]], all.bs,
               dir, paste0(out_dir, 'add_meth_left_', i), bin_dir)
  }
  left.gr <- unsplit(left.list, splits)

  right.list <- split(right.gr, splits)
  for(i in 1:length(right.list)) {
    right.list[[i]] <- condor_add_meth_cpg_cov(right.list[[i]], all.bs,
               dir, paste0(out_dir, 'add_meth_right_', i), bin_dir)
  }
  right.gr <- unsplit(right.list, splits)
} else {
  splits <- rep(1:(length(left.gr)/10000+1), each=10000)[1:length(left.gr)]
  left.list <- split(left.gr, splits)
  for(i in 1:length(left.list)) {
    left.list[[i]]$meth <- getMeth(cov4.bs, regions=left.list[[i]], type='raw', what='perRegion')
    left.list[[i]]$cov <- getCoverage(all.bs, regions=left.list[[i]], what='perRegionAverage')
    left.list[[i]]$cpgs <- sapply(getCoverage(all.bs, regions=left.list[[i]], type="Cov"), length)
    left.list[[i]]$cpgs_w_cov <- sapply(getCoverage(cov4.bs, regions=left.list[[i]], type="Cov"), length)
  }
  left.gr <- unsplit(left.list, splits)

  right.list <- split(right.gr, splits)
  for(i in 1:length(right.list)) {
    right.list[[i]]$meth <- getMeth(cov4.bs, regions=right.list[[i]], type='raw', what='perRegion')
    right.list[[i]]$cov <- getCoverage(all.bs, regions=right.list[[i]], what='perRegionAverage')
    right.list[[i]]$cpgs <- sapply(getCoverage(all.bs, regions=right.list[[i]], type="Cov"), length)
    right.list[[i]]$cpgs_w_cov <- sapply(getCoverage(cov4.bs, regions=right.list[[i]], type="Cov"), length)
  }
  right.gr <- unsplit(right.list, splits)
}


# Calculate ad values for each region in ad.gr
################################################
ad.gr$meth <- abs(left.gr$meth - right.gr$meth)
ad.gr$cov <- abs(left.gr$cov - right.gr$cov)
ad.gr$cpgs <- abs(left.gr$cpgs - right.gr$cpgs)
ad.gr$cpgs_w_cov <- abs(left.gr$cpgs_w_cov - right.gr$cpgs_w_cov)

# Write AD values out to bedgraph files
########################################
ad.df <- data.frame(chr=sub('chr', '', seqnames(ad.gr)),
  start=start(ad.gr), end=end(ad.gr), meth=ad.gr$meth, 
  cov=ad.gr$cov, cpgs=ad.gr$cpgs, cpgs_w_cov=ad.gr$cpgs_w_cov)

writeLines(paste0('track type=bedGraph name=', 'meth_', size,
  ' color=0,0,0 viewLimits=0:1 visibility=full'), 
  con=paste0(out_dir, 'ad_meth.bedgraph'))
write.table(ad.df[,c('chr', 'start', 'end', 'meth')],
  file=paste0(out_dir, 'ad_meth.bedgraph'), col.names=F,
  row.names=F, quote=F, append=T)

writeLines(paste0('track type=bedGraph name=', 'cov_', size,
  ' color=0,0,0 visibility=full'),
  con=paste0(out_dir, 'ad_cov.bedgraph'))
write.table(ad.df[,c('chr', 'start', 'end', 'cov')],
  file=paste0(out_dir, 'ad_cov.bedgraph'), col.names=F,
  row.names=F, quote=F, append=T)

writeLines(paste0('track type=bedGraph name=', 'cpgs_', size,
  ' color=0,0,0 visibility=full'),
  con=paste0(out_dir, 'ad_cpgs.bedgraph'))
write.table(ad.df[,c('chr', 'start', 'end', 'cpgs')],
  file=paste0(out_dir, 'ad_cpgs.bedgraph'), col.names=F,
  row.names=F, quote=F, append=T)

writeLines(paste0('track type=bedGraph name=', 'cpgs_w_cov_', size,
  ' color=0,0,0 visibility=full'),
  con=paste0(out_dir, 'ad_cpgs_w_cov.bedgraph'))
write.table(ad.df[,c('chr', 'start', 'end', 'cpgs_w_cov')],
  file=paste0(out_dir, 'ad_cpgs_w_cov.bedgraph'), col.names=F,
  row.names=F, quote=F, append=T)
