#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# Compare methylation, coverage and CpG counts of genes containing 
# LAVA elements to genes without LAVA elements in gibbon NomLeu1.0

# Usage: LAVA_gibbons.R \
#   /mnt/lustre1/users/lazar/APE_METH/POST_CRASH/ \
#   APE_METH_bin/ \
#   LAVA/Gibbon \
#   NomLeu1_0/seq_len.txt \
#   MAPPED_GOOD/Gibbon/DNA111101LC_62_HSA_normal_NoIndex_L006_1_val_1_trim.fq.gz/cpg10/ \
#   GIBBON_FEATURES/NomLeu1.0.70.genes.gtf
#   LAVA/LAVAs.txt
#   LAVA/LAVA_genes.txt

library(bsseq)

args <- commandArgs(TRUE)
dir <- args[1]
bindir <- args[2]
outdir <- args[3]
len_file <- args[4]
cpg_drive <- args[5]
gene_file <- args[6]
lava_file <- args[7]
lava_genes <- args[8]

source(paste0(bindir, 'R_meth_functions.R'))

# Functions used later
######################
get_p <- function(obs, center, perm) {
  diff <- abs(obs - center)
  p <- sum(perm > (center + diff) | perm < (center-diff))/length(perm)
  p
}


perm <- function(gr, num, all.bs, reps=1000, min_cpgs_w_cov=4) {
  # Function to get 1000 sets of features w/o LAVA
  # insertions and measure their mean and median
  # methylation, coverage and cpg counts

  # Get random gene sets
  rand.gr <- gr[sample(length(gr), num*reps, replace=T)]

  # Add methylation, coverage, etc. info to genes
  rand.gr <- add_meth_cpg_cov(rand.gr, all.bs, parallel=F, min_cov=4)

  perm.res <- data.frame(perm=1:reps, meth.mean=0, cov.mean=0, 
                         cpgs.mean=0, cpgs_w_cov.mean=0,
                         meth.median=0, cov.median=0,
                         cpgs.median=0, cpgs_w_cov.median=0)
  for(i in 1:reps) {
    perm.res$meth.mean[i] <- mean(rand.gr[(i*(num-1)+1):(i*num)]$meth, na.rm=T)
    perm.res$cov.mean[i] <- mean(rand.gr[(i*(num-1)+1):(i*num)]$cov, na.rm=T)
    perm.res$cpgs.mean[i] <- mean(rand.gr[(i*(num-1)+1):(i*num)]$cpgs, na.rm=T)
    perm.res$cpgs_w_cov.mean[i] <- mean(rand.gr[(i*(num-1)+1):(i*num)]$cpgs_w_cov, na.rm=T)
    perm.res$meth.median[i] <- median(rand.gr[(i*(num-1)+1):(i*num)]$meth, na.rm=T)
    perm.res$cov.median[i] <- median(rand.gr[(i*(num-1)+1):(i*num)]$cov, na.rm=T)
    perm.res$cpgs.median[i] <- median(rand.gr[(i*(num-1)+1):(i*num)]$cpgs, na.rm=T)
    perm.res$cpgs_w_cov.median[i] <- median(rand.gr[(i*(num-1)+1):(i*num)]$cpgs_w_cov, na.rm=T)
  }
  return(perm.res)
}

# Read in lengths of chromsomes and make seqinfo object
#######################################################
seqinfo <- make_seqinfo(len_file, strsplit(len_file, "/")[[1]][1])

# Read in CpG data and make BSeq object for all CpGs
############################################################
all.bs <- make_all_bs(cpg_drive, strsplit(len_file, "/")[[1]][1], seqinfo, 0)

# Subset CpGs with coverage of at least 4
#########################################
cov4.bs <- all.bs[getCoverage(all.bs) > 3]

# Read in genes
###############
gene.gr.list <- gtf2GRanges(gene_file, seqinfo, prom_size=1000)
gene.gr.list <- add_meth_cpg_cov(gene.gr.list, all.bs, parallel=F, min_cov=4)

# Read in LAVA elements and make GRanges object
###############################################
lava <- read.table(lava_file, stringsAsFactors=F, sep='\t', header=F)
names(lava) <- c('chr', 'start', 'end')
# Add 'chr'
lava$chr <- paste0('chr', lava$chr)
# Change to 1 based coordinates
lava$start <- lava$start + 1
# Make GRanges object
lava.gr <- makeGRangesFromDataFrame(lava,
               keep.extra.columns=T)
lava.gr <- lava.gr[seqnames(lava.gr) %in% seqlevels(seqinfo)]
seqlevels(lava.gr) <- seqlevels(seqinfo)
seqlengths(lava.gr) <- seqlengths(seqinfo)
lava.gr <- sort(lava.gr)

# Get methylation, coverage & CpG count for all LAVAs
#####################################################
lava.gr <- add_meth_cpg_cov(lava.gr, all.bs, parallel=F, min_cov=4)

# Overlap LAVAs w/ genes
########################
lava_in_gene <- list()
lava_in_gene$gene <- GenomicRanges::subsetByOverlaps(lava.gr, gene.gr.list$gene, ignore.strand=T)
lava_in_gene$exon <- GenomicRanges::subsetByOverlaps(lava.gr, gene.gr.list$exon, ignore.strand=T)
lava_in_gene$intron <- GenomicRanges::subsetByOverlaps(lava.gr, gene.gr.list$intron, ignore.strand=T)
lava_in_gene$promoter <- GenomicRanges::subsetByOverlaps(lava.gr, gene.gr.list$promoter, ignore.strand=T)

# Overlap genes w/ LAVAs
########################
gene_w_lava <- list()
gene_w_lava$gene <- GenomicRanges::subsetByOverlaps(gene.gr.list$gene, lava.gr, ignore.strand=T)
gene_w_lava$exon <- GenomicRanges::subsetByOverlaps(gene.gr.list$exon, lava.gr, ignore.strand=T)
gene_w_lava$intron <- GenomicRanges::subsetByOverlaps(gene.gr.list$intron, lava.gr, ignore.strand=T)
gene_w_lava$promoter <- GenomicRanges::subsetByOverlaps(gene.gr.list$promoter, lava.gr, ignore.strand=T)

# Get genes not in LAVAs
########################
gene_wo_lava <- list()
gene_wo_lava$gene <- gene.gr.list$gene[!(gene.gr.list$gene %in% gene_w_lava$gene)]
gene_wo_lava$exon <- gene.gr.list$exon[!(gene.gr.list$exon %in% gene_w_lava$exon)]
gene_wo_lava$intron <- gene.gr.list$intron[!(gene.gr.list$intron %in% gene_w_lava$intron)]
gene_wo_lava$promoter <- gene.gr.list$promoter[!(gene.gr.list$promoter %in% gene_w_lava$promoter)]

# Get mean and median methylation, coverage and CpG counts for genes 
# with and without LAVA insertions
####################################################################
res <- list()

res$w_lava.mean <- data.frame(type=names(gene.gr.list))
res$w_lava.mean$counts <- sapply(gene_w_lava, length)
res$w_lava.mean$meth <- sapply(gene_w_lava, function(x) mean(x$meth, na.rm=T))
res$w_lava.mean$cov <- sapply(gene_w_lava, function(x) mean(x$cov, na.rm=T))
res$w_lava.mean$cpgs <- sapply(gene_w_lava, function(x) mean(x$cpgs, na.rm=T))
res$w_lava.mean$cpgs_w_cov <- sapply(gene_w_lava, function(x) mean(x$cpgs_w_cov, na.rm=T))

res$wo_lava.mean <- data.frame(type=names(gene.gr.list))
res$wo_lava.mean$counts <- sapply(gene_wo_lava, length)
res$wo_lava.mean$meth <- sapply(gene_wo_lava, function(x) mean(x$meth, na.rm=T))
res$wo_lava.mean$cov <- sapply(gene_wo_lava, function(x) mean(x$cov, na.rm=T))
res$wo_lava.mean$cpgs <- sapply(gene_wo_lava, function(x) mean(x$cpgs, na.rm=T))
res$wo_lava.mean$cpgs_w_cov <- sapply(gene_wo_lava, function(x) mean(x$cpgs_w_cov, na.rm=T))

res$w_lava.median <- data.frame(type=names(gene.gr.list))
res$w_lava.median$counts <- sapply(gene_w_lava, length)
res$w_lava.median$meth <- sapply(gene_w_lava, function(x) median(x$meth, na.rm=T))
res$w_lava.median$cov <- sapply(gene_w_lava, function(x) median(x$cov, na.rm=T))
res$w_lava.median$cpgs <- sapply(gene_w_lava, function(x) median(x$cpgs, na.rm=T))
res$w_lava.median$cpgs_w_cov <- sapply(gene_w_lava, function(x) median(x$cpgs_w_cov, na.rm=T))

res$wo_lava.median <- data.frame(type=names(gene.gr.list))
res$wo_lava.median$counts <- sapply(gene_wo_lava, length)
res$wo_lava.median$meth <- sapply(gene_wo_lava, function(x) median(x$meth, na.rm=T))
res$wo_lava.median$cov <- sapply(gene_wo_lava, function(x) median(x$cov, na.rm=T))
res$wo_lava.median$cpgs <- sapply(gene_wo_lava, function(x) median(x$cpgs, na.rm=T))
res$wo_lava.median$cpgs_w_cov <- sapply(gene_wo_lava, function(x) median(x$cpgs_w_cov, na.rm=T))

# Get mean and median methylation, coverage & cpgs of genic sequence not counting LAVA
# instertions 
######################################################################################
types <- c(names(gene.gr.list), 'gene_w_intron_insert', 'exon_w_intron_insert', 
           'intron_w_intron_insert', 'prom_w_intron_insert')

res$minus_lava.mean <- data.frame(type=types)
res$minus_lava.median <- data.frame(type=types)

gene_id_w_intron_insert <- setdiff(setdiff(gene_w_lava$gene$gene_id, gene_w_lava$exon$gene_id),
                                   gene_w_lava$promoter$gene_id)
w_intron_insert.gr.list <- lapply(gene.gr.list, function(x) x[x$gene_id %in% gene_id_w_intron_insert])
w_intron_insert.gr.list$intron <- gene_w_lava$intron

minus_lava <- c(lapply(gene_w_lava, GenomicRanges::setdiff, lava.gr, ignore.strand=T),
                lapply(w_intron_insert.gr.list, GenomicRanges::setdiff, lava.gr, ignore.strand=T))
names(minus_lava) <- types

minus_lava <- lapply(minus_lava, add_meth_cpg_cov, all.bs, parallel=F, min_cov=4)

res$minus_lava.mean$counts <- sapply(minus_lava, length)
res$minus_lava.mean$meth <- sapply(minus_lava, function(x) mean(x$meth, na.rm=T))
res$minus_lava.mean$cov <- sapply(minus_lava, function(x) mean(x$cov, na.rm=T))
res$minus_lava.mean$cpgs <- sapply(minus_lava, function(x) mean(x$cpgs, na.rm=T))
res$minus_lava.mean$cpgs_w_cov <- sapply(minus_lava, function(x) mean(x$cpgs_w_cov, na.rm=T))

res$minus_lava.median$counts <- sapply(minus_lava, length)
res$minus_lava.median$meth <- sapply(minus_lava, function(x) median(x$meth, na.rm=T))
res$minus_lava.median$cov <- sapply(minus_lava, function(x) median(x$cov, na.rm=T))
res$minus_lava.median$cpgs <- sapply(minus_lava, function(x) median(x$cpgs, na.rm=T))
res$minus_lava.median$cpgs_w_cov <- sapply(minus_lava, function(x) median(x$cpgs_w_cov, na.rm=T))

# Run permutations
##################
perm.list <- list()
perm.list$gene <- perm(gene_wo_lava$gene, length(gene_w_lava$gene), all.bs, reps=1000, min_cov=4)
perm.list$exon <- perm(gene_wo_lava$exon, length(gene_w_lava$exon), all.bs, reps=1000, min_cov=4)
perm.list$intron <- perm(gene_wo_lava$intron, length(gene_w_lava$intron), all.bs, reps=1000, min_cov=4)
perm.list$promoter <- perm(gene_wo_lava$promoter, length(gene_w_lava$promoter), all.bs, reps=1000, min_cov=4)

# Save permutation stats and get p-values
#########################################
res$perm_stats <- data.frame(type=names(gene.gr.list))
res$p_value_lava <- data.frame(type=names(gene.gr.list))
res$p_value_minus_lava <- data.frame(type=types)

res$perm_stats$meth.mean <- sapply(perm.list, function(x) mean(x$meth.mean))
res$perm_stats$cov.mean <- sapply(perm.list, function(x) mean(x$cov.mean))
res$perm_stats$cpgs.mean <- sapply(perm.list, function(x) mean(x$cpgs.mean))
res$perm_stats$cpgs_w_cov.mean <- sapply(perm.list, function(x) mean(x$cpgs_w_cov.mean))
res$perm_stats$meth.sd <- sapply(perm.list, function(x) sd(x$meth.mean))
res$perm_stats$cov.sd <- sapply(perm.list, function(x) sd(x$cov.mean))
res$perm_stats$cpgs.sd <- sapply(perm.list, function(x) sd(x$cpgs.mean))
res$perm_stats$cpgs_w_cov.sd <- sapply(perm.list, function(x) sd(x$cpgs_w_cov.mean))
res$perm_stats$meth.median <- sapply(perm.list, function(x) median(x$meth.median))
res$perm_stats$cov.median <- sapply(perm.list, function(x) median(x$cov.median))
res$perm_stats$cpgs.median <- sapply(perm.list, function(x) median(x$cpgs.median))
res$perm_stats$cpgs_w_cov.median <- sapply(perm.list, function(x) median(x$cpgs_w_cov.median))

res$p_value_lava$meth.mean.p_value <- sapply(1:4, 
  function(i) get_p(res$w_lava.mean$meth[i], res$perm_stats$meth.mean[i], perm.list[[i]]$meth.mean))
res$p_value_lava$cov.mean.p_value <- sapply(1:4,
  function(i) get_p(res$w_lava.mean$cov[i], res$perm_stats$cov.mean[i], perm.list[[i]]$cov.mean))
res$p_value_lava$cpg.mean.p_value <- sapply(1:4,
  function(i) get_p(res$w_lava.mean$cpgs[i], res$perm_stats$cpgs.mean[i], perm.list[[i]]$cpgs.mean))
res$p_value_lava$cpg_w_cov.mean.p_value <- sapply(1:4,
  function(i) get_p(res$w_lava.mean$cpgs_w_cov[i], res$perm_stats$cpgs_w_cov.mean[i], perm.list[[i]]$cpgs_w_cov.mean))

res$p_value_lava$meth.median.p_value <- sapply(1:4, 
  function(i) get_p(res$w_lava.median$meth[i], res$perm_stats$meth.median[i], perm.list[[i]]$meth.median))
res$p_value_lava$cov.median.p_value <- sapply(1:4,
  function(i) get_p(res$w_lava.median$cov[i], res$perm_stats$cov.median[i], perm.list[[i]]$cov.median))
res$p_value_lava$cpg.median.p_value <- sapply(1:4,
  function(i) get_p(res$w_lava.median$cpgs[i], res$perm_stats$cpgs.median[i], perm.list[[i]]$cpgs.median))
res$p_value_lava$cpg_w_cov.median.p_value <- sapply(1:4,
  function(i) get_p(res$w_lava.median$cpgs_w_cov[i], res$perm_stats$cpgs_w_cov.median[i], perm.list[[i]]$cpgs_w_cov.median))

res$p_value_minus_lava$meth.mean.p_value <- c(sapply(1:4, 
  function(i) get_p(res$minus_lava.mean$meth[i], res$perm_stats$meth.mean[i], perm.list[[i]]$meth.mean)),
  sapply(1:4, function(i) get_p(res$minus_lava.mean$meth[i+4], res$perm_stats$meth.mean[i],  perm.list[[i]]$meth.mean)))
res$p_value_minus_lava$cov.mean.p_value <- c(sapply(1:4, 
  function(i) get_p(res$minus_lava.mean$cov[i], res$perm_stats$cov.mean[i], perm.list[[i]]$cov.mean)),
  sapply(1:4, function(i) get_p(res$minus_lava.mean$cov[i+4], res$perm_stats$cov.mean[i],  perm.list[[i]]$cov.mean)))
res$p_value_minus_lava$cpgs.mean.p_value <- c(sapply(1:4,
  function(i) get_p(res$minus_lava.mean$cpgs[i], res$perm_stats$cpgs.mean[i], perm.list[[i]]$cpgs.mean)),
  sapply(1:4, function(i) get_p(res$minus_lava.mean$cpgs[i+4], res$perm_stats$cpgs.mean[i],  perm.list[[i]]$cpgs.mean)))
res$p_value_minus_lava$cpgs_w_cov.mean.p_value <- c(sapply(1:4,
  function(i) get_p(res$minus_lava.mean$cpgs_w_cov[i], res$perm_stats$cpgs_w_cov.mean[i], perm.list[[i]]$cpgs_w_cov.mean)),
  sapply(1:4, function(i) get_p(res$minus_lava.mean$cpgs_w_cov[i+4], res$perm_stats$cpgs_w_cov.mean[i],  perm.list[[i]]$cpgs_w_cov.mean)))

res$p_value_minus_lava$meth.median.p_value <- c(sapply(1:4, 
  function(i) get_p(res$minus_lava.median$meth[i], res$perm_stats$meth.median[i], perm.list[[i]]$meth.median)),
  sapply(1:4, function(i) get_p(res$minus_lava.median$meth[i+4], res$perm_stats$meth.median[i],  perm.list[[i]]$meth.median)))
res$p_value_minus_lava$cov.median.p_value <- c(sapply(1:4, 
  function(i) get_p(res$minus_lava.median$cov[i], res$perm_stats$cov.median[i], perm.list[[i]]$cov.median)),
  sapply(1:4, function(i) get_p(res$minus_lava.median$cov[i+4], res$perm_stats$cov.median[i],  perm.list[[i]]$cov.median)))
res$p_value_minus_lava$cpgs.median.p_value <- c(sapply(1:4,
  function(i) get_p(res$minus_lava.median$cpgs[i], res$perm_stats$cpgs.median[i], perm.list[[i]]$cpgs.median)),
  sapply(1:4, function(i) get_p(res$minus_lava.median$cpgs[i+4], res$perm_stats$cpgs.median[i],  perm.list[[i]]$cpgs.median)))
res$p_value_minus_lava$cpgs_w_cov.median.p_value <- c(sapply(1:4,
  function(i) get_p(res$minus_lava.median$cpgs_w_cov[i], res$perm_stats$cpgs_w_cov.median[i], perm.list[[i]]$cpgs_w_cov.median)),
  sapply(1:4, function(i) get_p(res$minus_lava.median$cpgs_w_cov[i+4], res$perm_stats$cpgs_w_cov.median[i],  perm.list[[i]]$cpgs_w_cov.median)))

save(res, file='LAVA_perm_results.Rdata')
