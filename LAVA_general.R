#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# Compare methylation, coverage and CpG counts of genes in other species
# that contain LAVA elements in gibbons to genes without LAVA elements in 
# gibbon NomLeu1.0

# Usage: LAVA_general.R \
#   /mnt/lustre1/users/lazar/APE_METH/POST_CRASH/ \
#   APE_METH_bin/ \
#   LAVA/Human \
#   human_hg19_noIUPAC/lengths.txt \
#   MAPPED_GOOD/Human/combined_cpg/ \
#   RefSeqGenes/Human_hg19_OtherRefSeq.txt \
#   LAVA/Gibbon/gene_w_lava.Rdata

.libPaths("/home/users/lazar/R/x86_64-redhat-linux-gnu-library/3.1")
library(bsseq)

args <- commandArgs(TRUE)
dir <- args[1]
bindir <- args[2]
outdir <- args[3]
len_file <- args[4]
cpg_drive <- args[5]
gene_file <- args[6]
lava_file <- args[7]

source(paste0(bindir, 'R_meth_functions.R'))

# Functions used later
######################
get_p <- function(obs, center, perm) {
  diff <- abs(obs - center)
  p <- sum(perm > (center + diff) | perm < (center-diff), na.rm=T)/length(perm)
  p
}

perm <- function(gr, num, all.bs, reps=1000, min_cpgs_w_cov=4, parallel=T, cores=NA) {
  # Function to get 1000 sets of features w/o LAVA
  # insertions and measure their mean and median
  # methylation, coverage and cpg counts

  # Get random gene sets
  rand.gr <- gr[sample(length(gr), num*reps, replace=T)]

  # Add methylation, coverage, etc. info to genes
  rand.gr <- add_meth_cpg_cov(rand.gr, all.bs, parallel=T, min_cov=4)

  perm.res <- data.frame(perm=1:reps, meth.mean=0, cov.mean=0, 
                         cpgs.mean=0, cpgs_w_cov.mean=0,
                         meth.median=0, cov.median=0,
                         cpgs.median=0, cpgs_w_cov.median=0)

  if(parallel) {
    # Import packages
    library(doMC)

    # If the number of cores isn't specified use all available
    if(is.na(cores)) cores <- getOption("mc.cores", 2L)

    # Setup parallel backend
    registerDoMC(cores)

    # Run parallel operations
    perm.res$meth.mean <- foreach(i=1:reps) %dopar% mean(rand.gr$meth[(i*(num-1)+1):(i*num)], na.rm=T)
    perm.res$cov.mean <- foreach(i=1:reps) %dopar% mean(rand.gr$cov[(i*(num-1)+1):(i*num)], na.rm=T)
    perm.res$cpgs.mean <- foreach(i=1:reps) %dopar% mean(rand.gr$cpgs[(i*(num-1)+1):(i*num)], na.rm=T)
    perm.res$cpgs_w_cov.mean <- foreach(i=1:reps) %dopar% mean(rand.gr$cpgs_w_cov[(i*(num-1)+1):(i*num)], na.rm=T)
    perm.res$meth.median <- foreach(i=1:reps) %dopar% median(rand.gr$meth[(i*(num-1)+1):(i*num)], na.rm=T)
    perm.res$cov.median <- foreach(i=1:reps) %dopar% median(rand.gr$cov[(i*(num-1)+1):(i*num)], na.rm=T)
    perm.res$cpgs.median <- foreach(i=1:reps) %dopar% median(rand.gr$cpgs[(i*(num-1)+1):(i*num)], na.rm=T)
    perm.res$cpgs_w_cov.median <- foreach(i=1:reps) %dopar% median(rand.gr$cpgs_w_cov[(i*(num-1)+1):(i*num)], na.rm=T)

  } else {
    for(i in 1:reps) {
      perm.res$meth.mean[i] <- mean(rand.gr$meth[(i*(num-1)+1):(i*num)], na.rm=T)
      perm.res$cov.mean[i] <- mean(rand.gr$cov[(i*(num-1)+1):(i*num)], na.rm=T)
      perm.res$cpgs.mean[i] <- mean(rand.gr$cpgs[(i*(num-1)+1):(i*num)], na.rm=T)
      perm.res$cpgs_w_cov.mean[i] <- mean(rand.gr$cpgs_w_cov[(i*(num-1)+1):(i*num)], na.rm=T)
      perm.res$meth.median[i] <- median(rand.gr$meth[(i*(num-1)+1):(i*num)], na.rm=T)
      perm.res$cov.median[i] <- median(rand.gr$cov[(i*(num-1)+1):(i*num)], na.rm=T)
      perm.res$cpgs.median[i] <- median(rand.gr$cpgs[(i*(num-1)+1):(i*num)], na.rm=T)
      perm.res$cpgs_w_cov.median[i] <- median(rand.gr$cpgs_w_cov[(i*(num-1)+1):(i*num)], na.rm=T)
    }
  }
  return(perm.res)
}

# Read in lengths of chromsomes and make seqinfo object
#######################################################
seqinfo <- make_seqinfo(len_file, strsplit(len_file, "/")[[1]][1])

# Read in CpG data and make BSeq object for all CpGs
############################################################
# If the 'all_meth.Rdata' file is already there, just load it.
if(file.exists(paste0(outdir, 'all_meth.Rdata'))) {
  load(paste0(outdir, 'all_meth.Rdata'))
} else {
  all.bs <- make_all_bs(cpg_drive, strsplit(len_file, "/")[[1]][1], seqinfo, 0)
  save(all.bs, file=paste0(outdir, 'all_meth.Rdata'))
}

# Subset CpGs with coverage of at least 4
#########################################
cov4.bs <- all.bs[getCoverage(all.bs) > 3]

# Read in genes
###############
if(grepl('.gtf', gene_file)) {
  gene.gr.list <- gtf2GRanges(gene_file, seqinfo, prom_size=1000)
} else {
  gene.gr.list <- UCSC2GRanges(gene_file, seqinfo, prom_size=1000, filter=F)
}
# Only keep the first isoform 'name' is isoform_id 'symbol' is gene_id
gene.gr.list$gene <- gene.gr.list$gene[!duplicated(gene.gr.list$gene$symbol)]
gene.gr.list$promoter <- gene.gr.list$promoter[gene.gr.list$promoter$name %in% gene.gr.list$gene$name]
gene.gr.list$utr5 <- gene.gr.list$utr5[gene.gr.list$utr5$name %in% gene.gr.list$gene$name]
gene.gr.list$exon <- gene.gr.list$exon[gene.gr.list$exon$name %in% gene.gr.list$gene$name]
gene.gr.list$intron <- gene.gr.list$intron[gene.gr.list$intron$name %in% gene.gr.list$gene$name]
gene.gr.list$utr3 <- gene.gr.list$utr3[gene.gr.list$utr3$name %in% gene.gr.list$gene$name]

gene.gr.list <- lapply(gene.gr.list, add_meth_cpg_cov, all.bs, parallel=T, min_cov=4)
save(gene.gr.list, file=paste0(outdir, 'gene_gr_list.Rdata'))

# Get list of genes that have  LAVA insertions in gibbon
########################################################
load(lava_file)
gibbon_lavas <- gene_w_lava

# Get gene mapping between gibbon & human ############## Change this for other species ############
mapping <- read.table(paste0(dir, 'EnsemblGenes/hg19_nomLeu1_0_Ensembl_mapping.txt'),
           sep='\t', stringsAsFactors=F, header=T)

# Get gene elements with LAVA in gibbon
gene_w_lava <- list()
for(i in 1:length(gene.gr.list)) {
  type <- names(gene.gr.list)[i]
  gene_w_lava[[type]] <- gene.gr.list[[i]][gene.gr.list[[i]]$symbol %in% 
    unique(mapping$Ensembl.Gene.ID[mapping$Gibbon.Ensembl.Gene.ID %in% gibbon_lavas[[type]]$symbol])]
}
# The above gives all exons & introns in genes w/ insertions. 
# note, there isn't a mapping between these, so this could introduce errors

# Get a list of gene symbols with LAVAs in gibbon
gene_symbol_w_lava <- unique(mapping$Ensembl.Gene.ID[mapping$Gibbon.Ensembl.Gene.ID %in%
  unique(unlist(lapply(gibbon_lavas, function(x) x$symbol)))])

# Report the number of LAVA gene names found in the gene set
print(paste('There are ', length(unique(unlist(lapply(gibbon_lavas, function(x) x$symbol)))), 
      'genes with LAVA insertions in gibbon'))
print(paste('In this genome', length(gene_symbol_w_lava), 'have been identified'))

# Get GRanges list of genes w/ LAVAs
####################################
gene_w_lava_any_ins <- lapply(gene.gr.list, function(x) x[x$symbol %in% gene_symbol_w_lava])
names(gene_w_lava_any_ins) <- paste0(names(gene_w_lava_any_ins), '_w_gene_insert')
gene_w_lava <- c(gene_w_lava, gene_w_lava_any_ins)

save(gene_w_lava, file=paste0(outdir, 'gene_w_lava.Rdata'))

# Get genes not in LAVAs
########################
gene_wo_lava <- list()
gene_wo_lava$gene <- gene.gr.list$gene[!(gene.gr.list$gene$symbol %in% gene_w_lava$gene$symbol)]
gene_wo_lava$promoter <- gene.gr.list$promoter[!(gene.gr.list$promoter$symbol %in% gene_w_lava$promoter$symbol)]
gene_wo_lava$utr5 <- gene.gr.list$utr5[!(gene.gr.list$utr5$symbol %in% gene_w_lava$utr5$symbol)]
gene_wo_lava$exon <- gene.gr.list$exon[!(gene.gr.list$exon$symbol %in% gene_w_lava$exon$symbol)]
gene_wo_lava$intron <- gene.gr.list$intron[!(gene.gr.list$intron$symbol %in% gene_w_lava$intron$symbol)]
gene_wo_lava$utr3 <- gene.gr.list$utr3[!(gene.gr.list$utr3$symbol %in% gene_w_lava$utr3$symbol)]

save(gene_wo_lava, file=paste0(outdir, 'gene_wo_lava.Rdata'))

# Get mean and median methylation, coverage and CpG counts for genes 
# with and without LAVA insertions
####################################################################
res <- list()

res$w_lava.mean <- data.frame(type=names(gene_w_lava))
res$w_lava.mean$counts <- sapply(gene_w_lava, length)
res$w_lava.mean$meth <- sapply(gene_w_lava, function(x) mean(x$meth, na.rm=T))
res$w_lava.mean$cov <- sapply(gene_w_lava, function(x) mean(x$cov, na.rm=T))
res$w_lava.mean$cpgs <- sapply(gene_w_lava, function(x) mean(x$cpgs, na.rm=T))
res$w_lava.mean$cpgs_w_cov <- sapply(gene_w_lava, function(x) mean(x$cpgs_w_cov, na.rm=T))

res$wo_lava.mean <- data.frame(type=names(gene_wo_lava))
res$wo_lava.mean$counts <- sapply(gene_wo_lava, length)
res$wo_lava.mean$meth <- sapply(gene_wo_lava, function(x) mean(x$meth, na.rm=T))
res$wo_lava.mean$cov <- sapply(gene_wo_lava, function(x) mean(x$cov, na.rm=T))
res$wo_lava.mean$cpgs <- sapply(gene_wo_lava, function(x) mean(x$cpgs, na.rm=T))
res$wo_lava.mean$cpgs_w_cov <- sapply(gene_wo_lava, function(x) mean(x$cpgs_w_cov, na.rm=T))

res$w_lava.median <- data.frame(type=names(gene_w_lava))
res$w_lava.median$counts <- sapply(gene_w_lava, length)
res$w_lava.median$meth <- sapply(gene_w_lava, function(x) median(x$meth, na.rm=T))
res$w_lava.median$cov <- sapply(gene_w_lava, function(x) median(x$cov, na.rm=T))
res$w_lava.median$cpgs <- sapply(gene_w_lava, function(x) median(x$cpgs, na.rm=T))
res$w_lava.median$cpgs_w_cov <- sapply(gene_w_lava, function(x) median(x$cpgs_w_cov, na.rm=T))

res$wo_lava.median <- data.frame(type=names(gene_wo_lava))
res$wo_lava.median$counts <- sapply(gene_wo_lava, length)
res$wo_lava.median$meth <- sapply(gene_wo_lava, function(x) median(x$meth, na.rm=T))
res$wo_lava.median$cov <- sapply(gene_wo_lava, function(x) median(x$cov, na.rm=T))
res$wo_lava.median$cpgs <- sapply(gene_wo_lava, function(x) median(x$cpgs, na.rm=T))
res$wo_lava.median$cpgs_w_cov <- sapply(gene_wo_lava, function(x) median(x$cpgs_w_cov, na.rm=T))

# Run permutations
##################               *** FIX PROBLEM W/ INTRON COVERAGE **
perm.list <- list()
perm.list$gene <- perm(gene_wo_lava$gene, length(gene_w_lava$gene), all.bs, reps=1000, min_cpgs_w_cov=4, parallel=T, cores=24)
perm.list$promoter <- perm(gene_wo_lava$promoter, length(gene_w_lava$promoter), all.bs, reps=1000, min_cpgs_w_cov=4, parallel=T, cores=24)
perm.list$utr5 <- perm(gene_wo_lava$utr5, length(gene_w_lava$utr5), all.bs, reps=1000, min_cpgs_w_cov=4, parallel=T, cores=24)
perm.list$exon <- perm(gene_wo_lava$exon, length(gene_w_lava$exon), all.bs, reps=1000, min_cpgs_w_cov=4, parallel=T, cores=24)
# Split up introns 
tmp1 <- perm(gene_wo_lava$intron, length(gene_w_lava$intron), all.bs, reps=500, min_cpgs_w_cov=4, parallel=T, cores=24)
tmp2 <- perm(gene_wo_lava$intron, length(gene_w_lava$intron), all.bs, reps=500, min_cpgs_w_cov=4, parallel=T, cores=24)
perm.list$intron <- rbind(tmp1, tmp2)
perm.list$utr3 <- perm(gene_wo_lava$utr3, length(gene_w_lava$utr3), all.bs, reps=1000, min_cpgs_w_cov=4, parallel=T, cores=24)

perm.list$gene_w_gene_insert <- perm(gene_wo_lava$gene, length(gene_w_lava$gene_w_gene_insert), all.bs, reps=1000, min_cpgs_w_cov=4, parallel=T, cores=24)
perm.list$promoter_w_gene_insert <- perm(gene_wo_lava$promoter, length(gene_w_lava$promoter_w_gene_insert), all.bs, reps=1000, min_cpgs_w_cov=4, parallel=T, cores=24)
perm.list$utr5_w_gene_insert <- perm(gene_wo_lava$utr5, length(gene_w_lava$utr5_w_gene_insert), all.bs, reps=1000, min_cpgs_w_cov=4, parallel=T, cores=24)
perm.list$exon_w_gene_insert <- perm(gene_wo_lava$exon, length(gene_w_lava$exon_w_gene_insert), all.bs, reps=1000, min_cpgs_w_cov=4, parallel=T, cores=24)
# Split up introns
tmp1 <- perm(gene_wo_lava$intron, length(gene_w_lava$intron_w_gene_insert), all.bs, reps=500, min_cpgs_w_cov=4, parallel=T, cores=24)
tmp2 <- perm(gene_wo_lava$intron, length(gene_w_lava$intron_w_gene_insert), all.bs, reps=500, min_cpgs_w_cov=4, parallel=T, cores=24)
perm.list$intron_w_gene_insert <- rbind(tmp1, tmp2)
perm.list$utr3_w_gene_insert <- perm(gene_wo_lava$utr3, length(gene_w_lava$utr3_w_gene_insert), all.bs, reps=1000, min_cpgs_w_cov=4, parallel=T, cores=24)

save(res, perm.list, file=paste0(outdir, 'LAVA_perm_results.Rdata'))

# Get p-values
##############
res$perm_stats <- data.frame(type=names(perm.list))
res$p_value_lava <- data.frame(type=names(perm.list))

res$perm_stats$meth.mean <- sapply(perm.list, function(x) mean(unlist(x$meth.mean), na.rm=T))
res$perm_stats$cov.mean <- sapply(perm.list, function(x) mean(unlist(x$cov.mean), na.rm=T))
res$perm_stats$cpgs.mean <- sapply(perm.list, function(x) mean(unlist(x$cpgs.mean), na.rm=T))
res$perm_stats$cpgs_w_cov.mean <- sapply(perm.list, function(x) mean(unlist(x$cpgs_w_cov.mean), na.rm=T))
res$perm_stats$meth.sd <- sapply(perm.list, function(x) sd(unlist(x$meth.mean), na.rm=T))
res$perm_stats$cov.sd <- sapply(perm.list, function(x) sd(unlist(x$cov.mean), na.rm=T))
res$perm_stats$cpgs.sd <- sapply(perm.list, function(x) sd(unlist(x$cpgs.mean), na.rm=T))
res$perm_stats$cpgs_w_cov.sd <- sapply(perm.list, function(x) sd(unlist(x$cpgs_w_cov.mean), na.rm=T))
res$perm_stats$meth.median <- sapply(perm.list, function(x) median(unlist(x$meth.median), na.rm=T))
res$perm_stats$cov.median <- sapply(perm.list, function(x) median(unlist(x$cov.median), na.rm=T))
res$perm_stats$cpgs.median <- sapply(perm.list, function(x) median(unlist(x$cpgs.median), na.rm=T))
res$perm_stats$cpgs_w_cov.median <- sapply(perm.list, function(x) median(unlist(x$cpgs_w_cov.median), na.rm=T))
res$perm_stats$meth.ct <- sapply(perm.list, function(x) sum(!is.na(unlist(x$meth.mean))))
res$perm_stats$cov.ct <- sapply(perm.list, function(x) sum(!is.na(unlist(x$cov.mean))))
res$perm_stats$cpgs.ct <- sapply(perm.list, function(x) sum(!is.na(unlist(x$cpgs.mean))))
res$perm_stats$cpgs_w_cov.ct <- sapply(perm.list, function(x) sum(!is.na(unlist(x$cpgs_w_cov.mean))))

res$p_value_lava$meth.mean.p_value <- sapply(1:12, 
  function(i) get_p(res$w_lava.mean$meth[i], res$perm_stats$meth.mean[i], perm.list[[i]]$meth.mean))
res$p_value_lava$cov.mean.p_value <- sapply(1:12,
  function(i) get_p(res$w_lava.mean$cov[i], res$perm_stats$cov.mean[i], perm.list[[i]]$cov.mean))
res$p_value_lava$cpg.mean.p_value <- sapply(1:12,
  function(i) get_p(res$w_lava.mean$cpgs[i], res$perm_stats$cpgs.mean[i], perm.list[[i]]$cpgs.mean))
res$p_value_lava$cpg_w_cov.mean.p_value <- sapply(1:12,
  function(i) get_p(res$w_lava.mean$cpgs_w_cov[i], res$perm_stats$cpgs_w_cov.mean[i], perm.list[[i]]$cpgs_w_cov.mean))

res$p_value_lava$meth.median.p_value <- sapply(1:12, 
  function(i) get_p(res$w_lava.median$meth[i], res$perm_stats$meth.median[i], perm.list[[i]]$meth.median))
res$p_value_lava$cov.median.p_value <- sapply(1:12,
  function(i) get_p(res$w_lava.median$cov[i], res$perm_stats$cov.median[i], perm.list[[i]]$cov.median))
res$p_value_lava$cpg.median.p_value <- sapply(1:12,
  function(i) get_p(res$w_lava.median$cpgs[i], res$perm_stats$cpgs.median[i], perm.list[[i]]$cpgs.median))
res$p_value_lava$cpg_w_cov.median.p_value <- sapply(1:12,
  function(i) get_p(res$w_lava.median$cpgs_w_cov[i], res$perm_stats$cpgs_w_cov.median[i], perm.list[[i]]$cpgs_w_cov.median))

save(res, perm.list, file=paste0(outdir, 'LAVA_perm_results.Rdata'))
