#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# Compare methylation, coverage and CpG counts of genes containing 
# LAVA elements to genes without LAVA elements in gibbon NomLeu1.0

# Usage: LAVA_gibbons.R \
#   /mnt/lustre1/users/lazar/APE_METH/POST_CRASH/ \
#   APE_METH_bin/ \
#   LAVA/Gibbon/ \
#   NomLeu1_0/seq_len.txt \
#   MAPPED_GOOD/Gibbon/DNA111101LC_62_HSA_normal_NoIndex_L006_1_val_1_trim.fq.gz/cpg10/ \
#   EnsemblGenes/nomLeu1_0_Ensembl_genes.txt \
#   LAVA/Gibbon/LAVAs.txt

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
  p <- sum((perm > (center + diff) | perm < (center-diff)), na.rm=T)/length(perm)
  p
}

perm <- function(gr, num, all.bs, cov4.bs, reps=1000, parallel=T, cores=NA) {
  # Function to get 1000 sets of features w/o LAVA
  # insertions and measure their mean and median
  # methylation, coverage and cpg counts

  # Get random gene sets
  rand.gr <- gr[sample(length(gr), num*reps, replace=T)]

  # Add methylation, coverage, etc. info to genes
  rand.gr <- add_meth_cpg_cov(rand.gr, all.bs, cov4.bs, parallel=T)

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
      perm.res$meth.mean[i] <- mean(rand.gr[(i*(num-1)+1):(i*num)]$meth, na.rm=T)
      perm.res$cov.mean[i] <- mean(rand.gr[(i*(num-1)+1):(i*num)]$cov, na.rm=T)
      perm.res$cpgs.mean[i] <- mean(rand.gr[(i*(num-1)+1):(i*num)]$cpgs, na.rm=T)
      perm.res$cpgs_w_cov.mean[i] <- mean(rand.gr[(i*(num-1)+1):(i*num)]$cpgs_w_cov, na.rm=T)
      perm.res$meth.median[i] <- median(rand.gr[(i*(num-1)+1):(i*num)]$meth, na.rm=T)
      perm.res$cov.median[i] <- median(rand.gr[(i*(num-1)+1):(i*num)]$cov, na.rm=T)
      perm.res$cpgs.median[i] <- median(rand.gr[(i*(num-1)+1):(i*num)]$cpgs, na.rm=T)
      perm.res$cpgs_w_cov.median[i] <- median(rand.gr[(i*(num-1)+1):(i*num)]$cpgs_w_cov, na.rm=T)
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

gene.gr.list <- lapply(gene.gr.list, add_meth_cpg_cov, all.bs, cov4.bs, parallel=T)
save(gene.gr.list, file=paste0(outdir, 'gene_gr_list.Rdata'))

# Reorder if necessary
gene.gr.list <- list(promoter = gene.gr.list$promoter,
                     gene = gene.gr.list$gene,
                     utr5 = gene.gr.list$utr5,
                     exon = gene.gr.list$exon,
                     intron = gene.gr.list$intron,
                     utr3 = gene.gr.list$utr3)

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
#lava.gr <- add_meth_cpg_cov(lava.gr, all.bs, cov4.bs, parallel=T)

# Overlap LAVAs w/ genes
########################
#lava_in_gene <- list()
#lava_in_gene$promoter <- GenomicRanges::subsetByOverlaps(lava.gr, gene.gr.list$promoter, ignore.strand=T)
#lava_in_gene$gene <- GenomicRanges::subsetByOverlaps(lava.gr, gene.gr.list$gene, ignore.strand=T)
#lava_in_gene$utr5 <- GenomicRanges::subsetByOverlaps(lava.gr, gene.gr.list$utr5, ignore.strand=T)
#lava_in_gene$exon <- GenomicRanges::subsetByOverlaps(lava.gr, gene.gr.list$exon, ignore.strand=T)
#lava_in_gene$intron <- GenomicRanges::subsetByOverlaps(lava.gr, gene.gr.list$intron, ignore.strand=T)
#lava_in_gene$utr3 <- GenomicRanges::subsetByOverlaps(lava.gr, gene.gr.list$utr3, ignore.strand=T)

# Overlap genes w/ LAVAs
########################
gene_w_lava <- list()
gene_w_lava$promoter <- GenomicRanges::subsetByOverlaps(gene.gr.list$promoter, lava.gr, ignore.strand=T)
gene_w_lava$gene <- GenomicRanges::subsetByOverlaps(gene.gr.list$gene, lava.gr, ignore.strand=T)
gene_w_lava$utr5 <- GenomicRanges::subsetByOverlaps(gene.gr.list$utr5, lava.gr, ignore.strand=T)
gene_w_lava$exon <- GenomicRanges::subsetByOverlaps(gene.gr.list$exon, lava.gr, ignore.strand=T)
gene_w_lava$intron <- GenomicRanges::subsetByOverlaps(gene.gr.list$intron, lava.gr, ignore.strand=T)
gene_w_lava$utr3 <- GenomicRanges::subsetByOverlaps(gene.gr.list$utr3, lava.gr, ignore.strand=T)

gene_w_lava$promoter <- trim(gene_w_lava$promoter)

# Write out the list of gene names and symbols that contain LAVA elements
# record gene name, symbol and where the LAVA inserted (intron, exon, etc.)
lava_genes <- data.frame(rbind(mcols(gene_w_lava$gene)[,c('name', 'symbol')],
                               mcols(gene_w_lava$promoter)[,c('name', 'symbol')]))
lava_genes <- unique(lava_genes)
lava_genes$type <- ''
lava_genes$type[lava_genes$name %in% gene_w_lava$intron$name] <- 'intron'
lava_genes$type[lava_genes$name %in% gene_w_lava$utr5$name] <- 
  paste0(lava_genes$type[lava_genes$name %in% gene_w_lava$utr5$name], ',utr5')
lava_genes$type[lava_genes$name %in% gene_w_lava$utr3$name] <- 
  paste0(lava_genes$type[lava_genes$name %in% gene_w_lava$utr3$name], ',utr3')
lava_genes$type[lava_genes$name %in% gene_w_lava$promoter$name] <- 
  paste0(lava_genes$type[lava_genes$name %in% gene_w_lava$promoter$name], ',promoter')
lava_genes$type[lava_genes$name %in% gene_w_lava$exon$name] <- 
  paste0(lava_genes$type[lava_genes$name %in% gene_w_lava$exon$name], ',exon')
lava_genes$type <- gsub('^,', '', lava_genes$type)

write.table(lava_genes, file=paste0(outdir, 'genes_w_lava.txt'), 
            sep='\t', quote=F, row.names=F)

# Add list of genes with insertions anywhere to gene_w_lava
w_insert.gr.list <- lapply(gene.gr.list, function(x) x[x$name %in% unique(lava_genes$name)])
names(w_insert.gr.list) <- paste0(names(w_insert.gr.list), '_w_gene_insert')

# Combine lists of elements w/ insertions and elements of genes w/ insertions anywhere
gene_w_lava <- c(gene_w_lava, w_insert.gr.list)
save(gene_w_lava, file=paste0(outdir, 'gene_w_lava.Rdata'))


# Get genes not in LAVAs
########################
gene_wo_lava <- list()
gene_wo_lava$promoter <- gene.gr.list$promoter[!(gene.gr.list$promoter %in% gene_w_lava$promoter)]
gene_wo_lava$gene <- gene.gr.list$gene[!(gene.gr.list$gene %in% gene_w_lava$gene)]
gene_wo_lava$utr5 <- gene.gr.list$utr5[!(gene.gr.list$utr5 %in% gene_w_lava$utr5)]
gene_wo_lava$exon <- gene.gr.list$exon[!(gene.gr.list$exon %in% gene_w_lava$exon)]
gene_wo_lava$intron <- gene.gr.list$intron[!(gene.gr.list$intron %in% gene_w_lava$intron)]
gene_wo_lava$utr3 <- gene.gr.list$utr3[!(gene.gr.list$utr3 %in% gene_w_lava$utr3)]

gene_wo_lava$promoter <- trim(gene_wo_lava$promoter)
gene_wo_lava$gene <- trim(gene_wo_lava$gene)
gene_wo_lava$utr5 <- trim(gene_wo_lava$utr5)
gene_wo_lava$exon <- trim(gene_wo_lava$exon)

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

# Get mean and median methylation, coverage & cpgs of genic sequence not counting LAVA
# instertions 
######################################################################################
types <- c(names(gene.gr.list), 'promoter_w_gene_insert', 'gene_w_gene_insert', 
          'utr5_w_gene_insert', 'exon_w_gene_insert', 
          'intron_w_gene_insert', 'utr3_w_gene_insert')    

res$minus_lava.mean <- data.frame(type=types)
res$minus_lava.median <- data.frame(type=types)

# Split each list of gene elements by their name.
split.list <- lapply(gene_w_lava, function(x) split(x, x$name))

# Remove LAVA sequences. This splits genes, etc. in half usually.
# Set diff with '*' stranded elements doesn't perform as expected. Run it twice with 
# '+' and '-' for the lava elements
strand(lava.gr) <- '+' 
minus_lava <- lapply(split.list, function(x) lapply(x, GenomicRanges::setdiff, lava.gr))
strand(lava.gr) <- '-'
minus_lava <- lapply(minus_lava, function(x) lapply(x, GenomicRanges::setdiff, lava.gr))

# Remove empty ranges
for(i in 1:length(minus_lava)) {
  minus_lava[[i]] <- minus_lava[[i]][sapply(minus_lava$exon, function(x) length(width(x))>0)]
}

minus_lava <- lapply(minus_lava, GRangesList)
minus_lava <- lapply(minus_lava, unlist, use.names=T)

minus_lava <- lapply(minus_lava, add_meth_cpg_cov, all.bs, cov4.bs, parallel=F)

# Add name as metadata
for(i in 1:length(minus_lava)) {
  minus_lava[[i]]$name <- names(minus_lava[[i]])
}

# Summarize meth, cov, cpgs, by name
minus_lava_sum <- list()
for(i in 1:length(minus_lava)) {
  minus_lava_sum[[names(minus_lava)[i]]] <- mcols(minus_lava[[i]]) %>% 
    data.frame %>%
    tbl_df %>%
    group_by(name) %>%
    summarise(meth=sum(meth*cpgs_w_cov)/sum(cpgs_w_cov),
              cov=sum(cov*cpgs)/sum(cpgs),
              cpgs_w_cov=sum(cpgs_w_cov), 
              cpgs=sum(cpgs))
}

save(minus_lava, minus_lava_sum, file=paste0(outdir, 'minus_lava.Rdata'))

res$minus_lava.mean$counts <- sapply(minus_lava_sum, nrow)
res$minus_lava.mean$meth <- sapply(minus_lava_sum, function(x) mean(x$meth, na.rm=T))
res$minus_lava.mean$cov <- sapply(minus_lava_sum, function(x) mean(x$cov, na.rm=T))
res$minus_lava.mean$cpgs <- sapply(minus_lava_sum, function(x) mean(x$cpgs, na.rm=T))
res$minus_lava.mean$cpgs_w_cov <- sapply(minus_lava_sum, function(x) mean(x$cpgs_w_cov, na.rm=T))

res$minus_lava.median$counts <- sapply(minus_lava_sum, nrow)
res$minus_lava.median$meth <- sapply(minus_lava_sum, function(x) median(x$meth, na.rm=T))
res$minus_lava.median$cov <- sapply(minus_lava_sum, function(x) median(x$cov, na.rm=T))
res$minus_lava.median$cpgs <- sapply(minus_lava_sum, function(x) median(x$cpgs, na.rm=T))
res$minus_lava.median$cpgs_w_cov <- sapply(minus_lava_sum, function(x) median(x$cpgs_w_cov, na.rm=T))

# Run permutations
##################
perm.list <- list()
reps <- 1000
parallel <- T
cores <- 24
perm.list$promoter <- perm(gene_wo_lava$promoter, nrow(minus_lava_sum$promoter), all.bs, cov4.bs, reps=reps, parallel=parallel, cores=cores)
perm.list$gene <- perm(gene_wo_lava$gene, nrow(minus_lava_sum$gene), all.bs, cov4.bs, reps=reps, parallel=parallel, cores=cores)
perm.list$utr5 <- perm(gene_wo_lava$utr5, nrow(minus_lava_sum$utr5), all.bs, cov4.bs, reps=reps, parallel=parallel, cores=cores)
perm.list$exon <- perm(gene_wo_lava$exon, nrow(minus_lava_sum$exon), all.bs, cov4.bs, reps=reps, parallel=parallel, cores=cores)
perm.list$intron <- perm(gene_wo_lava$intron, nrow(minus_lava_sum$intron), all.bs, cov4.bs, reps=reps, parallel=parallel, cores=cores)
perm.list$utr3 <- perm(gene_wo_lava$utr3, nrow(minus_lava_sum$utr3), all.bs, cov4.bs, reps=reps, parallel=parallel, cores=cores)

perm.list$promoter_w_gene_insert <- perm(gene_wo_lava$promoter, 
  nrow(minus_lava_sum$promoter_w_gene_insert), all.bs, cov4.bs, reps=reps, parallel=parallel, cores=cores)
perm.list$gene_w_gene_insert <- perm(gene_wo_lava$gene, 
  nrow(minus_lava_sum$gene_w_gene_insert), all.bs, cov4.bs, reps=reps, parallel=parallel, cores=cores)
perm.list$utr5_w_gene_insert <- perm(gene_wo_lava$utr5, 
  nrow(minus_lava_sum$utr5_w_gene_insert), all.bs, cov4.bs, reps=reps, parallel=parallel, cores=cores)
perm.list$exon_w_gene_insert <- perm(gene_wo_lava$exon, 
  nrow(minus_lava_sum$exon_w_gene_insert), all.bs, cov4.bs, reps=reps, parallel=parallel, cores=cores)
perm.list$intron_w_gene_insert <- perm(gene_wo_lava$intron, 
  nrow(minus_lava_sum$intron_w_gene_insert), all.bs, cov4.bs, reps=reps, parallel=parallel, cores=cores)
perm.list$utr3_w_gene_insert <- perm(gene_wo_lava$utr3, 
  nrow(minus_lava_sum$utr3_w_gene_insert), all.bs, cov4.bs, reps=reps, parallel=parallel, cores=cores)

# Convert columns of perm list to vectors 
for(i in 1:length(perm.list)) {
  for (j in 1:ncol(perm.list[[i]])) {
    perm.list[[i]][,j] <- unlist(perm.list[[i]][j])
  }
}

save(res, perm.list, file=paste0(outdir, 'perm_results.Rdata'))

# Get p-values
###############
res$perm_stats <- data.frame(type=types)
res$p_value_lava <- data.frame(type=types)
res$p_value_minus_lava <- data.frame(type=types)

res$perm_stats$meth.mean <- sapply(perm.list, function(x) mean(x$meth.mean, na.rm=T))
res$perm_stats$cov.mean <- sapply(perm.list, function(x) mean(x$cov.mean, na.rm=T))
res$perm_stats$cpgs.mean <- sapply(perm.list, function(x) mean(x$cpgs.mean, na.rm=T))
res$perm_stats$cpgs_w_cov.mean <- sapply(perm.list, function(x) mean(x$cpgs_w_cov.mean, na.rm=T))
res$perm_stats$meth.sd <- sapply(perm.list, function(x) sd(x$meth.mean, na.rm=T))
res$perm_stats$cov.sd <- sapply(perm.list, function(x) sd(x$cov.mean, na.rm=T))
res$perm_stats$cpgs.sd <- sapply(perm.list, function(x) sd(x$cpgs.mean, na.rm=T))
res$perm_stats$cpgs_w_cov.sd <- sapply(perm.list, function(x) sd(x$cpgs_w_cov.mean, na.rm=T))
res$perm_stats$meth.median <- sapply(perm.list, function(x) median(x$meth.median, na.rm=T))
res$perm_stats$cov.median <- sapply(perm.list, function(x) median(x$cov.median, na.rm=T))
res$perm_stats$cpgs.median <- sapply(perm.list, function(x) median(x$cpgs.median, na.rm=T))
res$perm_stats$cpgs_w_cov.median <- sapply(perm.list, function(x) median(x$cpgs_w_cov.median, na.rm=T))
res$perm_stats$meth.ct <- sapply(perm.list, function(x) sum(!is.na(x$meth.median)))
res$perm_stats$cov.ct <- sapply(perm.list, function(x) sum(!is.na(x$cov.median)))
res$perm_stats$cpgs.ct <- sapply(perm.list, function(x) sum(!is.na(x$cpgs.median)))
res$perm_stats$cpgs_w_cov.ct <- sapply(perm.list, function(x) sum(!is.na(x$cpgs_w_cov.median)))

res$p_value_lava$meth.mean.p_value <- sapply(1:6, 
  function(i) get_p(res$w_lava.mean$meth[i], res$perm_stats$meth.mean[i], perm.list[[i]]$meth.mean))
res$p_value_lava$cov.mean.p_value <- sapply(1:6,
  function(i) get_p(res$w_lava.mean$cov[i], res$perm_stats$cov.mean[i], perm.list[[i]]$cov.mean))
res$p_value_lava$cpg.mean.p_value <- sapply(1:6,
  function(i) get_p(res$w_lava.mean$cpgs[i], res$perm_stats$cpgs.mean[i], perm.list[[i]]$cpgs.mean))
res$p_value_lava$cpg_w_cov.mean.p_value <- sapply(1:6,
  function(i) get_p(res$w_lava.mean$cpgs_w_cov[i], res$perm_stats$cpgs_w_cov.mean[i], perm.list[[i]]$cpgs_w_cov.mean))

res$p_value_lava$meth.median.p_value <- sapply(1:6, 
  function(i) get_p(res$w_lava.median$meth[i], res$perm_stats$meth.median[i], perm.list[[i]]$meth.median))
res$p_value_lava$cov.median.p_value <- sapply(1:6,
  function(i) get_p(res$w_lava.median$cov[i], res$perm_stats$cov.median[i], perm.list[[i]]$cov.median))
res$p_value_lava$cpg.median.p_value <- sapply(1:6,
  function(i) get_p(res$w_lava.median$cpgs[i], res$perm_stats$cpgs.median[i], perm.list[[i]]$cpgs.median))
res$p_value_lava$cpg_w_cov.median.p_value <- sapply(1:6,
  function(i) get_p(res$w_lava.median$cpgs_w_cov[i], res$perm_stats$cpgs_w_cov.median[i], perm.list[[i]]$cpgs_w_cov.median))

res$p_value_minus_lava$meth.mean.p_value <- c(sapply(1:6, 
  function(i) get_p(res$minus_lava.mean$meth[i], res$perm_stats$meth.mean[i], perm.list[[i]]$meth.mean)),
  sapply(1:6, function(i) get_p(res$minus_lava.mean$meth[i+6], res$perm_stats$meth.mean[i],  perm.list[[i]]$meth.mean)))
res$p_value_minus_lava$cov.mean.p_value <- c(sapply(1:6, 
  function(i) get_p(res$minus_lava.mean$cov[i], res$perm_stats$cov.mean[i], perm.list[[i]]$cov.mean)),
  sapply(1:6, function(i) get_p(res$minus_lava.mean$cov[i+6], res$perm_stats$cov.mean[i],  perm.list[[i]]$cov.mean)))
res$p_value_minus_lava$cpgs.mean.p_value <- c(sapply(1:6,
  function(i) get_p(res$minus_lava.mean$cpgs[i], res$perm_stats$cpgs.mean[i], perm.list[[i]]$cpgs.mean)),
  sapply(1:6, function(i) get_p(res$minus_lava.mean$cpgs[i+6], res$perm_stats$cpgs.mean[i],  perm.list[[i]]$cpgs.mean)))
res$p_value_minus_lava$cpgs_w_cov.mean.p_value <- c(sapply(1:6,
  function(i) get_p(res$minus_lava.mean$cpgs_w_cov[i], res$perm_stats$cpgs_w_cov.mean[i], perm.list[[i]]$cpgs_w_cov.mean)),
  sapply(1:6, function(i) get_p(res$minus_lava.mean$cpgs_w_cov[i+6], res$perm_stats$cpgs_w_cov.mean[i],  perm.list[[i]]$cpgs_w_cov.mean)))

res$p_value_minus_lava$meth.median.p_value <- c(sapply(1:6, 
  function(i) get_p(res$minus_lava.median$meth[i], res$perm_stats$meth.median[i], perm.list[[i]]$meth.median)),
  sapply(1:6, function(i) get_p(res$minus_lava.median$meth[i+6], res$perm_stats$meth.median[i],  perm.list[[i]]$meth.median)))
res$p_value_minus_lava$cov.median.p_value <- c(sapply(1:6, 
  function(i) get_p(res$minus_lava.median$cov[i], res$perm_stats$cov.median[i], perm.list[[i]]$cov.median)),
  sapply(1:6, function(i) get_p(res$minus_lava.median$cov[i+6], res$perm_stats$cov.median[i],  perm.list[[i]]$cov.median)))
res$p_value_minus_lava$cpgs.median.p_value <- c(sapply(1:6,
  function(i) get_p(res$minus_lava.median$cpgs[i], res$perm_stats$cpgs.median[i], perm.list[[i]]$cpgs.median)),
  sapply(1:6, function(i) get_p(res$minus_lava.median$cpgs[i+6], res$perm_stats$cpgs.median[i],  perm.list[[i]]$cpgs.median)))
res$p_value_minus_lava$cpgs_w_cov.median.p_value <- c(sapply(1:6,
  function(i) get_p(res$minus_lava.median$cpgs_w_cov[i], res$perm_stats$cpgs_w_cov.median[i], perm.list[[i]]$cpgs_w_cov.median)),
  sapply(1:6, function(i) get_p(res$minus_lava.median$cpgs_w_cov[i+6], res$perm_stats$cpgs_w_cov.median[i],  perm.list[[i]]$cpgs_w_cov.median)))

save(res, file=paste0(outdir, 'LAVA_perm_results.Rdata'))
