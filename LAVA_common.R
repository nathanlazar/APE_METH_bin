#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# Compare methylation, coverage and CpG counts of genes containing 
# LAVA elements in gibbon to other species. 
# Restrict analysis to genes shared in all genomes:
# gibbon, human, chimp, rhesus & orangutan

# Usage: LAVA_common.R \
#   /mnt/lustre1/users/lazar/APE_METH/POST_CRASH/ \
#   APE_METH_bin/ \
#   LAVA/All/ \
#   LAVA/Gibbon/LAVAs.txt

.libPaths("/home/users/lazar/R/x86_64-redhat-linux-gnu-library/3.1")
library(bsseq)

args <- commandArgs(TRUE)
dir <- args[1]
bindir <- args[2]
outdir <- args[3]
lava_file <- args[7]

source(paste0(bindir, 'R_meth_functions.R'))

# Functions used later
######################
get_p <- function(obs, center, perm) {
  diff <- abs(obs - center)
  p <- sum((perm > (center + diff) | perm < (center-diff)), na.rm=T)/length(perm)
  p
}

perm <- function(gr, num, all.bs, reps=1000, min_cpgs_w_cov=4) {
  # Function to get 1000 sets of features w/o LAVA
  # insertions and measure their mean and median
  # methylation, coverage and cpg counts

  # Get random gene sets
  rand.gr <- gr[sample(length(gr), num*reps, replace=T)]

  # Add methylation, coverage, etc. info to genes
  rand.gr <- add_meth_cpg_cov(rand.gr, all.bs, parallel=T, min_cov=min_cpgs_w_cov)

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

# Read in lengths of chromsomes and make seqinfo objects for each genome
########################################################################
gibbon.seqinfo <- make_seqinfo('NomLeu1_0/seq_len.txt', 'NomLeu1_0')
human.seqinfo <- make_seqinfo('human_hg19_noIUPAC/lengths.txt', 'hg19')
chimp.seqinfo <- make_seqinfo('chimp_panTro4/lengths.txt', 'panTro4')
rhesus.seqinfo <- make_seqinfo('rhesus_v2/lengths.txt', 'rheMac2')
orangutan.seqinfo <- make_seqinfo('orangutan_ponAbe2/lengths.txt', 'ponAbe2')

# Read in genes
###############
gibbon.gene.gr.list <- UCSC2GRanges('EnsemblGenes/nomLeu1_0_Ensembl_genes.txt', 
                                    gibbon.seqinfo, prom_size=1000, filter=F)
human.gene.gr.list <- UCSC2GRanges('EnsemblGenes/hg19_Ensembl_genes.txt', 
                                    human.seqinfo, prom_size=1000, filter=F)
chimp.gene.gr.list <- UCSC2GRanges('EnsempanTro4_Ensembl_genesblGenes/.txt', 
                                   chimp.seqinfo, prom_size=1000, filter=F)
rhesus.gene.gr.list <- UCSC2GRanges('EnsemblGenes/rheMac2_Ensembl_genes.txt', 
                                    rhesus.seqinfo, prom_size=1000, filter=F)
orangutan.gene.gr.list <- UCSC2GRanges('EnsemblGenes/ponAbe2_Ensembl_genes.txt', 
                                       orangutan.seqinfo, prom_size=1000, filter=F)

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

# Overlap gibbon genes w/ LAVAs
########################
gene_w_lava <- list()
gene_w_lava$promoter <- GenomicRanges::subsetByOverlaps(gene.gr.list$promoter, lava.gr, ignore.strand=T)
gene_w_lava$gene <- GenomicRanges::subsetByOverlaps(gene.gr.list$gene, lava.gr, ignore.strand=T)
gene_w_lava$utr5 <- GenomicRanges::subsetByOverlaps(gene.gr.list$utr5, lava.gr, ignore.strand=T)
gene_w_lava$exon <- GenomicRanges::subsetByOverlaps(gene.gr.list$exon, lava.gr, ignore.strand=T)
gene_w_lava$intron <- GenomicRanges::subsetByOverlaps(gene.gr.list$intron, lava.gr, ignore.strand=T)
gene_w_lava$utr3 <- GenomicRanges::subsetByOverlaps(gene.gr.list$utr3, lava.gr, ignore.strand=T)

gene_w_lava$promoter <- trim(gene_w_lava$promoter)



gene.gr.list <- lapply(gene.gr.list, add_meth_cpg_cov, all.bs, parallel=T, min_cov=4)
save(gene.gr.list, file=paste0(outdir, 'gene_gr_list.Rdata'))


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

# Subset these genes to just genes present in other genomes 
# (human, chimp, rhesus, orangutan)
gibbon.genes <- gene.gr.list
load('LAVA/Human/gene_gr_list.Rdata')
human.genes <- gene.gr.list
for(x in names(human.genes)) human.genes[[x]]$symbol <- toupper(human.genes[[x]]$symbol)
load('LAVA/Chimp/gene_gr_list.Rdata')
chimp.genes <- gene.gr.list
for(x in names(chimp.genes)) chimp.genes[[x]]$symbol <- toupper(chimp.genes[[x]]$symbol)
load('LAVA/Rhesus/gene_gr_list.Rdata')
rhesus.genes <- gene.gr.list
for(x in names(rhesus.genes)) rhesus.genes[[x]]$symbol <- toupper(rhesus.genes[[x]]$symbol)
load('LAVA/Orangutan/gene_gr_list.Rdata')
orangutan.genes <- gene.gr.list
for(x in names(orangutan.genes)) orangutan.genes[[x]]$symbol <- toupper(orangutan.genes[[x]]$symbol)
gene.gr.list <- gibbon.genes

lava_in_all <- list()
for(type in names(gibbon.genes)) {
  lava_in_all[[type]] <- lava_genes$symbol %>% toupper %>%
  base::intersect(human.genes[[type]]$symbol) %>%
  base::intersect(chimp.genes[[type]]$symbol) %>%
  base::intersect(rhesus.genes[[type]]$symbol) %>%
  base::intersect(orangutan.genes[[type]]$symbol)
}

gene_w_lava <- list()
gene_w_lava$promoter <- GenomicRanges::subsetByOverlaps(gene.gr.list$promoter, lava.gr, ignore.strand=T)
gene_w_lava$gene <- GenomicRanges::subsetByOverlaps(gene.gr.list$gene, lava.gr, ignore.strand=T)
gene_w_lava$utr5 <- GenomicRanges::subsetByOverlaps(gene.gr.list$utr5, lava.gr, ignore.strand=T)
gene_w_lava$exon <- GenomicRanges::subsetByOverlaps(gene.gr.list$exon, lava.gr, ignore.strand=T)
gene_w_lava$intron <- GenomicRanges::subsetByOverlaps(gene.gr.list$intron, lava.gr, ignore.strand=T)
gene_w_lava$utr3 <- GenomicRanges::subsetByOverlaps(gene.gr.list$utr3, lava.gr, ignore.strand=T)

gene_w_lava.all <- list()
for(type in names(gene_w_lava)) {
  gene_w_lava.all[[type]] <- unique(gene_w_lava[[type]][toupper(gene_w_lava[[type]]$symbol) %in% lava_in_all[[type]]])
}

##########*************HERE*********************##############

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
types <- c(names(gene.gr.list), 'prom_w_gene_insert', 'gene_w_gene_insert', 
          'utr5_w_gene_insert', 'exon_w_gene_insert', 
          'intron_w_gene_insert', 'utr3_w_gene_insert')    

res$minus_lava.mean <- data.frame(type=types)
res$minus_lava.median <- data.frame(type=types)

w_insert.gr.list <- lapply(gene.gr.list, function(x) x[x$name %in% unique(lava_genes$name)])

# Set diff with '*' stranded elements doesn't perform as expected. Run it twice with 
# '+' and '-' for the lava elements

strand(lava.gr) <- '+' 
minus_lava <- c(lapply(gene_w_lava, GenomicRanges::setdiff, lava.gr),
                lapply(w_insert.gr.list, GenomicRanges::setdiff, lava.gr))
strand(lava.gr) <- '-' 
minus_lava <- lapply(minus_lava, GenomicRanges::setdiff, lava.gr)
strand(lava.gr) <- '*'

names(minus_lava) <- types

minus_lava <- lapply(minus_lava, add_meth_cpg_cov, all.bs, parallel=T, min_cov=4)
save(minus_lava, file=paste0(outdir, 'minus_lava.Rdata'))

# Add back in the gene name and symbol for genes etc. that got split by LAVA elements ?
###*** Figure out how to get measurements for the whole gene even though it's now cut in two by ***###
###*** the lava coordinates... use group_by?                                                    ***###
# I decided not to do this. We're looking means across elements, so this shouldn't make that much 
# of a difference.

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
perm.list$promoter <- perm(gene_wo_lava$promoter, length(gene_w_lava$promoter), all.bs, reps=1000, min_cpgs_w_cov=4)
perm.list$gene <- perm(gene_wo_lava$gene, length(gene_w_lava$gene), all.bs, reps=1000, min_cpgs_w_cov=4)
perm.list$utr5 <- perm(gene_wo_lava$utr5, length(gene_w_lava$utr5), all.bs, reps=1000, min_cpgs_w_cov=4)
perm.list$exon <- perm(gene_wo_lava$exon, length(gene_w_lava$exon), all.bs, reps=1000, min_cpgs_w_cov=4)
perm.list$intron <- perm(gene_wo_lava$intron, length(gene_w_lava$intron), all.bs, reps=1000, min_cpgs_w_cov=4)
perm.list$utr3 <- perm(gene_wo_lava$utr3, length(gene_w_lava$utr3), all.bs, reps=1000, min_cpgs_w_cov=4)

save(res, perm.list, file=paste0(outdir, 'perm_results.Rdata'))

# Get p-values
###############
res$perm_stats <- data.frame(type=names(gene.gr.list))
res$p_value_lava <- data.frame(type=names(gene.gr.list))
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
