#!/usr/bin/env Rscript

#nathan dot lazar at gmail dot com

#R script to colate data on CpG methylation from cpg10 evidence files

# Usage:
# Rscript ./bin/full_R_analysis.R
#   <base directory>
#   <bin directory>
#   <output directory> 
#   <file of chromosome lengths>
#   <sorted and tabulated cpg evidence file>
#   <breakpoint region file>
#   <gene_file>
#   <repmask file>
#   <cpg_island_file>

# Example: #######**********Change These *******#########
# Rscript ./bin/full_R_analysis.R 
#   /mnt/lustre1/users/lazar/APE_METH/POST_CRASH/
#   APE_METH_bin/
#   PERM_ANALYSIS/
#   NomLeu1_0/seq_len.txt
#   MAPPED_GOOD/Gibbon/DNA111101LC_62_HSA_normal_NoIndex_L006_1_val_1_trim.fq.gz/cpg10/
#   BREAKPOINTS/Gibbon_Breakpointslist_1_27_15.txt
#   GIBBON_FEATURES/NomLeu1.0.70.genes.gtf
#   GIBBON_FEATURES/NomLeu1.0_repmask.txt
#   GIBBON_FEATURES/gibbon_cpgislands.gff

library(bsseq)
condor <- TRUE  # Set variable to use cluster

args <- commandArgs(TRUE)
dir <- args[1]
bindir <- args[2]
outdir <- args[3]
len_file <- args[4]
cpg_drive <- args[5]
bp_file <- args[6]
gene_file <- args[7]
rep_file <- args[8]
cpg_isl_file <- args[9]

#bac_file <- 'BAC_amp_on_NomLeu1.0.txt'

source(paste0(bindir, 'R_meth_functions.R'))
source(paste0(bindir, 'plot_sides.R'))

# Read in lengths of chromsomes and make seqinfo object
#######################################################
seqinfo <- make_seqinfo(len_file, strsplit(len_file, "/")[[1]][1])

# Read in CpG data and make BSeq object with all CpGs
############################################################
all.bs <- make_all_bs(cpg_drive, strsplit(len_file, "/")[[1]][1], seqinfo, 0)
cov4.bs <- all.bs[getCoverage(all.bs) >=4]

# Write bed files of CpG coverage and methylation
#################################################
make_tracks(all.bs, outdir, 'meth.bedgraph', 'cov.bedgraph')

#Measure methylation of CpGs with coverage of 4 or more
########################################################
all.cpg.meth <- mean(mcgetMeth(cov4.bs, type='raw', what='perBase'), na.rm=T)

############################
# Breakpoint region analysis
############################

# Read in Breakpoint region data and make GRanges object
bp.gr <- read_bp(bp_file, seqinfo)

# Make GRanges object of 10kb regions on each side of breaks
# unless bp_file has 'BAC' in it or outdir name has 'inside' 
# or the species isn't gibbon (hack)
if(grepl('BAC', bp_file) | grepl('inside', outdir) |
   seqinfo@genome[1] != 'NomLeu1_0') {
  bp.lr.gr <- bp.gr
} else {
  bp.lr.gr <- make_lr_bp(bp.gr, 10000)
}

# Write breakpoint regions to bed file
bp.bed <- data.frame(seqnames(bp.lr.gr), start(bp.lr.gr), end(bp.lr.gr))
write.table(bp.bed, file=paste0(outdir,'bp_regions.bed'), quote=F, 
            sep='\t', row.names=F, col.names=F)

# Add mean methylation, number of CpGs and coverage
# of sides of bp regions to GRanges object
if(condor) {
  bp.lr.gr <- condor_add_meth_cpg_cov(bp.lr.gr, all.bs,
                dir, paste0(outdir, 'bp_add/'), bindir)
} else {
  bp.lr.gr <- add_meth_cpg_cov(bp.lr.gr, all.bs, parallel=F, min_cov=4)
}

# Write out coverage of targeted regions
write.table(data.frame(chr=as.vector(seqnames(bp.lr.gr)), 
                       start=start(bp.lr.gr), end=end(bp.lr.gr), 
                       cov=bp.lr.gr$cov), 
            'amp.cov.txt', quote=F, sep='\t', row.names=F,
            col.names=T)

# Run permutation analysis to determine whether the breakpoint
# regions have lower methylation, coverage or CpG counts than are 
# seen in random regions
bp.lr.gr$size <- width(bp.lr.gr)
min.chr.size <- max(bp.lr.gr$size)+2000 # Choose only from chroms that
                                        # are big enough to 
                                        # contain regions
bp.permute <- par_permute(outdir, bindir, bp.lr.gr, bp.lr.gr, all.bs, 
                          n=1000, type='all', 
                          min.chr.size=min.chr.size, end.exclude=1000)

###############
# Gene analysis
###############

gene.gr.list <- gtf2GRanges(gene_file, seqinfo, prom_size=1000)

# Add methylation, number of CpGs and mean CpG coverage to each
# of these GRanges objects
gene.gr.list <- lapply(gene.gr.list, add_meth_cpg_cov, all.bs)

# Run permutation analyses
gene.permute <- par_permute(outdir, bindir, gene.gr.list$gene, 
                            bp.lr.gr, all.bs, n=1000, type='gene', 
                            min.chr.size=min.chr.size, end.exclude=1000)
exon.permute <- par_permute(outdir, bindir, gene.gr.list$exon, 
                            bp.lr.gr, all.bs, n=1000, type='exon', 
                            min.chr.size=min.chr.size, end.exclude=1000)
intron.permute <- par_permute(outdir, bindir, gene.gr.list$intron, 
                              bp.lr.gr, all.bs, n=1000, type='intron',
                              min.chr.size=min.chr.size, end.exclude=1000)
promoter.permute <- par_permute(outdir, bindir, gene.gr.list$promoter, 
                                bp.lr.gr, all.bs, n=1000, type='promoter',
                                min.chr.size=min.chr.size, end.exclude=1000)
gene.permute.list <- list(gene=gene.permute, exon=exon.permute, 
                          intron=intron.permute, promoter=promoter.permute)

#################
# Repeat analysis
#################

# Read in repmask data
rep.gr <- make_rep_gr(rep_file, seqinfo(all.bs))

# Add mean methylation, number of CpGs and mean CpG
# coverage to repeat GRanges object
# rep.gr <- add_meth_cpg_cov(rep.gr, all.bs)

rep.gr <- condor_add_meth_cpg_cov(rep.gr, all.bs, 
            bin, paste0(outdir, 'rep_add'), bindir)

rep.permute <- par_permute(outdir, bindir, rep.gr, bp.lr.gr, 
                           all.bs, n=1000,
                           type='repeat', min.chr.size=min.chr.size,
                           end.exclude=1000)

# Run permutation analysis on major repeat classes 
# LINE, SINE, DNA, LTR, Satellite

rep.gr.list <- list()
rep.gr.list$LINE <- rep.gr[grepl('LINE', rep.gr$rep_family)]
rep.gr.list$SINE <- rep.gr[grepl('SINE', rep.gr$rep_family)]
rep.gr.list$DNA <- rep.gr[grepl('DNA', rep.gr$rep_family)]
rep.gr.list$LTR <- rep.gr[grepl('LTR', rep.gr$rep_family)]
rep.gr.list$Satellite <- rep.gr[grepl('Satellite', rep.gr$rep_family)]

LINE.permute <- par_permute(outdir, bindir, rep.gr.list$LINE, bp.lr.gr,
                            all.bs, n=1000, type='LINE',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
SINE.permute <- par_permute(outdir, bindir, rep.gr.list$SINE, bp.lr.gr,
                            all.bs, n=1000, type='SINE',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
DNA.permute <- par_permute(outdir, bindir, rep.gr.list$DNA, bp.lr.gr,
                            all.bs, n=1000, type='DNA',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
LTR.permute <- par_permute(outdir, bindir, rep.gr.list$LTR, bp.lr.gr,
                            all.bs, n=1000, type='LTR',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
Satellite.permute <- par_permute(outdir, bindir, rep.gr.list$Satellite, 
                            bp.lr.gr, all.bs, n=1000, type='Satellite',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
rep.permute.list <- list(LINE=LINE.permute, SINE=SINE.permute,
                         DNA=DNA.permute, LTR=LTR.permute,
                         Satellite=Satellite.permute)

# Run permutation analysis on major classes of SINEs
# Alu, MIR, AluS, AluJ, AluY

SINE.gr.list <- list()
SINE.gr.list$Alu <- rep.gr[grepl('Alu', rep.gr$rep_class)]
SINE.gr.list$MIR <- rep.gr[grepl('MIR', rep.gr$rep_class)]
SINE.gr.list$AluS <- rep.gr[grepl('AluS', rep.gr$rep_class)]
SINE.gr.list$AluJ <- rep.gr[grepl('AluJ', rep.gr$rep_class)]
SINE.gr.list$AluY <- rep.gr[grepl('AluY', rep.gr$rep_class)] 

Alu.permute <-  par_permute(outdir, bindir, SINE.gr.list$Alu, bp.lr.gr,
                            all.bs, n=1000, type='Alu',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
MIR.permute <-  par_permute(outdir, bindir, SINE.gr.list$MIR, bp.lr.gr,
                            all.bs, n=1000, type='MIR',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
AluJ.permute <-  par_permute(outdir, bindir, SINE.gr.list$AluJ, bp.lr.gr,
                            all.bs, n=1000, type='AluJ',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
AluS.permute <-  par_permute(outdir, bindir, SINE.gr.list$AluS, bp.lr.gr,
                            all.bs, n=1000, type='AluS',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
AluY.permute <-  par_permute(outdir, bindir, SINE.gr.list$AluY, bp.lr.gr,
                            all.bs, n=1000, type='AluY',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
SINE.permute.list <- list(Alu=Alu.permute, MIR=MIR.permute,
                          AluJ=AluJ.permute, AluS=AluS.permute,
                          AluY=AluY.permute)

#####################
# CpG island analysis
#####################

# Read in CpG island data
cpg_island.gr <- gff2GRanges(cpg_isl_file, seqinfo)

cpg_shore.gr <- make_cpg_shore(cpg_island.gr, 1000)

# Add mean methylation, number of CpGs and mean CpG
# coverage to CpGisland GRanges object

cpg_island.gr <- condor_add_meth_cpg_cov(cpg_island.gr, all.bs,
            bin, paste0(outdir, 'cpg_island_add'), bindir)
cpg_shore.gr <- condor_add_meth_cpg_cov(cpg_shore.gr, all.bs,
            bin, paste0(outdir, 'cpg_shore_add'), bindir)

#cpg_island.gr <- add_meth_cpg_cov(cpg_island.gr, all.bs)
#cpg_shore.gr <- add_meth_cpg_cov(cpg_shore.gr, all.bs)

# Run permutation analysis
cpg_island.permute <- par_permute(outdir, bindir, cpg_island.gr, bp.lr.gr,
                                  all.bs, n=1000, type='CpGisl',
                                  min.chr.size=min.chr.size, end.exclude=1000)
cpg_shore.permute <- par_permute(outdir, bindir, cpg_shore.gr, bp.lr.gr,
                                 all.bs, n=1000, type='CpGshore',
                                 min.chr.size=min.chr.size, end.exclude=1000)

cpg.permute.list <- list(cpg_island=cpg_island.permute, 
                         cpg_shore=cpg_shore.permute)

#################################################
# Combine permutation results into one data frame
#################################################

per.results <- combine_per(bp.permute, gene.permute.list, rep.permute, 
  rep.permute.list, SINE.permute.list, cpg.permute.list)
save(per.results, file=paste0(outdir, 'per_results.Rdat'))

#################################################
# Plot methylation vs. coverage and cpg counts for 
# the sides of the breakpoints
#################################################

# Plots of methylation, coverage and cpgs
if(ncol(mcols(bp.lr.gr)) >= 6) {
  bp.lr.df <- data.frame(mcols(bp.lr.gr)[,-1])
  names(bp.lr.df)[7] <- 'methylation'
} else  {
  bp.lr.df <- data.frame(mcols(bp.lr.gr))
  names(bp.lr.df)[2] <- 'methylation'    ## This probably needs to be changed ##
}

png(file=paste0(outdir, "meth_cov_plot.png"), height=480, width=480)
  plot_meth_cov(bp.lr.df)
dev.off()
png(file=paste0(outdir, "meth_cpgs_plot.png"), height=480, width=480)
  plot_meth_cpgs(bp.lr.df)
dev.off()

#################################################
# Compare mean absolute differences in methylation, 
# coverage and cpg counts between the breakpoint 
# regions and randomly selected regions (both 
# adjacent and non-adjacent)
#################################################

# Run permutation analysis to see if the distribution of MAD values 
# for the set of breakpoints differs from what would be expected 
# by chance. This is done by obtaining 1,000 sets of adjacent regions
# with the same sizes as the BP regions and measuring the MAD
# Additionally, do the same with non-adjacent regions.

adj.sides <- par_permute_sides(paste0(dir,outdir), paste0(dir, bindir), bp.lr.gr, all.bs,
                               n=1000, adjacent=T, end.exclude=1000)

disj.sides <- par_permute_sides(paste0(dir,outdir), bindir, bp.lr.gr, all.bs,
                                n=1000, adjacent=F, end.exclude=1000)

# Use this to get a line on the plots showing the percentiles of MAD
# scores for adj and non-adj. individual regions

bp.lr.mad <- get_mad(bp.lr.gr)

source(paste0(bindir, 'plot_sides.R'))

png(file=paste0(outdir, "mad_meth.png"), height=480, width=480*1.62)
  plot_mad_meth(bp.lr.mad) +  
  geom_hline(data=adj.sides$region.percentiles[4:7,], aes(yintercept=mad.meth), linetype=2) +
  annotate("text", label=adj.sides$region.percentiles$pers[4:7], x=0, 
           y=adj.sides$region.percentiles$mad.meth[4:7], size=6, colour = "black", vjust=-0.2) +
  geom_hline(data=disj.sides$region.percentiles[4:7,], aes(yintercept=mad.meth), colour="blue", linetype=2) +
  annotate("text", label=disj.sides$region.percentiles$pers[4:7], x=.1, 
           y=disj.sides$region.percentiles$mad.meth[4:7], size=6, colour = "blue", vjust=-0.2)
dev.off()

png(file=paste0(outdir, "mad_cov.png"), height=480, width=480*1.62)
  plot_mad_cov(bp.lr.mad) +  
  geom_hline(data=adj.sides$region.percentiles[4:7,], aes(yintercept=mad.cov), linetype=2) +
  annotate("text", label=adj.sides$region.percentiles$pers[4:7], x=18 , 
           y=adj.sides$region.percentiles$mad.cov[4:7], size=6, colour = "black", vjust=-0.2) +
  geom_hline(data=disj.sides$region.percentiles[4:7,], aes(yintercept=mad.cov), colour="blue", linetype=2) +
  annotate("text", label=disj.sides$region.percentiles$pers[4:7], x=22 , 
           y=disj.sides$region.percentiles$mad.cov[4:7], size=6, colour = "blue", vjust=-0.2)
dev.off()

png(file=paste0(outdir, "mad_cpgs.png"), height=480, width=480*1.62)
  plot_mad_cpgs(bp.lr.mad) +  
  geom_hline(data=adj.sides$region.percentiles[4:7,], aes(yintercept=mad.cpgs), linetype=2) +
  annotate("text", label=adj.sides$region.percentiles$pers[4:7], x=-30 , 
           y=adj.sides$region.percentiles$mad.cpgs[4:7], size=6, colour = "black", vjust=-0.2) +
  geom_hline(data=disj.sides$region.percentiles[4:7,], aes(yintercept=mad.cpgs), colour="blue", linetype=2) +
  annotate("text", label=disj.sides$region.percentiles$pers[4:7], x=330 , 
           y=disj.sides$region.percentiles$mad.cpgs[4:7], size=6, colour = "blue", vjust=-0.2)
dev.off()

# Print out permutation group p-values
print("P-values for comparisons to groups of adjacent regions:")
print(paste('Methylation:', adj.sides$group.p_values$meth))
print(paste('Coverage:', adj.sides$group.p_values$cov))
print(paste('CpG counts:', adj.sides$group.p_values$cpgs))

print("P-values for comparisons to groups of disjoint regions:")
print(paste('Methylation:', disj.sides$group.p_values$meth))
print(paste('Coverage:', disj.sides$group.p_values$cov))
print(paste('CpG counts:', disj.sides$group.p_values$cpgs))

# Write out data
#********************************************************************
save.image(file=paste0(outdir, 'R.dat'))
#load(paste0(outdir, 'R.dat'))
#********************************************************************

# Look at differences between classes of breakpoints
####################################################

#bp_Class1.gr <- bp_region.gr[bp_region.gr$class=='Class_I']
#bp_Class2a.gr <- bp_region.gr[bp_region.gr$class=='Class_II-a']
#bp_Class2b.gr <- bp_region.gr[bp_region.gr$class=='Class_II-b']
