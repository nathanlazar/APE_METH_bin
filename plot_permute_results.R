#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# Script to plot the results of permutation analyses

library(ggplot2)
library(dplyr)
library(reshape2)

# Load data stored from permutation analyses

load('PERM_ANALYSIS/Gibbon/per_results.Rdat')
gibbon.per.results <- data.frame(per.results)
gibbon.per.results$species <- 'gibbon'
gibbon.per.results$name <- row.names(gibbon.per.results)
load('PERM_ANALYSIS/Human/per_results.Rdat')
human.per.results <- data.frame(per.results)
human.per.results$species <- 'human'
human.per.results$name <- row.names(human.per.results)
load('PERM_ANALYSIS/Chimp/per_results.Rdat')
chimp.per.results <- data.frame(per.results)
chimp.per.results$species <- 'chimp'
chimp.per.results$name <- row.names(chimp.per.results)
load('PERM_ANALYSIS/Rhesus/per_results.Rdat')
rhesus.per.results <- data.frame(per.results)
rhesus.per.results$species <- 'rhesus'
rhesus.per.results$name <- row.names(rhesus.per.results)
load('PERM_ANALYSIS/Orangutan/per_results.Rdat')
orangutan.per.results <- data.frame(per.results)
orangutan.per.results$species <- 'orangutan'
orangutan.per.results$name <- row.names(orangutan.per.results)

# P-values for permutation results on SINEs for chimp can't be trusted b/c 
# there weren't enough elements with coverage
chimp.per.results[chimp.per.results$meth.n < 500,
                  c('meth.p', 'cov.p', 'cpg.p')] <- NaN

per.results <- rbind(gibbon.per.results, human.per.results, 
                     chimp.per.results, orangutan.per.results,
                     rhesus.per.results)

# Calculate two-sided permutation p-values

# Replace zero values with 1/(number of repetitions run) and
# one values with 1-1/(number of repetitions)
per.results$meth.p[per.results$meth.p==0 & !is.na(per.results$meth.p)] <- 
  1/per.results$meth.n[per.results$meth.p==0 & !is.na(per.results$meth.p)]
per.results$meth.p[per.results$meth.p==1 & !is.na(per.results$meth.p)] <- 
  1 - 1/per.results$meth.n[per.results$meth.p==1 & !is.na(per.results$meth.p)]

per.results$cov.p[per.results$cov.p==0 & !is.na(per.results$cov.p)] <- 
  1/per.results$cov.n[per.results$cov.p==0 & !is.na(per.results$cov.p)]
per.results$cov.p[per.results$cov.p==1 & !is.na(per.results$cov.p)] <- 
  1 - 1/per.results$cov.n[per.results$cov.p==1 & !is.na(per.results$cov.p)]

per.results$cpg.p[per.results$cpg.p==0 & !is.na(per.results$cpg.p)] <- 
  1/per.results$cpg.n[per.results$cpg.p==0 & !is.na(per.results$cpg.p)]
per.results$cpg.p[per.results$cpg.p==1 & !is.na(per.results$cpg.p)] <- 
  1 - 1/per.results$cpg.n[per.results$cpg.p==1 & !is.na(per.results$cpg.p)]

per.results$per.p[per.results$per.p==0 & !is.na(per.results$per.p)] <- 
  1/per.results$per.n[per.results$per.p==0 & !is.na(per.results$per.p)]
per.results$per.p[per.results$per.p==1 & !is.na(per.results$per.p)] <- 
  1 - 1/per.results$per.n[per.results$per.p==1 & !is.na(per.results$per.p)]

# Reshape results for plotting with ggplot2
perm.df <- per.results %>%
  select(meth.p, cov.p, cpg.p, per.p, species, name) %>%
  melt(value.name="p") %>%
  filter(!(is.na(p)))

# Make two-sided p-values
p_val <- 2 * sapply(perm.df$p, function(x) min(1-x, x))

# Take negative log in base 10
perm.df$log_p <- -log10(p_val)

# Make these values negative if they represent a depletion in BP for the given species
perm.df$log_p[perm.df$p<0.5] <- -perm.df$log_p[perm.df$p<0.5]

# Rename variables for plotting
perm.df$variable <- as.vector(perm.df$variable)
perm.df$variable[perm.df$variable=='meth.p'] <- 'Methylation'
perm.df$variable[perm.df$variable=='cov.p'] <- 'Coverage'
perm.df$variable[perm.df$variable=='cpg.p'] <- 'CpG Count'
perm.df$variable[perm.df$variable=='per.p'] <- 'Cov %'
perm.df$variable <- factor(perm.df$variable,
                           levels=c('Methylation', 'CpG Count', 
                                    'Coverage', 'Cov %'))

plot_ps <- function(perm.df, features=c()) {
  # Subset and reorder features if given
  if(length(features)>0) {
    sub.df <- perm.df[perm.df$name %in% features,] 
    sub.df$name <- factor(sub.df$name, levels=features)
  } else {sub.df <- perm.df}
  
  # Order species
  sub.df$species <- factor(sub.df$species, levels=c('human', 'chimp', 
                                                    'orangutan', 'gibbon', 
                                                    'rhesus'))
  
  p <- ggplot(data=sub.df, aes(x=species, y=log_p, color=species)) +
    ylim(-3, 3) +
    #  facet_grid(. ~ species) +
    facet_grid(variable ~ name) +
    ggtitle("Negative log10 p-values") +
    geom_point(size=5) +  
    geom_segment(aes(x=species, y=0, yend=log_p, xend=species), size=1.5) +
    geom_hline(yintercept=c(log10(0.05),-log10(0.05)), linetype=2) +
    geom_hline(yintercept=0, size=1.5) +
    theme(axis.line.x=element_blank(), 
          axis.text.x=element_blank(),
          axis.title.x=element_blank(), 
          axis.ticks.x=element_blank(),
          strip.text.x=element_text(size = 16),
          strip.text.y=element_text(size = 16),
          plot.title=element_text(size=20, face="bold"),
          axis.title=element_text(size=16),
          axis.title.y=element_blank())
    # scale_fill_brewer(type='div', palette=5)
  p  
}

png(file='PERM_ANALYSIS/perm_p_plot_all.png', width=480*1.5)
plot_ps(perm.df)
dev.off()

png(file='PERM_ANALYSIS/perm_p_plot.png', width=480*1.5)
plot_ps(perm.df, c('all', 'reps', 'gene', 'cpg_isl', 'cpg_shore'))
dev.off()

png(file='PERM_ANALYSIS/perm_p_plot_genes.png', width=480*1.5)
plot_ps(perm.df, c('gene', 'exon', 'intron', 'promoter'))
dev.off()

png(file='PERM_ANALYSIS/perm_p_plot_reps.png', width=480*1.5)
plot_ps(perm.df, c('reps', 'SINE', 'LINE', 'DNA', 
                   'LTR', 'MIR', 'Satellite'))
dev.off()

png(file='PERM_ANALYSIS/perm_p_plot_SINE.png', width=480*1.5)
plot_ps(perm.df, c('SINE', 'Alu', 'AluJ', 'AluS', 'AluY'))
dev.off()
