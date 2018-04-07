#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com
# code to make mbias plot for the given mbias.tsv file

# Usage: mbias_plot.R mbias.tsv 

.libPaths("/home/users/lazar/R/x86_64-redhat-linux-gnu-library/3.1")

library(ggplot2)
library(reshape2)
library(dplyr)
library(grid)

args <- commandArgs(TRUE)

file <- args[1] 

mbias <- read.table(file, header=T, sep='\t', stringsAsFactors=F)
mbias$per <- with(mbias, C/(C+T+Other)*100)
mbias$read_len <- max(mbias$Offset)

mk_plot <- function(data) {
  p <- ggplot(data=data, aes(x=Offset, y=per)) +
    geom_line(size=1.25) + ylim(50,90) +
    ggtitle("Methylation bias plot") +
    labs(x = "Position in read", y= "Methylation level")
  p
}

p <- mk_plot(mbias)

png('mbias_plot.png', height=(480*1.5)/1.618, width=480*1.5)
p
dev.off()
