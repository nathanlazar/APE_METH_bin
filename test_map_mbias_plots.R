#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com
# code to make mbias plots
# for the test mapped 100,000 reads from all species

library(ggplot2)
library(reshape2)
library(dplyr)
library(grid)

dir <- 'TEST_MAP/'

files <- sort(list.files(path=dir, pattern = "test_mbias.tsv"))
names <- c('Chimp 1', 'Chimp 2', 'Gibbon 1', 'Orangutan 1', 'Orangutan 2', 
           'Orangutan 3', 'Gorilla 1', 'Gorilla 2', 'Gorilla 3', 'Human 1', 
           'Human 2', 'Rhesus 1', 'Gibbon 2', 'Gibbon 3')
trim5 <- rep(5, 14)
trim3 <- rep(5, 14)

mbias <- list()
mbias.df <- data.frame()
for(i in 1:length(files)) {
  file <- files[i]
  mbias[[i]] <- read.table(paste0(dir, file), header=T, sep='\t')
  mbias[[i]]$file <- sub('_test_mbias.tsv', '', file)
  mbias[[i]]$per <- with(mbias[[i]], C/(C+T+Other)*100)
  mbias[[i]]$name <- names[i]
  mbias[[i]]$read_len <- max(mbias[[i]]$Offset)
  mbias[[i]]$trim3 <- mbias[[i]]$read_len[1]-trim3[i]
  mbias[[i]]$trim5 <- trim5[i]
  mbias.df <- rbind(mbias.df, mbias[[i]])
}

mk_plot <- function(data) {
  p <- ggplot(data=data, aes(x=Offset, y=per, color=name)) +
    geom_line(size=1.25) + ylim(50,90) +
    geom_segment(aes(x=trim3, y=55, xend=trim3, yend=85),
                 size=1.2, linetype=2) +
    geom_segment(aes(x=trim5, y=55, xend=trim5, yend=85),
                 size=1.2, linetype=2) +
    ggtitle("Methylation bias of 100k reads") +
    labs(x = "Position in read", y= "Methylation level",
         color="Sequencing\nbatch") + 
    annotation_custom(grobTree(textGrob("Dotted lines show trimming", 
                                        x=0.1,  y=0.1, hjust=0,
                                        gp=gpar(fontsize=12)))) 
  # This complex annotation allows for relative location
  p
}

gibbon.df <- mbias.df %>% filter(grepl("Gibbon",name))
chimp.df <- mbias.df %>% filter(grepl("Chimp",name))
gorilla.df <- mbias.df %>% filter(grepl("Gorilla",name))
orangutan.df <- mbias.df %>% filter(grepl("Orangutan",name))
human.df <- mbias.df %>% filter(grepl("Human",name))
rhesus.df <- mbias.df %>% filter(grepl("Rhesus",name))

all.p <- mk_plot(mbias.df)
gibbon.p <- mk_plot(gibbon.df)
chimp.p <- mk_plot(chimp.df)
gorilla.p <- mk_plot(gorilla.df)
orangutan.p <- mk_plot(orangutan.df)
human.p <- mk_plot(human.df)
rhesus.p <- mk_plot(rhesus.df)

png('test_map_mbias.png', height=(480*1.5)/1.618, width=480*1.5)
all.p
dev.off()

pdf('test_map_mbias.pdf', height=(7*1.5)/1.618, width=7*1.5)
all.p
dev.off()
