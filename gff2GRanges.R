# Modified from Jiang (River) Li
# nathan dot lazar at gmail dot com

library(IRanges)
library(GenomicRanges)

# Returns GRanges object cpg_islands

gff2GRanges <- function(myfile="my.gff", seqinfo) {
  gff <- read.delim(myfile, header=FALSE)
  colnames(gff) <- c("chr", "source", "type", "start", "end", "score", "strand", "frame",      
                     "attributes")

  # Add in 'chr' if necessary
  if(sum(!grepl('chr', gff$chr))>0) gff$chr[!grepl('chr', gff$chr)] <- 
    paste0('chr', gff$chr[!grepl('chr', gff$chr)])

  len <- nrow(gff)

  Size <- rep('', len)                                      #get size from attributes column
  idx <- grepl('Size', gff$attributes)
  Size[idx] <- gsub(".*Size=(.*?);.*", "\\1", gff$attributes[idx])

  PercentCG <- rep('', len)                                 #get Percent CG
  idx <- grepl('PercentCG', gff$attributes)
  PercentCG[idx] <- 
    gsub(".*PercentCG=(.*?),.*", "\\1", gff$attributes[idx])

  ObsExp <- rep('', len)                                    #get ObsExp
  idx <- grepl('ObsExp', gff$attributes)
  ObsExp[idx] <- 
    gsub(".*ObsExp=(.*?)", "\\1", gff$attributes[idx])

  all.gr<-GRanges(seqnames=gff$chr,
                  ranges=IRanges(gff$start,gff$end),
                  strand=rep('*', len),
                  source=gff$source,
	          type=gff$type,
                  Size=as.numeric(Size),
                  PercentCG=as.numeric(PercentCG),
                  ObsExp=as.numeric(ObsExp))

  all.gr <- all.gr[seqnames(all.gr) %in% seqlevels(seqinfo)]
  seqlevels(all.gr) <- seqlevels(seqinfo)
  seqlengths(all.gr) <- seqlengths(seqinfo)

  all.gr
}