# Modified from Jiang (River) Li
# nathan dot lazar at gmail dot com

library(IRanges)
library(GenomicRanges)

# Returns GRangesList object of genes, exons, introns, 
# promoters of length <prom_size>, 3' UTR, 5' UTR
# from information in gtf file

gtf2GRanges <- function(myfile="my.gff", seqinfo, prom_size=1000) {
  gtf <- read.delim(myfile, header=FALSE)
  colnames(gtf) <- c("chr", "source", "type", "start", "end", "score", "strand", "frame",      
                     "attributes")
  gtf$chr <- paste0('chr', gtf$chr)

  len <- nrow(gtf)

  gene_id <- rep('', len)                                      #get gene_id from attributes column
  idx <- grepl('gene_id', gtf$attributes)
  gene_id[idx] <- 
    gsub(".*gene_id (.*?);.*", "\\1", gtf$attributes[idx])

  transcript_id <- rep('', len)                                #get transcript_id
  idx <- grepl('transcript_id', gtf$attributes)
  transcript_id[idx] <- 
    gsub(".*transcript_id (.*?);.*", "\\1", gtf$attributes[idx])

  exon_number <- rep('', len)                                  #get exon_number
  idx <- grepl('exon_number', gtf$attributes)
  exon_number[idx] <- 
    gsub(".*exon_number (.*?);.*", "\\1", gtf$attributes[idx])

  gene_name <- rep('', len)                                    #get gene_name
  idx <- grepl('gene_name', gtf$attributes)
  gene_name[idx] <- 
    gsub(".*gene_name (.*?);.*", "\\1", gtf$attributes[idx])

  gene_biotype <- rep('', len)                                 #get gene_biotype
  idx <- grepl('gene_biotype', gtf$attributes)
  gene_biotype[idx] <- 
    gsub(".*gene_biotype (.*?);.*", "\\1", gtf$attributes[idx])

  transcript_name <- rep('', len)                              #get transcript_name
  idx <- grepl('transcript_name', gtf$attributes)
  transcript_name[idx] <- 
    gsub(".*transcript_name (.*?);.*", "\\1", gtf$attributes[idx])

  exon_id <- rep('', len)                                      #get exon_id
  idx <- grepl('exon_id', gtf$attributes)
  exon_id[idx] <- 
    gsub(".*exon_id (.*?);.*", "\\1", gtf$attributes[idx])

  all.gr<-GRanges(seqnames=gtf$chr,
                  ranges=IRanges(gtf$start,gtf$end),
                  strand=gtf$strand,
                  source=gtf$source,
	          type=gtf$type,
                  gene_id=gene_id,
                  transcript_id=transcript_id,
                  exon_number=as.numeric(exon_number),
	          gene_name=gene_name,
	          gene_biotype=gene_biotype,
	          transcript_name=transcript_name,
	          exon_id=exon_id)

  seqlevels(all.gr) <- seqlevels(seqinfo)
  seqlengths(all.gr) <- seqlengths(seqinfo)

  # Make gene GRanges object
  ##########################
  idx <- !grepl('pseudo', all.gr$source)
  gene.gr.list <- split(all.gr[idx], all.gr$gene_id[idx])
  gene.gr <- unlist(range(gene.gr.list))
  gene.gr$gene_id <- names(gene.gr.list)

  gene.tree <- GIntervalTree(gene.gr)     #used to subsetByOverlaps efficiently

  # Make exon GRanges object
  ##########################
  exon.gr <- all.gr[all.gr$type=='exon' & !grepl('pseudo', all.gr$source)]

  # Make intron GRanges object
  ############################
  intron.gr <- setdiff(gene.gr, exon.gr)
  intron.gr <- intron.gr[start(intron.gr) < end(intron.gr)]

  # Make promoter GRanges object
  ##############################
  promoter.gr <- promoters(gene.gr, upstream=prom_size,
                           downstream=0)

  # Make 3' UTR and 5' UTR GRanges objects
  ########################################
  # I had some difficulty with these. left out for now

  gr <- list()
  gr$gene <- gene.gr
  gr$exon <- exon.gr
  gr$intron <- intron.gr
  gr$promoter <- promoter.gr
#  gr$utr3 <- utr3.gr
#  gr$utr5 <- utr5.gr
  gr
}