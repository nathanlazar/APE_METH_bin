# nathan dot lazar at gmail dot com

# Convert genes in the UCSC format below to a GRanges list object
#bin    name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount 
#exonStarts    exonEnds       score   name2   cdsStartStat    cdsEndStat      exonFrames

library(IRanges)
library(GenomicRanges)

UCSC2GRanges <- function(myfile='my.txt', seqinfo, prom_size=1000, filter=F) {
  UCSC <- read.delim(myfile, sep='\t', header=T, stringsAsFactors=F)

  # Add chr to chrom names if necessary
  if(sum(!grepl('chr', UCSC$chrom))>0) {
    UCSC$chrom[!grepl('chr', UCSC$chrom)] <- paste0('chr', UCSC$chrom[!grepl('chr', UCSC$chrom)])
  }

  if(filter) {
  # Remove all but the longest transcript for each overlapping set of genes
  # The original gibbon RefSeq gene set had 296,965 genes,
  # filtering to keep only largest reduces this to 19,271
    UCSC <- filter_genes(UCSC)
  }

  gene.gr <- GRanges(seqnames=UCSC$chrom,
                     ranges=IRanges(UCSC$txStart, UCSC$txEnd),
                     strand=UCSC$strand,
                     name=UCSC$name,
                     symbol=UCSC$name2)

  utr5.gr <- GRanges(seqnames=UCSC$chrom,
                     ranges=IRanges(UCSC$txStart, UCSC$cdsStart),
                     strand=UCSC$strand,
                     name=UCSC$name,
                     symbol=UCSC$name2)

  utr3.gr <- GRanges(seqnames=UCSC$chrom,
                     ranges=IRanges(UCSC$cdsEnd, UCSC$txEnd),
                     strand=UCSC$strand,
                     name=UCSC$name,
                     symbol=UCSC$name2)

  exon.gr <- GRanges(seqnames=rep(UCSC$chrom, UCSC$exonCount),
                     ranges=IRanges(as.numeric(unlist(strsplit(UCSC$exonStarts, split=','))),
                                    as.numeric(unlist(strsplit(UCSC$exonEnds, split=',')))),
                     strand=rep(UCSC$strand, UCSC$exonCount),
                     name=rep(UCSC$name, UCSC$exonCount),
                     symbol=rep(UCSC$name2, UCSC$exonCount),
                     exon_number=unlist(lapply(UCSC$exonCount, function(x) 1:x)))

  intron.gr <- GRanges(seqnames=rep(UCSC$chrom, UCSC$exonCount-1),
                       ranges=IRanges(as.numeric(unlist(sapply(strsplit(UCSC$exonEnds, split=','), function(x) x[-length(x)])))+1,
                                      as.numeric(unlist(sapply(strsplit(UCSC$exonStarts, split=','), '[', -1)))-1),
                       strand=rep(UCSC$strand, UCSC$exonCount-1),
                       name=rep(UCSC$name, UCSC$exonCount-1),
                       symbol=rep(UCSC$name2, UCSC$exonCount-1),
                       intron_number=unlist(lapply((UCSC$exonCount[UCSC$exonCount>1]-1), function(x) 1:x)))

  promoter.gr <- promoters(gene.gr, upstream=prom_size,  downstream=0)

  gr <- list()
  gr$promoter <- promoter.gr
  gr$gene <- gene.gr
  gr$utr5 <- utr5.gr
  gr$exon <- exon.gr
  gr$intron <- intron.gr
  gr$utr3 <- utr3.gr

  # Remove genes that aren't on chromosomes in provided seqinfo
  for(i in 1:length(gr)) {
    gr[[i]] <- gr[[i]][seqnames(gr[[i]]) %in% seqlevels(seqinfo)]
    seqlevels(gr[[i]], force=T) <- seqlevels(seqinfo)
    seqlengths(gr[[i]]) <- seqlengths(seqinfo)
  }

  # Remove ranges with zero base pairs
  gr$utr5 <- gr$utr5[width(gr$utr5) != 1]
  gr$utr3 <- gr$utr3[width(gr$utr3) != 1]
  gr$intron <- gr$intron[width(gr$intron) != 1]
  gr$promoter <- gr$promoter[width(gr$promoter) != 1]

  # Trim any ranges that are out of range
  gr$promoter <- trim(gr$promoter)
  gr$utr5 <- trim(gr$utr5)
  gr$gene <- trim(gr$gene)
  gr$exon <- trim(gr$exon)
  gr$intron <- trim(gr$intron)
  gr$utr3 <- trim(gr$utr3)

  gr
}