# nathan dot lazar at gmail dot com

read_cpg_island <- function(cpg_island_file, seqinfo) {
################################################
# Reads in cpg_island gff file and makes GRanges
# object
################################################

  cpg_isl.df <- read.table(cpg_island_file, header=F,
                           stringsAsFactors=F)
  cpg_isl.df <- cpg_isl.df[,c(1,4,5,9)]
  names(cpg_isl.df) <- c('chr', 'start', 'end', 'info')

  cpg_isl.gr <- makeGRangesFromDataFrame(cpg_isl.df,
                  keep.extra.columns=T)

  seqlevels(cpg_isl.gr) <- seqlevels(seqinfo)
  seqlengths(cpg_isl.gr) <- seqlengths(seqinfo)

  cpg_isl.gr
}
