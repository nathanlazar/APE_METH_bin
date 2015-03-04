# nathan dot lazar at gmail dot com

read_bp <- function(bp_file, seqinfo) {
########################################################
# Read in Breakpoint region data and make GRanges object
########################################################
  bps <- read.table(bp_file, sep="\t", header=T,
    blank.lines.skip=T, strip.white=T)
  if(ncol(bps) > 7)
    bps <- bps[,1:7]
  if(ncol(bps)==5)
    names(bps)=c('BP_name', 'chr', 'start', 'end', 'size')
  if(ncol(bps)==7)
    names(bps)=c('BP_name', 'chr', 'start', 'end', 'size', 'notes', 'class')
  bps$chr <- paste0('chr', bps$chr)

  bp.gr <- makeGRangesFromDataFrame(bps, keep.extra.columns=T)
  seqlevels(bp.gr) <- seqlevels(seqinfo)
  seqlengths(bp.gr) <- seqlengths(seqinfo)
  genome(bp.gr) <- genome(seqinfo)

  sort(bp.gr)
}

