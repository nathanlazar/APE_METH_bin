# nathan dot lazar at gmail dot com

read_bp <- function(bp_file, seqinfo) {
########################################################
# Read in Breakpoint region data and make GRanges object
########################################################
  bps <- read.table(bp_file, sep="\t", header=F,
    blank.lines.skip=T, strip.white=T, fill=T)

  # Drop NA columns
  bps <- bps[,apply(bps, 2, function(x) sum(!is.na(x))) !=0]

  if(ncol(bps) > 7)
    bps <- bps[,1:7]
  if(ncol(bps)==7) {
    names(bps)=c('chr', 'start', 'end', 'BP_name', 'size', 'notes', 'class')
    bps$chr <- paste0('chr', bps$chr)
  } else {
    if(ncol(bps)==4)
      names(bps)=c('chr', 'start', 'end', 'BP_name')
    if(ncol(bps) %in% 5:6) {
      bps <- bps[,1:4]
      names(bps)=c('chr', 'start', 'end', 'BP_name')
    }
  }

  # Replace 2A and 2B with 2a and 2b
  bps$BP_name<- sub('A', 'a', bps$BP_name)
  bps$BP_name <- sub('B', 'b', bps$BP_name)

  bp.gr <- makeGRangesFromDataFrame(bps, keep.extra.columns=T)
  seqlevels(bp.gr) <- seqlevels(seqinfo)
  seqlengths(bp.gr) <- seqlengths(seqinfo)
  genome(bp.gr) <- genome(seqinfo)

  bp.gr <- trim(bp.gr)
  sort(bp.gr)
}

