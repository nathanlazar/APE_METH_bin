# nathan dot lazar at gmail dot com

make_lr_bp <- function(bp.gr, span) {
##############################################################
# Make GRanges object of <span> regions on each side of breaks
##############################################################

  bp.classIIb.gr <- bp.gr[bp.gr$class=="II-B"]
  bp.classIIb.gr$side <- 'neither'
  bp.gr <- bp.gr[!bp.gr$class=="II-B"]

  bp.left.gr <- bp.gr[!grepl('Start', bp.gr$notes)]
  end(bp.left.gr) <- start(bp.left.gr)
  start(bp.left.gr) <- start(bp.left.gr) - span
  bp.left.gr$size <- end(bp.left.gr) - start(bp.left.gr)
  bp.left.gr$side <- 'left'

  bp.right.gr <- bp.gr[!grepl('End', bp.gr$notes)]
  start(bp.right.gr) <- end(bp.right.gr)
  end(bp.right.gr) <- end(bp.right.gr) + span
  bp.right.gr$size <- end(bp.right.gr) - start(bp.right.gr)
  bp.right.gr$side <- 'right'

  bp.lr.gr <- c(bp.left.gr, bp.right.gr, bp.classIIb.gr)
  bp.lr.gr <- sort(bp.lr.gr)

  # Trim ranges to be within chromosomes
  bp.lr.gr <- trim(bp.lr.gr)

  # Combine ranges that overlap
  tmp <- reduce(bp.lr.gr, with.mapping=T, min.gapwidth=0)

  # Remove ranges with zero width
  tmp <- tmp[width(tmp) > 1]

  # Keep metadata for first of overlapping regions
  rows <- unlist(lapply(tmp$mapping, '[', 1))
  mcols(tmp) <- cbind(mcols(tmp), mcols(bp.lr.gr)[rows,])

  # If regions were combined, concatenate metadata BP_name,
  # notes, class, side with a '-'
  cats <- which(lapply(tmp$mapping, length) > 1)
  if(length(cats) > 0) {
    tmp.names <- as.character(tmp$BP_name)
    tmp.notes <- as.character(tmp$notes)
    tmp.side <- as.character(tmp$side)
    for(rng in cats) {
      tmp.names[rng] <-  paste(as.character(
        bp.lr.gr$BP_name[unlist(tmp[rng]$mapping)]), collapse='-')
      tmp.notes[rng] <-  paste(as.character(
        bp.lr.gr$notes[unlist(tmp[rng]$mapping)]), collapse='-')
      tmp.side[rng] <-  paste(as.character(
        bp.lr.gr$side[unlist(tmp[rng]$mapping)]), collapse='-')
    }
    tmp$BP_name <- as.factor(tmp.names)
    tmp$notes <- as.factor(tmp.notes)
    tmp$side <- as.factor(tmp.side)
  }
  tmp$size <- width(tmp)
  return(tmp)
}
