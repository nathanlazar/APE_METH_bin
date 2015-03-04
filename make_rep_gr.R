# nathan dot lazar at gmail dot com

make_rep_gr <- function(rep_file, seqinfo) {
##################################################
# Reads repeat data from repeat masker output file
# <rep_file> and makes a GRanges object
##################################################
  reps <- read.table(rep_file, skip=3)
  names(reps) <- c('SW_score', 'perc_div', 'per_del',
                   'per_ins', 'chr', 'start',
                   'end', 'rep_left', 'matching_rep',
                   'rep_class', 'rep_family', 'rep_start',
                   'rep_end', 'rep_left', 'ID')
  #Remove 5 repeats on "random" chromosomes
  reps <- reps[grepl('chr', reps$chr)==F, ]

  reps$chr <- paste0('chr', reps$chr)

  reps.gr <- makeGRangesFromDataFrame(reps,
               keep.extra.columns=T)
  seqlevels(reps.gr) <- seqlevels(seqinfo)
  seqlengths(reps.gr) <- seqlengths(seqinfo)

  sort(reps.gr)
}

