# nathan dot lazar at gmail dot com

make_rep_gr <- function(rep_file, seqinfo) {
##################################################
# Reads repeat data from repeat masker output file
# <rep_file> and makes a GRanges object
##################################################

  peek <- readLines(rep_file, 1)
  if(grepl('perc', peek)) {
    reps <- read.table(rep_file, skip=3, stringsAsFactors=F)
    names(reps) <- c('SW_score', 'perc_div', 'per_del',
                     'per_ins', 'chr', 'start',
                     'end', 'rep_left', 'matching_rep',
                     'rep_class', 'rep_family', 'rep_start',
                     'rep_end', 'rep_left', 'ID')
  } else {
    peek <- strsplit(peek, split='\t')
    if(length(peek[[1]])==17) {
      reps <- read.table(rep_file, stringsAsFactors=F, sep='\t')
      names(reps) <- c('dunno1', 'dunno2', 'dunno3', 'dunno4', 'dunno5',
                       'chr', 'start', 'end', 'dunno6', 'strand',
                       'rep_class', 'rep_family', 'rep_name', 'rep_start',
                       'rep_end', 'dunno7', 'dunno8')
    } else {
      if(length(peek[[1]]==9)) {
        reps <- read.table(rep_file, stringsAsFactors=F)
        names(reps) <- c('chr', 'rep_class', 'rep_family', 'start', 'end', 
                         'dunno1', 'strand', 'dunno2', 'notes')
        reps$rep_class <- sub('RPMSK_', '', reps$rep_class)
      }
    }
  }

  if(grepl('human', rep_file)) { 
    #Remove repeats on 'random', 'hap', 'Un'
    reps <- reps[!grepl('random', reps$chr),]
    reps <- reps[!grepl('hap', reps$chr),] 
    reps <- reps[!grepl('Un', reps$chr),]
    reps <- reps[!grepl('chrM', reps$chr),]
  }

  # Add "chr" if necessary
  reps$chr[!grepl('chr', reps$chr)] <- paste0('chr', reps$chr[!grepl('chr', reps$chr)])

  # Change to 1 based coordinates
  reps$start <- reps$start+1

  reps.gr <- makeGRangesFromDataFrame(reps,
               keep.extra.columns=T)
  reps.gr <- reps.gr[seqnames(reps.gr) %in% seqlevels(seqinfo)]
  seqlevels(reps.gr) <- seqlevels(seqinfo)
  seqlengths(reps.gr) <- seqlengths(seqinfo)

  sort(reps.gr)
}

