# nathan dot lazar at gmail dot com

make_seqinfo <- function(len_file, species='') {
#########################################
# Read in lengths and make seqinfo object
#########################################
  lengths <- read.table(len_file, sep="\t")

  # Order
  lengths <- lengths[with(lengths, order(V1)), ]

  # Fix chr naming
  lengths$V1 <- sub('Chr', 'chr', lengths$V1)

  if(sum(grepl('chr', as.character(lengths$V1)))==0) {
    seqinf <- Seqinfo(paste0('chr', as.character(lengths[,1])),
                      lengths[,2], NA, species)
  } else {
    seqinf <- Seqinfo(as.character(lengths[,1]), lengths[,2], NA, species)
  }
  return(seqinf)
}
