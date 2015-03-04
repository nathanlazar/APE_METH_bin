# nathan dot lazar at gmail dot com

make_seqinfo <- function(len_file, species='') {
#########################################
# Read in lengths and make seqinfo object
#########################################
  lengths <- read.table(len_file, sep="\t")
  #order
  lengths <- lengths[with(lengths, order(V1)), ]
  Seqinfo(paste0('chr', as.character(lengths[,1])),
    lengths[,2], NA, species)
}
