# nathan dot lazar at gmail dot com

get_region <- function(size, breaks, lengths, end.exclude) {
#############################################
# Choose a random region from the genome
# Chromsomes are chosen by length and regions
# exclude 1kb at the start and end of each
# chromosome
#############################################
  idx <- min(which(runif(1) < breaks))
  #Find random start (excluding first & last end.exclude)
  len <- lengths[idx]
  start <- round(runif(1, end.exclude, len-size-end.exclude))
  end <- start + size
  c(names(lengths)[idx], start, end)
}
