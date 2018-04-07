# Function to filter a gene data frame (from UCSC table browser) to keep only the longest gene 
# amongst overlapping genes/isoforms

# nathan dot lazar at gmail dot com

library(dplyr)

filter_genes <- function(genes.df) {

  # Sort genes.df
  genes.df <- genes.df[order(genes.df$chrom, genes.df$txStart, genes.df$txEnd),]

  # Add in length
  genes.df$len <- genes.df$txEnd-genes.df$txStart

  # Vector to indicate whether to keep or remove each gene in the list
  keep <- rep(1, nrow(genes.df))

  # Get the first row
  best <- genes.df[1,]
  best.ind <- 1

  # Loop over rows of the data frame deciding whether to keep each
  for(i in 2:nrow(genes.df)) {
    current <- genes.df[i,]

    # If on a new chrom set new gene to 'best'
    if(best$chrom != current$chrom) {
      best <- current
      best.ind <- i
      next
    }

    if(current$txStart <= best$txEnd) {  # overlaps with best
      if(current$len > best$len) {       # longer than best
        keep[best.ind] <- 0  # discard best
        best <- current      # make current the new best
        best.ind <- i        # update the index of the best
      } else {
        keep[i] <- 0         # New isoform isn't longer
      }
    } else {                 # Doesn't overlap, i.e. new gene
      best <- current
      best.ind <- i
    }
  }
  filt.df <- data.frame(genes.df[keep==1,])
  return(filt.df)
}