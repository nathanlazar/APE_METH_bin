#nathan dot lazar at gmail dot com

# Combines two tallied and sorted CpG evidence directories

# Usage: Rscript combine_cpg.R <drive1> <drive2> <out_dir>

library(dplyr)
args <- commandArgs(TRUE)

files1 <- list.files(path=args[1])
files2 <- list.files(path=args[2])
files <- unique(c(files1, files2))

empty.df <- data.frame(ref=character(0), off=numeric(0), strand=character(0),
                           Mstr=character(0), Mcy=character(0), Ustr=character(0),
                           Ucy=numeric(0), filt_cycle=numeric(0),
                           filt_readlen=numeric(0), filt_allele=numeric(0),
                           filt_mapq=numeric(0), filt_baseq=numeric(0))

lapply(files, function(x) {
  if(file.exists(paste0(args[1], "/" ,x)))  {
     f1 <- read.table(paste0(args[1], "/" ,x), 
                      header=T, sep='\t', stringsAsFactors=F, 
                      check.names=F, comment.char="", quote="")
     f1[is.na(f1)] <- ''
  } else {f1 <- empty.df}
  if(file.exists(paste0(args[2], "/" ,x)))  {
     f2 <- read.table(paste0(args[2], "/" ,x), 
                      header=T, sep='\t', stringsAsFactors=F, 
                      check.names=F, comment.char="", quote="")
     f2[is.na(f2)] <- ''
  } else {f2 <- empty.df}
  f <- merge(f1,f2, all=T)
  if(nrow(f) > 0) {
    f.sum <- f %>% 
      group_by(ref, off, strand) %>%
      summarize(Mstr=paste(Mstr, collapse=""), 
                Mcy=sum(Mcy),
                Ustr=paste(Ustr, collapse=""), 
                Ucy=sum(Ucy),
                filt_cycle=sum(filt_cycle),
                filt_readlen=sum(filt_readlen),
                filt_allele=sum(filt_allele),
                filt_mapq=sum(filt_mapq),
                filt_baseq=sum(filt_baseq))
    write.table(f.sum, paste0(args[3], "/", x), sep="\t", quote=F, row.names=F, col.names=T)
  } else {
    write.table(f, paste0(args[3], "/", x), sep="\t", quote=F, row.names=F, col.names=T)
  }
})
