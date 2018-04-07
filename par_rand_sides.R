# nathan dot lazar at gmail dot com

par_rand_sides <- function(all.bs, bp.lr.mad, breaks, lengths, end.exclude, 
                           adjacent, reps=1000, min_cov=4) {
# Finds random regions in a given genome and gets the
# mean absolute differnce in methylation coverage and cpg count 

# If adjacent=T then the regions are next to each other, 
# If adjacent=F we compare random pairs of regions
  
  # Subset all.bs by coverage >= min_cov
  cov_min.bs <- all.bs[getCoverage(all.bs) >= min_cov]

  if('class' %in% names(bp.lr.mad)) {
    bp.lr.mad.size <- bp.lr.mad %>%
      select(s_name, size, side, class) %>%
      dcast(s_name + class ~ side, value.var='size', fun.aggregate=sum) 
  } else {
    bp.lr.mad.size <- bp.lr.mad %>%
      select(s_name, size, side) %>%
      dcast(s_name ~ side, value.var='size', fun.aggregate=sum) 
  }

  if(adjacent) {
    num <- nrow(bp.lr.mad.size)
    sizes <- bp.lr.mad.size$L + bp.lr.mad.size$R
    nn <- reps*num
  } else {
    num <- nrow(bp.lr.mad)
    sizes <- bp.lr.mad$size
    nn <- reps*num
  }

  rand_regions <- data.frame(chr=rep('', nn),
                             start=rep(0, nn),
                             end=rep(0, nn),
                             stringsAsFactors=F)

  for(i in 1:reps) {
    # Choose a set of regions randomly with probability
    # proportional to chromosome lengths
    regions <- mclapply(sizes, get_region, breaks, lengths, end.exclude)
    for(j in 1:num) {
      rand_regions[(i-1)*num + j,] <- regions[[j]]
    }
  }

  if(adjacent) {
    rand_regions.l <- rand_regions
    rand_regions.l$start <- as.numeric(rand_regions$start)
    rand_regions.l$end <- rand_regions.l$start + bp.lr.mad.size$L
    rand_regions.r <- rand_regions
    rand_regions.r$start <- rand_regions.l$end + 1
    rand_regions.r$end <- as.numeric(rand_regions$end)
    rand.l.gr <- makeGRangesFromDataFrame(rand_regions.l)
    rand.r.gr <- makeGRangesFromDataFrame(rand_regions.r)
  } else {
    rand_regions$start <- as.numeric(rand_regions$start)
    rand_regions$end <- as.numeric(rand_regions$end)
    rand.l.gr <- makeGRangesFromDataFrame(rand_regions[bp.lr.mad$side=='L',])
    rand.r.gr <- makeGRangesFromDataFrame(rand_regions[bp.lr.mad$side=='R',])
    num <- length(rand.r.gr)/reps
  }

  # Get methylation and coverage for all random regions at once
  meth.l <- mcgetMeth(cov_min.bs, regions=rand.l.gr, type='raw', what='perBase')
  meth.r <- mcgetMeth(cov_min.bs, regions=rand.r.gr, type='raw', what='perBase')

  cov.l <- mcgetCoverage(all.bs, regions=rand.l.gr, what='perBase')
  cov.r <- mcgetCoverage(all.bs, regions=rand.r.gr, what='perBase')

  # Get distribution of the absolute difference in methylation for all regions
  region.ad.meth <- abs(sapply(meth.l, mean, na.rm=T) -
                        sapply(meth.r, mean, na.rm=T))
  region.ad.cov <- abs(sapply(cov.l, mean, na.rm=T) -
                        sapply(cov.r, mean, na.rm=T))
  region.ad.cpgs <- abs(sapply(cov.l, length) -
                        sapply(cov.r, length))
    
  # Get group level MAD distributions
  # Take the mean methylation of each region, subtract two sides,
  # take absolute value, take the mean of each group

  group.mad.meth <- foreach(i=1:reps) %dopar% {
    mean(abs(sapply(meth.l[((i-1)*num+1):(i*num)], mean, na.rm=T) -
             sapply(meth.r[((i-1)*num+1):(i*num)], mean, na.rm=T)),na.rm=T)
  }

  group.mad.cov <- foreach(i=1:reps) %dopar% {
    mean(abs(sapply(cov.l[((i-1)*num+1):(i*num)], mean, na.rm=T) - 
             sapply(cov.r[((i-1)*num+1):(i*num)], mean, na.rm=T)),na.rm=T)
  }

  group.mad.cpg <- foreach(i=1:reps) %dopar% {
    mean(abs(sapply(cov.l[((i-1)*num+1):(i*num)], length) - 
             sapply(cov.r[((i-1)*num+1):(i*num)], length)))
  }

  rand <- data.frame(group.mad.meth=unlist(group.mad.meth),
                     group.mad.cov=unlist(group.mad.cov),
                     group.mad.cpgs=unlist(group.mad.cpg))

  return(list(group.rand=rand, region.ad.meth=region.ad.meth, 
              region.ad.cov=region.ad.cov, region.ad.cpgs=region.ad.cpgs))
}
