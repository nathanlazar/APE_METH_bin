# nathan dot lazar at gmail dot com

par_rand <- function(all.bs, feat.gr, breaks, sizes, lengths, end.exclude, type, reps) {
# Finds random regions in a given genome,
# overlaps with features (if given) and gets the
# mean methylation in features, mean coverage of 
# CpGs in features, total count of CpGs in features
# and percentage of the random regions covered by
# the features

  num <- length(sizes)
  tot.size <- sum(sizes)
  nn <- reps*num

  rand_regions <- data.frame(chr=rep('', nn),
                             start=rep(0, nn),
                             end=rep(0, nn),
                             stringsAsFactors=F)

  for(i in 1:reps) {
    # Choose a set of regions randomly with probability
    # proportional to chromosome lengths
    regions <- mclapply(sizes, get_region, breaks, lengths,
                      end.exclude)
    for(j in 1:num) {
      rand_regions[(i-1)*num + j,] <- regions[[j]]
    }
  }

  rand_regions$start <- as.numeric(rand_regions$start)
  rand_regions$end <- as.numeric(rand_regions$end)
  rand.gr <- makeGRangesFromDataFrame(rand_regions)

  rand <- data.frame(mean.meth=rep(0,reps),
                     mean.cov=rep(0,reps),
                     cpgs.per.kb=rep(0,reps),
                     per.cov=rep(0,reps))

  if(type=='all') {
    # Get methylation and coverage for all random regions at once
    # to minimize overhead
    meth <- mcgetMeth(all.bs, regions=rand.gr, type='raw', what='perBase')
    cov <- mcgetCoverage(all.bs, regions=rand.gr, what='perBase')

    # Get means in groups by the number of regions in each group
    rand$mean.meth <- foreach(i=1:reps) %dopar%
      mean(unlist(meth[((i-1)*num+1):(i*num)]), na.rm=T)
    rand$mean.cov <- foreach(i=1:reps) %dopar% 
      mean(unlist(cov[((i-1)*num+1):(i*num)]), na.rm=T)
    rand$cpgs.per.kb <- foreach(i=1:reps) %dopar% {
      length(unlist(cov[((i-1)*num+1):(i*num)]))/tot.size*1000}
    
  } else {
    
    feat.in.rand <- foreach(i=1:reps) %dopar%
      subsetByOverlaps(feat.gr, rand.gr[((i-1)*num+1):(i*num)])
    rand$mean.meth <- foreach(i=1:reps) %dopar%
      weighted.mean(feat.in.rand[[i]]$meth, feat.in.rand[[i]]$cpgs, na.rm=T)
    rand$mean.cov <- foreach(i=1:reps) %dopar%
      weighted.mean(feat.in.rand[[i]]$cov, feat.in.rand[[i]]$cpgs, na.rm=T)
    rand$cpgs.per.kb <- foreach(i=1:reps) %dopar% {
      sum(feat.in.rand[[i]]$cpgs)/sum(width(feat.in.rand[[i]]))*1000}
    rand$per.cov <- foreach(i=1:reps) %dopar% {
      sum(width(intersect(rand.gr[((i-1)*num+1):(i*num)], feat.in.rand[[i]],
                          ignore.strand=T)))/tot.size}
  }
  rand
}
