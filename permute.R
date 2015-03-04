# nathan dot lazar at gmail dot com

permute <- function(feat.gr, bp.gr, all.bs, n=1000,
                    type=c('gene', 'exon', 'intron', 'promoter', '3UTR',
                           '5UTR', 'CpGisl', 'CpGshore', 'repeat', 'LINE',
                           'SINE', 'DNA', 'LTR', 'SINE', 'Alu', 'AluS',
                           'AluJ', 'AluY', 'MIR', 'all'),
                    min.chr.size=12000, end.exclude=1000) {
# feat.gr is a GRanges object of features of the given type with
# metadata columns meth, cpgs and cov telling the average methylation,
# the number of cpgs in each range and the average coverage of those cpgs.
# bp.gr is a GRanges object of BP regions
#
# Random regions of the same size and number as in bp.gr are chosen
# from the genome (excluding 1kb on the end of scaffolds).
# Methylation of features in these regions is calculated as the average
# of the feature methylation weighted by how many CpGs there are in each
# feature. This is equivalent to an average of all CpG methylation for
# CpGs in the given set of features.
#
# The percentage of the regions covered by the type of feature is also
# recorded and compared to the BP regions.
#
# If type=='all' then all CpGs in the random regions are measured
# Prints out p-values and returns dataframe of mean methylation,
# mean coverage and area covered (if type is a feature)


  # Make a list for sampling according to chrom lengths
  # only sample from chroms at least <min.chr.size>
  lengths <- seqlengths(bp.gr)[seqlengths(bp.gr) >=
                               min.chr.size]
  breaks <- cumsum(as.numeric(lengths-(end.exclude*2))) /
              sum(as.numeric(lengths-(end.exclude*2)))
  names(breaks) <- names(lengths)

  sizes <- bp.gr$size
  num <- length(sizes)

  nn <- n*num
  rand_regions <- data.frame(chr=rep('', nn),
                             start=rep(0, nn),
                             end=rep(0, nn),
                             stringsAsFactors=F)

  for(i in 1:n) {
    # Choose a set of regions randomly with probability
    # proportional to chromosome lengths
    regions <- lapply(sizes, get_region, breaks, lengths, 
                      end.exclude)
    for(j in 1:num) {
      rand_regions[(i-1)*num + j,] <- regions[[j]]
    }
  }

  rand_regions$start <- as.numeric(rand_regions$start)
  rand_regions$end <- as.numeric(rand_regions$end)
  rand.gr <- makeGRangesFromDataFrame(rand_regions)

  rand <- data.frame(mean.meth=rep(0,n),
                     mean.cov=rep(0,n),
                     tot.cpgs=rep(0,n),
                     per.cov=rep(0,n))


  if(type=='all') {
    meth <- mcgetMeth(all.bs, regions=rand.gr, type='raw', what='perBase')
    cov <- mcgetCoverage(all.bs, regions=rand.gr, what='perBase')

    # Get means in groups by the number of regions in each group
    for(i in 1:n) {
      rand$mean.meth[i] <- mean(unlist(meth[((i-1)*num+1):(i*num)]), na.rm=T)
      rand$mean.cov[i] <- mean(unlist(cov[((i-1)*num+1):(i*num)]), na.rm=T)
      rand$tot.cpgs[i] <- length(unlist(cov[((i-1)*num+1):(i*num)]))
    }

    # See how many of these permutations have methylation as low as
    # the breakpoint regions
    bp.w.av.meth <- weighted.mean(bp.gr$meth, bp.gr$cpgs)
    bp.w.av.cov <- weighted.mean(bp.gr$cov, bp.gr$cpgs)

    ######Report p-values############
    n <- length(!is.na(rand$mean.cov))
    cat('Permutation p-values (random < observed):\n')
    cat('Methylation:\t',  sum(rand$mean.meth < bp.w.av.meth, na.rm=T)/n, '\n')
    cat('Coverage: \t', sum(rand$mean.cov < bp.w.av.cov, na.rm=T)/n, '\n')
    cat('CpG count: \t', sum(rand$tot.cpgs < sum(bp.gr$cpgs))/n, '\n')

  } else {

    tot.size <- sum(sizes)
    feat.in.bp <- subsetByOverlaps(feat.gr, bp.gr)

    bp.w.av.meth <- weighted.mean(feat.in.bp$meth, feat.in.bp$cpgs,
                                  na.rm=T)
    bp.w.av.cov <- weighted.mean(feat.in.bp$cov, feat.in.bp$cpgs)
    bp.cpgs <- sum(feat.in.bp$cpgs)
    overlap <- intersect(bp.gr, feat.in.bp, ignore.strand=T)
    bp.per.cov <- sum(width(overlap))/tot.size

    for(i in 1:n) {
      feat.in.rand <- subsetByOverlaps(feat.gr, rand.gr[((i-1)*num+1):(i*num)])

      rand$mean.meth[i] <- weighted.mean(feat.in.rand$meth, feat.in.rand$cpgs,
      			                 na.rm=T)
      rand$mean.cov[i] <- weighted.mean(feat.in.rand$cov, feat.in.rand$cpgs,
      		                        na.rm=T)
      rand$tot.cpgs[i] <- sum(feat.in.rand$cpgs)

      #get percentage of regions covered by features
      rand.lap <- intersect(rand.gr[((i-1)*num+1):(i*num)], feat.in.rand,
                            ignore.strand=T)
      rand$per.cov[i] <- sum(width(rand.lap))/tot.size
    }

    ######Report p-values############
    n <- length(!is.na(rand$mean.cov))
    cat('Permutation p-values (random < observed):\n')
    cat('Methylation:\t',  sum(rand$mean.meth < bp.w.av.meth, na.rm=T)/n, '\n')
    cat('Coverage: \t', sum(rand$mean.cov < bp.w.av.cov, na.rm=T)/n, '\n')
    cat('CpG count: \t', sum(rand$tot.cpgs < sum(bp.gr$cpgs, na.rm=T))/n, '\n')
    cat('Percent region covered: \t', sum(rand$per.cov < bp.per.cov, na.rm=T)/n, '\n')
  }
  ####return dataframe of permutation values#########
  rand
}
