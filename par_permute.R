# nathan dot lazar at gmail dot com

par_permute <- function(wdir, bindir, feat.gr, bp.lr.gr, all.bs, n=1000,
                        type=c('gene', 'exon', 'intron', 'promoter', '3UTR',
                               '5UTR', 'CpGisl', 'CpGshore', 'repeat', 'LINE',
                               'SINE', 'DNA', 'LTR', 'SINE', 'Alu', 'AluS',
                               'AluJ', 'AluY', 'MIR', 'all'),
                        min.chr.size=12000, end.exclude=1000) {
# This function performs the same function as permute, but utilizes parallel 
# processing through HTCondor
#
# feat.gr is a GRanges object of features of the given type with
# metadata columns meth, cpgs and cov telling the average methylation,
# the number of cpgs in each range and the average coverage of those cpgs.
# bp.lr.gr is a GRanges object of BP regions
#
# Random regions of the same size and number as in bp.lr.gr are chosen
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
  lengths <- seqlengths(bp.lr.gr)[seqlengths(bp.lr.gr) >=
                               min.chr.size]
  breaks <- cumsum(as.numeric(lengths-(end.exclude*2))) /
              sum(as.numeric(lengths-(end.exclude*2)))
  names(breaks) <- names(lengths)
  sizes <- bp.lr.gr$size

  # Make directory to store output
  wdir <- paste0(wdir, type)
  if(!file.exists(wdir)) {
    dir.create(wdir)
    dir.create(paste0(wdir, '/logs')) 
  }

  # Save the objects feat.gr, bp.lr.gr and all.bs, breaks and lengths
  #  to a file that can be read by all workers
  save(feat.gr, bp.lr.gr, all.bs, breaks, sizes, lengths, 
       file=paste0(wdir, '/par_permute.dat'))

  # Make condor submit script
  cores <- 16
  jobs <- 63
  make_per_submit(wdir, paste0(bindir, '/wrap_par_rand.R'),
    c('$(dir)/par_permute.dat', type, '1000', '$$(Cpus)'),
    wdir, cores, '2 GB', '2 GB', jobs, 
    paste0(wdir, '/condor.submit'))

  # Run condor script
  system(paste0('condor_submit ', wdir, '/condor.submit'))

  # If writing over output files, wait a minute so condor can create empty files
  if(file.exists(paste0(wdir, '/permute.', as.character(jobs-1), '.txt')))
    Sys.sleep(60) 

  # Wait for these to be done (there's probably a better way to do this)
  written.files <- rep(0,jobs)
  while(sum(written.files) < jobs) {
    for(i in 1:jobs) {
      if(written.files[i] < 1) {
        f <- paste0(wdir, '/permute.', as.character(i-1), '.txt')
        if (length(readLines(f)) > 0)
          written.files[i] <- 1
      }
    }
    Sys.sleep(5) #check every 5 seconds
  }

  # Read in files written by HTCondor and combine into one data.frame
  rand <- data.frame()
  for(i in 0:(jobs-1)) {
    file <- paste0(wdir, '/permute.', i, '.txt')
    results <- read.table(file, header=T)
    rand <- rbind(rand, data.frame(results))
  }

  tot.size <- sum(sizes)
  results <- list()

  if(type=='all') {

    # See how many of these permutations have methylation as low as
    # the breakpoint regions
    results$bp.count <- length(bp.lr.gr)
    results$bp.w.av.meth <- weighted.mean(bp.lr.gr$meth, bp.lr.gr$cpgs)
    results$bp.w.av.cov <- weighted.mean(bp.lr.gr$cov, bp.lr.gr$cpgs)
    results$bp.cpgs.per.kb <- sum(bp.lr.gr$cpgs)/tot.size*1000

    ######Report p-values############
    n <- length(!is.nan(rand$mean.cov))
    cat(type, 'Permutation p-values (random < observed):\n')
    results$meth.p <- sum(rand$mean.meth < results$bp.w.av.meth, na.rm=T)/n
    cat('Methylation:\t',  results$meth.p, '\n')
    results$cov.p <- sum(rand$mean.cov < results$bp.w.av.cov, na.rm=T)/n
    cat('Coverage: \t', results$cov.p, '\n')
    results$cpg.p <-  sum(rand$cpgs.per.kb < results$bp.cpgs.per.kb)/n
    cat('CpGs per Kb: \t', results$cpg.p, '\n')

  } else {

    feat.in.bp <- subsetByOverlaps(feat.gr, bp.lr.gr)

    results$bp.count <- length(feat.in.bp)
    results$bp.w.av.meth <- weighted.mean(feat.in.bp$meth, feat.in.bp$cpgs,
                                          na.rm=T)
    results$bp.w.av.cov <- weighted.mean(feat.in.bp$cov, feat.in.bp$cpgs)
    results$bp.cpgs <- sum(feat.in.bp$cpgs)
    results$bp.cpgs.per.kb <- results$bp.cpgs/sum(width(feat.in.bp))*1000
    overlap <- intersect(bp.lr.gr, feat.in.bp, ignore.strand=T)
    results$bp.per.cov <- sum(width(overlap))/tot.size

    ######Report p-values############
    cat('For ', results$bp.count, 'features of type ', type, 'the permutation p-values are: (random < observed)\n')
    results$meth.n <- sum(!is.nan(rand$mean.meth))
    results$meth.p <- sum(rand$mean.meth < 
                          results$bp.w.av.meth, na.rm=T)/results$meth.n
    cat('Methylation ( n=', results$meth.n, ') :\t', results$meth.p, '\n')

    results$cov.n <- sum(!is.nan(rand$mean.cov))
    results$cov.p <- sum(rand$mean.cov < 
                         results$bp.w.av.cov, na.rm=T)/results$cov.n
    cat('Coverage ( n=', results$cov.n, ') :\t', results$cov.p, '\n')

    results$cpg.n <- sum(!is.nan(rand$cpgs.per.kb))
    results$cpg.p <- sum(rand$cpgs.per.kb < 
                         results$bp.cpgs.per.kb, na.rm=T)/results$cpg.n
    cat('CpGs per Kb ( n=', results$cpg.n, ') :\t', results$cpg.p, '\n')

    results$per.n <- sum(!is.nan(rand$per.cov))
    results$per.p <- sum(rand$per.cov <= 
                         results$bp.per.cov, na.rm=T)/results$per.n
    cat('% region covered: ( n=', results$per.n, ') :\t', results$per.p, '\n')
  }

  # Return list containing p-values and dataframe of permutation values
  results$rand <- rand
  results
}
