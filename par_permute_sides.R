# nathan dot lazar at gmail dot com

par_permute_sides <- function(wdir, bindir, bp.lr.gr, all.bs, n=1000,
                              adjacent=T, end.exclude=1000) {
# This function utilizes parallel processing through HTCondor to
# find 1,000 sets of regions with the same size distribution of those
# in bp.lr.gr and make comparisons between the two sides of the breakpoints

# bp.lr.gr is a GRanges object of BP regions

# If adjacent=T then find regions of size (left + right) look at 
# differences between the two sides

# If adjacent=F then find randomly located regions and compare
# differences between them.

# Returns the mean absolute difference (MAD) in methylation, coverage and
# CpG counts for each pair of regions (adjacent or non)

# Random regions of the same size and number as in bp.lr.gr are chosen
# from the genome (excluding 1kb on the end of scaffolds).
# Methylation in these regions is calculated as the average
# of the methylation weighted by how many CpGs there are in each
# feature. This is equivalent to an average of all CpG methylation.

  bp.lr.mad <- get_mad(bp.lr.gr)

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
  } else {
    num <- nrow(bp.lr.mad)
    sizes <- bp.lr.mad$size
  }

  min.chr.size <- max(sizes) + 2*end.exclude

  # Make a list for sampling according to chrom lengths
  # only sample from chroms at least <min.chr.size>
  lengths <- seqlengths(bp.lr.gr)[seqlengths(bp.lr.gr) >=
                               min.chr.size]
  breaks <- cumsum(as.numeric(lengths-(end.exclude*2))) /
              sum(as.numeric(lengths-(end.exclude*2)))
  names(breaks) <- names(lengths)

  # Make directory to store output
  if(adjacent) {type <- 'adj'} else {type='disj'} 
  wdir <- paste0(wdir, type)
  if(!file.exists(wdir)) {
    dir.create(wdir)
    dir.create(paste0(wdir, '/logs')) 
  }

  # Save the objects all.bs, bp.lr.mad, breaks and lengths
  #  to a file that can be read by all workers
  save(all.bs, bp.lr.mad, breaks, lengths, file=paste0(wdir, '/par_permute_sides.dat'))

  adjacent <- as.numeric(adjacent)

  # Make condor submit script
  cores <- 16
  jobs <- ceiling(n/cores)
  make_per_submit(wdir, paste0(bindir, '/wrap_par_rand_sides.R'),
    c('$(dir)/par_permute_sides.dat', '1000', adjacent, '$$(Cpus)', paste0(wdir,'/par_rand_sides.$(Process)')),
    wdir, cores, '10 GB', '5 GB', jobs, 
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

  # Clean up temporary output files
  for(i in 1:jobs) system(paste0('rm ', wdir, '/permute.', as.character(i-1), '.txt'))

  # Read in files written by HTCondor and combine into one data.frame
  all.rand <- data.frame(group.mad.meth=numeric(0),
                         group.mad.cov=numeric(0),
                         group.mad.cpgs=numeric(0))
  all.reg.ad.meth <- c()
  all.reg.ad.cov <- c()
  all.reg.ad.cpgs <- c()
  for(i in 0:(jobs-1)) {
    file <- paste0(wdir, '/par_rand_sides.', i)
    load(file)                                   # loads the list 'rand'
    all.rand <- rbind(all.rand, rand$group.rand)
    all.reg.ad.meth <- c(all.reg.ad.meth, rand$region.ad.meth)
    all.reg.ad.cov <- c(all.reg.ad.cov, rand$region.ad.cov)
    all.reg.ad.cpgs <- c(all.reg.ad.cpgs, rand$region.ad.cpgs)
  }
  all.rand <- as.data.frame(apply(all.rand, 2, sort))
  all.reg.ad.meth <- sort(all.reg.ad.meth[!is.na(all.reg.ad.meth)])
  all.reg.ad.cov <- sort(all.reg.ad.cov[!is.na(all.reg.ad.cov)])
  all.reg.ad.cpgs <- sort(all.reg.ad.cpgs[!is.na(all.reg.ad.cpgs)])

  # Store results
  results <- list()
  results$all.rand <- all.rand
  results$all.reg.ad.meth <- all.reg.ad.meth
  results$all.reg.ad.cov <- all.reg.ad.cov
  results$all.reg.ad.cpgs <- all.reg.ad.cpgs
  results$bp.count <- length(unique(bp.lr.mad$s_name))

  # Get percentiles for region level permutation
  pers <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
  meth.ranks <- sapply(pers, function(x) round(x * length(all.reg.ad.meth)))
  meth.ranks[meth.ranks==0] <- 1
  cov.ranks <- sapply(pers, function(x) round(x * length(all.reg.ad.cov)))
  cov.ranks[cov.ranks==0] <- 1
  cpgs.ranks <- sapply(pers, function(x) round(x * length(all.reg.ad.cpgs)))
  cpgs.ranks[cpgs.ranks==0] <- 1
  results$region.percentiles <- data.frame(mad.meth=all.reg.ad.meth[meth.ranks],
                                           mad.cov=all.reg.ad.cov[cov.ranks],
                                           mad.cpgs=all.reg.ad.cpgs[cpgs.ranks])
  results$region.percentiles$pers <- pers

  # Get percentiles and p-values for group level permutation
  ranks <- sapply(pers, function(x) round(x * nrow(all.rand)))
  ranks[ranks==0] <- 1
  results$group.percentiles <- as.data.frame(apply(all.rand, 2, "[", ranks))
  names(results$group.percentiles) <- names(all.rand)
  results$group.percentiles$pers <- pers

  results$group.p_values <- list(meth=sum(all.rand$group.mad.meth < mean(bp.lr.mad$mad_meth, na.rm=T), na.rm=T)/nrow(all.rand),
    cov=sum(all.rand$group.mad.cov < mean(bp.lr.mad$mad_cov, na.rm=T), na.rm=T)/nrow(all.rand),
    cpgs=sum(all.rand$group.mad.cpgs < mean(bp.lr.mad$mad_cpgs, na.rm=T), na.rm=T)/nrow(all.rand))

  return(results)
}
