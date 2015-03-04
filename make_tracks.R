# nathan dot lazar at gmail dot com

# Make bed files to view coverage and methylation on UCSC (or other) browser

make_tracks <- function(all.bs, out_dir='', meth_file='meth.bedgraph', 
                        cov_file='cov.bedgraph') {

  meth_name <- strsplit(meth_file, split='.', fixed=T)[[1]][1]
  cov_name <- strsplit(cov_file, split='.', fixed=T)[[1]][1]

  chrs <- gsub('chr', '', as.vector(seqnames(all.bs)))

  df <- data.frame(chr=chrs, start=start(all.bs)-1,
                   end=end(all.bs), meth=getMeth(all.bs, type='raw'),
                   cov=getCoverage(all.bs))
  names(df)[4:5] <- c('meth', 'cov')

  # Write out methylation bedGraph header
  writeLines(paste0('track type=bedGraph name=', meth_name,
    ' description=', meth_name,
    ' visibility=display_mode color=0,0,255 altColor=0,0,255 ',
    ' priority=20 autoScale=off graphType=bar',
    ' viewLimits=0:1 yLineOnOff=off windowingFunction=mean+whiskers',
    ' smoothingWindow=off'), paste0(outdir,meth_file))

  # Write out methylation data
  write.table(df[,c('chr', 'start', 'end', 'meth')], 
    file=paste0(outdir, meth_file), append=T, quote=F, sep='\t', 
    row.names=F, col.names=F)  

  # Write out coverage bedGraph header
  writeLines(paste0('track type=bedGraph name=', cov_name,
    ' description=', cov_name,
    ' visibility=display_mode color=0,0,255 altColor=0,0,255 ',
    ' priority=20 autoScale=on graphType=bar',
    ' yLineOnOff=off windowingFunction=mean+whiskers',
    ' smoothingWindow=off'), paste0(outdir, cov_file))

  # Write out coverage data
  write.table(df[,c('chr', 'start', 'end', 'cov')], 
    file=paste0(outdir, cov_file), append=T, quote=F, sep='\t', 
    row.names=F, col.names=F)

  # Compress files
  system(paste0('gzip ', outdir, meth_file))
  system(paste0('gzip ', outdir, cov_file))

  1
}