# nathan dot lazar at gmail dot com

add_meth_cpg_cov <- function(gr, all.bs, parallel=T) {
########################################################
# Add metadata column of mean methylation per region,
# number of CpGs per region and mean coverage per region
########################################################
  if(parallel) {
    gr$meth <- mcgetMeth(all.bs, regions=gr, type='raw',
                         what='perRegion')
    cov <- mcgetCoverage(all.bs, regions=gr,
                         what='perBase')
  } else {
    gr$meth <- getMeth(all.bs, regions=gr, type='raw',
                       what='perRegion')
    cov <- getCoverage(all.bs, regions=gr,
                       what='perBase')
  }
  gr$cpgs <- unlist(lapply(cov, length))
  gr$cov <- unlist(lapply(cov, mean, na.rm=T))
  gr
}
