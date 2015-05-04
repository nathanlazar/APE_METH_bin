# nathan dot lazar at gmail dot com

add_meth_cpg_cov <- function(gr, all.bs, parallel=T, min_cov=4) {
########################################################
# Add metadata column of mean methylation per region,
# for CpGs w/ coverage >= min_cov, total number of CpGs per region 
# and mean coverage per region for all CpGs
########################################################
  cov_min.bs <- all.bs[getCoverage(all.bs) >= min_cov]
  if(parallel) {
    meth <- mcgetMeth(cov_min.bs, regions=gr, type='raw',
                         what='perBase')
    cov <- mcgetCoverage(all.bs, regions=gr,
                         what='perBase')
  } else {
    meth <- getMeth(cov_min.bs, regions=gr, type='raw',
                       what='perBase')
    cov <- getCoverage(all.bs, regions=gr,
                       what='perBase')
  }
  gr$cpgs_w_cov <- unlist(lapply(meth, length))
  gr$meth <- unlist(lapply(meth, mean, na.rm=T))
  gr$cpgs <- unlist(lapply(cov, length))
  gr$cov <- unlist(lapply(cov, mean, na.rm=T))
  gr
}
