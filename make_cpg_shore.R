# nathan dot lazar at gmail dot com

make_cpg_shore <- function(cpg_island.gr, shore_size) {
########################################################
# Makes GRanges object of shores up to size <shore_size>
# that are not overlapping other CpG islands
########################################################

  upstream.gr <- promoters(cpg_island.gr, 
                           upstream=shore_size)
  cpg_isl.rev.gr <- cpg_island.gr
  strand(cpg_isl.rev.gr) <- '-'
  downstream.gr <- promoters(cpg_isl.rev.gr, 
                             upstream=shore_size)
  strand(downstream.gr) <- '*'

  rough.gr <- c(upstream.gr, downstream.gr)

  out_isl.gr <- gaps(cpg_island.gr)
  
  cpg_shore.gr <- GenomicRanges::intersect(out_isl.gr, rough.gr)
  cpg_shore.gr
}
