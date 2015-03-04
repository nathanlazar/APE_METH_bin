# nathan dot lazar at gmail dot com

# Combines results from all permutation analyses into one 
# data frame

combine_per <- function(bp.permute, gene.permute.list, 
  rep.permute, rep.permute.list, SINE.permute.list, cpg.permute.list) {

  df <- matrix(0, nrow=18 , ncol=9, 
               dimnames=list(c('all', 'gene', 'exon', 'intron',
                               'promoter', 'reps', 'LINE', 
                               'SINE', 'DNA', 'LTR', 'Satellite',
                               'Alu', 'MIR', 'AluJ', 'AluS', 
                               'AluY', 'cpg_isl', 'cpg_shore'),
                             c('count', 'meth.p', 'cov.p', 'cpg.p', 'per.p',
                               'meth.n', 'cov.n', 'cpg.n', 'per.n')))

  df[1,] <- c(unlist(bp.permute[c('bp.count', 'meth.p', 'cov.p', 'cpg.p')]),NA,rep(1008,4))
  for(i in 1:length(gene.permute.list))
    df[i+1,] <- unlist(gene.permute.list[[i]][c('bp.count', 'meth.p', 'cov.p', 'cpg.p',
                                                'per.p', 'meth.n', 'cov.n', 'cpg.n', 'per.n')])
  df[6,] <- unlist(rep.permute[c('bp.count', 'meth.p', 'cov.p', 'cpg.p',
                                                'per.p', 'meth.n', 'cov.n', 'cpg.n', 'per.n')])
  for(i in 1:length(rep.permute.list))
    df[i+6,] <- unlist(rep.permute.list[[i]][c('bp.count', 'meth.p', 'cov.p', 'cpg.p',
                                                'per.p', 'meth.n', 'cov.n', 'cpg.n', 'per.n')])
  for(i in 1:length(SINE.permute.list))
    df[i+11,] <- unlist(SINE.permute.list[[i]][c('bp.count', 'meth.p', 'cov.p', 'cpg.p',
                                                'per.p', 'meth.n', 'cov.n', 'cpg.n', 'per.n')])
  for(i in 1:length(cpg.permute.list))
    df[i+16,] <- unlist(cpg.permute.list[[i]][c('bp.count', 'meth.p', 'cov.p', 'cpg.p',
                                                'per.p', 'meth.n', 'cov.n', 'cpg.n', 'per.n')])
  df
}