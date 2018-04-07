# nathan dot lazar at gmail dot com

# R script to colate group level p-values from permutations looking 
# at the sides of breakpoints in gibbon, human, chimp, orangutan & rhesus

# Usage: colate_side_permute_data.R PERMUTE_ANALYSIS perm_group_sum.txt

# Output is a table of p-values for analyses run written to the second argument

args <- commandArgs(TRUE)
permdir <- args[1]
outfile <- args[2]

gibbon <- readLines(paste0(permdir, '/Gibbon/full_analysis.out'))
human <- readLines(paste0(permdir, '/Human/full_analysis.out'))
chimp <- readLines(paste0(permdir, '/Chimp/full_analysis.out'))
orangutan <- readLines(paste0(permdir, '/Orangutan/full_analysis.out'))
rhesus <- readLines(paste0(permdir, '/Rhesus/full_analysis.out'))

sp_list <- list(gibbon=gibbon, human=human,
                chimp=chimp, orangutan=orangutan,
                rhesus=rhesus)

get_p_vals <- function(species) {
  idx <- which(grepl('\"Methylation:', species))
  meth.p_val <- as.numeric(gsub('\"', '', 
    gsub('\\[1] \"Methylation: ', '', species[idx])))
  meth.p_val <- 1-meth.p_val

  idx <- which(grepl('\"CpG counts: ', species))
  cpg.p_val <- as.numeric(gsub('\"', '', 
    gsub('\\[1] \"CpG counts: ', '', species[idx])))
  cpg.p_val <- 1-cpg.p_val

  idx <- which(grepl('\"Coverage: ', species))
  cov.p_val <- as.numeric(gsub('\"', '', 
    gsub('\\[1] \"Coverage: ', '', species[idx])))
  cov.p_val <- 1-cov.p_val

  c(meth.p_val, cpg.p_val, cov.p_val)
}

species <- c('gibbon', 'human', 'chimp', 'orangutan', 'rhesus')

p_vals <- data.frame(species=species,
                     meth.adj=0, meth.disj=0,
                     cpg.adj=0, cpg.disj=0,
                     cov.adj=0, cov.disj=0)

for(sp in species) p_vals[p_vals$species==sp,-1] <- get_p_vals(sp_list[[sp]])

write.table(file=outfile, p_vals, sep='\t', quote=F, row.names=F, col.names=T)