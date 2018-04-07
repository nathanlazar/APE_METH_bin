# nathan dot lazar at gmail dot com

get_mad <- function(bp.lr.gr) {
  # Gets dataframe to calculate mean absolute difference between sides of breakpoints.
  
#  if(ncol(mcols(bp.lr.gr)) > 5) {
#    bp.lr.df <- data.frame(mcols(bp.lr.gr)[,-1])
#    names(bp.lr.df)[7] <- 'methylation'
#  } else {
#    bp.lr.df <- data.frame(mcols(bp.lr.gr))
#    names(bp.lr.df)[2] <- 'methylation'
#    bp.lr.df$size <- width(bp.lr.gr)
#  }

  bp.lr.df <- data.frame(mcols(bp.lr.gr))
  names(bp.lr.df)[names(bp.lr.df)=='meth'] <- 'methylation'

  # Add in info on left right for class II-B breakpoints
  bp.lr.df <- bp.lr.df %>%
    arrange(BP_name) %>%
    mutate(s_name = sub("(NLE[0-9|X]+_[0-9]+[a|b]*)_[L|R]*_[l|r]*", "\\1", BP_name)) %>%
    mutate(side = sub("NLE[0-9|X]+_[0-9]+[a|b]*_([L|R]*)_[l|r]", "\\1", BP_name))

#  # For class II-B mark "a" or "A" as left and "b" or "B" as right
#  bp.lr.df$side[bp.lr.df$a_or_b %in% c("a", "A", "a_A")] <- "left"
#  bp.lr.df$side[bp.lr.df$a_or_b %in% c("b", "B", "a_B")] <- "right"

  bp.lr.mad <- bp.lr.df %>%
    group_by(s_name) %>%
    summarise(mad_meth=mad(methylation),
              mad_cov=mad(cov),
              mad_cpgs=mad(cpgs)) %>%
    filter(!is.na(mad_meth))

  bp.lr.mad <- merge(bp.lr.mad, bp.lr.df)

  return(bp.lr.mad)
}
