# nathan dot lazar at gmail dot com

# Some plotting functions for comparing the two sides of breakpoints

library(GenomicRanges)
library(ggplot2)
library(dplyr)

mad <- function(list) abs(list[1] - list[2])

get_mad <- function(bp.lr.gr) {
  bp.lr.df <- data.frame(mcols(bp.lr.gr)[,-1])
  names(bp.lr.df)[6] <- 'methylation'
  
  # Add in info on left right for class II-B breakpoints
  bp.lr.df <- bp.lr.df %>%
    arrange(BP_name) %>%
    mutate(s_name = sub("(NLE[0-9|X]+[a|b]*_[0-9]+)[a|b]*", "\\1", BP_name)) %>%
    mutate(a_or_b = sub("NLE[0-9|X]+[a|b]*_[0-9]+([a|b]*)", "\\1", BP_name)) 
  
  # For class II-B mark "a" as left and "b" as right
  bp.lr.df$side[bp.lr.df$a_or_b=="a"] <- "left"
  bp.lr.df$side[bp.lr.df$a_or_b=="b"] <- "right"
  
  bp.lr.mad <- bp.lr.df %>%
    group_by(s_name) %>%
    summarise(mad_meth=mad(methylation),
              mad_cov=mad(cov),
              mad_cpgs=mad(cpgs)) %>%
    filter(!is.na(mad_meth)) 
  
  bp.lr.mad <- merge(bp.lr.mad, bp.lr.df) %>%
    select(-a_or_b)
  
  return(bp.lr.mad)
}

get_mad_perm <- function(regions.gr) {
  # Get the mean absolute difference of methylation, coverage and cpgs in
  # a set of randomly selected regions.
  bp.lr.df <- data.frame(mcols(bp.lr.gr)[,-1])
  names(bp.lr.df)[6] <- 'methylation'
  
  # Add in info on left right for class II-B breakpoints
  bp.lr.df <- bp.lr.df %>%
    arrange(BP_name) %>%
    mutate(s_name = sub("(NLE[0-9|X]+[a|b]*_[0-9]+)[a|b]*", "\\1", BP_name)) %>%
    mutate(a_or_b = sub("NLE[0-9|X]+[a|b]*_[0-9]+([a|b]*)", "\\1", BP_name)) 
  
  # For class II-B mark "a" as left and "b" as right
  bp.lr.df$side[bp.lr.df$a_or_b=="a"] <- "left"
  bp.lr.df$side[bp.lr.df$a_or_b=="b"] <- "right"
  
  bp.lr.mad <- bp.lr.df %>%
    group_by(s_name) %>%
    summarise(mad_meth=mad(methylation),
              mad_cov=mad(cov),
              mad_cpgs=mad(cpgs)) %>%
    filter(!is.na(mad_meth)) 
  
  bp.lr.mad <- merge(bp.lr.mad, bp.lr.df) %>%
    select(-a_or_b)
  
  return(bp.lr.mad)
}


# Function to get linear model parameters and format for printing on plot
lm_eqn = function(df.x, df.y){
  m = lm(df.y ~ df.x);
  eq <- substitute(italic(y) == a + b %.% italic(x)*
                     ","~~italic(r)^2~"="~r2*
                     ","~~italic(p)~"="~p_value,
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p_value = format(summary(m)$coef[2,4], digits=2)))
  as.character(as.expression(eq))
}

# Plot methylation by cpgs
plot_meth_cpgs <- function(bp.lr.df) {
  p <- ggplot(data=bp.lr.df, aes(x=cpgs, y=methylation)) +
#   aes(x=cpgs, y=methylation, color=side, shape=class)
    geom_point(size=5) +
#    facet_wrap(. ~ class) +
    geom_smooth(method=lm, formula=y~x) +
    geom_text(aes(x = 75, y = .2, label=lm_eqn(bp.lr.df$cpgs, bp.lr.df$methylation)), parse = TRUE)
  return(p)
}

# Plot methyation by coverage
plot_meth_cov <- function(bp.lr.df) {
  p <- ggplot(data=bp.lr.df, aes(x=cov, y=methylation)) +
#   aes(x=cov, y=methylation, color=side, shape=class)
    geom_point(size=5) +
#    facet_grid(. ~ class) +
    geom_smooth(method=lm) +
    geom_text(aes(x = 15, y = .25, label=lm_eqn(bp.lr.df$cov, bp.lr.df$methylation)), parse = TRUE)
  return(p)
}


plot_mad_meth <- function(bp.lr.mad) {
  p <- ggplot(data=bp.lr.mad, aes(x=methylation, y=mad_meth, color=class)) +
    #   aes(x=cov, y=methylation, color=side, shape=class)
    geom_point(size=5)
  
  # Get data frame used to connect dots form the same bp with lines
  bp.lr.mad2 <- bp.lr.mad %>%
    select(s_name, mad_meth, side, methylation, class) %>%
    dcast(s_name + mad_meth + class ~ side, value.var='methylation', fun.aggregate=max) %>%
    filter(left != -Inf & right != -Inf)
    
  p <- p + geom_segment(data=bp.lr.mad2, aes(x=left, y=mad_meth, yend=mad_meth, xend=right, color=class))
  return(p)
}

plot_mad_cov <- function(bp.lr.mad) {
  p <- ggplot(data=bp.lr.mad, aes(x=cov, y=mad_cov, color=class)) +
    geom_point(size=5)
  
  # Get data frame used to connect dots form the same bp with lines
  bp.lr.mad2 <- bp.lr.mad %>%
    select(s_name, mad_cov, side, cov, class) %>%
    dcast(s_name + mad_cov + class ~ side, value.var='cov', fun.aggregate=max) %>%
    filter(left != -Inf & right != -Inf)
    
  p <- p + geom_segment(data=bp.lr.mad2, aes(x=left, y=mad_cov, yend=mad_cov, xend=right, color=class))
  return(p)
}

plot_mad_cpgs <- function(bp.lr.mad) {
  p <- ggplot(data=bp.lr.mad, aes(x=cpgs, y=mad_cpgs, color=class)) +
    geom_point(size=5)
  
  # Get data frame used to connect dots form the same bp with lines
  bp.lr.mad2 <- bp.lr.mad %>%
    select(s_name, mad_cpgs, side, cpgs, class) %>%
    dcast(s_name + mad_cpgs + class ~ side, value.var='cpgs', fun.aggregate=max) %>%
    filter(left != -Inf & right != -Inf)
    
  p <- p + geom_segment(data=bp.lr.mad2, aes(x=left, y=mad_cpgs, yend=mad_cpgs, xend=right, color=class))
  return(p)
}


