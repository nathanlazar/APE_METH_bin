# nathan dot lazar at gmail dot com

# Some plotting functions for comparing the two sides of breakpoints

library(GenomicRanges)
if (!is.element('ggplot2', installed.packages()[,1]))
        install.packages('ggplot2', dependencies = TRUE, 
                         repos='http://cran.us.r-project.org')
library(ggplot2)
library(dplyr)

mad <- function(list) abs(list[1] - list[2])

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
    geom_point(size=5) +
    geom_smooth(method=lm, formula=y~x) +
    geom_text(mapping=aes(x = 180, y = .2, label=lm_eqn(bp.lr.df$cpgs, bp.lr.df$methylation)), 
              parse = TRUE, hjust=0.5)
  return(p)
}

# Plot methyation by coverage
plot_meth_cov <- function(bp.lr.df) {
  p <- ggplot(data=bp.lr.df, aes(x=cov, y=methylation)) +
    geom_point(size=5) +
    geom_smooth(method=lm) +
    geom_text(aes(x = 12, y = .25, label=lm_eqn(bp.lr.df$cov, bp.lr.df$methylation)), 
              parse = TRUE, hjust=0.5)
  return(p)
}

plot_mad_meth <- function(bp.lr.mad) {
  if('class' %in% names(bp.lr.mad)) {
    p <- ggplot(data=bp.lr.mad, aes(x=methylation, y=mad_meth, color=class)) +
      #   aes(x=cov, y=methylation, color=side, shape=class)
      geom_point(size=5)
    bp.lr.mad2 <- bp.lr.mad %>%
      select(s_name, mad_meth, side, methylation, class) %>%
      dcast(s_name + mad_meth + class ~ side, value.var='methylation', fun.aggregate=max) %>%
      filter(L != -Inf & R != -Inf)
    # Get data frame used to connect dots form the same bp with lines
    p <- p + geom_segment(data=bp.lr.mad2, aes(x=L, y=mad_meth, yend=mad_meth, xend=R, color=class))
  } else {
    p <- ggplot(data=bp.lr.mad, aes(x=methylation, y=mad_meth)) + geom_point(size=5)
    bp.lr.mad2 <- bp.lr.mad %>%
      select(s_name, mad_meth, side, methylation) %>%
      dcast(s_name + mad_meth ~ side, value.var='methylation', fun.aggregate=max) %>%
      filter(L != -Inf & R != -Inf)
    p <- p + geom_segment(data=bp.lr.mad2, aes(x=L, y=mad_meth, yend=mad_meth, xend=R))
  }
  return(p)
}

plot_mad_cov <- function(bp.lr.mad) {
  if('class' %in% names(bp.lr.mad)) {
    p <- ggplot(data=bp.lr.mad, aes(x=cov, y=mad_cov, color=class)) +
      geom_point(size=5) 
    # Get data frame used to connect dots form the same bp with lines
    bp.lr.mad2 <- bp.lr.mad %>%
      select(s_name, mad_cov, side, cov, class) %>%
      dcast(s_name + mad_cov + class ~ side, value.var='cov', fun.aggregate=max) %>%
      filter(L != -Inf & R != -Inf)
    p <- p + geom_segment(data=bp.lr.mad2, aes(x=L, y=mad_cov, yend=mad_cov, xend=R, color=class))
    return(p)
  } else {
    p <- ggplot(data=bp.lr.mad, aes(x=cov, y=mad_cov)) + geom_point(size=5) 
    # Get data frame used to connect dots form the same bp with lines
    bp.lr.mad2 <- bp.lr.mad %>%
      select(s_name, mad_cov, side, cov) %>%
      dcast(s_name + mad_cov ~ side, value.var='cov', fun.aggregate=max) %>%
      filter(L != -Inf & R != -Inf)
    p <- p + geom_segment(data=bp.lr.mad2, aes(x=L, y=mad_cov, yend=mad_cov, xend=R))
    return(p)
  }
}

plot_mad_cpgs <- function(bp.lr.mad) {
  if('class' %in% names(bp.lr.mad)) {
    p <- ggplot(data=bp.lr.mad, aes(x=cpgs, y=mad_cpgs, color=class)) +
      geom_point(size=5)
    # Get data frame used to connect dots form the same bp with lines
    bp.lr.mad2 <- bp.lr.mad %>%
      select(s_name, mad_cpgs, side, cpgs, class) %>%
      dcast(s_name + mad_cpgs + class ~ side, value.var='cpgs', fun.aggregate=max) %>%
      filter(L != -Inf & R != -Inf)
    p <- p + geom_segment(data=bp.lr.mad2, aes(x=L, y=mad_cpgs, yend=mad_cpgs, xend=R, color=class))
    return(p)
  } else {
    p <- ggplot(data=bp.lr.mad, aes(x=cpgs, y=mad_cpgs)) + geom_point(size=5)  
    # Get data frame used to connect dots form the same bp with lines
    bp.lr.mad2 <- bp.lr.mad %>%
      select(s_name, mad_cpgs, side, cpgs) %>%
      dcast(s_name + mad_cpgs ~ side, value.var='cpgs', fun.aggregate=max) %>%
      filter(L != -Inf & R != -Inf)  
    p <- p + geom_segment(data=bp.lr.mad2, aes(x=L, y=mad_cpgs, yend=mad_cpgs, xend=R))
    return(p)
  }
}
