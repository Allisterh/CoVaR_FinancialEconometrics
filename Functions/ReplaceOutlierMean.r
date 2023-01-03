replace_outlier_with_mean <- function(x) {
  replace(x, x %in% boxplot.stats(x)$out, mean(x))  
}