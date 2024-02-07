#' deepKin function
#' This returns a list of deepKin results based on three principles, including the suggested threshold for each degree of relatedness
#'
#' @param n The population size
#' @param me The effective number of markers
#' @param alpha Type I error rate
#' @param beta Type II error rate
#' @param max.degree User specified deepest interesting degree
#'
#' @return a list
#' @export
#'
#' @examples deepKin(n = 1000, me = 300, alpha = 0.05, beta = 0.1, max.degree = 5)
deepKin <- function(n, me, alpha = 0.05, beta = 0.1, max.degree = 5){

  npairs = n*(n-1)/2
  df1 = data.frame(degree = 0:max.degree,
                   Theta  = (1/2)^(0:max.degree),
                   Me.min = format(me.min(theta = (1/2)^(0:max.degree), alpha = alpha/npairs, beta = beta), digits = 4))

  theta.min = theta.min(me = me, alpha = alpha/npairs, beta = beta)
  delta = log(theta.min, base = 1/2)

  thrd = (1/2)^seq(0.5, 10.5, 1)
  thrd = thrd[which(thrd>theta.min)]
  thrd = format(c(1,thrd), digits = 2)
  range = paste0("(", thrd[-1], ",", thrd[-length(thrd)],"]")

  df2 = data.frame(degree = 0:(length(range)-1),
                   Range = range)
  df2 = rbind(df2, c("delta", format(theta.min, digits = 2)))

  df3 = data.frame(degree = 0:max.degree,
                   Theta = (1/2)^(0:max.degree),
                   Power.max = format(power.max(me, theta = (1/2)^(0:max.degree), alpha = 0.05/npairs), digits = 4))

  deepKin = list(n = n,
                 me = me,
                 npairs = npairs,
                 me.min = df1,
                 theta.min = theta.min,
                 delta = delta,
                 threshold = df2,
                 power.max = df3)


  return(deepKin)
}
