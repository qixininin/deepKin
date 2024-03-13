#' theta.min function
#' This function calculates the deepest relatedness supported by the data, based on me and siginificant level
#'
#' @param me The effective number of markers
#' @param alpha Significant level
#'
#' @return theta.min deepest theta relatedness score
#' @export
#'
#' @examples theta.min(me = 5000, alpha = 0.05)
theta.min <- function(me, alpha){
  # theta.min = sqrt( 2/me ) * (qnorm(1-alpha) + qnorm(1-beta)) # approximation
  # theta.min = (qnorm(1-alpha) + qnorm(1-beta))/(sqrt( me/2 ) + qnorm(1-beta))
  theta.min = qnorm(alpha, mean = 0, sd = sqrt(2/me), lower.tail = F)
  return(theta.min)
}
