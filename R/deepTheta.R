#' deepTheta function
#' This function is related to the second principle in deepKin analysis.
#' This function calculates the deepest relatedness supported by the data, based on me, type I and type II error rates
#'
#' @param me The effective number of markers
#' @param alpha Type I error rate
#' @param beta Type II error rate
#'
#' @return deepTheta relatedness score
#' @export
#'
#' @examples deepTheta(me = 5000, alpha = 0.05, beta = 0.1)
deepTheta <- function(me, alpha, beta){
  deepTheta = sqrt( 2/me ) * (qnorm(1-alpha) + qnorm(1-beta))
  return(deepTheta)
}
