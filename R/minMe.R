#' minMe function
#' This function is related to the first principle in deepKin analysis
#' This function calculates the minimum size of m_e (the effective number of markers) based on
#' target relatedness, type I and type II error rates.
#'
#' @param theta Target relatedness
#' @param alpha Type I error rate
#' @param beta Type II error rate
#'
#' @return minMe the minimum size of m_e
#' @export
#' @importFrom stats qnorm
#'
#' @examples minMe(theta = 0.5, alpha = 0.05, beta = 0.1)
minMe <- function(theta, alpha, beta){
  minMe = ( 2/theta/theta ) * (qnorm(1-alpha) + qnorm(1-beta)*(1-theta))^2
  return(minMe)
}
