#' logalpha function
#' This function is related to the third principle in deepKin analysis
#' This function calculates the type I error rate based on me, target relatedness, and type II error rate.
#'
#' @param me The effective number of markers
#' @param theta Target relatedness
#' @param beta Type II error rate
#'
#' @return log10 alpha
#' @export
#' @importFrom stats pnorm
#'
#' @examples logalpha(me = 5000, theta = 0.5, beta = 0.1)
logalpha <- function(me, theta, beta){
  # logalpha = pnorm(sqrt(me/2)*theta-qnorm(1-beta), lower.tail = FALSE, log.p = TRUE) / log(10) # approximation
  logalpha = pnorm((sqrt(me/2)+qnorm(1-beta))*theta-qnorm(1-beta), lower.tail = FALSE, log.p = TRUE) / log(10)
  return(logalpha)
}
