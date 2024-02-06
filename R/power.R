#' power function
#'
#' @param me The effective number of markers
#' @param theta Target relatedness
#' @param alpha Type I error rate, usually sets to 0.05/N after Bonferroni correction
#'
#' @return power
#' @export
#'
#' @examples deepKin_power(me = 5000, theta = 0.5, alpha = 0.05/100000)
deepKin_power <- function(me, theta, alpha)
{
  power = pnorm((sqrt(me/2)*theta-qnorm(1-alpha))/(1-theta), lower.tail = TRUE)
  return(power)
}
