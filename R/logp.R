#' logp function
#'
#' @param x test theta vector
#' @param null.theta the theta of null hypothesis, null.theta can not been 0,
#'        because it is one-side test for null.theta = 0, we considered absolute value in this function
#' @param me the number of effective markers
#'
#' @return logp
#' @export
#'
#' @examples logp(x = seq(0.3,0.4,0.01), null.theta = 0.25, me = 1000)
logp <- function(x, null.theta, me)
{
  logp = pnorm(abs(x-null.theta) /sqrt(2*(1-null.theta)^2/me) , lower.tail = F, log.p = T) / log(10)
  return(logp)
}
