#' deepKin.classification
#'
#' @param theta.x input a vector of theta.x to be tested
#' @param me the number of effective markers
#' @param alpha significant level
#'
#' @return a data frame with theta, logpt, logpt1, class, quality
#' @export
#'
#' @examples \dontrun{qc = kingx[which(kingx$king>theta.min), "king"]
#'                    qc.rst = deepKin.classification(theta.x = qc, me = me, alpha = alpha)}
deepKin.classification <- function(theta.x, me, alpha)
{
  ## Settle degree vector into intervals (t,t+1)
  degree.x = log(theta.x, base = 0.5)
  degree.t = floor(degree.x)
  theta.t = (1/2)^degree.t
  theta.t1 = (1/2)^(degree.t+1)

  ## p values for two null distribution
  logp.t  = logp(x = theta.x, null.theta = theta.t, me = me)
  logp.t1 = logp(x = theta.x, null.theta = theta.t1, me = me)

  ## Classification
  class = ifelse(logp.t>logp.t1, degree.t, degree.t+1)
  # ## Classification between 0 and 1 degree
  # ## If p.t1 is significant,     --> 0 degree
  # ## If p.t1 is not significant, --> 1 degree
  # class[which(degree.t==0)] = ifelse(logp.t1[which(degree.t==0)]<log10(alpha), 0,1)
  ## using 0.5 degree as the threshold between 0 and 1
  class[which(degree.t==0)] = ifelse(degree.x[which(degree.t==0)]<0.5, 0, 1)

  ## See if both p.t and p.t1 are significant
  quality = ifelse(logp.t<log10(alpha) & logp.t1<log10(alpha), "0", "1")
  quality[which(degree.t==0)] = "1"

  return(data.frame(theta = theta.x,
                    logpt = logp.t,
                    logpt1 = logp.t1,
                    class = class,
                    quality = quality))

}
