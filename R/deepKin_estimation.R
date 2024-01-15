#' deepKin_estimation function
#'
#' @param grm.diag an n vector
#' @param grm.tri when xcohort = F, grm.tri is an n\*(n-1) vector; when xcohort = T, it is an n1\*n2 vector
#' @param xcohort xcohort = T for cross-cohort GRM extraction; xcohort = F for one cohort GRM;
#'                Make sure that when you use xcohort = T, the individual order is organized just the same order cohort by cohort.
#' @param me the effective number of markers
#' @param pop_size1 population size for the first cohort when xcohort = T
#' @param pop_size2 population size for the second cohort when xcohort = T
#'
#' @return deepKin estimates and p-values
#' @export
#'
#' @examples \dontrun{KINGX = deepKin_estimation(grm.diag = grm.diag, grm.tri = grm.tri,
#'                    xcohort = F, me = me)}
deepKin_estimation <- function(grm.diag, grm.tri, xcohort = F, me, pop_size1, pop_size2){

  if(!xcohort){
    n = length(grm.diag)
    a = 0
    KINGX = matrix(NA, n*(n-1)/2, 1)
    for(i in 2:n){
      for(j in 1:(i-1)){
        a = a + 1
        KINGX[a,1] = 1-(grm.diag[i]+grm.diag[j])/2+grm.tri[a]
      }
    }
  } else {
    if(is.null(pop_size1)|is.null(pop_size2)){
      stop("Error: pop_size1 and pop_size2 is not specified in a deepkin estimation")
    }
    n = pop_size1 + pop_size2
    grm.diag1 = grm.diag[1:pop_size1]
    grm.diag2 = grm.diag[(pop_size1+1):n]

    a = 0
    KINGX = matrix(NA, pop_size1*pop_size2, 1)
    for(j in 1:pop_size2){
      for(i in 1:pop_size1){
        a = a + 1
        KINGX[a,1] = 1-(grm.diag1[i]+grm.diag2[j])/2+grm.tri[a]
      }
    }
  }

  KINGX = data.frame(KINGX)
  colnames(KINGX) = c("king")
  KINGX$king = as.numeric(KINGX$king)

  ## Calculate p-value
  KINGX$var =  ( 2 * (1-KINGX$king)^2 ) / me
  KINGX$t = KINGX$king / sqrt(KINGX$var)
  KINGX$`-logp` = -pnorm(KINGX$t, lower.tail = F, log.p = T) / log(10)

  return(KINGX[,c("king","-logp")])
}
