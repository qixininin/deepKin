#' deepKin.estimation function
#'
#' @param bfileprefix plink bfile prefix
#' @param plink_path path to plink executable file
#' @param xcohort xcohort = T for cross-cohort GRM extraction; xcohort = F for one cohort GRM;
#'                Make sure that when you use xcohort = T, the individual order is organized just the same order cohort by cohort.
#' @param me the effective number of markers
#' @param pop_size1 population size for the first cohort when xcohort = T
#' @param pop_size2 population size for the second cohort when xcohort = T
#'
#' @return a dataframe with two columns
#'         $theta: deepKin estimates
#'         $minuslogp: deepKin -log10(p-values)
#' @export
#'
#' @examples \dontrun{KINGX = deepKin.estimation(grm.diag = grm.diag, grm.tri = grm.tri,
#'                    xcohort = F, me = me)}
deepKin.estimation <- function(bfileprefix, plink_path, me, xcohort = F, pop_size1, pop_size2){

  if(xcohort){
    if(is.null(pop_size1) | is.null(pop_size2)){
      stop("Error deepKin.estimation(): pop_size1 or pop_size2 is not provided.")
    } else {
      cat("*** Cross-cohort deepKin estimation START.\n")
      cat(paste("    Cohort1 size:", pop_size1,"\n"))
      cat(paste("    Cohort2 size:", pop_size2,"\n"))
    }
  } else {
    if(is.null(pop_size1)){
      stop("Error deepKin.estimation(): pop_size1  is not provided.")
    } else {
      cat("*** Within-cohort deepKin estimation START.\n")
      cat(paste("    Cohort size:", pop_size1,"\n"))
    }
  }

  ## Perform plink GRM
  if(!file.exists(paste0(bfileprefix,".rel.bin"))){
    system(paste0(plink_path, " --silent --bfile ", bfileprefix, " --make-rel triangle bin4 --out ", bfileprefix))
  }
  cat("*** plink GRM calculation DONE. \n")

  ## Extract plink GRM
  grm.rst = extract.plink.grm(bfileprefix, xcohort = xcohort, pop_size1 = pop_size1, pop_size2 = pop_size2)
  grm.diag = grm.rst$diag
  grm.tri  = grm.rst$tri

  cat("*** plink GRM extraction DONE. \n")

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
  colnames(KINGX) = c("theta")
  KINGX$theta = as.numeric(KINGX$theta)

  ## Calculate p-value
  KINGX$minuslogp = -pnorm(KINGX$theta/sqrt(2/me), mean = 0, lower.tail = F, log.p = T) / log(10)

  cat("*** deepKin estimation DONE. \n")

  return(KINGX[,c("theta","minuslogp")])
}
