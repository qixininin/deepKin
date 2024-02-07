#' extract.plink.grm function
#'
#' @param bfileprefix plink bfile prefix
#' @param xcohort xcohort = T for cross-cohort GRM extraction; xcohort = F for one cohort GRM;
#'                Make sure that when you use xcohort = T, the individual order is organized just the same order cohort by cohort.
#' @param pop_size1 population size; population size for the first cohort when xcohort = T
#' @param pop_size2 population size for the second cohort when xcohort = T
#'
#' @return a list $diag the diagnal elements, $tri the triangle elements
#' @export
#'
#' @examples \dontrun{grm.rst = extract.plink.grm(bfileprefix, xcohort = F, pop_size1 = n)
#'                    grm.diag = grm.rst$diag
#'                    grm.tri  = grm.rst$tri}
extract.plink.grm <- function(bfileprefix, xcohort = F, pop_size1, pop_size2){

  if(!xcohort){  ## single cohort grm calculation

    if(!file.exists(paste0(bfileprefix,".rel.bin"))){
      stop(paste0("Error:", bfileprefix,".rel.bin and .rel.id does not exist!"))
    }
    if(is.null(pop_size1)){
      stop("Error: pop_size1 is not specified in the within-cohort grm extraction")
    }

    grm = readBin(paste0(bfileprefix, ".rel.bin"), what="numeric", n=pop_size1*(pop_size1+1)/2, size=4)
    grm.diag = grm[cumsum(1:pop_size1)]   # n vector
    grm.tri  = grm[-cumsum(1:pop_size1)]  # n*(n-1)/2 vector

    return(list(diag = grm.diag, tri = grm.tri))

  } else { ## cross-cohort grm calculation

    if(!file.exists(paste0(bfileprefix,".rel.bin"))){
      stop(paste0("Error:", bfileprefix,".rel.bin and .rel.id does not exist!"))
    }

    if(is.null(pop_size1) | is.null(pop_size2)){
      stop("Error: pop_size1 and pop_size2 is not specified in a cross-cohort grm extraction")
    }
    n = pop_size1 + pop_size2
    grm = readBin(paste0(bfileprefix,".rel.bin"), what="numeric", n=n*(n+1)/2, size=4)
    # base index for each column, totaling n+1 column
    base = c(1, cumsum(1:n) + 1)
    # fetch index for diag and n1-by-n2 individuals
    index.diag = base[-1]-1
    index.tri = c()
    for(y in (pop_size1+1):n)
    {
      index.tri = c(index.tri, base[y]:(base[y]+pop_size1-1))
    }

    return(list(diag = grm[index.diag], tri = grm[index.tri]))

  }

}
