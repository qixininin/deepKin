#' extractPlinkGRM function
#'
#' @param bfileprefix plink bfile prefix
#' @param xcohort xcohort = T for cross-cohort GRM extraction; xcohort = F for one cohort GRM;
#'                Make sure that when you use xcohort = T, the individual order is organized just the same order cohort by cohort.
#' @param pop_size population size; population size for the first cohort when xcohort = T
#' @param pop_size2 population size for the second cohort when xcohort = T
#'
#' @return a list $diag the diagnal elements, $tri the triangle elements
#' @export
#'
#' @examples \dontrun{grm.rst = extractPlinkGRM(bfileprefix, xcohort = F, pop_size = n)
#'                    grm.diag = grm.rst$diag
#'                    grm.tri  = grm.rst$tri}
extractPlinkGRM <- function(bfileprefix, xcohort = F, pop_size, pop_size2){

  if(!xcohort){  ## single cohort grm calculation

    if(!file.exists(paste0(bfileprefix,".rel.bin"))){
      stop(paste0("Error:", bfileprefix,".rel.bin and .rel.id does not exist!"))
    }
    if(is.null(pop_size)){
      stop("Error: pop_size is not specified in the within-cohort grm extraction")
    }

    grm = readBin(paste0(bfileprefix, ".rel.bin"), what="numeric", n=pop_size*(pop_size+1)/2, size=4)
    grm.diag = grm[cumsum(1:pop_size)]   # n vector
    grm.tri  = grm[-cumsum(1:pop_size)]  # n*(n-1)/2 vector

    return(list(diag = grm.diag, tri = grm.tri))

  } else { ## cross-cohort grm calculation

    if(!file.exists(paste0(bfileprefix,".rel.bin"))){
      stop(paste0("Error:", bfileprefix,".rel.bin and .rel.id does not exist!"))
    }

    if(is.null(pop_size) | is.null(pop_size2)){
      stop("Error: pop_size and pop_size2 is not specified in a cross-cohort grm extraction")
    }
    n = pop_size + pop_size2
    grm = readBin(paste0(bfileprefix,".rel.bin"), what="numeric", n=n*(n+1)/2, size=4)
    G = matrix(NA, n, n)
    G[upper.tri(G, diag = T)] = grm

    grm.diag = as.matrix(diag(G))
    grm.tri  = G[1:pop_size,(pop_size+1):n]   # n1*n2 matrix

    return(list(diag = grm.diag, tri = grm.tri))
  }

}
