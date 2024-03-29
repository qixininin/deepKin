#' extract.indi.id function
#' This function will return individual ID based on index in deepKin results and individual IDs
#' @param id individual IDs
#' @param index index in deepKin results
#' @param xcohort xcohort = T for cross-cohort GRM extraction; xcohort = F for one cohort GRM;
#'                Make sure that when you use xcohort = T, the individual order is organized just the same order cohort by cohort.
#' @param pop_size1 population size for the first cohort when xcohort = T
#' @param pop_size2 population size for the second cohort when xcohort = T
#'
#' @return two columns, the first one is ID1 and the second column is ID2
#' @export
#'
#' @examples \dontrun{deepkin.qc.id = extract.indi.id(id = id, index = index, xcohort = F)}
extract.indi.id <- function(id, index, xcohort, pop_size1 = NULL, pop_size2 = NULL)
{
  pop_size = length(id)

  if(!xcohort){

    base = c(1, cumsum(1:(pop_size-1)) + 1)
    row = findInterval(index, base)
    col = index-base[row]+1
    row = row + 1

    return(data.frame(ID1=id[row],ID2=id[col]))
  } else {
    id1 = id[1:pop_size1]
    id2 = id[(pop_size1+1):pop_size]
    row = index %% pop_size1
    col = index %/% pop_size1 + 1
    col[which(row==0)]=index[which(row==0)]/pop_size1
    row[which(row==0)]=pop_size1

    return(data.frame(ID1=id1[row],ID2=id2[col]))
  }
}
