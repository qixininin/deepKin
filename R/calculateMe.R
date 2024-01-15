#' calculateMe function
#'
#' @param bfileprefix bfileprefix of your plink bfiles.
#' @param method
#' method = "LB": a approximation approach. LB method is applied using "GEAR". The default interation time is 100.
#' method = "GRM": based on 1/var(grm.off.diag); GRM method is applied using "plink"
#' We highly suggest using method = "LB" for biobank-scale data.
#' @param plink_path plink path
#' @param gear_path gear path
#' @param freq_path frequency file (.freq) path
#' @param pop_size population size
#' @return me
#' @export
#' @importFrom stats var
#' @importFrom utils read.table
#'
#' @examples \dontrun{calculateMe(bfileprefix = "1KG-EUR.qc.inter", method = "LB")}

calculateMe <- function(bfileprefix, method, plink_path = NULL, gear_path = NULL, freq_path = NULL, pop_size = NULL){
  switch(method,
         "GRM" = {
           if(is.null(plink_path)) { stop( "Error in calculateMe(): No plink_path was specified! ") }

           ## plink1.9
           # if(!file.exists(paste0(bfileprefix,".grm.gz"))){
           #   system(paste0(plink_path, " --silent --bfile ", bfileprefix, " --make-grm-gz --out ", bfileprefix)) # plink
           # }
           # grm = read.table(gzfile(paste0(bfileprefix,".grm.gz")), as.is = T)
           # me = 1/var(grm[grm[,1]!=grm[,2], 4], na.rm = TRUE)

           ## plink2.0
           if(!file.exists(paste0(bfileprefix,".rel.bin"))){
             system(paste0(plink_path, " --silent --bfile ", bfileprefix, " --make-rel triangle bin4 --out ", bfileprefix))
           }
           grm = readBin(paste0(bfileprefix, ".rel.bin"), what="numeric", n=pop_size*(pop_size+1)/2, size=4)
           me = 1/var(grm[-cumsum(1:pop_size)])
         },
         "LB"  = {
           if(is.null(gear_path)) { stop( "Error in calculateMe(): No plink_path was specified! ") }

           if(!file.exists(paste0(bfileprefix,".it.me"))){
             system(paste0(gear_path, " --me --bfile ",bfileprefix," --iter 100 --out ",bfileprefix))
           }
           itme = read.table(paste0(bfileprefix,".it.me"), header = T)
           me = itme[nrow(itme),"Me"]
         }
  )
  return(me)
}
