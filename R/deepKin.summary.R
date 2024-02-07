#' deepKin.summary function
#' This returns the summary report of deepKin three principles, including the suggested threshold for each degree of relatedness
#'
#' @param deepKin the output of deepKin()
#'
#' @return a string
#' @export
#' @importFrom readr format_delim
#'
#' @examples \dontrun{deepKin.summary(deepKin)}
deepKin.summary <- function(deepKin){

  full_summary <- paste0("                        DeepKin Summary                        ", "\n",
                         "  -------------------------------------------------------------", "\n",
                         "  Data summary:", "\n",
                         "  -------------------------------------------------------------", "\n",
                         "  **The number of samples (n):            ", deepKin$n,  "\n",
                         "  **The number of effective markers (me): ", deepKin$me, "\n",
                         "  **The number of comparisons (N):        ", deepKin$npairs, "\n","\n",



                         "  -------------------------------------------------------------", "\n",
                         "  Principle I: the minimum number of me", "\n",
                         "  -------------------------------------------------------------", "\n",
                         "  ", readr::format_delim(deepKin$me.min, "\t", eol = "\n  "), "\n",


                         "  -------------------------------------------------------------", "\n",
                         "  Principle II: the deepest degree of relatedness supported by data", "\n",
                         "  -------------------------------------------------------------", "\n",
                         "  **Deepest theta (theta.min):    ", round(deepKin$theta.min,4),  "\n",
                         "  **Deepest degree (delta): ", round(deepKin$delta,4),  "\n", "\n",
                         "  **Suggested thresholds: ", "\n",
                         "  ", readr::format_delim(deepKin$threshold, "\t", eol = "\n  "), "\n", "\n",


                         "  -------------------------------------------------------------", "\n",
                         "  Principle III: the maximum power for each degree", "\n",
                         "  -------------------------------------------------------------", "\n",

                         "  ", readr::format_delim(deepKin$power.max, "\t", eol = "\n  "), "\n", "\n"
                        )

  return(full_summary)
}
