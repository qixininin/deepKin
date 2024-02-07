

deepKin.summary <- function(n, me, alpha = 0.05, beta = 0.1){


  full_summary <- paste("DeepKin Summary", "\n",
                        "------------------------------------", "\n",
                        "Data summary:", "\n",
                        "The number of samples:          ", n,  "\n",
                        "The number of effective markers:", me, "\n",
                        "The number of comparisons:      ", n*(n-1)/2, "\n",
                        "------------------------------------", "\n",
                        "Principle I: the minimum number of me", "\n",
                        )




  df = data.frame(degree = 0:4,
                  me.min = me.min(theta = (1/2)^(0:4), alpha = 0.05/npairs, beta = 0.1))




  return(full_summary)
}
