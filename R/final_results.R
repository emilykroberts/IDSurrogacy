#' final results function
#'
#' @description tabulate results from mcmc 
#'
#' @param params_matrix parameter matrix of results
#' @param write logical value if results should be written to a file
#'
#' @return table
#'
#' @examples 
#' example(final_results(params_matrix = params_matrix, write = T))
final_results = function(params_matrix, write){
  ## save results
  param = params_matrix$params
  
  n = params_matrix$args$n
  burnin = params_matrix$args$burnin
  sim = params_matrix$args$SIM
  
  params = data.frame(c12_0 = numeric(1), c12_0SE = numeric(1), 
                      c13_0 = numeric(1), c13_0SE = numeric(1), 
                      c23_0 = numeric(1), c23_0SE = numeric(1), 
                      c12_1 = numeric(1), c12_1SE = numeric(1), 
                      c13_1 = numeric(1), c13_1SE = numeric(1), 
                      c23_1 = numeric(1), c23_1SE = numeric(1), 
                      theta23_0 = numeric(1), theta23_0SE = numeric(1), 
                      theta23_1 = numeric(1), theta23_1SE = numeric(1), 
                      scale12_0 = numeric(1), scale12_0SE = numeric(1), 
                      scale13_0 = numeric(1), scale13_0SE = numeric(1), 
                      scale23_0 = numeric(1), scale23_0SE = numeric(1), 
                      scale12_1 = numeric(1), scale12_1SE = numeric(1), 
                      scale13_1 = numeric(1), scale13_1SE = numeric(1), 
                      scale23_1 = numeric(1), scale23_1SE = numeric(1),
                      shape12_0 = numeric(1), shape12_0SE = numeric(1),
                      shape13_0 = numeric(1), shape13_0SE = numeric(1),
                      shape23_0 = numeric(1), shape23_0SE = numeric(1),
                      shape12_1 = numeric(1), shape12_1SE = numeric(1),
                      shape13_1 = numeric(1), shape13_1SE = numeric(1),
                      shape23_1 = numeric(1), shape23_1SE = numeric(1),
                      beta12_0 = numeric(1), beta12_0SE = numeric(1), 
                      beta13_0 = numeric(1), beta13_0SE = numeric(1), 
                      beta23_0 = numeric(1), beta23_0SE = numeric(1), 
                      beta12_1 = numeric(1), beta12_1SE = numeric(1), 
                      beta13_1 = numeric(1), beta13_1SE = numeric(1), 
                      beta23_1 = numeric(1), beta23_1SE = numeric(1)
  )
  
  params[1] = mean(params_matrix$params$c12_0[burnin : (SIM - 1)], na.rm = T)
  params[2] = sqrt(var(params_matrix$params$c12_0[burnin : (SIM - 1)], na.rm = T))
  params[3] = mean(params_matrix$params$c13_0[burnin : (SIM - 1)], na.rm = T)
  params[4] = sqrt(var(params_matrix$params$c13_0[burnin : (SIM - 1)], na.rm = T))
  params[5] = mean(params_matrix$params$c23_0[burnin : (SIM - 1)], na.rm = T)
  params[6] = sqrt(var(params_matrix$params$c23_0[burnin : (SIM - 1)], na.rm = T))
  params[7] = mean(params_matrix$params$c12_1[burnin : (SIM - 1)], na.rm = T)
  params[8] = sqrt(var(params_matrix$params$c12_1[burnin : (SIM - 1)], na.rm = T))
  params[9] = mean(params_matrix$params$c13_1[burnin : (SIM - 1)], na.rm = T)
  params[10] = sqrt(var(params_matrix$params$c13_1[burnin : (SIM - 1)], na.rm = T))
  params[11] = mean(params_matrix$params$c23_1[burnin : (SIM - 1)], na.rm = T)
  params[12] = sqrt(var(params_matrix$params$c23_1[burnin : (SIM - 1)], na.rm = T))
  params[13] = mean(params_matrix$params$theta23_0[burnin : (SIM - 1)], na.rm = T)
  params[14] = sqrt(var(params_matrix$params$theta23_0[burnin : (SIM - 1)], na.rm = T))
  params[15] = mean(params_matrix$params$theta23_1[burnin : (SIM - 1)], na.rm = T)
  params[16] = sqrt(var(params_matrix$params$theta23_1[burnin : (SIM - 1)], na.rm = T))
  params[17] = mean(params_matrix$params$scale12_0[burnin : (SIM - 1)], na.rm = T)
  params[18] = sqrt(var(params_matrix$params$scale12_0[burnin : (SIM - 1)], na.rm = T))
  params[19] = mean(params_matrix$params$scale13_0[burnin : (SIM - 1)], na.rm = T)
  params[20] = sqrt(var(params_matrix$params$scale13_0[burnin : (SIM - 1)], na.rm = T))
  params[21] = mean(params_matrix$params$scale23_0[burnin : (SIM - 1)], na.rm = T)
  params[22] = sqrt(var(params_matrix$params$scale23_0[burnin : (SIM - 1)], na.rm = T))
  params[23] = mean(params_matrix$params$scale12_1[burnin : (SIM - 1)], na.rm = T)
  params[24] = sqrt(var(params_matrix$params$scale12_1[burnin : (SIM - 1)], na.rm = T))
  params[25] = mean(params_matrix$params$scale13_1[burnin : (SIM - 1)], na.rm = T)
  params[26] = sqrt(var(params_matrix$params$scale13_1[burnin : (SIM - 1)], na.rm = T))
  params[27] = mean(params_matrix$params$scale23_1[burnin : (SIM - 1)], na.rm = T)
  params[28] = sqrt(var(params_matrix$params$scale23_1[burnin : (SIM - 1)], na.rm = T))
  
  params[29] = mean(params_matrix$params$shape12_0[burnin : (SIM - 1)], na.rm = T)
  params[30] = sqrt(var(params_matrix$params$shape12_0[burnin : (SIM - 1)], na.rm = T))
  params[31] = mean(params_matrix$params$shape13_0[burnin : (SIM - 1)], na.rm = T)
  params[32] = sqrt(var(params_matrix$params$shape13_0[burnin : (SIM - 1)], na.rm = T))
  params[33] = mean(params_matrix$params$shape23_0[burnin : (SIM - 1)], na.rm = T)
  params[34] = sqrt(var(params_matrix$params$shape23_0[burnin : (SIM - 1)], na.rm = T))
  params[35] = mean(params_matrix$params$shape12_1[burnin : (SIM - 1)], na.rm = T)
  params[36] = sqrt(var(params_matrix$params$shape12_1[burnin : (SIM - 1)], na.rm = T))
  params[37] = mean(params_matrix$params$shape13_1[burnin : (SIM - 1)], na.rm = T)
  params[38] = sqrt(var(params_matrix$params$shape13_1[burnin : (SIM - 1)], na.rm = T))
  params[39] = mean(params_matrix$params$shape23_1[burnin : (SIM - 1)], na.rm = T)
  params[40] = sqrt(var(params_matrix$params$shape23_1[burnin : (SIM - 1)], na.rm = T))
  
  params[41] = mean(params_matrix$params$beta12_0[burnin : (SIM - 1)], na.rm = T)
  params[42] = sqrt(var(params_matrix$params$beta12_0[burnin : (SIM - 1)], na.rm = T))
  params[43] = mean(params_matrix$params$beta13_0[burnin : (SIM - 1)], na.rm = T)
  params[44] = sqrt(var(params_matrix$params$beta13_0[burnin : (SIM - 1)], na.rm = T))
  params[45] = mean(params_matrix$params$beta23_0[burnin : (SIM - 1)], na.rm = T)
  params[46] = sqrt(var(params_matrix$params$beta23_0[burnin : (SIM - 1)], na.rm = T))
  params[47] = mean(params_matrix$params$beta12_1[burnin : (SIM - 1)], na.rm = T)
  params[48] = sqrt(var(params_matrix$params$beta12_1[burnin : (SIM - 1)], na.rm = T))
  params[49] = mean(params_matrix$params$beta13_1[burnin : (SIM - 1)], na.rm = T)
  params[50] = sqrt(var(params_matrix$params$beta13_1[burnin : (SIM - 1)], na.rm = T))
  params[51] = mean(params_matrix$params$beta23_1[burnin : (SIM - 1)], na.rm = T)
  params[52] = sqrt(var(params_matrix$params$beta23_1[burnin : (SIM - 1)], na.rm = T))
  
  params$int = mean(params_matrix$params$int, na.rm = T)
  params$intse = sd(params_matrix$params$int, na.rm = T)
  params$slope = mean(params_matrix$params$slope, na.rm = T)
  params$slopese = sd(params_matrix$params$slope, na.rm = T)
  
  params = round(params, 3)
  
  print(params[c(T, F)])
  
  if(write){ fname = paste('params', ",", scenario, ",", array_id, '.txt', sep="")
  write.table(params, file = fname, sep = "\t", row.names = F, col.names = T) }
  
}
