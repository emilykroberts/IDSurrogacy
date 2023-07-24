#' plot results from mcmc
#'
#' @description plot the results in a cep curve
#'
#' @param params_matrix matrix of parameters from mcmc
#' @param write logical value if results should be written to a file
#' @param write logical value if results should be written to a file

#' 
#' @return cep curve
#'
#' @examples 
#' example(plot_results(params_matrix = params_matrix, write = F))
plot_results = function(params_matrix, write, fignum){
  
  if(is.null(fignum))
  fignum = 'X'
    
  param = params_matrix$params
  n = params_matrix$args$n
  burnin = params_matrix$args$burnin
  sim = params_matrix$args$SIM
  
  saveCEPx_plot = params_matrix$saveCEPx
  saveCEPy_plot = params_matrix$saveCEPy
  
  dat = data.frame(matrix(data = NA, nrow = n, ncol = 1))
  dat$X = rowMeans(saveCEPx_plot[, burnin : ncol(saveCEPx_plot)], na.rm = T)
  dat$Y = rowMeans(saveCEPy_plot[, burnin : ncol(saveCEPy_plot)], na.rm = T)
  
  slope = params_matrix$params$slope
  int = params_matrix$params$int
  
  d = ggplot(dat, aes(X, Y)) + ylim(c(-1, 1)) +
    ggtitle("Illness-Death CEP Curve" ) + xlab("Delta S_i")+ 
    ylab(TeX("$\\Delta T_i = P(T_i(1) > \\tau_T| \\omega_{.i}^1) - P(T_i(0) > \\tau_T | \\omega_{.i}^0)$$ ")) +
    xlab(TeX("$\\Delta S_i = log \\frac{\\Lambda_{12}^{0} (\\tau_S| \\omega_{12i}^0)}{\\Lambda_{12}^{1} (\\tau_S| \\omega_{12i}^1)}$$ ")) + 
    geom_hline(yintercept = 0, linetype = "solid") + geom_vline(xintercept = 0, linetype = "solid") + 
    coord_cartesian(ylim=c(-1, 1)) + theme_bw(base_size = 20) + geom_point()
  
  for(i in burnin : sim) d = d + geom_abline(slope = slope[i], intercept = int[i], col = 'gray' )
  d3 = d + geom_point(aes(color = 1)) + theme(legend.position = "none") + geom_vline(xintercept = mean((dat$X), na.rm = T), linetype = "dashed", col = 'red') +
    geom_hline(yintercept = mean((dat$Y), na.rm = T), linetype = "dashed", col = 'red') 
  
  print(d3)
  # if(write) ggsave(paste0("estimatedCEP", ",", scenario, ",", array_id, ".jpeg"), d3)
  if(write) ggsave(paste0("Figure", fignum, ".jpeg"), d3)

  
}
