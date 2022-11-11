#' plot results from mcmc
#'
#' @description plot the results in a cep curve
#'
#' @param results from mcmc
#'
#' @return cep curve
#'
#' @examples example cep from simulations
#' 
plot_results = function(params_matrix, write){
 
 param = params_matrix$params
 n = params_matrix$args$n
 burnin = params_matrix$args$burnin
 sim = params_matrix$args$SIM
 
 saveCEPx = params_matrix$saveCEPx
 saveCEPy = params_matrix$saveCEPy
 
 dat = data.frame(matrix(data = NA, nrow = n, ncol = 1))
 dat$X = rowMeans(saveCEPx[,burnin:ncol(saveCEPx)], na.rm = T)
 dat$Y = rowMeans(saveCEPy[,burnin:ncol(saveCEPx)], na.rm = T)

 slope = params_matrix$params$slope
 int = params_matrix$params$int

 theme_set( theme_classic(base_size = 16))

 d = ggplot(dat, aes(X, Y), col(c(rep(0, n/2), rep(1, n/2)))) + ylim(c(-1, 1)) +
 ggtitle("Illness-Death CEP Curve" ) + xlab("Delta S_i")+ 
 ylab(TeX("$\\Delta T_i = P(T_i(1) > \\tau_T| \\omega_{.i}^1) - P(T_i(0) > \\tau_T | \\omega_{.i}^0)$$ ")) +
 xlab(TeX("$\\Delta S_i = log \\frac{\\Lambda_{12}^{0} (\\tau_S| x_i,\\omega_{12i}^0)}{\\Lambda_{12}^{1} (\\tau_S| x_i,\\omega_{12i}^1)}$$ ")) + 
 geom_hline(yintercept = 0, linetype = "solid") + geom_vline(xintercept = 0, linetype = "solid") + 
 coord_cartesian(ylim=c(-1,1)) + theme_bw() + geom_point()#aes(color = factor(c(rep(0, n/2), rep(1, n/2))))) 

 for(i in burnin:sim) d = d + geom_abline(slope = slope[i], intercept = int[i] )
 d3 = d + geom_point(aes(color = 1)) + theme(legend.position = "none")

 print(d3)
 if(write) ggsave(paste0("estimatedCEP", ",", scenario, ",", array_id, ".jpeg"), d3)

}
