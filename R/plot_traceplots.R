#' traceplots
#'
#' @description creates traceplots
#'
#' @param simulation results
#'
#' @return traceplots
#'
#' @examples example traceplots
#' 
plot_traceplots = function(params_matrix, variable){
 param = params_matrix$params
 
 plot(eval(parse(text = paste0("param$", variable))), ylab = "Parameter Draw",
 xlab = "MCMC Iteration", main = paste("Traceplot of Parameter", variable))
}
