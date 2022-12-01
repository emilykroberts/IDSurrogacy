#' traceplots
#'
#' @description creates traceplots
#'
#' @param 
#' params_matrix: matrix of parameters from mcmc
#' variable: variable of which to plot
#'
#' @return traceplots
#'
#' @examples
#' example(plot_traceplots(params_matrix = params_matrix, variable = "int"))
plot_traceplots = function(params_matrix, variable){
 param = params_matrix$params
 
 plot(eval(parse(text = paste0("param$", variable))), ylab = "Parameter Draw",
 xlab = "MCMC Iteration", main = paste("Traceplot of Parameter", variable))
}
