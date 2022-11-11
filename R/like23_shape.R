#' likelihood function
#'
#' @description evaluate likelihood ...
#'
#' @param weibull model parameters
#'
#' @return likelihood
#'
#' @examples likelihood for given parameter set
#' 
like23_shape = function(shape23_0, dat0, c23_0, scale23_0, theta23_0, omega13_z0){
 b = ((shape23_0 - 1) * (sum(log(dat0$y23 ^ dat0$s23), na.rm = T)) - 
 scale23_0 / shape23_0 * sum(dat0$y23 ^ shape23_0 * exp(c23_0 * omega13_z0 + 
 theta23_0 * dat0$y12
 ), na.rm = T))
 return(b)
}
