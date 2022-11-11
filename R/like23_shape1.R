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
like23_shape1 = function(shape23_1, dat1, c23_1, scale23_1, theta23_1, omega13_z1){
 b = ((shape23_1 - 1) * (sum(log(dat1$y23 ^ dat1$s23), na.rm = T)) - 
 scale23_1 / shape23_1 * sum(dat1$y23 ^ shape23_1 * exp(c23_1 * omega13_z1 + 
 theta23_1 * dat1$y12
 ), na.rm = T))

return(b)
}
