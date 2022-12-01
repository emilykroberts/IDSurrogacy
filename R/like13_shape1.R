#' likelihood function
#'
#' @description evaluate likelihood ...
#'
#' @param weibull model parameters
#'
#' @return likelihood
#'
#' @examples
#' example(like13_shape1(shape13_1 = 1, dat1 = 1, c13_1 = 1, scale13_1 = 1, omega13_z1 = 1))
like13_shape1 = function(shape13_1, dat1, c13_1, scale13_1, omega13_z1){

  b = (((1/shape13_1) - 1) * (sum(log(dat1$y13 ^ dat1$s13), na.rm = T)) - 
 scale13_1 * (shape13_1) * sum(dat1$y13 ^ (1/shape13_1) * exp(c13_1 * omega13_z1
 ), na.rm = T))
 
 return(b)
 }
