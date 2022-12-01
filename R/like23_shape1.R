#' likelihood function
#'
#' @description evaluate likelihood ...
#'
#' @param 
#' omega13_z1: frailty term
#' shape23_1: shape of Weibull distribution
#' scale23_1: scale of Weibull distribution
#' c23_1: frailty coefficient
#' theta23_1: covariate coefficient
#' dat1: data
#'
#' @return likelihood
#'
#' @examples 
#' example(like23_shape1(shape23_1 = 1, dat1 = 1, c23_1 = 1, scale23_1 = 1, theta23_1 = 1, omega13_z1 = 1))
like23_shape1 = function(shape23_1, dat1, c23_1, scale23_1, theta23_1, omega13_z1){

b = (((1/shape23_1) - 1) * (sum(log(dat1$y23 ^ dat1$s23), na.rm = T)) -
 scale23_1 * shape23_1 * sum(dat1$y23 ^ (1/shape23_1) * exp(c23_1 * omega13_z1
 ), na.rm = T))
 
 
return(b)

}
