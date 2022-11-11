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
like12_shape1 = function(shape12_1, scale12_1, dat1, c12_1, omega12_z1){

 a =
 log ( (prod(dat1$y12 ^ dat1$s12, na.rm = T))^(shape12_1 - 1) *
 exp( -scale12_1/shape12_1 * sum(dat1$y12 ^ (shape12_1) * exp(c12_1 * omega12_z1 ), na.rm = T)))
 
 
return(a)
}
