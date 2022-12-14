#' likelihood function
#'
#' @description evaluate likelihood ...
#'
#' @param omega12_z1 frailty term
#' @param shape12_1 shape of Weibull distribution
#' @param scale12_1 scale of Weibull distribution
#' @param c12_1 frailty coefficient
#' @param dat1 data
#'
#' @return likelihood
#'
#' @examples 
#' example(like12_shape1(shape12_1 = 1, scale12_1 = 1, dat1 = 1,
#' c12_1 = 1, omega12_z1 = 1))
like12_shape1 = function(shape12_1, scale12_1, dat1, c12_1, omega12_z1){

b = ((1/shape12_1 - 1) * (sum(log(dat1$y12 ^ dat1$s12), na.rm = T)) - 
 scale12_1 * shape12_1 * sum(dat1$y12 ^ (1/shape12_1) * exp(c12_1 * omega12_z1
 ), na.rm = T))

return(b)
}
