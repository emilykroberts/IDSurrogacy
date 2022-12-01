#' likelihood function
#'
#' @description evaluate likelihood ...
#'
#' @param weibull model parameters
#' omega12_z0: frailty term
#' shape12_0: shape of Weibull distribution
#' scale12_0: scale of Weibull distribution
#' c12_0: frailty coefficient
#' dat0: data
#'
#' @return likelihood
#'
#' @examples 
#' example(like12_shape(shape12_0 = 1, scale12_0 = 1, dat0 = 1, c12_0 = 1, omega12_z0 = 1))
like12_shape = function(shape12_0, scale12_0, dat0, c12_0, omega12_z0){

b = ((1/shape12_0 - 1) * (sum(log(dat0$y12 ^ dat0$s12), na.rm = T))- 
 scale12_0 * shape12_0 * sum(dat0$y12 ^ (1/shape12_0)* exp(c12_0 * omega12_z0)
 , na.rm = T))

return(b)
}
