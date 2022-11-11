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
like12_omega1_i = function(omega12_z1, i, scale12_1, dat1, shape12_1, c12_1){

 e = sum((c12_1 * omega12_z1)* dat1$s12[i], na.rm = T) *
 ((scale12_1* sum(dat1$y12[i] ^ (shape12_1) *exp(c12_1 * omega12_z1), na.rm = T)) )
 
 
 return(e)
}

