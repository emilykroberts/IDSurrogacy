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
like12_omega_i = function(omega12_z0, i, scale12_0, z, dat0, shape12_0, c12_0){

 e = ((1/shape12_0 - 1) * (sum(log(dat0$y12[i] ^ dat0$s12[i]), na.rm = T))- 
 scale12_0 * shape12_0 * sum(dat0$y12[i] ^ (1/shape12_0)* exp(c12_0 * omega12_z0)
 , na.rm = T))
 
 return(e)
}