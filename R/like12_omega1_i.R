#' likelihood function
#'
#' @description evaluate likelihood ...
#'
#' @param weibull model parameters
#'
#' @return likelihood
#'
#' @examples 
#' example(like12_omega1_i(omega12_z1 = 1, i = 1, scale12_1 = 1, dat1 = dat1, shape12_1 = 1, c12_1 = 1))
like12_omega1_i = function(omega12_z1, i, scale12_1, dat1, shape12_1, c12_1){

 e = ((1/shape12_1 - 1) * (sum(log(dat1$y12[i] ^ dat1$s12[i]), na.rm = T)) - 
 scale12_1 * shape12_1 * sum(dat1$y12[i] ^ (1/shape12_1) * exp(c12_1 * omega12_z1
 ), na.rm = T))


 return(e)
}

