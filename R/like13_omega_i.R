#' likelihood function for frailty
#'
#' @description evaluate likelihood for frailty for 1-3 transition
#'
#' @param omega13_z0 frailty term
#' @param i index for individual in dataset
#' @param dat0 dataset
#' @param shape13_0 shape for Weibull distribution
#' @param scale13_0 scale for Weibull distribution
#' @param c13_0 frailty coefficient
#' @param shape23_0 shape for Weibull distribution
#' @param scale23_0 scale for Weibull distribution
#' @param c23_0 frailty coefficient
#' @param theta23_0 theta parameter for 2-3 transition
#'
#' @return likelihood
#'
#' @examples 
#' example(like13_omega_i(omega13_z0 = 1, i = 1, scale13_0 = 1, 
#' dat1 = dat1, shape13_0 = 1, c13_0 = 1, scale23_0 = 1, shape23_0 = 1, 
#' c23_0 = 1, theta23_0 = 0))
like13_omega_i = function(omega13_z0, i, scale13_0, dat0, shape13_0, c13_0, 
                          scale23_0, shape23_0, c23_0, theta23_0){
  
  c =  (((1 / shape13_0) - 1) * (sum(log(dat0$y13[i] ^ dat0$s13), na.rm = T)) - 
          scale13_0 * shape13_0 * sum(dat0$y13[i] ^ (1/shape13_0) * exp(c13_0 * omega13_z0
          ), na.rm = T))
  
  d = ((1 / shape23_0 - 1) * (sum(log(dat0$y23[i] ^ dat0$s23[i]), na.rm = T)) -
         scale23_0 * shape23_0 * sum(dat0$y23[i] ^ (1 / shape23_0) * exp(c23_0 * omega13_z0
         ), na.rm = T))
  
  e = sum(c, d, na.rm = T)
  
  return(e)
  
}
