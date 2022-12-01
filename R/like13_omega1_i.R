#' likelihood function for frailty
#'
#' @description evaluate likelihood for frailty for 1-3 transition
#'
#' @param omega13_z1 frailty term
#' @param i index for individual in dataset
#' @param dat1 dataset
#' @param shape13_1 shape for Weibull distribution
#' @param scale13_1 scale for Weibull distribution
#' @param c13_1 frailty coefficient
#' @param shape23_1 shape for Weibull distribution
#' @param scale23_1 scale for Weibull distribution
#' @param c23_1 frailty coefficient
#' @param theta23_1 theta parameter for 2-3 transition
#'
#' @return likelihood
#'
#' @examples 
#' example(like13_omega1_i(omega13_z1 = 1, i = 1, scale13_1 = 1, 
#'  dat1 = dat1, shape13_1 = 1, c13_1 = 1, scale23_1 = 1, shape23_1 = 1, 
#' c23_1 = 1, theta23_1 = 0))
like13_omega1_i = function(omega13_z1, i, scale13_1, dat1, shape13_1, c13_1, 
 scale23_1, shape23_1, c23_1, theta23_1){

c = (((1/shape23_1) - 1) * (sum(log(dat1$y23[i] ^ dat1$s23[i]), na.rm = T)) -
 scale23_1 * shape23_1 * sum(dat1$y23[i] ^ (1/shape23_1) * exp(c23_1 * omega13_z1
 ), na.rm = T))
 
 d = (((1/shape13_1) - 1) * (sum(log(dat1$y13[i] ^ dat1$s13[i]), na.rm = T)) - 
 scale13_1 * (shape13_1) * sum(dat1$y13[i] ^ (1/shape13_1) * exp(c13_1 * omega13_z1
 ), na.rm = T))
 
 
d = sum(c, d, na.rm = T)

return(d)
}
