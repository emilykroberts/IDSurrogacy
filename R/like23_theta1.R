#' likelihood function
#'
#' @description evaluate likelihood ...
#'
#' @param omega13_z1 frailty term
#' @param shape23_1 shape of Weibull distribution
#' @param scale23_1 scale of Weibull distribution
#' @param c23_1 frailty coefficient
#' @param theta23_1 covariate coefficient
#' @param dat1 data
#'
#' @return likelihood
#'
#' @examples 
#' example(like23_theta1(theta23_1 = 1, dat1 = 1, c23_1 = 1, 
#' scale23_1 = 1, shape23_1 = 1, omega13_z1 = 1))
like23_theta1 = function(theta23_1, dat1, c23_1, scale23_1, shape23_1, omega13_z1){

 c = (c(c23_1, theta23_1) * (dat1$s23 %*% (cbind(omega13_z1, dat1$y12))) -
 scale23_1 * shape23_1 * sum(dat1$y23 ^ (1/shape23_1) * exp(c23_1 * omega13_z1 + 
 theta23_1 * dat1$y12), na.rm = T))[2]
 return(c)
}
