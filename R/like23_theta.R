#' likelihood function
#'
#' @description evaluate likelihood ...
#'
#' @param omega13_z0 frailty term
#' @param shape23_0 shape of Weibull distribution
#' @param scale23_0 scale of Weibull distribution
#' @param c23_0 frailty coefficient
#' @param theta23_0 covariate coefficient
#' @param dat0 data
#'
#' @return likelihood
#'
#' @examples 
#' example(like23_theta(theta23_0 = 1, dat0 = 1, c23_0 = 1, 
#' shape23_0 = 1, scale23_0 = 1, omega13_z0 = 1))
like23_theta = function(theta23_0, dat0, c23_0, shape23_0, scale23_0, omega13_z0){
  a = (c(c23_0, theta23_0) * (dat0$s23 %*% (cbind(omega13_z0, dat0$y12))) -
         scale23_0 / shape23_0 * sum(dat0$y23 ^ shape23_0 * exp(c23_0 * omega13_z0 + 
                                                                  theta23_0 * dat0$y12), na.rm = T))[2]
  
  return(a)
}
