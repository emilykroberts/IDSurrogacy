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
like13_omega1_i = function(omega13_z1, i, scale13_1, dat1, shape13_1, c13_1, 
 scale23_1, shape23_1, c23_1, theta23_1){

c = (c(c23_1, theta23_1) * (dat1$s23[i] %*% (cbind(omega13_z1, dat1$y12[i]))) -
 scale23_1/shape23_1 * sum(dat1$y23[i] ^ shape23_1 * exp(c23_1 * omega13_z1 +
 theta23_1 * dat1$y12[i]), na.rm = T))*dat1$s23[i] 
e = sum((c13_1 * omega13_z1)* dat1$s13[i], na.rm = T) *
 ((scale13_1* sum(dat1$y13[i] ^ (shape13_1) *exp(c13_1 * omega13_z1), na.rm = T)) )

d = sum(c, e, na.rm = T)

return(d)
}
