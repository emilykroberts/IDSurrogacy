like23_theta1 = function(theta23_1, dat1, c23_1, scale23_1, shape23_1, omega13_z1){
 
  c =   (c(c23_1, theta23_1) * (dat1$s23 %*% (cbind(omega13_z1, dat1$y12))) -
           scale23_1/shape23_1 * sum(dat1$y23 ^ shape23_1  * exp(c23_1 * omega13_z1 + 
                                                                                    theta23_1 * dat1$y12), na.rm = T))[2]
  return(c)
}
