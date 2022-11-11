like13_omega_i = function(omega13_z0, i, scale13_0, dat0, shape13_0, c13_0, 
                          scale23_0, shape23_0, c23_0, theta23_0){
 
 c  =   sum((c13_0 * omega13_z0)* dat0$s13[i], na.rm = T) *
            ((scale13_0* sum(dat0$y13[i] ^ (shape13_0) *exp(c13_0 * omega13_z0), na.rm = T)) )
 
 d  =    (c(c23_0,  theta23_0) * (dat0$s23[i] %*% (cbind(omega13_z0, dat0$y12[i]))) -
            scale23_0/shape23_0 * sum(dat0$y23[i] ^ shape23_0  * exp(c23_0 * omega13_z0 +
                           theta23_0 * dat0$y12[i]), na.rm = T))*dat0$s23[i]  
e = sum(c, d, na.rm = T)
 
 return(e)
  
}
