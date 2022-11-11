like12_shape = function(shape12_0, scale12_0, dat0, c12_0, omega12_z0){
 b =   ((shape12_0 - 1) * (sum(log(dat0$y12 ^ dat0$s12), na.rm = T))- 
          scale12_0/shape12_0 * sum(dat0$y12 ^ (shape12_0)* exp(c12_0 * omega12_z0)
                                        , na.rm = T))
 
return(b)
 }
