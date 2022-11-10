library(MASS); library(mvtnorm); library(MCMCpack); library(LaplacesDemon)
library(ggplot2); library(survival); library(survminer); library(frailtyEM); library(frailtypack)
library(truncnorm); library(latex2exp); library(kableExtra); library(Rfast)
library(xtable); library(ggforce); library(wesanderson); library(cowplot)

array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(is.na(array_id)) array_id = 1
{n = 600
rhost = rhos = .5; rhot = .5
SIM = 1000
holdtheta = F
holdscale12 = F
holdscale23 = F
holdscale13 = F
holdfrail12 = F
holdfrail13 = F
equalfrail = T
independent = T
diffscale1323 = T
holdshape = F
holdc = T
effecttheta = F

tau_s = 1
tau_t = 2
frailtysd = 1
scenario = 2

proposalsd = 0.02
proposalsdtheta = 0.04
proposalsdfrail = 0.005

write = T
plotwrite = as.numeric(array_id == 10)
}

sim_data = function(n, array_id, scenario, effecttheta, frailtysd, rhot, rhos, rhost, effectsize){

 {
  if(scenario == 1){ effect12 = F; effect13 = F; effect23 = F}
  if(scenario == 2){ effect12 = T; effect13 = F; effect23 = F}
  if(scenario == 3){ effect12 = T; effect13 = F; effect23 = T}
  if(scenario == 4){ effect12 = T; effect13 = T; effect23 = F}
  if(scenario == 5){ effect12 = T; effect13 = T; effect23 = T}
  if(scenario == 6){ effect12 = F; effect13 = T; effect23 = T}
  if(scenario == 7){ effect12 = F; effect13 = F; effect23 = T}
  if(scenario == 8){ effect12 = F; effect13 = T; effect23 = F}
 }
 
 beta12_1 = 0; beta13_1 = 0; beta23_1 = 0; beta12_0 = 0; beta13_0 = 0; beta23_0 = 0
 theta23_1 = 0; theta23_0 = 0

 if(effecttheta){theta23_1 = -0.1; theta23_0 = -0.1
 }

 c13_0 = c13_1 = 1
 c12_1 = c12_0 = 1
 c23_1 = c23_0 = 1
 
 scale12_1 = 1 
 shape12_1 = 1 
 scale13_1 = 0.5 
 if(!diffscale1323){scale13_1 = 1} 
 shape13_1 = 1 
 scale23_1 = 1 
 shape23_1 = 1 
 shape12_0 = shape12_1
 scale12_0 = scale12_1 
 scale13_0 = scale13_1 
 shape13_0 = shape13_1 
 scale23_0 = scale23_1 
 shape23_0 = shape23_1 
 
 if(effect12){ scale12_1 = effectsize }
 if(effect13){ scale13_1 = scale13_1 * (effectsize) }
 if(effect23){ scale23_1 = effectsize }
 
 # generate data
 R=matrix(rep(1,6*6),6,6); 
 R[1,2] = R[2,1] = rhos
 R[3,1] = R[1,3] = 0.6
 R[1,4] = R[4,1] = 0
 R[1,5] = R[5,1] = 0
 R[1,6] = R[6,1] = 0
 
 R[2,3] = R[3,2] = 0
 R[2,4] = R[4,2] = 0
 R[2,5] = R[5,2] = 0
 R[2,6] = R[6,2] = 0
 
 R[3,4] = R[4,3] = rhot
 R[3,5] = R[5,3] = 0
 R[3,6] = R[6,3] = 0
 
 R[4,5] = R[5,4] = 0
 R[4,6] = R[6,4] = 0
 
 R[5,6] = R[6,5] = rhost
 
 if(independent){ R[1,2] = R[2,1] = rhos
 R[3,1] = R[1,3] = 0
 R[1,4] = R[4,1] = 0
 R[1,5] = R[5,1] = 0
 R[1,6] = R[6,1] = 0
 
 R[2,3] = R[3,2] = 0
 R[2,4] = R[4,2] = 0
 R[2,5] = R[5,2] = 0
 R[2,6] = R[6,2] = 0
 
 R[3,4] = R[4,3] = rhot
 R[3,5] = R[5,3] = 0
 R[3,6] = R[6,3] = 0
 
 R[4,5] = R[5,4] = 0
 R[4,6] = R[6,4] = 0
 
 R[5,6] = R[6,5] = rhost
 }
 
 S = diag(c(frailtysd, frailtysd, frailtysd, frailtysd, frailtysd, frailtysd))
 o = mvtnorm::rmvnorm(n, c(0,0,0,0,0,0), S%*%R%*%S)
 omega12true0 = o[,1]; omega12true1 = o[,2]
 omega13true0 = o[,3]; omega13true1 = o[,4]
 omega23true0 = o[,5]; omega23true1 = o[,6]
 
 o = mvtnorm::rmvnorm(n, c(0,0,0,0,0,0), S%*%R%*%S)
 
 if(equalfrail){omega23true0 = omega13true0; omega23true1 = omega13true1}
 
 U1 = runif(n); U2 = runif(n); U3 = runif(n)
 U4 = runif(n); U5 = runif(n); U6 = runif(n)
 xtrue = x = rbinom(n, 1, 0.5)
 x_0 = x[1:(n/2)]; x_1 = x[1:(n/2)]
 
 ST = status = ST_raw = cbind(
  (1/scale12_0*(-log(1-U1))/exp((c12_0*omega12true0 + beta12_0 * x)))^shape12_0,
  (1/scale23_0*(-log(1-U2))/exp((c23_0*omega23true0 + beta23_0 * x)))^shape23_0,
  (1/scale13_0*(-log(1-U3))/exp((c13_0*omega13true0 + beta13_0 * x)))^shape13_0,
  (1/scale12_1*(-log(1-U4))/exp((c12_1*omega12true1 + beta12_1 * x)))^shape12_1,
  (1/scale23_1*(-log(1-U5))/exp((c23_1*omega23true1 + beta23_1 * x)))^shape23_1,
  (1/scale13_1*(-log(1-U6))/exp((c13_1*omega13true1 + beta13_1 * x)))^shape13_1
 ); status[,1:6] = 1
 
 
 ST[,2] =  (1/scale23_0*(-log(1-U2))/exp((c23_0*omega23true0 + theta23_0 * ST[,1])))^shape23_0
 ST[,5] =  (1/scale23_1*(-log(1-U5))/exp((c23_1*omega23true1 + theta23_1 * ST[,4])))^shape23_1
 
 ST[is.infinite(ST)] = NA
 ST[(ST) == 0]= NA


if(T){
 status[ST[,6] < ST[,4], 4] = 0
 status[ST[,3] < ST[,1], 1] = 0
 
 ST[ST[,3] < ST[,1], 1] = ST[ST[,3] < ST[,1], 3]
 ST[ST[,6] < ST[,4], 4] = ST[ST[,6] < ST[,4], 6]

 status[ST[,6] <= ST[,4], 5] = 0
 status[ST[,3] <= ST[,1], 2] = 0
 ST[ST[,3] <= ST[,1], 2] = NA
 ST[ST[,6] <= ST[,4], 5] = NA
 
 status[ST[,3] > ST[,1], 3] = 0 
 status[ST[,6] > ST[,4], 6] = 0
 ST[ST[,3] > ST[,1], 3] = ST[ST[,3] > ST[,1], 1]
 ST[ST[,6] > ST[,4], 6] = ST[ST[,6] > ST[,4], 4]
}

trt = c(rep(0, n), rep(1, n))

dat0 = data.frame(y12 = c(ST[1:(n/2),1]), s12 = c(status[1:(n/2),1]),
         y13 = c(ST[1:(n/2),3]), s13 = c(status[1:(n/2),3]),
         y23 = c(ST[1:(n/2),2]), s23 = c(status[1:(n/2),2])
)

dat1 = data.frame(y12 = c(ST[1:(n/2),4]), s12 = c(status[1:(n/2),4]),
         y13 = c(ST[1:(n/2),6]), s13 = c(status[1:(n/2),6]),
         y23 = c(ST[1:(n/2),5]), s23 = c(status[1:(n/2),5])
)

o12save0 = omega12_z0 = omega12true0[1:(n/2)]
o13save0 = omega13_z0 = omega13true0[1:(n/2)]
o23save0 = omega23_z0 = omega23true0[1:(n/2)]

o12save1 = omega12_z1 = omega12true1[1:(n/2)]
o13save1 = omega13_z1 = omega13true1[1:(n/2)]
o23save1 = omega23_z1 = omega23true1[1:(n/2)]

params_list = data.frame(scale12_0 = scale12_0, scale13_0 = scale13_0, scale23_0 = scale23_0,
          scale12_1 = scale12_1, scale13_1 = scale13_1, scale23_1 = scale23_1,
          shape12_0 = shape12_0, shape13_0 = shape13_0, shape23_0 = shape23_0,
          shape12_1 = shape12_1, shape13_1 = shape13_1, shape23_1 = shape23_1,
          theta23_0 = theta23_0, theta23_1 = theta23_1,
          beta12_0 = beta12_0, beta13_0 = beta13_0, beta23_0 = beta23_0, 
          beta12_1 = beta12_1, beta13_1 = beta13_1, beta23_1 = beta23_1,
          c12_0 = c12_0, c13_0 = c13_0, c23_0 = c23_0, c12_1 = c12_1, c13_1 = c13_1, c23_1 = c23_1
          )
 
return(list( dat0 = dat0, dat1 = dat1, o12save0 = omega12true0, o12save1 = omega12true1, o13save0 = omega13true0, 
       o13save1 = omega13true1, o23save0 = omega23true0, o23save1 = omega23true1, params = params_list))
}

set.seed(1 + array_id)

dat = sim_data(n = n, array_id = array_id, scenario = scenario, effecttheta = effecttheta,
        rhos = rhos, rhot = rhot, rhost = rhost, frailtysd = frailtysd, effectsize = 0.61)

dat0 = dat$dat0
dat1 = dat$dat1
omega12true0 = o12save0 = dat$o12save0
omega12true1 = o12save1 = dat$o12save1
omega13true0 = o13save0 = dat$o13save0
omega13true1 = o13save1 = dat$o13save1
omega23true0 = o23save0 = dat$o23save0
omega23true1 = o23save1 = dat$o23save1

true_params = dat$params

## functions
{
Lambda23_frailty=function(x,xdata1,v_predict, omega2, shape23, scale23, theta23, c23, beta23_1) {
 return(-pweibull(x, scale= scale23, shape = shape23, lower = F, log = T)*exp(beta23_1*xdata1+theta23*v_predict + c23*omega2))}

Lambda12_frailty_tx=function(x,xdata1, omega1, shape12, scale12, c12, beta12_1, betatx12) {
 return( -pweibull(x, scale= scale12, shape= shape12, lower = F, log = T)*exp(beta12_1*xdata1 + c12*omega1 + betatx12))}

Lambda13_frailty=function(x,xdata1, omega2, scale13, shape13, c13, beta13_1) {
 return(-pweibull(x, scale= scale13, shape = shape13, lower = F, log = T)*exp(beta13_1*xdata1 + c13*omega2))}

lambda23_frailty=function(x,xdata1,v_predict, omega2, shape23, scale23, c23, theta23, beta23_1) {
 return(shape23*(1/scale23)^(shape23)*(x)^(shape23-1)*exp(beta23_1*xdata1+theta23*v_predict + c23*omega2))}


lambda12_frailty_tx=function(x,xdata1, omega1, scale12, shape12, c12, beta12_1, betatx12) {
 return(shape12*(scale12)^(shape12)*(x)^(shape12-1)*exp(beta12_1*xdata1 + c12*omega1 + betatx12))}

lambda13_frailty=function(x,xdata1, omega2, scale13, shape13, c13, beta13_1) {
 return(shape13*(1/scale13)^(shape13)*(x)^(shape13-1)*exp(beta13_1*xdata1 + c13*omega2))}

Lambda23=function(x,xdata1,v_predict, omega2, shape23, scale23, theta23, c23, beta23_1) {
 return(-pweibull(x, scale= scale23, shape = shape23, lower = F, log = T))}

Lambda12=function(x,xdata1, omega1, shape12, scale12, c12, beta12_1) {
 return( -pweibull(x, scale= scale12, shape = shape12, lower = F, log = T))}

Lambda13=function(x,xdata1, omega2, scale13, shape13, c13, beta13_1) {
 return(-pweibull(x, scale= scale13, shape = shape13, lower = F, log = T))}

lambda23=function(x,xdata1,v_predict, omega2, shape23, scale23, c23, theta23, beta23_1) {
 return(shape23*scale23*(x*scale23)^(shape23-1))}

lambda12=function(x,xdata1, omega1, scale12, c12, beta12_1) {
 return(shape12*scale12*(x*scale12)^(shape12-1))}

lambda13=function(x,xdata1, omega2, scale13, c13, beta13_1) {
 return(shape13*scale13*(x*scale13)^(shape13-1))}

Lambda13_frailty_lk_tx=function(x,xdata1, omega2, scale13, shape13, c13, beta13_1, betatx13) {
 return(-pweibull(x, scale= 1/scale13, shape = shape13, lower = F, log = T)*exp(beta13_1*xdata1 + c13*omega2 + betatx13))}

Lambda23_frailty_lk_tx=function(x,xdata1,v_predict, omega3, shape23, scale23, theta23, c23, beta23_1, betatx23) {
 return(-pweibull(x, scale= 1/scale23, shape = shape23, lower = F, log = T)*exp((beta23_1*xdata1+c23*omega3+theta23*v_predict+ betatx23)))}

Lambda12_frailty_lk_tx=function(x,xdata1, omega1, shape12, scale12, c12, beta12_1, betatx12) {
 return( -pweibull(x, scale= 1/scale12, shape= shape12, lower = F, log = T)*exp(beta12_1*xdata1 + c12*omega1 + betatx12)
 )}

Lambda13_frailty_lk=function(x,xdata1, omega2, scale13, shape13, c13, beta13_1) {
 return(-pweibull(x, scale= 1/scale13, shape = shape13, lower = F, log = T)*exp(beta13_1*xdata1 + c13*omega2))}

Lambda23_frailty_lk=function(x,xdata1,v_predict, omega2, shape23, scale23, omega1, theta23, c23, beta23_1) {
 return(-pweibull(x, scale= 1/scale23, shape = shape23, lower = F, log = T)*exp(-(beta23_1*xdata1+c23*omega2+theta23*v_predict)))}

Lambda12_frailty_lk=function(x,xdata1, omega1, shape12, scale12, c12, beta12_1) {
 return( -pweibull(x, scale= 1/scale12, shape= shape12, lower = F, log = T)*exp(beta12_1*xdata1 + c12*omega1))}

lambda12_frailty=function(x,xdata1, omega1, scale12, shape12, c12, beta12_1) {
 return(shape12*(scale12)^(shape12)*(x)^(shape12-1)*exp(beta12_1*xdata1 + c12*omega1))}


like12_shape1 = function(shape12_1, scale12_1, dat1, c12_1, omega12_z1){
 b =    ((shape12_1 - 1) * (sum(log(dat1$y12 ^ dat1$s12), na.rm = T)) - 
        scale12_1/shape12_1  * sum(dat1$y12 ^ (shape12_1) * exp(c12_1 * omega12_z1
    ), na.rm = T))

 a =
  log ( (prod(dat1$y12 ^ dat1$s12, na.rm = T))^(shape12_1 - 1) *
       exp( -scale12_1/shape12_1  * sum(dat1$y12 ^ (shape12_1) * exp(c12_1 * omega12_z1 ), na.rm = T)))
 
 
return(a)
}

like12_shape = function(shape12_0, scale12_0, dat0, c12_0, omega12_z0){
 b =  ((shape12_0 - 1) * (sum(log(dat0$y12 ^ dat0$s12), na.rm = T))- 
     scale12_0/shape12_0 * sum(dat0$y12 ^ (shape12_0)* exp(c12_0 * omega12_z0)
                    , na.rm = T))

return(b)
 }


like23_theta = function(theta23_0, dat0, c23_0, shape23_0, scale23_0, omega13_z0){
 a =  (c(c23_0, theta23_0) * (dat0$s23 %*% (cbind(omega13_z0, dat0$y12))) -
     scale23_0/shape23_0 * sum(dat0$y23 ^ shape23_0 * exp(c23_0 * omega13_z0 + 
                                       theta23_0 * dat0$y12), na.rm = T))[2]

return(a)
}

like23_shape = function(shape23_0, dat0, c23_0, scale23_0, theta23_0, omega13_z0){
 b = ((shape23_0 - 1) * (sum(log(dat0$y23 ^ dat0$s23), na.rm = T)) - 
     scale23_0 / shape23_0 * sum(dat0$y23 ^ shape23_0 * exp(c23_0 * omega13_z0 + 
                                     theta23_0 * dat0$y12
     ), na.rm = T))
 return(b)
 }


like23_shape1 = function(shape23_1, dat1, c23_1, scale23_1, theta23_1, omega13_z1){
 b =  ((shape23_1 - 1) * (sum(log(dat1$y23 ^ dat1$s23), na.rm = T)) - 
      scale23_1 / shape23_1 * sum(dat1$y23 ^ shape23_1 * exp(c23_1 * omega13_z1 + 
                                     theta23_1 * dat1$y12
      ), na.rm = T))
 a = log ( (prod(dat1$y23 ^ dat1$s23, na.rm = T))^(shape23_1 - 1) *
      exp( -scale23_1/shape23_1  * sum(dat1$y23 ^ (shape23_1) * exp(c23_1 * omega13_z1 + theta23_1 * dat1$y12), na.rm = T)))
return(b)
 }

like13_shape = function(shape13_0, dat0, c13_0, scale13_0, omega13_z0){
b =  ((shape13_0 - 1) * (sum(log(dat0$y13 ^ dat0$s13), na.rm = T)) - 
  scale13_0 / shape13_0 * sum(dat0$y13 ^ shape13_0 * exp(c13_0 * omega13_z0
   ), na.rm = T))


return(b)
}

like13_shape1 = function(shape13_1, dat1, c13_1, scale13_1, omega13_z1){
 b =  ((shape13_1 - 1) * (sum(log(dat1$y13 ^ dat1$s13), na.rm = T)) - 
      scale13_1 / shape13_1 * sum(dat1$y13 ^ shape13_1 * exp(c13_1 * omega13_z1
      ), na.rm = T))
 
 return(b)
 }


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

like13_omega_i = function(omega13_z0, i, scale13_0, dat0, shape13_0, c13_0, 
             scale23_0, shape23_0, c23_0, theta23_0){

 c =  sum((c13_0 * omega13_z0)* dat0$s13[i], na.rm = T) *
      ((scale13_0* sum(dat0$y13[i] ^ (shape13_0) *exp(c13_0 * omega13_z0), na.rm = T)) )

 d =  (c(c23_0, theta23_0) * (dat0$s23[i] %*% (cbind(omega13_z0, dat0$y12[i]))) -
      scale23_0/shape23_0 * sum(dat0$y23[i] ^ shape23_0 * exp(c23_0 * omega13_z0 +
              theta23_0 * dat0$y12[i]), na.rm = T))*dat0$s23[i] 
e = sum(c, d, na.rm = T)
 
 return(e)
 
}

like23_theta1 = function(theta23_1, dat1, c23_1, scale23_1, shape23_1, omega13_z1){

 c =  (c(c23_1, theta23_1) * (dat1$s23 %*% (cbind(omega13_z1, dat1$y12))) -
      scale23_1/shape23_1 * sum(dat1$y23 ^ shape23_1 * exp(c23_1 * omega13_z1 + 
                                          theta23_1 * dat1$y12), na.rm = T))[2]
 return(c)
}


like12_omega_i = function(omega12_z0, i, scale12_0, z, dat0, shape12_0, c12_0){

 e = sum((c12_0 * omega12_z0)* dat0$s12[i], na.rm = T) *
  ((scale12_0* sum(dat0$y12[i] ^ (shape12_0) *exp(c12_0 * omega12_z0), na.rm = T)) )
 return(e)
}



like12_omega1_i = function(omega12_z1, i, scale12_1, dat1, shape12_1, c12_1){

   e = sum((c12_1 * omega12_z1)* dat1$s12[i], na.rm = T) *
   ((scale12_1* sum(dat1$y12[i] ^ (shape12_1) *exp(c12_1 * omega12_z1), na.rm = T)) )
   
   
   return(e)
}

}

true_cep = function(dat0, dat1, write, params_list, plotwrite){
 
 params_list = attach(params_list)
 Fsave1 = Fsave = Fw_0 = Fw = Fw_1 = NULL
 x = rep(0, n)
 
 intfunction = function(j, i, t){
  exp(
   -Lambda13_frailty_lk_tx(x = t, xdata = x[i], omega2 = omega13true0[i], scale13 = scale13_0, 
              shape13 = shape13_0,  c13 = c13_0, beta13_1 = beta13_0, betatx13 = 0) - 
    Lambda12_frailty_lk_tx(x = t, xdata = x[i], 
              omega1 = omega12true0[i],  scale12 = scale12_0, shape12 = shape12_0, c12 = c12_0, beta12_1 = beta12_0, betatx12 = 0)) *
   lambda12_frailty_tx(t, xdata = x[i], omega1 = omega12true0[i], scale12 = scale12_0, shape12 = shape12_0,
            c12 = c12_0, beta12_1 = beta12_0, betatx12 = 0)* 
   exp(-Lambda23_frailty_lk_tx(x = tau_t - t, xdata = x[i], omega3 = omega23true0[i], scale23 = scale23_0, 
                shape23 = shape23_0, betatx23 = 0,
                c23 = c23_0, theta23 = theta23_0, v_predict = t, beta23_1 = beta23_0))
 }
 
 
 
 for(i in 1:n){
  y =tryCatch({ integrate(f = intfunction, lower = 0, upper = tau_t, i = i, j = j)$val
  },
  error = function(error_condition) {
   NA}
  ) 
  
  L12 = Lambda12_frailty_lk_tx(x = tau_t, xdata = x[i], shape12 = shape12_0, betatx12 = 0,
               omega1 = omega12true0[i],  scale12 = scale12_0, c12 = c12_0, beta12_1 = beta12_0)
  L13 = Lambda13_frailty_lk_tx(x = tau_t, xdata = x[i], omega2 = omega13true0[i], scale13 = scale13_0, shape13 = shape13_0,
               c13 = c13_0, beta13_1 = beta13_0, betatx13 = 0)
  
  Fw = (exp(-L13 - L12) + y)
  Fw_0 = c(Fw_0, Fw)
 }
 
 intfunction = function(j, i, t){
  exp(
   -Lambda13_frailty_lk_tx(x = t, xdata = x[i], omega2 = omega13true1[i], scale13 = scale13_1, 
              shape13 = shape13_1,  c13 = c13_1, beta13_1 = beta13_1, betatx13 = 0) - 
    Lambda12_frailty_lk_tx(x = t, xdata = x[i], 
              omega1 = omega12true1[i],  scale12 = scale12_1, shape12 = shape12_1, c12 = c12_1, beta12_1 = beta12_1, betatx12 = 0)) *
   lambda12_frailty_tx(t, xdata = x[i], omega1 = omega12true1[i], scale12 = scale12_1, shape12 = shape12_1,
            c12 = c12_1, beta12_1 = beta12_1, betatx12 = 0)* 
   exp(-Lambda23_frailty_lk_tx(x = tau_t - t, xdata = x[i], omega3 = omega23true1[i], scale23 = scale23_1, 
                shape23 = shape23_1, betatx23 = 0,
                c23 = c23_1, theta23 = theta23_1, v_predict = t, beta23_1 = beta23_1))
 }
 
 ## for z = 1
 for(i in 1:n){
  y =tryCatch({ integrate(f = intfunction, lower = 0, upper = tau_t, i = i, j = j)$val
  },
  error = function(error_condition) {  NA}
  ) 
  
  L12 = Lambda12_frailty_lk_tx(x = tau_t, xdata = x[i], shape12 = shape12_1, betatx12 = 0,
               omega1 = omega12true1[i],  scale12 = scale12_1, c12 = c12_1, beta12_1 = beta12_1)
  L13 = Lambda13_frailty_lk_tx(x = tau_t, xdata = x[i], omega2 = omega13true1[i], scale13 = scale13_1, shape13 = shape13_1,
               c13 = c13_1, beta13_1 = beta13_1, betatx13 = 0)
  
  Fw = (exp(-L13 - L12) + y)
  if(Fw>1) Fw = NA
  Fw_1 = c(Fw_1, Fw)
 }
 
 z = 1
 
 s0cumulative = -pweibull(q = tau_s, scale= 1/scale12_0, shape= shape12_0, lower = F, log = T)*exp(beta12_0*x + c12_0*omega12true0)
 s1cumulative = -pweibull(q = tau_s, scale= 1/scale12_1, shape= shape12_1, lower = F, log = T)*exp(beta12_1*x + c12_1*omega12true1 + 0)
 
 
 dat = data.frame(cbind(s0cumulative/s1cumulative, Fw_1 - Fw_0))
 dat$X = s0cumulative/s1cumulative
 dat$Y = c(Fw_1 - Fw_0)
 dat$x = x
 
 reg_true = reg = lm(formula = Y ~ log(X), data=dat) 
 reg_true
 
 
 fname = paste('cep','.n',n,array_id,'.txt',sep="")
 res = cbind(summary(reg)$coef[1,1], summary(reg)$coef[2,1])
 print(res)
 if(write) write.table(res, file=fname, sep="\t", row.names=F, col.names=T)

 d2 = ggplot(dat, aes(log(X), Y, alpha = 0.01)) + geom_point(alpha = 0.5) + theme_classic() + 
  ggtitle("Illness-Death CEP Curve with True Values" ) +
  theme(legend.position = "none") +
  xlab("Delta S_i")+ ylim(-1,1)  +
  ylab(TeX("$\\Delta T_i = P(T_i(1) > \\tau_T|x_i, \\omega_{.i}^1) - P(T_i(0) > \\tau_T |x_i, \\omega_{.i}^0)$$ ")) +
  xlab(TeX("$\\Delta S_i = log \\frac{\\Lambda_{12}^{0} (\\tau_S| x_i,\\omega_{12i}^0)}{\\Lambda_{12}^{1} (\\tau_S| x_i,\\omega_{12i}^1)}$$ ")) + 
  geom_hline(yintercept = 0, linetype = "solid") + geom_vline(xintercept = 0, linetype = "solid") + 
  geom_hline(yintercept = mean(dat$Y, na.rm = T), linetype = "dashed", col = 'red') + geom_vline(xintercept = mean(log(dat$X), na.rm = T), linetype = "dashed", col = 'red') + 
  
  theme(text = element_text(size = 16))+
  geom_smooth(method='lm', level = 0.95) + coord_cartesian(xlim = c(-5,5));     
 
 
 d2 = d2 + theme(axis.title.x = element_text(
  margin = margin(t = 5, r = 0, b = 50, l = 0,
          unit = "mm")))
 line_1 = expr(paste("Intercept ", !!c(round(summary(reg)$coef[1,1], 3)), " Slope ", !!c(round(summary(reg)$coef[2,1], 3)), " R\U00B2 = ", !!round(summary(reg_true)$r.squared, 3),
            sep=""))
 
 line_2 = expr(paste("Generating Parameters ", alpha[12]^0, "=", !!shape12_0, ", ",
            alpha[13]^0, "=", !!shape13_0, ", ",
            alpha[23]^0, "=", !!shape23_0, ", ",
            alpha[12]^1, "=", !!shape12_1, ", ",
            alpha[13]^1, "=", !!shape13_1, ", ",
            alpha[23]^1, "=", !!shape23_1, 
            sep=""))
 line_3 = expr(paste(gamma[12]^0, "=", !!scale12_0, ", ",
            gamma[13]^0, "=", !!scale13_0, ", ",
            gamma[23]^0, "=", !!scale23_0, ", ",
            gamma[12]^1, "=", !!scale12_1, ", ",
            gamma[13]^1, "=", !!scale13_1, ", ",
            gamma[23]^1, "=", !!scale23_1, 
            sep=""))
 line_4 = expr(paste(kappa[12]^0, "=", !!c12_0, ", ",
            kappa[13]^0, "=", !!c13_0, ", ",
            kappa[23]^0, "=", !!c23_0, ", ",
            kappa[12]^1, "=", !!c12_1, ", ",
            kappa[13]^1, "=", !!c13_1, ", ",
            kappa[23]^1, "=", !!c23_1, 
            sep=""))
 
 line_5 = expr(paste(  theta[23]^0, "=", !!theta23_0, ", ",
              theta[23]^1, "=", !!theta23_1, ", ",
              rho[s], "=", !!rhos, ", ",
              rho[t], "=", !!rhot, ", ",
              tau[s], "=", !!tau_s, ", ",
              tau[t], "=", !!tau_t,
              sep=""))
 
 d3 = ggdraw(d2)+
  draw_label(line_1, x = 0.55, y = 0.215) + 
  draw_label(line_2, x = 0.55, y= .175)+
  draw_label(line_3, x = 0.55, y= .135)+
  draw_label(line_4, x = 0.55, y = 0.095)+
  draw_label(line_5, x = 0.55, y = 0.055)
 print(d3)
 if(plotwrite) ggsave(paste0("CEPtrue", ",", scenario, ",", array_id,  ".jpeg"), d3)
 
}

true_cep(dat0 = dat0, dat1 = dat1, write = write, params_list = true_params, plotwrite = plotwrite)


params_list = true_params

run_sim = function(SIM, rhos, rhot, frailtysd, params_list, dat0, dat1, n){
 
 burnin = SIM * .3
 
 ## save parameters
 {accept1 = accept2 = accept3 = accept4 = accept5 = accept6 = accept7 = accept8 = accept9 = accept10 = accept11 = accept12 = 0
 holdint = holdslope = holdscale12_0 = holdshape12_0 = holdscale13_0 = holdshape13_0 = holdscale23_0 = holdshape23_0 = rep(NA, SIM)
 holdscale12_1 = holdshape12_1 = holdscale13_1 = holdshape13_1 = holdscale23_1 = holdshape23_1 = rep(NA, SIM)
 holdbeta12_1 = holdc12_1 = holdbeta13_1 = holdc13_1 = holdbeta23_1 = holdc23_1 = rep(NA, SIM)
 holdtheta23_1 = holdtheta23_0 = holdbeta12_0 = holdc12_0 = holdbeta13_0 = holdc13_0 = holdbeta23_0 = holdc23_0 = rep(NA, SIM)
 holdfrailsd23_1 = holdfrailsd23_0 = holdfrailmean23_0 = holdfrailmean23_1 = rep(NA, SIM)
 holdfrailsd13_1 = holdfrailsd13_0 = holdfrailsd12_1 = holdfrailsd12_0 = rep(NA, SIM)
 holdfrailmean13_1 = holdfrailmean13_0 = holdfrailmean12_1 = holdfrailmean12_0 = rep(NA, SIM)
 holdomega12_0 = matrix(NA, nrow = 10, ncol = SIM)
 holdomega23_0 = holdomega23_1 = holdomega13_0 = matrix(NA, nrow = 10, ncol = SIM)
 holdomega12_1 = matrix(NA, nrow = 10, ncol = SIM)
 holdomega13_1 = matrix(NA, nrow = 10, ncol = SIM)
 
 holdscale12_0[1] = 1
 holdshape12_0[1] = 1
 holdshape13_0[1] = 1
 holdscale23_0[1] = 1
 holdshape23_0[1] = 1
 holdshape12_1[1] = 1
 holdshape23_1[1] = 1
 holdshape13_1[1] = 1
 holdscale23_1[1] = 1
 holdscale12_1[1] = 1

 holdc12_0[1] = holdc13_0[1] = holdc12_1[1] = holdc13_1[1] = holdc23_1[1] = holdc23_0[1] = 1
 holdscale13_0[1] = 1
 holdscale13_1[1] = 1
 holdtheta23_1[1] = 0
 holdtheta23_0[1] = 0
 
 holdscale12_0[2] = holdshape12_0[2] = holdshape12_1[2] = holdscale12_1[2] = 1
 holdbeta23_1[1] = holdbeta13_1[1] = holdbeta12_0[1] = holdbeta12_1[1] = holdbeta13_0[1] = holdbeta23_0[1] = 0
 
 saveCEPx = matrix(data = NA, nrow = n, ncol = SIM )
 saveCEPy = matrix(data = NA, nrow = n, ncol = SIM )
 acceptfrail0 = acceptfrail1 = matrix(data = 0, nrow = 100, ncol = SIM)
 }
 
 dat0$y23[dat0$y23 == 0] = NA
 dat1$y23[dat1$y23 == 0] = NA
 
 if(T){
  
 f0 = emfrail(formula = Surv(c(dat0$y12), c(dat0$s12)) ~ cluster(rep(1:(n/2), 1)), data = dat0)
 #f0_2 = emfrail(formula = Surv(c(dat0$y23, dat0$y13), c(dat0$s23, dat0$s13)) ~ cluster(rep(1:(n/2), 2)), data = dat0)
 f0_2 = emfrail(formula = Surv(c(dat0$y23, dat0$y13), c(dat0$s23, dat0$s13)) ~ cluster(rep(1:(n/2), 2)) + 
          c(dat0$y12, rep(1, n/2)), data = dat0)
 
 f0_23 = emfrail(formula = Surv(c(dat0$y23), c(dat0$s23)) ~ cluster(rep(1:(n/2), 1)) + 
          c(dat0$y12), data = dat0,
         distribution = emfrail_dist(dist = "stable"))
 f0_23 = f0_2
 
 f1 = emfrail(formula = Surv(c(dat1$y12), c(dat1$s12)) ~ cluster(rep(1:(n/2), 1)), data = dat1)
 #f1_2 = emfrail(formula = Surv(c(dat1$y23, dat1$y13), c(dat1$s23, dat1$s13)) ~ cluster(rep(1:(n/2), 2)), data = dat1)
 f1_2 = emfrail(formula = Surv(c(dat1$y23, dat1$y13), c(dat1$s23, dat1$s13)) ~ cluster(rep(1:(n/2), 2)) + 
          c(dat1$y12, rep(1, n/2)), data = dat1)
 
 f1_23 = emfrail(formula = Surv(c(dat1$y23), c(dat1$s23)) ~ cluster(rep(1:(n/2), 1)) + 
          c(dat1$y12), data = dat1,
         distribution = emfrail_dist(dist = "stable"))
  f1_23 = f1_2
 
 factor = frailtysd/sd(log(f0$frail))
 o12save0[1:(n/2)] = omega12_z0 = log(f0$frail)#*factor
 o13save0[1:(n/2)] = omega13_z0 = log(f0_2$frail)#*factor
 o23save0[1:(n/2)] = omega23_z0 = log(f0_23$frail)#*factor
 
 #o23save0[1:(n/2)] = omega23_z0 = log(f0_23$frail)#*factor

 factor = frailtysd/sd(f1$frail)
 o12save1[1:(n/2)] = omega12_z1 = log(f1$frail)#*factor
 o13save1[1:(n/2)] = omega13_z1 = log(f1_2$frail)#*factor
 o23save1[1:(n/2)] = omega23_z1 = log(f1_23$frail)#*factor
 }
 
 ### option 2
 if(T){
  # dat1penal = data.frame(time = c(dat1$y23, dat1$y13, dat1$y12),
  #            status = c(dat1$s23, dat1$s13, dat1$s12),
  #            id = rep(1:(n/2), 3))
  
  dat1penal = data.frame(time = c(dat1$y12),
              status = c(dat1$s12),
              id = rep(1:(n/2), 1))
  
  dat1penal = dat1penal[complete.cases(dat1penal),]
  
  penalfrail1 = frailtyPenal(Surv(time, status)~cluster(id), data=dat1penal,n.knots=8,kappa= 10000,
               RandDist = "LogN")
  omega12_z1 = o12save1[1:(n/2)] = (penalfrail1$frailty.pred)
  
  dat1penal = data.frame(time = c(dat1$y23, dat1$y13), status = c(dat1$s23, dat1$s13), y12 = c(dat1$y12, rep(0, n/2)), 
              id = rep(1:(n/2), 2))
  dat1penal = dat1penal[complete.cases(dat1penal),]
  
  penalfrail1 = frailtyPenal(Surv(time, status)~cluster(id) + y12,
                data=dat1penal,n.knots=8,kappa= 10000,
                RandDist = "LogN")
  
  omega23_z1 = o23save1[1:(n/2)] = omega13_z1 = o13save1[1:(n/2)] = (penalfrail1$frailty.pred)

  dat0penal = data.frame(time = c(dat0$y12),status = c(dat0$s12),  id = rep(1:(n/2), 1))

  dat0penal = dat0penal[complete.cases(dat0penal),]

  penalfrail0 = frailtyPenal(Surv(time, status)~cluster(id),
                data=dat0penal,n.knots=8,kappa= 10000,
                RandDist = "LogN")

  omega12_z0 = o12save0[1:(n/2)] = (penalfrail0$frailty.pred)
  
  dat0penal = data.frame(time = c(dat0$y23, dat0$y13), status = c(dat0$s23, dat0$s13), y12 = c(dat0$y12, rep(0, n/2)), 
              id = rep(1:(n/2), 2))
  
  dat0penal = dat0penal[complete.cases(dat0penal),]

  penalfrail0 = frailtyPenal(Surv(time, status)~cluster(id) + y12,
                data=dat0penal,n.knots=8,kappa= 10000,
                RandDist = "LogN")
  
  omega23_z0 = o23save0[1:(n/2)] = omega13_z0 = o13save0[1:(n/2)] = (penalfrail0$frailty.pred)
 }
 
 if(F){
 factor1 = frailtysd/sd(o12save0[1:(n/2)])
 factor1 = 1
 o12save0[1:(n/2)] = omega12_z0 = omega12_z0*factor1
 
 factor2 = frailtysd/sd(o13save0[1:(n/2)])
 factor2 = 1
 
 #o23save0[1:(n/2)] = omega23_z0 =
  o13save0[1:(n/2)] = omega13_z0 = omega13_z0*factor2
 
 factor3 = frailtysd/sd(o12save1[1:(n/2)])
 factor3 = 1
 
 o12save1[1:(n/2)] = omega12_z1 = omega12_z1*factor3
 factor4 = frailtysd/sd(o13save1[1:(n/2)])
 factor4 = 1
 
 #o23save1[1:(n/2)] = omega23_z1 =
  o13save1[1:(n/2)] = omega13_z1 = omega13_z1*factor4
 }
 
 if(F){
  o12save0[1:(n/2)] = omega12_z0 = rnorm(n/2, mean = 0, sd = frailtysd)
  

  o23save0[1:(n/2)] = omega23_z0 =
  o13save0[1:(n/2)] = omega13_z0 = rnorm(n/2, mean = 0, sd = frailtysd)
  
  
  o12save1[1:(n/2)] = omega12_z1 = rnorm(n/2, mean = 0, sd = frailtysd)

  
  o23save1[1:(n/2)] = omega23_z1 =
  o13save1[1:(n/2)] = omega13_z1 = rnorm(n/2, mean = 0, sd = frailtysd)
 }
 # plot(o12save0[1:(n/2)], omega12true0[1:(n/2)])
 # plot(o13save0[1:(n/2)], omega13true0[1:(n/2)])
 # plot(o12save1[1:(n/2)], omega12true1[1:(n/2)])
 # plot(o13save1[1:(n/2)], omega13true1[1:(n/2)])
 # plot(o23save1[1:(n/2)], omega23true1[1:(n/2)])
 # plot(o23save1[1:(n/2)], omega23true1[1:(n/2)])
 
if(F){
 rgamma(1, shape = 0.01 + sum(dat1$s23, na.rm = T),
       scale = (0.01 + 1/1 * sum(dat1$y23^ 1 * exp(1 * omega23_z1), na.rm = T))^(-1))
 
 rgamma(1, shape = 0.01 + sum(dat1$s13, na.rm = T),
     scale = (0.01 + 1/1 *sum(dat1$y13 ^ 1 * exp(1 * omega13_z1), na.rm = T))^(-1))
 
 rgamma(1, shape = 0.01 + sum(dat1$s12, na.rm = T),
     scale = (0.01 + 1/1 *sum(dat1$y12 ^ 1 * exp(1 * omega12_z1), na.rm = T))^(-1))

 rgamma(1, shape = 0.01 + sum(dat0$s23, na.rm = T), scale = (0.01 + 1/1 * sum(dat0$y23^ 1 * 
              exp(1 * omega23_z0), na.rm = T))^(-1))
 
 rgamma(1, shape = 0.01 + sum(dat0$s13, na.rm = T),
     scale = (0.01 + 1/1 *sum(dat0$y13 ^ 1 * exp(1 * omega13_z0), na.rm = T))^(-1))
 
 rgamma(1, shape = 0.01 + sum(dat0$s12, na.rm = T),
     scale = (0.01 + 1/1 *sum(dat0$y12 ^ 1 * exp(1 * omega12_z0), na.rm = T))^(-1)) 

 }
 z = 2
 
 mod23_0 = summary(survreg(Surv(dat0$y23, dat0$s23) ~ dat0$y12, dist = "exponential"))
 mod23_0 = summary(survreg(Surv(dat0$y23, dat0$s23) ~ dat0$y12 + omega23true0[1:(n/2)], dist = "exponential"))
 mod23_0 = summary(survreg(Surv(dat0$y23, dat0$s23) ~ dat0$y12 + omega23_z0, dist = "exponential"))

 mod23_1 = summary(survreg(Surv(dat1$y23, dat1$s23) ~ dat1$y12 + omega23true1[1:(n/2)], dist = "exponential"))
 mod23_1 = summary(survreg(Surv(dat1$y23, dat1$s23) ~ dat1$y12 + omega23_z1, dist = "exponential"))
 
 x_0 = x_1 = rep(0, n/2) # no covariates 
 
 if(holdtheta){
  holdtheta23_0[1] = theta23_0 = params_list$theta23_0
  holdtheta23_1[1] = theta23_1 = params_list$theta23_1
 }
 
 if(holdscale12){
  scale12_0 = params_list$scale12_0
  scale12_1 = params_list$scale12_1
 }
 
 if(holdshape){
  shape12_0 = params_list$shape12_0
  shape12_1 = params_list$shape12_1
  shape13_0 = params_list$shape13_0
  shape13_1 = params_list$shape13_1
  shape23_0 = params_list$shape23_0
  shape23_1 = params_list$shape23_1
 }
 
# if(holdfrail13)  {
#   o23save1[1:(n/2)] = omega23_z1 =
#    o13save1[1:(n/2)] = omega13_z1 = omega13true1[1:(n/2)]
#   o23save0[1:(n/2)] = omega23_z0 =
#    o13save0[1:(n/2)] = omega13_z0 = omega13true0[1:(n/2)]
#  }
# 
#  
#  if(holdfrail12){
#   o12save1[1:(n/2)] = omega12_z1 = omega12true1[1:(n/2)]
#   o12save0[1:(n/2)] = omega12_z0 = omega12true0[1:(n/2)]
#  }
 
while(z < SIM){ 

 zseed = round(runif(1, 1, 100)); set.seed(zseed + z + array_id)
 
 omega12_star = rnorm(n/2, o12save0[1:(n/2)], sd = proposalsdfrail)
 omega23_star = rnorm(n/2, o23save0[1:(n/2)], sd = proposalsdfrail)
 omega13_star = rnorm(n/2, o13save0[1:(n/2)], sd = proposalsdfrail)

 omega12_z0 = o12save0[1:(n/2)]
 omega13_z0 = o13save0[1:(n/2)]
 omega23_z0 = o23save0[1:(n/2)]
 
 save = savestar = NULL
 sap = sap2 = z1 = p = rep(NA, n/2)
 
 for(i in 1:(n/2)){
  
  sap[i] =  like12_omega_i(omega12_star[i], i = i, scale12_0 = holdscale12_0[z],
               dat0 = dat0, shape12_0 = holdshape12_0[z], c12_0 = holdc12_0[z])[1]
  sap2[i] =  like12_omega_i(omega12_z0[i], i = i, scale12_0 = holdscale12_0[z],
                dat0 = dat0, shape12_0 = holdshape12_0[z], c12_0 = holdc12_0[z])[1]
  
 p[i] = exp( sap[i] - sap2[i] +
                   sum( log(pnorm(omega12_star, 0, sd = frailtysd)))-
                   sum(log(pnorm(omega12_z0, 0, sd = frailtysd)) ))
  
  if(is.na(p[i])) p[i] = 0; if((p[i]<0)) p[i] = 0; if((p[i]>1)) p[i] = 1
  z1[i] = rbinom(1, 1, prob = p[i])
  if(holdfrail12) z1[i] = 0
  
 }
 
 if(sum(z1) >= 1) accept1 = accept1 + 1

 omega12_star = omega12_z0 = omega_star = z1 * omega12_star + (1-z1) * omega12_z0
 holdfrailsd12_0[z] = sd(omega12_star)
 holdfrailmean12_0[z] = mean(omega12_star)
 
 omega_cf = rnorm(n/2, omega_star, sd = frailtysd * sqrt(frailtysd - rhos^2))
 omega_cf = rnorm(n/2, rhos*omega_star, sd = frailtysd * sqrt(frailtysd - rhos^2))
 
 o12save0 = c(omega_star, omega_cf)
 
 holdomega12_0[,z] = omega12_z0[1:10]
 if(holdc) c12_0_star = 1
 beta12_0_star = rnorm(1, holdbeta12_0[z-1], sd = proposalsd)
 shape12_0_star = rtruncnorm(n = 1, a = 0, b = 3, mean = holdshape12_0[z-1], sd = proposalsd)
 beta12_0_star = 0
 if(holdshape) shape12_0_star = 1
 
 scale12_0_star = rgamma(1, shape = 0.01 + sum(dat0$s12, na.rm = T),
             scale = (0.01 + 1/shape12_0_star *
                   sum(dat0$y12 ^ shape12_0_star * 
                      exp(c12_0_star * omega12_star), na.rm = T))^(-1))
 
 
 if(holdscale12) scale12_0_star = scale12_0
 
 c12prior_p = pnorm(c12_0_star, 0, sd = 3)
 beta12prior_p = pnorm(beta12_0_star, 0, sd = 3)
 shape12prior_p = pnorm(shape12_0_star, 0, sd = 3)
 shape12prior_p = ptruncnorm(shape12_0_star, a = 0, b = 3, mean = 0, sd = 3)
 
 c12prior_c = pnorm(holdc12_0[z-1], 0, sd = 3)
 beta12prior_c = pnorm(holdbeta12_0[z-1], 0, sd = 3)
 shape12prior_c = pnorm(holdshape12_0[z-1], 0, sd = 3)
 shape12prior_c = ptruncnorm(holdshape12_0[z-1], a = 0, b = 3, mean = 0, sd = 3)

  p = exp(like12_shape(shape12_0 = shape12_0_star, scale12_0 = scale12_0_star, omega12_z0 = omega12_z0,
      dat0 = dat0, c12_0 = c12_0_star) + log(shape12prior_p) - 
       (like12_shape(shape12_0 = holdshape12_0[z-1], scale12_0 = scale12_0_star, omega12_z0 = omega12_z0,
       dat0 = dat0, c12_0 = c12_0_star)+ log(shape12prior_c))
      )

 #  
 # like12_shape = function(shape12_0, scale12_0){
 #  b =  ((shape12_0 - 1) * (sum(log(dat0$y12 ^ dat0$s12), na.rm = T))- 
 #      scale12_0/shape12_0 * sum(dat0$y12 ^ (shape12_0)* exp(c12_0_star * omega12_z0+
 #                                  holdbeta12_0[z-1] * x_0), na.rm = T))
 #  
 #  return(b)
 # }
 # 
 # 
 # p = exp(     like12_shape(shape12_0 = shape12_0_star,scale12_0 = scale12_0_star) + log(shape12prior_p) - 
 #          (like12_shape(shape12_0 = holdshape12_0[z-1], scale12_0 = scale12_0_star) + log(shape12prior_c))
 # )
 #  
 #  
  # p = exp(     like12_shape(shape12_0_star, scale12_0_star) + log(shape12prior_p) +
  #          (like12_shape(holdshape12_0[z-1], scale12_0 = scale12_0_star)+ log(shape12prior_c)) 
  # )
  # 
  # 
  if(is.nan(p)) p = 0;if(is.na(p)) p = 0; if((p<0)) p = 0; if((p>1)) p = 1
  if(scale12_0_star<0 ) p = 0
  z13 = z1 = rbinom(1, 1, prob = p)
 
 if(z1 == 1) accept7 = accept7 + 1
 scale12_0 = z1*scale12_0_star + (1-z1) * holdscale12_0[z - 1]
 c12_0 = c12_0_star
 shape12_0 = z1*shape12_0_star + (1-z1) * holdshape12_0[z - 1]
 beta12_0 = z13*beta12_0_star + (1-z13) * holdbeta12_0[z - 1]
 holdbeta12_0[z] = beta12_0
 holdscale12_0[z] = scale12_0
 holdc12_0[z] = c12_0
 holdshape12_0[z] = shape12_0
 
 i= 1

 z1 = p = rep(NA, (n/2))
 sap = sap2 = rep(NA, (n/2))
 
 for(i in 1:(n/2)){
  sap[i] =  like13_omega_i(omega13_z0 = omega13_star[i], i = i, 
               dat0 = dat0, c13_0 = holdc13_0[z-1],
               shape13_0 = holdshape13_0[z-1], scale13_0 = holdscale13_0[z-1], c23_0 = holdc23_0[z-1],
               shape23_0 = holdshape23_0[z-1], scale23_0 = holdscale23_0[z-1], theta23_0 = holdtheta23_0[z-1])[1]
  sap2[i] =  like13_omega_i(omega13_z0[i], i = i, 
                dat0 = dat0, c13_0 = holdc13_0[z-1],
                shape13_0 = holdshape13_0[z-1], scale13_0 = holdscale13_0[z-1], c23_0 = holdc23_0[z-1],
                shape23_0 = holdshape23_0[z-1], scale23_0 = holdscale23_0[z-1], theta23_0 = holdtheta23_0[z-1])[1]
  
  p[i] = exp( sap[i] - sap2[i] +
          sum( log(pnorm(omega13_star, 0, sd = frailtysd)))-
          sum(log(pnorm(omega13_z0, 0, sd = frailtysd)) ))
  
  if(is.na(p[i])) p[i] = 0; if((p[i]<0)) p[i] = 0; if((p[i]>1)) p[i] = 1
  z1[i] = rbinom(1, 1, prob = p[i])
  if(holdfrail13) z1[i] = 0
  
 }
 if(sum(z1) >= 1) accept3 = accept3 + 1
 acceptfrail1[,z] = z1[1:100]
 
 omega13_star = omega13_z0 = omega_star = z1 * omega13_star + (1-z1) * omega13_z0
 omega_cf = rnorm(n/2, omega_star, sd = 1 * sqrt(frailtysd - rhot^2))
 omega_cf = rnorm(n/2, rhot * omega_star, sd = frailtysd * sqrt(frailtysd - rhot^2))
 
 o13save0 = c(omega_star, omega_cf)
 
 holdfrailsd13_0[z] = sd(omega13_star)
 holdfrailmean13_0[z] = mean(omega13_star)

 omega23_z0 = omega23_star = omega13_z0
 
 holdomega23_0[,z] = omega23_z0[1:10]
 holdfrailsd23_0[z] = sd(omega23_star)
 holdfrailmean23_0[z] = mean(omega23_star)

 o23save0 = o13save0 
 c13_0_star= rtruncnorm(1, a = -3, b = 3, mean = holdc13_0[z-1], sd = proposalsd)
 shape13_0_star = rtruncnorm(1, a = 0, b = 3, mean = holdshape13_0[z-1], sd = proposalsd)
 if(holdshape) shape13_0_star = 1
 
 if(holdc) c13_0_star = 1
 beta13_0_star = rnorm(1, holdbeta13_0[z-1], sd = proposalsd)
 beta13_0_star = 0
 scale13_0_star = 1/rinvgamma(1, shape = 0.01 + sum(dat0$s13),
                scale = 0.01 + 1/shape13_0_star * sum(dat0$y13^shape13_0_star * exp(c13_0_star *
              omega13_z0 + holdbeta13_0[z-1] * x_0)))

 beta23_0_star = rnorm(1, holdbeta23_0[z-1], sd = proposalsd)
 beta23_0_star = 0
 shape23_0_star = rtruncnorm(1, a = 0, b = 3, mean = holdshape23_0[z-1], sd = proposalsd)
 
 if(holdshape) shape23_0_star = 1
 
 theta23_0_star = rnorm(1, -mod23_0$coefficients[2], sd = proposalsdtheta)
 theta23_0_star = rnorm(1, holdtheta23_0[z-1], sd = proposalsd)
 
 if(holdtheta) theta23_0_star = theta23_0
 
 c23_0_star = rtruncnorm(1, a=-3, b= 3, mean = holdc23_0[z-1], sd = proposalsd)
 if(holdc) c23_0_star = 1
 
 scale23_0_star = rgamma(1, shape = 0.01 + sum(dat0$s23, na.rm = T),
             scale = (0.01 + 1/shape23_0_star *
                   sum(dat0$y23 ^ shape23_0_star * 
                      exp(c23_0_star * omega23_z0 +
                         theta23_0_star * dat0$y12), na.rm = T))^(-1))
 

 if(holdscale23) scale23_0_star = scale23_0

 c13prior_p = pnorm(c13_0_star, 0, sd = 3)
 beta13prior_p = pnorm(beta13_0_star, 0, sd = 3)
 shape13prior_p = pnorm(shape13_0_star, 0, sd = 3)
 shape13prior_p = ptruncnorm(shape13_0_star, a = 0, b = 3, mean = 0, sd = 3)
 
 c13prior_c = pnorm(holdc13_0[z-1], 0, sd = 3)
 beta13prior_c = pnorm(holdbeta13_0[z-1], 0, sd = 3)
 shape13prior_c = pnorm(holdshape13_0[z-1], 0, sd = 3)
 shape13prior_c = ptruncnorm(holdshape13_0[z-1], a = 0, b = 3, mean = 0, sd = 3)
 
 c23prior_p = pnorm(c23_0_star, 0, sd = 3)
 beta23prior_p = pnorm(beta23_0_star, 0, sd = 3)
 theta23prior_p = pnorm(theta23_0_star, 0, sd = 3)
 shape23prior_p = pnorm(shape23_0_star, 0, sd = 3)
 shape23prior_p = ptruncnorm(shape23_0_star, a = 0, b = 3, mean = 0, sd = 3)
 
 c23prior_c = pnorm(holdc23_0[z-1], 0, sd = 3)
 beta23prior_c = pnorm(holdbeta23_0[z-1], 0, sd = 3)
 theta23prior_c = pnorm(holdtheta23_0[z-1], 0, sd = 3)
 shape23prior_c = pnorm(holdshape23_0[z-1], 0, sd = 3)
 shape23prior_c = ptruncnorm(holdshape23_0[z-1], a = 0, b = 3, mean = 0, sd = 3)
 
 theta23prop_p = pnorm(theta23_0_star, -mod23_0$coefficients[2], sd = proposalsd)
 theta23prop_c = pnorm(holdtheta23_0[z-1], -mod23_0$coefficients[2], sd = proposalsd)
 
 theta23prop_p = pnorm(theta23_0_star, 0, sd = proposalsd)
 theta23prop_c = pnorm(holdtheta23_0[z-1], 0, sd = proposalsd)
 
  p = exp(
       like13_shape(dat0 = dat0, c13_0 = c13_0_star,
              shape13_0 = shape13_0_star, scale13_0 = scale13_0_star, omega13_z0 = omega13_z0) +
        log(shape13prior_p)-
       (like13_shape(dat0 = dat0, c13_0 = holdc13_0[z-1],
              shape13_0 = holdshape13_0[z-1], scale13_0 = scale13_0_star, omega13_z0 = omega13_z0)+
         log(shape13prior_c)) )[1]
  
  if(T){
  p2 = exp(  like23_theta(theta23 = theta23_0_star, dat0 = dat0,c23_0 = holdc23_0[z-1], omega13_z0 = omega13_z0,
               shape23_0 = holdshape23_0[z-1], scale23_0 = holdscale23_0[z-1])+ log(shape13prior_c) + 
    log(theta23prior_p) + log(theta23prop_p) -
         (like23_theta(theta23_0 = holdtheta23_0[z-1],dat0 = dat0, omega13_z0 = omega13_z0, c23_0 = holdc23_0[z-1],
                shape23_0 = holdshape23_0[z-1], scale23_0 = holdscale23_0[z-1])+ log(theta23prior_c) + log(theta23prop_c))
        + like23_shape(shape23_0 = shape23_0_star, dat0 = dat0, omega13_z0 = omega13_z0, c23_0 = holdc23_0[z-1],
                scale23_0 = holdscale23_0[z-1], theta23_0 = holdtheta23_0[z-1]) + log(shape23prior_p)-
         (like23_shape(shape23_0 = holdshape23_0[z-1], dat0 = dat0, c23_0 = holdc23_0[z-1], omega13_z0 = omega13_z0, 
                scale23_0 = holdscale23_0[z-1], theta23_0 = holdtheta23_0[z-1]) + log(shape23prior_c)))
  }
  
  if(F){
  like23_theta = function(theta23_0_star){
 
   a =  (c(c23_0_star, theta23_0_star) * (dat0$s23 %*% (cbind(omega13_z0, dat0$y12))) -
        scale23_0_star/shape23_0_star * sum(dat0$y23 ^ 1/scale23_0_star * exp(c23_0_star * omega13_z0 + 
                                            theta23_0_star * dat0$y12), na.rm = T))[2]

   return(a)
  }
  
  like23_shape = function(shape23_0){
   a = log ( (prod(dat0$y23 ^ dat0$s23, na.rm = T))^(shape23_0 - 1) *
          exp( -scale23_0_star/shape23_0  * sum(dat0$y23 ^ (shape23_0) * exp(c23_0_star * omega13_z0 + theta23_0_star * dat0$y12), na.rm = T)))
   return(a)
  }
  
  p2 = exp( 
         like23_theta(theta23_0_star) + log(theta23prior_p) + log(theta23prop_p) -
         (like23_theta(holdtheta23_0[z-1])+ log(theta23prior_c) + log(theta23prop_c))
        + like23_shape(shape23_0_star) + log(shape23prior_p)-
         (like23_shape(holdshape23_0[z-1])+ log(shape23prior_c)))
  
  
  }
  
  if(is.nan(p)) p = 0;if(is.na(p)) p = 0; if((p<0)) p = 0; if((p>1)) p = 1
  if(is.nan(p2)) p2 = 0;if(is.na(p2)) p2 = 0; if((p2<0)) p2 = 0; if((p2>1)) p2 = 1
  
  if(scale23_0_star<0) {p2 = 0}
  if(scale13_0_star<0) {p = 0}
  z14 = z1 = z12 = z13 = rbinom(1, 1, prob = p)
  z15 = rbinom(1, 1, prob = p2)
 
 if(z1 == 1) accept4 = accept4 + 1
 if(z15 == 1) accept5 = accept5 + 1
 

  holdomega13_0[,z] = omega13_z0[1:10]
 
 scale13_0 = z1*scale13_0_star + (1-z1) * holdscale13_0[z - 1]
 c13_0 = z1*c13_0_star + (1-z1) * holdc13_0[z - 1]
 beta13_0 = z13*beta13_0_star + (1-z13) * holdbeta13_0[z - 1]
 shape13_0 = z12*shape13_0_star + (1-z12) * holdshape13_0[z - 1]
 
 holdbeta13_0[z] = beta13_0
 holdscale13_0[z] = scale13_0
 holdc13_0[z] = c13_0
 holdshape13_0[z] = shape13_0
 
 scale23_0 = z15*scale23_0_star + (1-z15) * holdscale23_0[z - 1]
 theta23_0 = z15*theta23_0_star + (1-z15) * holdtheta23_0[z - 1]
 c23_0 = z15*c23_0_star + (1-z15) * holdc23_0[z - 1]
 beta23_0 = z15*beta23_0_star + (1-z15) * holdbeta23_0[z - 1]
 shape23_0 = z15*shape23_0_star + (1-z15) * holdshape23_0[z - 1]
 
 holdbeta23_0[z] = beta23_0
 holdscale23_0[z] = scale23_0
 holdc23_0[z] = c23_0
 holdshape23_0[z] = shape23_0
 holdtheta23_0[z] = theta23_0

 omega12_star = rnorm(n/2, o12save1[1:(n/2)], sd = proposalsdfrail)
 omega23_star = rnorm(n/2, o23save1[1:(n/2)], sd = proposalsdfrail)
 omega13_star = rnorm(n/2, o13save1[1:(n/2)], sd = proposalsdfrail)

omega12_z1 = o12save1[1:(n/2)]
omega13_z1 = o13save1[1:(n/2)]
omega23_z1 = o23save1[1:(n/2)]

save = savestar = NULL
z1 = p = rep(NA, n/2)

for(i in 1:(n/2)){
 
 sap[i] =  like12_omega1_i(omega12_star[i], i = i, scale12_1 = holdscale12_1[z-1], dat1 = dat1, shape12_1 = holdshape12_1[z-1],
               c12_1 = holdc12_1[z-1])[1]
 sap2[i] =  like12_omega1_i(omega12_z1[i], i = i, scale12_1 = holdscale12_1[z-1], dat1 = dat1, shape12_1 = holdshape12_1[z-1],
               c12_1 = holdc12_1[z-1])[1]
 p[i] = exp( sap[i] - sap2[i] +
         sum( log(pnorm(omega12_star, 0, sd = frailtysd)))-
         sum(log(pnorm(omega12_z1, 0, sd = frailtysd)) ))
 
 if(is.na(p[i])) p[i] = 0; if((p[i]<0)) p[i] = 0; if((p[i]>1)) p[i] = 1
 z1[i] = rbinom(1, 1, prob = p[i])
 if(holdfrail12) z1[i] = 0
 
}

 if(sum(z1) >= 1) accept6 = accept6 + 1

 omega12_star = omega12_z1 = omega_star = z1 * omega12_star + (1-z1) * omega12_z1
 holdfrailsd12_1[z] = sd(omega12_star)
 holdfrailmean12_1[z] = mean(omega12_star)
 
 omega_cf = rnorm(n/2, omega_star, sd = frailtysd * sqrt(frailtysd - rhos^2))
 omega_cf = rnorm(n/2, rhos * omega_star, sd = frailtysd * sqrt(frailtysd - rhos^2))
 
 o12save1 = c(omega_star, omega_cf)
 
 holdomega12_1[,z] = omega12_z1[1:10]
 if(holdc) c12_1_star = 1
 beta12_1_star = rnorm(1, holdbeta12_1[z-1], sd = proposalsd)
 shape12_1_star = rtruncnorm(n = 1, a = 0, b = 3, mean = holdshape12_1[z-1], sd = proposalsd)
 beta12_1_star = 0

 if(holdshape) shape12_1_star = shape12_1
 
  scale12_1_star = rgamma(1, shape = 0.01 + sum(dat1$s12, na.rm = T),
  scale = (0.01 + 1/shape12_1_star *
       sum(dat1$y12 ^ shape12_1_star * 
          exp(c12_1_star * omega12_star), na.rm = T))^(-1))
  
  

  
 if(holdscale12) scale12_1_star = scale12_1
  
 c12prior_p = pnorm(c12_1_star, 0, sd = 3)
 beta12prior_p = pnorm(beta12_1_star, 0, sd = 3)
 shape12prior_p = pnorm(shape12_1_star, 0, sd = 3)
 shape12prior_p = ptruncnorm(shape12_1_star, a = 0, b = 3, mean = 0, sd = 3)
 
 c12prior_c = pnorm(holdc12_1[z-1], 0, sd = 3)
 beta12prior_c = pnorm(holdbeta12_1[z-1], 0, sd = 3)
 shape12prior_c = pnorm(holdshape12_1[z-1], 0, sd = 3)
 shape12prior_c = ptruncnorm(holdshape12_1[z-1], a = 0, b = 3, mean = 0, sd = 3)
 
  # p = exp(     like12_shape1(shape12_1 = shape12_1_star, scale12_1 = scale12_1_star,
  #                dat1 = dat1, c12_1 = c12_1_star, omega12_z1 = omega12_z1
  #                ) + log(shape12prior_p) +
  #      (like12_shape1(shape12_1 = holdshape12_1[z-1], scale12_1 = scale12_1_star,
  #             dat1 = dat1, c12_1 = c12_1_star, omega12_z1 = omega12_z1)+ log(shape12prior_c)) 
  #    )
  # 
  like12_shape1 = function(shape12_1, scale12_1){
    b =  ((shape12_1 - 1) * (sum(log(dat1$y12 ^ dat1$s12), na.rm = T))- 
        scale12_1/shape12_1 * sum(dat1$y12 ^ (shape12_1)* exp(c12_1_star * omega12_z1+
       holdbeta12_1[z-1] * x_1), na.rm = T))
   
   return(b)
  }
  
  
  p = exp(     like12_shape1(shape12_1 = shape12_1_star,scale12_1 = scale12_1_star) + log(shape12prior_p) - 
            (like12_shape1(shape12_1 = holdshape12_1[z-1], scale12_1 = scale12_1_star) + log(shape12prior_c))
           )
  
  
 if(is.nan(p)) p = 0;if(is.na(p)) p = 0; if((p<0)) p = 0; if((p>1)) p = 1
 if(scale12_1_star<0 ) p = 0
  z12 = z13 = z1 = rbinom(1, 1, prob = p)

 if(z1 == 1) accept7 = accept7 + 1
 scale12_1 = z1*scale12_1_star + (1-z1) * holdscale12_1[z - 1]

 c12_1 = c12_1_star
 shape12_1 = z1*shape12_1_star + (1-z1) * holdshape12_1[z - 1]
 beta12_1 = z13*beta12_1_star + (1-z13) * holdbeta12_1[z - 1]
 
 holdbeta12_1[z] = beta12_1
 holdscale12_1[z] = scale12_1
 holdc12_1[z] = c12_1
 holdshape12_1[z] = shape12_1

 z1 = p = rep(NA, (n/2))
 sap = sap2 = rep(NA, (n/2))
 
 for(i in 1:(n/2)){
  sap[i] =  like13_omega1_i(omega13_star[i], i = i, scale13_1 = holdscale13_1[z-1], dat1 = dat1, shape13_1 = holdshape13_1[z-1],
                c13_1 = holdc13_1[z-1], shape23_1 = holdshape23_1[z-1], scale23_1 = holdscale23_1[z-1], theta23_1 = holdtheta23_1[z-1],
                c23_1 = holdc23_1[z-1])[1]
  sap2[i] =  like13_omega1_i(omega13_z1[i], i = i, scale13_1 = holdscale13_1[z-1], dat1 = dat1, shape13_1 = holdshape13_1[z-1],
                c13_1 = holdc13_1[z-1], shape23_1 = holdshape23_1[z-1], scale23_1 = holdscale23_1[z-1], theta23_1 = holdtheta23_1[z-1],
                c23_1 = holdc23_1[z-1])[1]
  
 p[i] = exp( sap[i] - sap2[i] +
          sum( log(pnorm(omega13_star, 0, sd = frailtysd)))-
          sum(log(pnorm(omega13_z1, 0, sd = frailtysd)) ))
 
  
  if(is.na(p[i])) p[i] = 0; if((p[i]<0)) p[i] = 0; if((p[i]>1)) p[i] = 1
  z1[i] = rbinom(1, 1, prob = p[i])
  if(holdfrail13) z1[i] = 0
  
 }
 
 omega23_z1 = omega13_star = omega13_z1 = omega_star = z1 * omega13_star + (1-z1) * omega13_z1
 omega_cf = rnorm(n/2, omega_star, sd = frailtysd * sqrt(frailtysd - rhos^2))
 omega_cf = rnorm(n/2, rhos * omega_star, sd = frailtysd * sqrt(frailtysd - rhos^2))
 
 o13save1 = c(omega_star, omega_cf)

 if(sum(z1[1]) >= 1) accept8 = accept8 + 1
 
  holdfrailsd13_1[z] = sd(omega13_star)
  holdfrailmean13_1[z] = mean(omega13_star)
  
  z1 = p = rep(NA, (n/2))
  sap = sap2 = rep(NA, (n/2))
 
  o23save1 = o13save1
  
  holdomega23_0[,z] = omega23_z1[1:10]
  holdfrailsd23_1[z] = sd(omega23_star)
  holdfrailmean23_1[z] = mean(omega23_star)
  
  holdomega23_1[,z] = omega23_z1[1:10]
  
 holdomega13_1[,z] = omega13_z1[1:10]
 c13_1_star = rnorm(1, holdc13_1[z-1], sd = proposalsd)

 if(holdc) c13_1_star = 1
 beta13_1_star = rnorm(1, holdbeta13_1[z-1], sd = proposalsd)
 beta13_1_star = 0
 shape13_1_star = rtruncnorm(1, a = 0, b = 3, mean = holdshape13_1[z-1], sd = proposalsd)
 if(holdshape) shape13_1_star = 1

 scale13_1_star = rgamma(1, shape = 0.01 + sum(dat1$s13, na.rm = T),
             scale = (0.01 + 1/shape13_1_star *
                   sum(dat1$y13 ^ shape13_1_star * 
                      exp(c13_1_star * omega13_z1), na.rm = T))^(-1))
 
 if(holdscale13) scale13_1_star = scale13_1

 beta23_1_star = 0
 shape23_1_star = rnorm(1, holdshape23_1[z-1], sd = proposalsd)
 
 shape23_1_star = rtruncnorm(1, a = 0, b = 3, mean = holdshape23_1[z-1], sd = proposalsd)
 theta23_1_star = rnorm(1, -mod23_1$coefficients[2], sd = proposalsdtheta)
 theta23_1_star = rnorm(1, holdtheta23_1[z-1], sd = proposalsd)
 
 if(holdtheta){ theta23_1_star = theta23_1 }

 c23_1_star = rtruncnorm(1, a=-3, b= 3, mean = holdc23_1[z-1], sd = proposalsd)
 
 if(holdshape) shape23_1_star = 1
 if(holdc) c23_1_star = 1
 
 scale23_1_star = rgamma(1, shape = 0.01 + sum(dat1$s23, na.rm = T),
             scale = (0.01 + 1/shape23_1_star *
                   sum(dat1$y23 ^ shape23_1_star * 
                      exp(c23_1_star * omega23_z1 +
                         theta23_1_star * dat1$y12), na.rm = T))^(-1))
 
 scale23_1_star = 1/rinvgamma(1, shape = 0.01 + sum(dat1$s23, na.rm = T),
                scale = (0.01 + 1/shape23_1_star *
                     sum(dat1$y23 ^ shape23_1_star * 
                        exp(c23_1_star * omega23_z1 +
                           theta23_1_star * dat1$y12), na.rm = T)))
 
 if(holdscale23) scale23_1_star = scale23_1

 c13prior_p = pnorm(c13_1_star, 0, sd = 3)
 beta13prior_p = pnorm(beta13_1_star, 0, sd = 3)
 shape13prior_p = pnorm(shape13_1_star, 0, sd = 3)
 shape13prior_p = ptruncnorm(shape13_1_star, a = 0, b = 3, mean = 0, sd = 3)
 
 c13prior_c = pnorm(holdc13_1[z-1], 0, sd = 3)
 beta13prior_c = pnorm(holdbeta13_1[z-1], 0, sd = 3)
 shape13prior_c = pnorm(holdshape13_1[z-1], 0, sd = 3)
 shape13prior_c = ptruncnorm(holdshape13_1[z-1], a = 0, b = 3, mean = 0, sd = 3)
 
 c23prior_p = pnorm(c23_1_star, 0, sd = 3)
 beta23prior_p = pnorm(beta23_1_star, 0, sd = 3)
 theta23prior_p = pnorm(theta23_1_star, 0, sd = 3)
 shape23prior_p = pnorm(shape23_1_star, 0, sd = 3)
 shape23prior_p = ptruncnorm(shape23_1_star, a = 0, b = 3, mean = 0, sd = 3)
 
 c23prior_c = pnorm(holdc23_1[z-1], 0, sd = 3)
 beta23prior_c = pnorm(holdbeta23_1[z-1], 0, sd = 3)
 shape23prior_c = pnorm(holdshape23_1[z-1], 0, sd = 3)
 shape23prior_c = ptruncnorm(holdshape23_1[z-1], a = 0, b = 3, mean = 0, sd = 3)
 theta23prior_c = pnorm(holdtheta23_1[z-1], 0, sd = 3)
 
 theta23prop_p = pnorm(theta23_1_star, -mod23_1$coefficients[2], sd = proposalsd)
 theta23prop_c = pnorm(holdtheta23_1[z-1], -mod23_1$coefficients[2], sd = proposalsd)
 
 theta23prop_p = pnorm(theta23_1_star, 0, sd = proposalsd)
 theta23prop_c = pnorm(holdtheta23_1[z-1], 0, sd = proposalsd)
 
 
  p = exp(like13_shape1(shape13_1 = shape13_1_star, scale13_1 = scale13_1_star, dat1 = dat1,
             omega13_z1 = omega13_z1, c13_1 = c13_1_star) + log(shape13prior_p)-
       (like13_shape1(shape13_1 = holdshape13_1[z-1], scale13_1 = scale13_1_star, dat1 = dat1,
               omega13_z1 = omega13_z1, c13_1 = c13_1_star)+ log(shape13prior_c))
       # like13_beta1(beta13_1_star) + log(beta13prior_p)-
      # (like13_beta1(holdbeta13_1[z-1])+ log(beta13prior_c)) +
       # like13_1(c13_1_star) + log(c13prior_p)-
      # (like13_1(holdc13_1[z-1])+ log(c13prior_c) )
     )
  
   #like23_1(c23_1_star) + log(c23prior_p)-
       #  (like23_1(holdc23_1[z-1])+ log(c23prior_c))+
       #  like23_beta1(beta23_1_star) + log(beta23prior_p)-
       #  (like23_beta1(holdbeta23_1[z-1])+ log(beta23prior_c))+
  # p2 = exp(   like23_theta1(theta23_1 = theta23_1_star, scale23_1 = scale23_1_star, dat1 = dat1, shape23_1 = shape23_1_star,
  #               omega13_z1 = omega13_z1, c23_1 = c23_1_star) + log(theta23prior_p) + log(theta23prop_p) -
  #        (like23_theta1(theta23_1 = holdtheta23_1[z-1], scale23_1 = scale23_1_star, dat1 = dat1, shape23_1 = shape23_1_star,
  #               omega13_z1 = omega13_z1, c23_1 = c23_1_star) + log(theta23prior_c) + log(theta23prop_c))
  #       + like23_shape1(shape23_1 = shape23_1_star, omega13_z1 = omega13_z1, theta23_1 = theta23_1_star, 
  #               scale23_1 = scale23_1_star, dat1 = dat1, c23_1 = c23_1_star) +
  #        log(shape23prior_p)-
  #        (like23_shape1(shape23_1 = holdshape23_1[z-1], omega13_z1 = omega13_z1, scale23_1 = scale23_1_star, 
  #               theta23_1 = theta23_1_star, dat1 = dat1, c23_1 = c23_1_star) +
  #         log(shape23prior_c))
  # )
  # 
  p2 = exp(  like23_theta1(theta23 = theta23_1_star, dat1 = dat1,c23_1 = holdc23_1[z-1], omega13_z1 = omega13_z1,
               shape23_1 = holdshape23_1[z-1], scale23_1 = holdscale23_1[z-1])+ log(shape13prior_c) +
         log(theta23prior_p) + log(theta23prop_p) -
         (like23_theta1(theta23_1 = holdtheta23_1[z-1],dat1 = dat1, omega13_z1 = omega13_z1, c23_1 = holdc23_1[z-1],
                shape23_1 = holdshape23_1[z-1], scale23_1 = holdscale23_1[z-1])+ log(theta23prior_c) + log(theta23prop_c))
        + like23_shape1(shape23_1 = shape23_1_star, dat1 = dat1, omega13_z1 = omega13_z1, c23_1 = holdc23_1[z-1],
                scale23_1 = holdscale23_1[z-1], theta23_1 = holdtheta23_1[z-1]) + log(shape23prior_p)-
         (like23_shape1(shape23_1 = holdshape23_1[z-1], dat1 = dat1, c23_1 = holdc23_1[z-1], omega13_z1 = omega13_z1,
                scale23_1 = holdscale23_1[z-1], theta23_1 = holdtheta23_1[z-1]) + log(shape23prior_c)))

  # 
  # like23_theta1 = function(theta23_1_star){
  # 
  #  c =  (c(c23_1_star, theta23_1_star) * (dat1$s23 %*% (cbind(omega13_z1, dat1$y12))) -
  #      scale23_1_star/shape23_1_star * sum(dat1$y23 ^ 1/scale23_1_star * exp(c23_1_star * omega13_z1 + 
  #                                           theta23_1_star * dat1$y12), na.rm = T))[2]
  #  return(c)
  # }
  # 
  # 
  # 
  # like23_shape1 = function(shape23_1){
  # 
  #  a = log ( (prod(dat1$y23 ^ dat1$s23, na.rm = T))^(shape23_1 - 1) *
  #         exp( -scale23_1_star/shape23_1  * sum(dat1$y23 ^ (shape23_1) * exp(c23_1_star * omega13_z1 + theta23_1_star * dat1$y12), na.rm = T)))
  #  return(a)
  # }
  # 
  # p2 = exp( (
  #        like23_theta1(theta23_1_star) + log(theta23prior_p) + log(theta23prop_p) -
  #        (like23_theta1(holdtheta23_1[z-1])+ log(theta23prior_c) + log(theta23prop_c))
  #       + like23_shape1(shape23_1_star) + log(shape23prior_p)-
  #        (like23_shape1(holdshape23_1[z-1])+ log(shape23prior_c))))
  # 
  # 
  # 
  if(is.nan(p)) p = 0;if(is.na(p)) p = 0; if((p<0)) p = 0; if((p>1)) p = 1
  if(is.nan(p2)) p2 = 0;if(is.na(p2)) p2 = 0; if((p2<0)) p2 = 0; if((p2>1)) p2 = 1
  
  if(scale23_1_star<0) {p2 = 0; p = 0}
  if(scale13_1_star<0) {p = 0}
  z12 = z14 = z1 = z13 = rbinom(1, 1, prob = p)
  z15 = rbinom(1, 1, prob = p2)
  
 if(z12 == 1) accept9 = accept9 + 1
 if(z15 == 1) accept10 = accept10 + 1
 
 scale13_1 = z1*scale13_1_star + (1-z1) * holdscale13_1[z - 1]

 c13_1 = c13_1_star
 beta13_1 = z1*beta13_1_star + (1-z1) * holdbeta13_1[z - 1]
 shape13_1 = z1*shape13_1_star + (1-z1) * holdshape13_1[z - 1]
 
 holdbeta13_1[z] = beta13_1 
 holdscale13_1[z] = scale13_1
 holdc13_1[z] = c13_1
 holdshape13_1[z] = shape13_1

 holdscale23_1[z] = scale23_1 = z15*scale23_1_star + (1-z15) * holdscale23_1[z - 1]
 holdtheta23_1[z] = z15*theta23_1_star + (1-z15) * holdtheta23_1[z - 1]
 holdc23_1[z] = c23_1 = z15*c23_1_star + (1-z15) * holdc23_1[z - 1]
 holdshape23_1[z] = shape23_1 = z15*shape23_1_star + (1-z15) * holdshape23_1[z - 1]
 holdbeta23_1[z] = beta23_1 = z15*beta23_1_star + (1-z15) * holdbeta23_1[z - 1]
 
 
 o12save1flip = c(o12save1[(n/2 + 1): n], o12save1[1:(n/2)])
 o13save1flip = c(o13save1[(n/2 + 1): n], o13save1[1:(n/2)])

 # o12save0flip = c(o12save0[(n/2 + 1): n], o12save0[1:(n/2)])
 # o13save0flip = c(o13save0[(n/2 + 1): n], o13save0[1:(n/2)])
 # 
 o12save0flip = o12save0
 o13save0flip = o13save0
 
 Fw_0 = Fw = Fw_1 = NULL
 j = 1
 xtrue = c(x_0, x_1)
 
 intfunction = function(j, i, t){
  exp(
   -Lambda13_frailty_lk(x = t, xdata = xtrue[i], omega2 = o13save0flip[i], scale13 = holdscale13_0[z], 
              shape13 = holdshape13_0[z],  c13 = holdc13_0[z], beta13_1 = holdbeta13_0[z]) - 
    Lambda12_frailty_lk(x = t, xdata = xtrue[i], 
              omega1 = o12save0flip[i],  scale12 = holdscale12_0[z], shape12 = holdshape12_0[z], c12 = holdc12_0[z], beta12_1 = holdbeta12_0[z])) *
   lambda12_frailty(t, xdata = xtrue[i], omega1 = o12save0[i], scale12 = holdscale12_0[z], shape12 = holdshape12_0[z],
            c12 = holdc12_0[z], beta12_1 = holdbeta12_0[z]) * 
   exp(-Lambda23_frailty_lk(x = tau_t - t, xdata = xtrue[i], omega1 = o12save0flip[i], omega2 = o13save0flip[i], scale23 = holdscale23_0[z], 
                shape23 = holdshape23_0[z],  c23 = holdc23_0[z], theta23 = holdtheta23_0[z], v_predict = t, beta23_1 = holdbeta23_0[z]
   ))
 }
 for(i in 1:n){
  y =tryCatch({ integrate(f = intfunction, lower = 0, upper = tau_t, i = i, j = j)$val
  },
  error = function(error_condition) {
   NA}
  ) 
  
  L12 = Lambda12_frailty_lk(x = tau_t, xdata = xtrue[i], shape12 = holdshape12_0[z],
               omega1 = o12save0flip[i],  scale12 =holdscale12_0[z], c12 = holdc12_0[z], beta12_1 = holdbeta12_0[z])
  L13 = Lambda13_frailty_lk(x = tau_t, xdata = xtrue[i], omega2 = o13save0flip[i], scale13 = holdscale13_0[z], shape13 = holdshape13_0[z],
               c13 = holdc13_0[z], beta13_1 = holdbeta13_0[z])
  
  Fw = (exp(-L13 - L12) + y)

  Fw_0 = c(Fw_0, Fw)
 }
 
 intfunction = function(j, i, t){
  exp(
   -Lambda13_frailty_lk(x = t, xdata = xtrue[i], omega2 = o13save1flip[i], scale13 = holdscale13_1[z], 
              shape13 = holdshape13_1[z],  c13 = holdc13_1[z], beta13_1 = holdbeta13_1[z]) - 
    Lambda12_frailty_lk(x = t, xdata = xtrue[i], 
              omega1 = o12save1flip[i],  scale12 = holdscale12_1[z], shape12 = holdshape12_1[z], c12 = holdc12_1[z], beta12_1 = holdbeta12_1[z])) *
   lambda12_frailty(t, xdata = xtrue[i], omega1 = o12save1flip[i], scale12 = holdscale12_1[z], shape12 = holdshape12_1[z],
            c12 = holdc12_1[z], beta12_1 = holdbeta12_1[z])* 
   exp(-Lambda23_frailty_lk(x = tau_t - t, xdata = xtrue[i], omega1 = o12save1flip[i], omega2 = o13save1flip[i], scale23 = holdscale23_1[z], 
                shape23 = holdshape23_1[z], c23 = holdc23_1[z], theta23 = holdtheta23_1[z], v_predict = t, beta23_1 = holdbeta23_1[z]))
 }

  for(i in 1:n){
  y =tryCatch({ integrate(f = intfunction, lower = 0, upper = tau_t, i = i, j = j)$val
  },
  error = function(error_condition) {
   NA}
  ) 
  L12 =Lambda12_frailty_lk(x = tau_t, xdata = xtrue[i], shape12 = holdshape12_1[z],
               omega1 = o12save1flip[i],  scale12 =holdscale12_1[z], c12 = holdc12_1[z], beta12_1 = holdbeta12_1[z])
  L13 = Lambda13_frailty_lk(x = tau_t, xdata = xtrue[i], omega2 = o13save1flip[i], scale13 = holdscale13_1[z], 
               shape13 = holdshape13_1[z],  c13 = holdc13_1[z], beta13_1 = holdbeta13_1[z])
  
  
  Fw = (exp(-L13 - L12) + y)
  
  Fw_1 = c(Fw_1, Fw)
  
 }
 

 s0cumulative = -pweibull(q = tau_s, scale= 1/holdscale12_0[z], shape= holdshape12_0[z], lower = F, log = T)*exp(holdc12_0[z]*o12save0flip)
 s1cumulative = -pweibull(q = tau_s, scale= 1/holdscale12_1[z], shape= holdshape12_1[z], lower = F, log = T)*exp(holdc12_1[z]*o12save1flip)
 
 dat = data.frame(cbind(s0cumulative/s1cumulative, Fw_1 - Fw_0))
 dat$X = log(s0cumulative/s1cumulative)
 dat$Y = c(Fw_1 - Fw_0) 

 reg=lm(formula = Y ~ X, data=dat)
 holdslope[z] = (summary(reg)$coef[2,1])
 holdint[z] = (summary(reg)$coef[1,1])
 
 saveCEPx[, z] = dat$X
 saveCEPy[, z] = dat$Y

 if(z %% 100 == 0) print(z)
 z = z + 1
 
}


 params = data.frame(int = holdint, slope = holdslope,
          scale12_0 = holdscale12_0, scale13_0 = holdscale13_0, scale23_0 = holdscale23_0,
          scale12_1 = holdscale12_1, scale13_1 = holdscale13_1, scale23_1 = holdscale23_1,
          shape12_0 = holdshape12_0, shape13_0 = holdshape13_0, shape23_0 = holdshape23_0,
          shape12_1 = holdshape12_1, shape13_1 = holdshape13_1, shape23_1 = holdshape23_1,
          c12_0 = holdc12_0, c13_0 = holdc13_0, c23_0 = holdc23_0,
          c12_1 = holdc12_1, c13_1 = holdc13_1, c23_1 = holdc23_1,
          beta12_0 = holdbeta12_0, beta13_0 = holdbeta13_0, beta23_0 = holdbeta23_0,
          beta12_1 = holdbeta12_1, beta13_1 = holdbeta13_1, beta23_1 = holdbeta23_1,
          theta23_0 = holdtheta23_0, theta23_1 = holdtheta23_1
)

colMeans(params, na.rm = T)

result = list(params = params, saveCEPx = saveCEPx, saveCEPy = saveCEPy, 
       accept.rate = accept1/(SIM - burnin), 
       args = list(SIM = SIM, burnin = burnin, n = n))


return(result)

}

set.seed(1 + array_id)

params_res = run_sim(SIM = SIM, rhos = rhos, rhot = rhot, frailtysd = frailtysd, params_list = true_params, dat0 = dat0, dat1 = dat1, n = n
         )


pseudo_data = sim_data(n = 760, rhot = 0.5, rhos = 0.5, rhost = 0.5, scenario = 2, effecttheta = T, frailtysd = 0.5, effectsize = 0.3)
pseudo_data0 = pseudo_data$dat0
pseudo_data1 = pseudo_data$dat1

prepdata = function(n){
 
 rtog = read.csv("~/Dropbox (University of Michigan)/rtog9601/NCT00002874-D1-Dataset.csv")
 rtog = rtog[rtog$include_in_analysis == 1,]
 ST = (cbind(rtog$metastatic_prostate_cancer_years, rtog$survival_years - rtog$metastatic_prostate_cancer_years,
       rtog$survival_years))
 status = cbind(as.numeric(rtog$metastatic_prostate_cancer == 1), as.numeric(rtog$metastatic_prostate_cancer == 1 &
                                        rtog$survival == 1), as.numeric(rtog$survival == 1 & rtog$metastatic_prostate_cancer == 2))
 
 ST = data.frame(cbind(ST, status))
 names(ST) = c("y12", "y23", "y13", "s12", "s23", "s13")
 
 ST$trt = rtog$rx - 1
 dat0 = data.frame(ST[ST$"trt" == 0,])
 dat1 = data.frame(ST[ST$"trt" == 1,])
 
 dat0 = dat0[1:(n/2), ]
 dat1 = dat1[1:(n/2), ]
 
 # run prentice criteria
 time23 = rowSums(cbind(ST[,1], ST[,2]), na.rm = T)
 time23[status[,1] == 0] = NA
 Ttime = apply(cbind(ST[,3], time23), 1, min, na.rm = T) 
 Tstat = apply(cbind(status[,3], status[,1]), 1, max, na.rm = T) 
 
 if(F){
 # 0 = censored, 1 = outcome of interest, 2 = competing event
 a = cmprsk::cuminc(rtog$metastatic_prostate_cancer_years, fstatus = rtog$metastatic_prostate_cancer, 
           cencode = 0)
 plot(a)

 ci_fit = 
  cmprsk::cuminc(
   ftime = rtog$metastatic_prostate_cancer_years, 
   fstatus = rtog$metastatic_prostate_cancer, 
   group = rtog$rx - 1,
   cencode = 0
  )

 
 # cumulative incidence of T
 a = cmprsk::cuminc(rtog$survival_years, fstatus = rtog$survival, group = rtog$rx, cencode = 0)
 plot(a)
 
 # cumulative incidence of S
 mstatcumul = mstate::Cuminc(time = rtog$metastatic_prostate_cancer_years, status = rtog$metastatic_prostate_cancer , data = rtog)
 
 plot(mstatcumul, conf.int = .95, xlab = "Time (years)",use.ggplot = TRUE, main = "Cumulative Incidence for Surrogate Endpoint S", lty = 1:4)
 
 ci = mstate::Cuminc(time = rtog$metastatic_prostate_cancer_years, status = rtog$metastatic_prostate_cancer , data = rtog, group = rtog$rx , failcodes = c(1,2))
 
 ci1 = ci[ci$group == 1,]
 ci2 = ci[ci$group == 2,]
 
 plot(c(0, ci1$time, 15), c(0, ci1$CI.1, max(ci1$CI.1)),
    type = 's', xlim = c(0, 14), ylim = c(0,
                       1), lty = 'dashed', ylab = "Probability", xlab = "Years from Randomization",
    main = "Estimates based on the cumulative incidence functions\nBlack dashed line indicates z = 1, solid red line indicates z = 0",
    lwd = 2); lines(c(0, ci2$time, 15), c(0, ci2$CI.2, max(ci2$CI.2)),
            type = 's', col = 'red')
 }
 
 t_mod = coxph(Surv(Ttime, Tstat) ~ trt, dat = ST)
 ts_mod = coxph(Surv(Ttime, Tstat) ~ y12 + trt, dat = ST)
 s_mod = coxph(Surv(y12, s12) ~ trt, dat = ST)
 
 as.numeric(summary(s_mod)$coefficients[,"Pr(>|z|)"] < 0.05)
 summary(t_mod)
 summary(ts_mod)
 
 exp(-coef(ts_mod)["trt"])
 exp(-coef(t_mod)["trt"])
 ST$id = 1:(nrow(ST))

 ST$Ttime = Ttime
 ST$Tstat = Tstat
 
 new = tmerge(data1 = ST, data2 = ST, id = id, tstop = rtog$survival_years)
 new2 = tmerge(data1 = new, data2 = ST, id = id, ev = event(y12))
 
 ts_mod2 = coxph(Surv(Ttime, Tstat) ~ y12 + trt, dat = new2, cluster = id)
 exp(-coef(ts_mod2)["trt"])
 
 # KM plots
 
 rtog$metastatic_prostate_cancer_cause = as.numeric(rtog$metastatic_prostate_cancer == 1)
 rtog$metastatic_prostate_cancer_years_comp = as.numeric(rtog$metastatic_prostate_cancer_years != 0)

 theme_set( theme_classic(base_size = 16))

 a = ggsurvplot(survfit(Surv(rtog$metastatic_prostate_cancer_years, rtog$metastatic_prostate_cancer_cause) ~ rtog$rx), data = rtog, risk.table = F, ggtheme = theme_classic2(base_size=20), legend.title = "Treatment", legend.labs=c("With antiandrogen\ntherapy","Without antiandrogen\ntherapy ")) + ggtitle("KM Curve of Intermediate Outcome S\n(Individuals are Censored at T)") + ylab("Freedom from S") + xlab("Time (years)")
 a

 b = ggsurvplot(survfit(Surv(rtog$survival_years, rtog$survival) ~ rtog$rx), data = rtog, risk.table = F, ggtheme = theme_classic2(base_size=20), legend.title = "Treatment", legend.labs=c("With antiandrogen\ntherapy","Without antiandrogen\ntherapy ")) + ggtitle("KM Curve of True Outcome T") + ylab("Freedom from T") + xlab("Time (years)")
 b

 rtogs = rtog[rtog$metastatic_prostate_cancer == 1,]
 d = ggsurvplot(survfit(Surv((rtogs$survival_years - rtogs$metastatic_prostate_cancer_years), rtogs$survival) ~ rtogs$rx), data = rtogs, risk.table = F, ggtheme = theme_classic2(base_size=20), legend.title = "Treatment", legend.labs=c("With antiandrogen\ntherapy","Without antiandrogen\ntherapy ")) + ggtitle("KM Curve of Time between S to T\nFor Those who Experienced S") + ylab("Freedom from T Survival (Post S)") + xlab("Time (years)")
 d
 
 if(F){
  ggsave(file = "KM_S.jpeg", a$plot, width = 8, height = 8)
  ggsave(file = "KM_T.jpeg", b$plot, width = 8, height = 8)
  ggsave(file = "KM_ST.jpeg", d$plot, width = 8, height = 8)
 }
 
 return(list(dat0, dat1, n))
}

dat = prepdata(n = 376*2)

params_res_data = run_sim(SIM = SIM, rhos = rhos, rhot = rhot, frailtysd = frailtysd, params_list = true_params,
             dat0 = dat[[1]], dat1 = dat[[2]], n = 376*2)

plot_traceplots = function(params_matrix, variable){
 param = params_matrix$params
 
 plot(eval(parse(text = paste0("param$", variable))), ylab = "Parameter Draw",
    xlab = "MCMC Iteration", main = paste("Traceplot of Parameter", variable))
}


plot_traceplots(params_matrix = params_res, variable = "int")
plot_traceplots(params_matrix = params_res, variable = "slope")
plot_traceplots(params_matrix = params_res, variable = "scale23_1")
plot_traceplots(params_matrix = params_res, variable = "scale13_1")

plot_traceplots(params_matrix = params_res, variable = "scale23_0")
plot_traceplots(params_matrix = params_res, variable = "scale13_0")

plot_traceplots(params_matrix = params_res, variable = "scale12_0")
plot_traceplots(params_matrix = params_res, variable = "scale12_1")

plot_traceplots(params_matrix = params_res, variable = "theta23_0")
plot_traceplots(params_matrix = params_res, variable = "theta23_1")


plot_traceplots(params_matrix = params_res_data, variable = "int")
plot_traceplots(params_matrix = params_res_data, variable = "slope")


final_results = function(params_matrix, write){
 ## save results
 param = params_matrix$params
 
 n = params_matrix$args$n
 burnin = params_matrix$args$burnin
 sim = params_matrix$args$SIM
 
params=data.frame(c12_0 = numeric(1), c12_0SE = numeric(1), 
          c13_0 = numeric(1), c13_0SE = numeric(1), 
          c23_0 = numeric(1), c23_0SE = numeric(1), 
          c12_1 = numeric(1), c12_1SE = numeric(1), 
          c13_1 = numeric(1), c13_1SE = numeric(1), 
          c23_1 = numeric(1), c23_1SE = numeric(1), 
          theta23_0 = numeric(1), theta23_0SE = numeric(1), 
          theta23_1 = numeric(1), theta23_1SE = numeric(1), 
          scale12_0 = numeric(1), scale12_0SE = numeric(1), 
          scale13_0 = numeric(1), scale13_0SE = numeric(1), 
          scale23_0 = numeric(1), scale23_0SE = numeric(1), 
          scale12_1 = numeric(1), scale12_1SE = numeric(1), 
          scale13_1 = numeric(1), scale13_1SE = numeric(1), 
          scale23_1 = numeric(1), scale23_1SE = numeric(1),
          
          shape12_0 = numeric(1), shape12_0SE = numeric(1),
          shape13_0 = numeric(1), shape13_0SE = numeric(1),
          shape23_0 = numeric(1), shape23_0SE = numeric(1),
          shape12_1 = numeric(1), shape12_1SE = numeric(1),
          shape13_1 = numeric(1), shape13_1SE = numeric(1),
          shape23_1 = numeric(1), shape23_1SE = numeric(1),
          beta12_0 = numeric(1), beta12_0SE = numeric(1), 
          beta13_0 = numeric(1), beta13_0SE = numeric(1), 
          beta23_0 = numeric(1), beta23_0SE = numeric(1), 
          beta12_1 = numeric(1), beta12_1SE = numeric(1), 
          beta13_1 = numeric(1), beta13_1SE = numeric(1), 
          beta23_1 = numeric(1), beta23_1SE = numeric(1)
          # accept1 = numeric(1), accept2= numeric(1),
          # accept3= numeric(1), accept4= numeric(1),
          # accept5= numeric(1), accept6= numeric(1),
          # accept7= numeric(1), accept8= numeric(1),
          # accept9= numeric(1), accept10 = numeric(1),
          # accept10 = numeric(1), accept12= numeric(1)
)

 params[1]=mean(params_matrix$params$c12_0[burnin:(SIM-1)], na.rm = T)
 params[2]= sqrt(var(params_matrix$params$c12_0[burnin:(SIM-1)], na.rm = T))
 params[3]=mean(params_matrix$params$c13_0[burnin:(SIM-1)], na.rm = T)
 params[4]= sqrt(var(params_matrix$params$c13_0[burnin:(SIM-1)], na.rm = T))
 params[5]=mean(params_matrix$params$c23_0[burnin:(SIM-1)], na.rm = T)
 params[6]= sqrt(var(params_matrix$params$c23_0[burnin:(SIM-1)], na.rm = T))
 params[7]=mean(params_matrix$params$c12_1[burnin:(SIM-1)], na.rm = T)
 params[8]= sqrt(var(params_matrix$params$c12_1[burnin:(SIM-1)], na.rm = T))
 params[9]=mean(params_matrix$params$c13_1[burnin:(SIM-1)], na.rm = T)
 params[10]= sqrt(var(params_matrix$params$c13_1[burnin:(SIM-1)], na.rm = T))
 params[11]=mean(params_matrix$params$c23_1[burnin:(SIM-1)], na.rm = T)
 params[12]= sqrt(var(params_matrix$params$c23_1[burnin:(SIM-1)], na.rm = T))
 params[13]=mean(params_matrix$params$theta23_0[burnin:(SIM-1)], na.rm = T)
 params[14]= sqrt(var(params_matrix$params$theta23_0[burnin:(SIM-1)], na.rm = T))
 params[15]=mean(params_matrix$params$theta23_1[burnin:(SIM-1)], na.rm = T)
 params[16]= sqrt(var(params_matrix$params$theta23_1[burnin:(SIM-1)], na.rm = T))
 params[17]=mean(params_matrix$params$scale12_0[burnin:(SIM-1)], na.rm = T)
 params[18]= sqrt(var(params_matrix$params$scale12_0[burnin:(SIM-1)], na.rm = T))
 params[19]=mean(params_matrix$params$scale13_0[burnin:(SIM-1)], na.rm = T)
 params[20]= sqrt(var(params_matrix$params$scale13_0[burnin:(SIM-1)], na.rm = T))
 params[21]=mean(params_matrix$params$scale23_0[burnin:(SIM-1)], na.rm = T)
 params[22]= sqrt(var(params_matrix$params$scale23_0[burnin:(SIM-1)], na.rm = T))
 params[23]=mean(params_matrix$params$scale12_1[burnin:(SIM-1)], na.rm = T)
 params[24]= sqrt(var(params_matrix$params$scale12_1[burnin:(SIM-1)], na.rm = T))
 params[25]=mean(params_matrix$params$scale13_1[burnin:(SIM-1)], na.rm = T)
 params[26]= sqrt(var(params_matrix$params$scale13_1[burnin:(SIM-1)], na.rm = T))
 params[27]=mean(params_matrix$params$scale23_1[burnin:(SIM-1)], na.rm = T)
 params[28]= sqrt(var(params_matrix$params$scale23_1[burnin:(SIM-1)], na.rm = T))
 
 params[29]=mean(params_matrix$params$shape12_0[burnin:(SIM-1)], na.rm = T)
 params[30]= sqrt(var(params_matrix$params$shape12_0[burnin:(SIM-1)], na.rm = T))
 params[31]=mean(params_matrix$params$shape13_0[burnin:(SIM-1)], na.rm = T)
 params[32]= sqrt(var(params_matrix$params$shape13_0[burnin:(SIM-1)], na.rm = T))
 params[33]=mean(params_matrix$params$shape23_0[burnin:(SIM-1)], na.rm = T)
 params[34]= sqrt(var(params_matrix$params$shape23_0[burnin:(SIM-1)], na.rm = T))
 params[35]=mean(params_matrix$params$shape12_1[burnin:(SIM-1)], na.rm = T)
 params[36]= sqrt(var(params_matrix$params$shape12_1[burnin:(SIM-1)], na.rm = T))
 params[37]=mean(params_matrix$params$shape13_1[burnin:(SIM-1)], na.rm = T)
 params[38]= sqrt(var(params_matrix$params$shape13_1[burnin:(SIM-1)], na.rm = T))
 params[39]=mean(params_matrix$params$shape23_1[burnin:(SIM-1)], na.rm = T)
 params[40]= sqrt(var(params_matrix$params$shape23_1[burnin:(SIM-1)], na.rm = T))
 
 params[41]=mean(params_matrix$params$beta12_0[burnin:(SIM-1)], na.rm = T)
 params[42]= sqrt(var(params_matrix$params$beta12_0[burnin:(SIM-1)], na.rm = T))
 params[43]=mean(params_matrix$params$beta13_0[burnin:(SIM-1)], na.rm = T)
 params[44]= sqrt(var(params_matrix$params$beta13_0[burnin:(SIM-1)], na.rm = T))
 params[45]=mean(params_matrix$params$beta23_0[burnin:(SIM-1)], na.rm = T)
 params[46]= sqrt(var(params_matrix$params$beta23_0[burnin:(SIM-1)], na.rm = T))
 params[47]=mean(params_matrix$params$beta12_1[burnin:(SIM-1)], na.rm = T)
 params[48]= sqrt(var(params_matrix$params$beta12_1[burnin:(SIM-1)], na.rm = T))
 params[49]=mean(params_matrix$params$beta13_1[burnin:(SIM-1)], na.rm = T)
 params[50]= sqrt(var(params_matrix$params$beta13_1[burnin:(SIM-1)], na.rm = T))
 params[51]=mean(params_matrix$params$beta23_1[burnin:(SIM-1)], na.rm = T)
 params[52]= sqrt(var(params_matrix$params$beta23_1[burnin:(SIM-1)], na.rm = T))
 # 
 # params[53]=accept1/z
 # params[54]=accept2/z
 # params[55]=accept3/z
 # params[56]=accept4/z
 # params[57]=accept5/z
 # params[58]=accept6/z
 # params[59]=accept7/z
 # params[60]=accept8/z
 # params[61]=accept9/z
 # params[62]=accept10/z
 # params[63]=accept11/z
 # params[64]=accept12/z
 # 
 params$int = mean(params_matrix$params$int, na.rm = T)
 params$intse = sd(params_matrix$params$int, na.rm = T)
 params$slope = mean(params_matrix$params$slope, na.rm = T)
 params$slopese = sd(params_matrix$params$slope, na.rm = T)

 print(params[c(T,F)])

 if(write){fname = paste('params',",", scenario, ",", array_id,'.txt',sep="")
 write.table(params, file=fname, sep="\t", row.names=F, col.names=T)}

}



plot_results = function(params_matrix, write){
 
 param = params_matrix$params
 n = params_matrix$args$n
 burnin = params_matrix$args$burnin
 sim = params_matrix$args$SIM
 
 if(sim - burnin > 1000) burnin = sim - 1000
 
 saveCEPx = params_matrix$saveCEPx
 saveCEPy = params_matrix$saveCEPy
 
 dat = data.frame(matrix(data = NA, nrow = n, ncol = 1))
 dat$X = rowMeans(saveCEPx[,burnin:ncol(saveCEPx)], na.rm = T)
 dat$Y = rowMeans(saveCEPy[,burnin:ncol(saveCEPx)], na.rm = T)

 slope = params_matrix$params$slope
 int = params_matrix$params$int

 theme_set( theme_classic(base_size = 26))

 d = ggplot(dat, aes(X, Y), col(c(rep(0, n/2), rep(1, n/2)))) + ylim(c(-1, 1)) +
 ggtitle("Illness-Death CEP Curve" ) + xlab("Delta S_i")+ 
  geom_hline(yintercept = mean(dat$Y, na.rm = T), linetype = "dashed", col = 'red') + 
 ylab(TeX("$\\Delta T_i = P(T_i(1) > \\tau_T| \\omega_{.i}^1) - P(T_i(0) > \\tau_T | \\omega_{.i}^0)$$ ")) +
 xlab(TeX("$\\Delta S_i = log \\frac{\\Lambda_{12}^{0} (\\tau_S| \\omega_{12i}^0)}{\\Lambda_{12}^{1} (\\tau_S| \\omega_{12i}^1)}$$ ")) + 
 geom_hline(yintercept = 0, linetype = "solid") + geom_vline(xintercept = 0, linetype = "solid") + 
 coord_cartesian(ylim=c(-1,1)) + theme_bw((base_size = 18)) + geom_point() 

 for(i in burnin:sim) d = d + geom_abline(slope = slope[i], intercept = int[i], alpha = 0.1, col = 'gray')
 d = d + geom_vline(xintercept = mean((dat$X), na.rm = T), linetype = "dashed", col = 'red') +
  geom_hline(yintercept = mean((dat$Y), na.rm = T), linetype = "dashed", col = 'red') 
 d3 = d + geom_point(aes(color = 1)) + theme(legend.position = "none")

 print(d3)
 if(write) ggsave(paste0("estimatedCEP", ",", scenario, ",", array_id, ".jpeg"), d3)

}


final_results(params_matrix = params_res, write = write)

final_results(params_matrix = params_res_data, write = write)

plot_results(params_matrix = params_res, write = plotwrite)

plot_results(params_matrix = params_res_data, write = plotwrite)
