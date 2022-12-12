#' simulate data for simulations
#'
#' @description simulate data
#'
#' @param n sample size
#' @param array_id ID of simulation (used for running in parallel)
#' @param scenario data setting for which treatment effects exist between transitions
#' @param effecttheta logical if theta23 should be zero
#' @param rhos correlation term
#' @param rhot correlation term
#' @param rhost correlation term
#' @param frailtysd standard deviation of generated frailties
#' @param diffscale1323 logical value if baseline hazards for 1-3 and 2-3 should be equal
#'
#' @return datasets dat0, dat1 and true frailties for use in other functions
#'
#' @examples 
#' equalfrail = TRUE
#' independent = TRUE
#' frailtysd = 0.5
#' example(sim_data(n = 600, array_id = 1, scenario = 2, effecttheta = FALSE,
#' rhos = 0.5, rhot = 0.5, rhost = 0.5, frailtysd = frailtysd, 
#' diffscale1323 = TRUE))
     
sim_data = function(n, array_id, scenario, effecttheta, frailtysd, rhot, rhos, rhost, diffscale1323){
	
	independent = FALSE
	equalfrail = F
  
  {
    if(scenario == 1){ effect12 = F; effect13 = F ; effect23 = F}
    if(scenario == 2){ effect12 = T; effect13 = F ; effect23 = F}
    if(scenario == 3){ effect12 = T; effect13 = F ; effect23 = T}
    if(scenario == 4){ effect12 = T; effect13 = T ; effect23 = F}
    if(scenario == 5){ effect12 = T; effect13 = T ; effect23 = T}
    if(scenario == 6){ effect12 = F; effect13 = T ; effect23 = T}
    if(scenario == 7){ effect12 = F; effect13 = F ; effect23 = T}
    if(scenario == 8){ effect12 = F; effect13 = T ; effect23 = F}
  }
  
  beta12_1<-0; beta13_1<-0; beta23_1<-0; beta12_0<-0; beta13_0<-0; beta23_0<-0
  theta23_1 = 0
  theta23_0 = 0
  
  if(effecttheta){theta23_1 = -.1
  theta23_0 = -.1
  }
  
  c13_0 = c13_1 = 1
  c12_1 = c12_0 = 1
  c23_1 = c23_0 = 1
  
  scale12_1<-1 #Weibull scale for latent illness time
  shape12_1<-1 #Weibull shape for latent illness time
  scale13_1<-0.5 #Weibull scale for latent life time
  
  if(!diffscale1323){scale13_1_star = scale13_1<-1} #Weibull scale for latent life time
  shape13_1<-1 #Weibull shape for latent life time
  scale23_1<-1 #Weibull scale for latent wait time
  shape23_1<-1 #Weibull shape for latent wait time'
  
  shape12_0<-shape12_1
  scale12_0<-scale12_1 #Weibull scale for latent illness time
  if(effect12){ scale12_1 = 0.61 }
  
  scale13_0 <-scale13_1 #Weibull scale for latent life time
  
  if(effect13){  scale13_1 = scale13_1 * (.61) }
  shape13_0<-shape13_1 #Weibull shape for latent life time
  scale23_0<-scale23_1 #Weibull scale for latent wait time
  shape23_0<-shape23_1 #Weibull shape for latent wait time
  if(effect23){  scale23_1 = 0.61 }
  
  # generate data
  R<-matrix(rep(1,6*6),6,6); 
  
  if(independent){  R[1,2] = R[2,1] = rhos
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
  eigen(R)
  
  S = diag(c(frailtysd, frailtysd, frailtysd, frailtysd, frailtysd, frailtysd))
  o = mvtnormrmvnorm(n, c(0,0,0,0,0,0), S%*%R%*%S)
  omega12true0 = o[,1]; omega12true1 = o[,2]
  
  omega13true0 = o[,3]
  omega23true0 = o[,5] 
  omega13true1 = o[,4]
  omega23true1 = o[,6]
  
  o = mvtnormrmvnorm(n, c(0,0,0,0,0,0), S%*%R%*%S)
  
  if(equalfrail){omega23true0 = omega13true0 ; omega23true1 = omega13true1}
  #omega13true0 = omega23true0 = omega12true0; omega13true1 = omega23true1 = omega12true1
  
  U1 = runif(n); U2 = runif(n); U3 = runif(n)
  U4 = runif(n); U5 = runif(n); U6 = runif(n)
  
  # U1 = U1/2
  # U2 = U2/2
  # U3 = U3/2
  xtrue = x = rbinom(n, 1, 0.5)
  x_0 = x[1(n/2)]; x_1 = x[1(n/2)]
  
  # columns S12, S23, S13 # 1 -> S, S -> T, 1 -> T
  ST = status = ST_raw = cbind(
    (1/scale12_0*(-log(1-U1))/exp((c12_0*omega12true0 + beta12_0 * x)))^shape12_0,
    (1/scale23_0*(-log(1-U2))/exp((c23_0*omega23true0 + beta23_0 * x)))^shape23_0,
    (1/scale13_0*(-log(1-U3))/exp((c13_0*omega13true0 + beta13_0 * x)))^shape13_0,
    (1/scale12_1*(-log(1-U4))/exp((c12_1*omega12true1 + beta12_1 * x)))^shape12_1,
    (1/scale23_1*(-log(1-U5))/exp((c23_1*omega23true1 + beta23_1 * x)))^shape23_1,
    (1/scale13_1*(-log(1-U6))/exp((c13_1*omega13true1 + beta13_1 * x)))^shape13_1
  ); status[,16] = 1
  
  
  ST[,2] =   (1/scale23_0*(-log(1-U2))/exp((c23_0*omega23true0 + theta23_0 * ST[,1])))^shape23_0
  ST[,5] =   (1/scale23_1*(-log(1-U5))/exp((c23_1*omega23true1 + theta23_1 *  ST[,4])))^shape23_1
  
  rgamma(1, shape = 0.01 + sum(ST[,2], na.rm = T),
         scale = (0.01 + 1/1 * sum(status[,2]^ 1 * exp(1 * omega23true0), na.rm = T))^(-1))
  
  
  rgamma(1, shape = 0.01 + sum(ST[,2], na.rm = T),
         scale = (0.01 + 1/1 * sum(status[,2]^ 1 * exp(0), na.rm = T))^(-1))
  
  
  if(T){
    ST[is.infinite(ST)] = NA
  
  status[ST[,3] < ST[,1], 12] = 0
  status[ST[,6] < ST[,4], 45] = 0 

  ST[ST[,3] < ST[,1], 1] =  ST[ST[,3] < ST[,1], 3]
  ST[ST[,6] < ST[,4], 4] = ST[ST[,6] < ST[,4], 6]
  
  ST[ST[,3] <= ST[,1], 2] = NA
  ST[ST[,6] <= ST[,4], 5] = NA
  
  status[ST[,3] > ST[,1], 3] = 0 # censor baseline -> S -> T
  status[ST[,6] > ST[,4], 6] = 0 
  
  ST[ST[,3] > ST[,1], 3] = ST[ST[,3] > ST[,1], 1]
  ST[ST[,6] > ST[,4], 6] = ST[ST[,6] > ST[,4], 4]
  
  ST[which(status[,1] == 0), 2] = NA
  ST[which(status[,4] == 0), 5] = NA
  
  status[which(status[,1] == 0), 2] = 0
  status[which(status[,4] == 0), 5] = 0

  }
  trt = c(rep(0, n), rep(1, n))
  
  dat0 = data.frame(y12 = c(ST[1(n/2),1]), s12 = c(status[1(n/2),1]),
                    y13 = c(ST[1(n/2),3]), s13 = c(status[1(n/2),3]),
                    y23 = c(ST[1(n/2),2]), s23 = c(status[1(n/2),2])
  )
  
  dat1 = data.frame(y12 = c(ST[1(n/2),4]), s12 = c(status[1(n/2),4]),
                    y13 = c(ST[1(n/2),6]), s13 = c(status[1(n/2),6]),
                    y23 = c(ST[1(n/2),5]), s23 = c(status[1(n/2),5])
  )
  
  o12save0 = omega12_z0 = omega12true0[1(n/2)]
  o13save0 = omega13_z0 = omega13true0[1(n/2)]
  o23save0 = omega23_z0 = omega23true0[1(n/2)]
  
  o12save1 = omega12_z1 = omega12true1[1(n/2)]
  o13save1 = omega13_z1 = omega13true1[1(n/2)]
  o23save1 = omega23_z1 = omega23true1[1(n/2)]
  
  params_list = data.frame(scale12_0 = scale12_0, scale13_0 = scale13_0, scale23_0 = scale23_0,
                           scale12_1 = scale12_1, scale13_1 = scale13_1, scale23_1 = scale23_1,
                           shape12_0 = shape12_0, shape13_0 = shape13_0, shape23_0 = shape23_0,
                           shape12_1 = shape12_1, shape13_1 = shape13_1, shape23_1 = shape23_1,
                           theta23_0 = theta23_0, theta23_1 = theta23_1,
                           beta12_0 = beta12_0,  beta13_0 = beta13_0,  beta23_0 = beta23_0, 
                           beta12_1 = beta12_1,  beta13_1 = beta13_1,  beta23_1 = beta23_1,
                           c12_0 = c12_0, c13_0 = c13_0, c23_0 = c23_0, c12_1 = c12_1, c13_1 = c13_1, c23_1 = c23_1
  )
  
  return(list( dat0 = dat0, dat1 = dat1, o12save0 = omega12true0, o12save1 = omega12true1, o13save0 = omega13true0, 
               o13save1 = omega13true1, o23save0 = omega23true0, o23save1 = omega23true1, params = params_list))
}

