#' run the simulation
#'
#' @description run mcmc simulation
#'
#' @param SIM number of iterations of mcmc to run
#' @param rhos correlation parameter
#' @param rhot correlation parameter
#' @param params_list list of true generating parameters
#' @param frailtysd standard deviation of frailty term
#'
#' @return simulation results
#'
#' @examples 
#' example(run_sim(SIM = 3000, rhos = 0.5, rhot = 0.5, 
#' frailtysd = 1, params_list = params_list))

run_sim = function(SIM, rhos, rhot, frailtysd, params_list){
  
  burnin = SIM * .3
  
  ## allocate memory to save parameter draws
  accept1 = accept2 = accept3 = accept4 = accept5 = accept6 = accept7 = accept8 = accept9 = accept10 = accept11 = accept12 = 0
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
  
  holdscale13_0[1] = 
  holdscale13_1[1] = 
  holdscale12_0[1] = 
  holdshape12_0[1] = 
  holdshape13_0[1] = 
  holdscale23_0[1] = 
  holdshape23_0[1] = 
  holdshape12_1[1] = 
  holdshape23_1[1] = 
  holdshape13_1[1] = 
  holdscale23_1[1] = 
  holdscale12_1[1] = 1

  holdc12_0[1] = 
  holdc13_0[1] =  
  holdc12_1[1] = 
  holdc13_1[1] = 
  holdc23_1[1] = 
  holdc23_0[1] = 1

  holdtheta23_1[1] = 0
  holdtheta23_0[1] = 0
  
  holdbeta23_1[1] = 
  holdbeta13_1[1] = 
  holdbeta12_0[1] = 
  holdbeta12_1[1] =  
  holdbeta13_0[1] = 
  holdbeta23_0[1] = 0
  
  saveCEPx = matrix(data = NA, nrow = n, ncol = SIM)
  saveCEPy = matrix(data = NA, nrow = n, ncol = SIM)
  
# estimate starting values for frailty terms
  if(T){
    dat1penal = data.frame(time = c(dat1$y12),  status = c(dat1$s12),
                           id = rep(1  (n/2), 1), death = c(dat1$s13))
    dat1penal = dat1penal[complete.cases(dat1penal),]

    penalfrail1 = frailtyPenal(Surv(time, status)~ cluster(id) + terminal(death),
                           jointGeneral = TRUE, data=dat1penal, n.knots = 4, kappa = c(1))
    
    omega12_z1 = o12save1[1  (n/2)] = log(penalfrail1$frailty.pred)

    dat1penal = data.frame(time = c(dat1$y23, dat1$y13),  status = c(dat1$s23, dat1$s13), y12 = c(dat1$y12, rep(0, n/2)), 
                           id = rep(1  (n/2), 2), death = c(rep(0, n/2)), dat1$s13)
    dat1penal = dat1penal[complete.cases(dat1penal),]


    dat1penal = data.frame(time = c(dat1$y23, dat1$y13),  status = c(dat1$s23, dat0$s13), y12 = c(dat1$y12, rep(1, n/2)), 
                           id = rep(1  (n/2), 2))
    dat1penal = dat1penal[complete.cases(dat1penal),]
    
    penalfrail1 = frailtyPenal(Surv(time, status)~cluster(id) + y12 ,
                               data=dat1penal, hazard = "Weibull",
                               RandDist = "LogN")
    
    omega23_z1 = o23save1[1  (n/2)] = omega13_z1 = o13save1[1  (n/2)] = (penalfrail1$frailty.pred)

    dat0penal = data.frame(time = c(dat0$y12),  status = c(dat0$s12),
                           id = rep(1  (n/2), 1), death = c(dat0$s13))

    dat0penal = dat0penal[complete.cases(dat0penal),]
    
    penalfrail0 = frailtyPenal(Surv(time, status)~ cluster(id) + terminal(death),
                               jointGeneral = TRUE, data=dat0penal, n.knots = 4, kappa = c(1))
    
    omega12_z0 = o12save0[1  (n/2)] = log(penalfrail0$frailty.pred)
      
    dat0penal = data.frame(time = c(dat0$y23, dat0$y13),  status = c(dat0$s23, dat0$s13), y12 = c(dat0$y12, rep(0, n/2)), 
                           id = rep(1  (n/2), 2))
    dat0penal = dat0penal[complete.cases(dat0penal),]

    penalfrail0 = frailtyPenal(Surv(time, status)~cluster(id) + y12 ,
                               data=dat0penal, hazard = "Weibull",
                               RandDist = "LogN")
    
    omega23_z0 = o23save0[1  (n/2)] = omega13_z0 = o13save0[1  (n/2)] = (penalfrail0$frailty.pred)
    
      }
  
  mod23_0 = summary(survreg(Surv(dat0$y23, dat0$s23) ~ dat0$y12 + omega23_z0, dist = "exponential"))
  mod23_1 = summary(survreg(Surv(dat1$y23, dat1$s23) ~ dat1$y12 + omega23_z1, dist = "exponential"))
  
  z = 2
  x_0 = x_1 = rep(0, n/2) # no covariates 
  
  	if(holdtheta){ # options for holding certain values fixed during simulation
     holdtheta23_0[1] = theta23_0 = params_list$theta23_0
     holdtheta23_1[1] = theta23_1 = params_list$theta23_1
  	}
  
  	if(holdscale12){
     scale12_0 = params_list$scale12_0
     scale12_1 = params_list$scale12_1
  	}
  
    if(holdscale23){
      scale23_0 = params_list$scale23_0
      scale23_1 = params_list$scale23_1
    }
    
    if(holdscale13){
      scale13_0 = params_list$scale13_0
      scale13_1 = params_list$scale13_1
    }
    
  	if(holdshape){
     shape12_0 = params_list$shape12_0
     shape12_1 = params_list$shape12_1
     shape13_0 = params_list$shape13_0
     shape13_1 = params_list$shape13_1
     shape23_0 = params_list$shape23_0
     shape23_1 = params_list$shape23_1
  	}
  
  while(z < SIM){ 

  zseed = round(runif(1, 1, 100)) ;  set.seed(zseed + z + array_id)
  
  omega12_star = rnorm(n/2, o12save0[1  (n/2)], sd = proposalsdfrail)
  omega13_star = rnorm(n/2, o13save0[1  (n/2)], sd = proposalsdfrail)

   if(z == 2){
     omega12_star = rnorm(n/2, o12save0[1  (n/2)], sd = frailtysd)
     omega13_star = rnorm(n/2, o13save0[1  (n/2)], sd = frailtysd)
   }
  
  omega12_z0 = o12save0[1  (n/2)]
  omega13_z0 = o13save0[1  (n/2)]
  omega23_z0 = o23save0[1  (n/2)]
  
  save = savestar = NULL
  sap = sap2 = z1 = p = rep(NA, n/2)
  
  for(i in 1:  (n/2)){
    
    sap[i] = like12_omega_i(omega12_star[i], i = i, scale12_0 = holdscale12_0[z - 1],
                              dat0 = dat0, shape12_0 = holdshape12_0[z - 1], c12_0 = holdc12_0[z - 1])
    sap2[i] = like12_omega_i(omega12_z0[i], i = i, scale12_0 = holdscale12_0[z - 1],
                               dat0 = dat0, shape12_0 = holdshape12_0[z - 1], c12_0 = holdc12_0[z - 1])

 	p[i] = exp(sap[i] - sap2[i] + sum(log(pnorm(omega12_star, 0, sd = frailtysd)))-
                                  sum(log(pnorm(omega12_z0, 0, sd = frailtysd))))
    
    if(is.na(p[i])) p[i] = 0; if((p[i] < 0)) p[i] = 0; if((p[i] > 1)) p[i] = 1
    z1[i] = rbinom(1, 1, prob = p[i])
    if(holdfrail12) z1[i] = 0
  }
  
  if(sum(z1) >= 1) accept1 = accept1 + 1

  omega12_star = omega12_z0 = omega_star = z1 * omega12_star + (1-z1) * omega12_z0
  holdfrailsd12_0[z] = sd(omega12_star)
  holdfrailmean12_0[z] = mean(omega12_star)
  
  omega_cf = rnorm(n/2, omega_star, sd = frailtysd * sqrt(frailtysd - rhos^2))
  
  o12save0 = c(omega_star, omega_cf)
  
  holdomega12_0[,z] = omega12_z0[1 : 10]
  if(holdc) c12_0_star = 1
  beta12_0_star = rnorm(1, holdbeta12_0[z - 1], sd = proposalsd)
  shape12_0_star = rtruncnorm(n = 1, a = 0, b = 3, mean = holdshape12_0[z - 1], sd = proposalsd)
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
  
  c12prior_c = pnorm(holdc12_0[z - 1], 0, sd = 3)
  beta12prior_c = pnorm(holdbeta12_0[z - 1], 0, sd = 3)
  shape12prior_c = pnorm(holdshape12_0[z - 1], 0, sd = 3)
  shape12prior_c = ptruncnorm(holdshape12_0[z - 1], a = 0, b = 3, mean = 0, sd = 3)

    p = exp(like12_shape(shape12_0 = shape12_0_star, scale12_0 = scale12_0_star, omega12_z0 = omega12_z0,
            dat0 = dat0, c12_0 = c12_0_star) + log(shape12prior_p) - 
            (like12_shape(shape12_0 = holdshape12_0[z - 1], scale12_0 = scale12_0_star, omega12_z0 = omega12_z0,
            dat0 = dat0, c12_0 = c12_0_star)+ log(shape12prior_c)))

    if(is.nan(p)) p = 0; if(is.na(p)) p = 0; if((p < 0)) p = 0; if((p > 1)) p = 1
    if(scale12_0_star < 0) p = 0
    z1 = rbinom(1, 1, prob = p)
  
  if(z1 == 1) accept2 = accept2 + 1
    
  c12_0 = c12_0_star
  shape12_0 = z1 * shape12_0_star + (1-z1) * holdshape12_0[z - 1]
  beta12_0 = z1 * beta12_0_star + (1-z1) * holdbeta12_0[z - 1]
  holdbeta12_0[z] = beta12_0
  
  if(holdscale12) z1 = 1
  scale12_0 = z1 * scale12_0_star + (1-z1) * holdscale12_0[z - 1]
  
  holdscale12_0[z] = scale12_0
  holdc12_0[z] = c12_0
  holdshape12_0[z] = shape12_0

  z1 = p = rep(NA, (n/2))
  sap = sap2 = rep(NA, (n/2))
  
  for(i in 1:  (n/2)){
    sap[i] =   like13_omega_i(omega13_z0 = omega13_star[i], i = i, 
                              dat0 = dat0, c13_0 = holdc13_0[z - 1],
                              shape13_0 = holdshape13_0[z - 1], scale13_0 = holdscale13_0[z - 1], c23_0 = holdc23_0[z - 1],
                              shape23_0 = holdshape23_0[z - 1], scale23_0 = holdscale23_0[z - 1], theta23_0 = holdtheta23_0[z - 1])[1]
    sap2[i] =   like13_omega_i(omega13_z0[i], i = i, 
                               dat0 = dat0, c13_0 = holdc13_0[z - 1],
                               shape13_0 = holdshape13_0[z - 1], scale13_0 = holdscale13_0[z - 1], c23_0 = holdc23_0[z - 1],
                               shape23_0 = holdshape23_0[z - 1], scale23_0 = holdscale23_0[z - 1], theta23_0 = holdtheta23_0[z - 1])[1]
    
  	p[i]  = exp(sap[i] - sap2[i] + sum(log(pnorm(omega13_star, 0, sd = frailtysd)))-
                   sum(log(pnorm(omega13_z0, 0, sd = frailtysd)) ))
    
    if(is.na(p[i])) p[i] = 0; if((p[i] < 0)) p[i] = 0;  if((p[i] > 1)) p[i] = 1
    z1[i] = rbinom(1, 1, prob = p[i])
    if(holdfrail13) z1[i] = 0
  }
  
  if(sum(z1) >= 1) accept3 = accept3 + 1
  
  omega13_star = omega13_z0 = omega_star = z1 * omega13_star + (1-z1) * omega13_z0
  omega_cf = rnorm(n/2, omega_star, sd = 1 * sqrt(frailtysd - rhot^2))
  
  o13save0 = c(omega_star, omega_cf)
  
  holdfrailsd13_0[z] = sd(omega13_star)
  holdfrailmean13_0[z] = mean(omega13_star)

  omega23_z0 = omega23_star = omega13_z0
  
  holdomega23_0[,z] = omega23_z0[1 : 10]
  holdfrailsd23_0[z] = sd(omega23_star)
  holdfrailmean23_0[z] = mean(omega23_star)

  o23save0 = o13save0 
  
  c13_0_star =  rnorm(1, holdc13_0[z - 1], sd = proposalsd)
  c13_0_star= rtruncnorm(1, a = 0, b = 3, mean = holdc13_0[z - 1], sd = proposalsd)
  shape13_0_star = rtruncnorm(1, a = 0, b = 3, mean = holdshape13_0[z - 1], sd = proposalsd)
  if(holdshape) shape13_0_star = 1
  
  if(holdc) c13_0_star =  1
  beta13_0_star = rnorm(1, holdbeta13_0[z - 1], sd = proposalsd)
  beta13_0_star = 0
  scale13_0_star = rgamma(1, shape = 0.01 + sum(dat0$s13, na.rm = T),
                           scale = (0.01 + 1/shape13_0_star * sum(dat0$y13^shape13_0_star * exp(c13_0_star *
                           omega13_z0), na.rm = T))^(-1))
  
  if(holdscale13) scale13_0_star = scale13_0
  beta23_0_star = rnorm(1, holdbeta23_0[z - 1], sd = proposalsd)
  beta23_0_star = 0
  shape23_0_star = rnorm(1, holdshape23_0[z - 1], sd = proposalsd)
  shape23_0_star = rtruncnorm(1, a = 0, b = 3, mean = holdshape23_0[z - 1], sd = proposalsd)
  
  if(holdshape) shape23_0_star = 1
  
  theta23_0_star = rnorm(1, -mod23_0$coefficients[2], sd = proposalsd)
  theta23_0_star = rnorm(1, holdtheta23_0[z - 1], sd = proposalsd)
  
  if(holdtheta) theta23_0_star = theta23_0
  
  c23_0_star = rnorm(1, holdc23_0[z - 1], sd = proposalsd)
  c23_0_star = rtruncnorm(1, a=-3, b= 3, mean = holdc23_0[z - 1], sd = proposalsd)
  if(holdc23) c23_0_star = 1
  
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
  
  c13prior_c = pnorm(holdc13_0[z - 1], 0, sd = 3)
  beta13prior_c = pnorm(holdbeta13_0[z - 1], 0, sd = 3)
  shape13prior_c = pnorm(holdshape13_0[z - 1], 0, sd = 3)
  shape13prior_c = ptruncnorm(holdshape13_0[z - 1], a = 0, b = 3, mean = 0, sd = 3)
  
  c23prior_p = pnorm(c23_0_star, 0, sd = 3)
  beta23prior_p = pnorm(beta23_0_star, 0, sd = 3)
  theta23prior_p = pnorm(theta23_0_star, 0, sd = 3)
  shape23prior_p = pnorm(shape23_0_star, 0, sd = 3)
  shape23prior_p = ptruncnorm(shape23_0_star, a = 0, b = 3, mean = 0, sd = 3)
  
  c23prior_c = pnorm(holdc23_0[z - 1], 0, sd = 3)
  beta23prior_c = pnorm(holdbeta23_0[z - 1], 0, sd = 3)
  theta23prior_c = pnorm(holdtheta23_0[z - 1], 0, sd = 3)
  shape23prior_c = pnorm(holdshape23_0[z - 1], 0, sd = 3)
  shape23prior_c = ptruncnorm(holdshape23_0[z - 1], a = 0, b = 3, mean = 0, sd = 3)
  
  theta23prop_p = pnorm(theta23_0_star, -mod23_0$coefficients[2], sd = proposalsd)
  theta23prop_c = pnorm(holdtheta23_0[z - 1], -mod23_0$coefficients[2], sd = proposalsd)
  
  theta23prop_p = pnorm(theta23_0_star, 0, sd = proposalsd)
  theta23prop_c = pnorm(holdtheta23_0[z - 1], 0, sd = proposalsd)
  
   p = exp(like13_shape(dat0 = dat0, c13_0 = c13_0_star, shape13_0 = shape13_0_star, scale13_0 = scale13_0_star, omega13_z0 = omega13_z0) +
               log(shape13prior_p)-(like13_shape(dat0 = dat0, c13_0 = holdc13_0[z - 1],
               shape13_0 = holdshape13_0[z - 1], scale13_0 = scale13_0_star, omega13_z0 = omega13_z0)+ log(shape13prior_c)) )
    
   p2 = exp( like23_theta(theta23_0 = theta23_0_star, dat0 = dat0,c23_0 = holdc23_0[z - 1], omega13_z0 = omega13_z0,
               	  shape23_0 = holdshape23_0[z - 1], scale23_0 = holdscale23_0[z - 1]) + log(theta23prior_p) + log(theta23prop_p) -
                  (like23_theta(theta23_0 = holdtheta23_0[z - 1],dat0 = dat0, omega13_z0 = omega13_z0, c23_0 = holdc23_0[z - 1],
                  shape23_0 = holdshape23_0[z - 1], scale23_0 = holdscale23_0[z - 1])+ log(theta23prior_c) + log(theta23prop_c)) +
                  like23_shape(shape23_0 = shape23_0_star, dat0 = dat0, omega13_z0 = omega13_z0, c23_0 = holdc23_0[z - 1],
                  scale23_0 = holdscale23_0[z - 1], theta23_0 = holdtheta23_0[z - 1]) + log(shape23prior_p)-
                  (like23_shape(shape23_0 = holdshape23_0[z - 1], dat0 = dat0, c23_0 = holdc23_0[z - 1], omega13_z0 = omega13_z0, 
                  scale23_0 = holdscale23_0[z - 1], theta23_0 = holdtheta23_0[z - 1]) + log(shape23prior_c)))

   if(is.nan(p)) p = 0; if(is.na(p)) p = 0; if((p < 0)) p = 0;  if((p > 1)) p = 1
   if(is.nan(p2)) p2 = 0; if(is.na(p2)) p2 = 0; if((p2<0)) p2 = 0;  if((p2>1)) p2 = 1
    
   if(scale23_0_star < 0) {p2 = 0}
   if(scale13_0_star < 0) {p = 0}
   z1 = rbinom(1, 1, prob = p)
   z15 = rbinom(1, 1, prob = p2)
  
  if(z1 == 1) accept4 = accept4 + 1
  if(z15 == 1) accept5 = accept5 + 1

  holdomega13_0[,z] = omega13_z0[1 : 10]
  
  scale13_0 = z1 * scale13_0_star + (1-z1) *  holdscale13_0[z - 1]
  c13_0 = z1 * c13_0_star + (1-z1) *  holdc13_0[z - 1]
  beta13_0 = z1 * beta13_0_star + (1-z1) *  holdbeta13_0[z - 1]
  shape13_0 = z1 * shape13_0_star + (1-z1) *  holdshape13_0[z - 1]
  
  holdbeta13_0[z] = beta13_0
  holdscale13_0[z] = scale13_0
  holdc13_0[z] = c13_0
  holdshape13_0[z] = shape13_0
  
  scale23_0 = z15 * scale23_0_star + (1-z15) *  holdscale23_0[z - 1]
  theta23_0 = z15*theta23_0_star + (1-z15) *  holdtheta23_0[z - 1]
  c23_0 = z15 * c23_0_star + (1-z15) *  holdc23_0[z - 1]
  beta23_0 = z15 * beta23_0_star + (1-z15) *  holdbeta23_0[z - 1]
  shape23_0 = z15 * shape23_0_star + (1-z15) *  holdshape23_0[z - 1]
  
  holdbeta23_0[z] = beta23_0
  holdscale23_0[z] = scale23_0
  holdc23_0[z] = c23_0
  holdshape23_0[z] = shape23_0
  holdtheta23_0[z] = theta23_0

  omega12_star = rnorm(n/2, o12save1[1  (n/2)], sd = proposalsdfrail)
  omega13_star = rnorm(n/2, o13save1[1  (n/2)], sd = proposalsdfrail)

   if(z == 2){
     omega12_star = rnorm(n/2, o12save1[1  (n/2)], sd = frailtysd)
     omega13_star = rnorm(n/2, o13save1[1  (n/2)], sd = frailtysd)
   }
   
  omega12_z1 = o12save1[1  (n/2)]
  omega13_z1 = o13save1[1  (n/2)]
  omega23_z1 = o23save1[1  (n/2)]

  save = savestar = NULL
  z1 = p = rep(NA, n/2)

  for(i in 1:  (n/2)){
  
  sap[i] = like12_omega1_i(omega12_star[i], i = i, scale12_1 = holdscale12_1[z - 1], dat1 = dat1, shape12_1 = holdshape12_1[z - 1],
                             c12_1 = holdc12_1[z - 1])[1]
  sap2[i] = like12_omega1_i(omega12_z1[i], i = i, scale12_1 = holdscale12_1[z - 1], dat1 = dat1, shape12_1 = holdshape12_1[z - 1],
                              c12_1 = holdc12_1[z - 1])[1]
  
  p[i] = exp(sap[i] - sap2[i] + sum(log(pnorm(omega12_star, 0, sd = frailtysd)))-
                 sum(log(pnorm(omega12_z1, 0, sd = frailtysd)) ))
  
  if(is.na(p[i])) p[i] = 0; if((p[i] < 0)) p[i] = 0; if((p[i] > 1)) p[i] = 1
  z1[i] = rbinom(1, 1, prob = p[i])
  if(holdfrail12) z1[i] = 0
  }

  if(sum(z1) >= 1) accept6 = accept6 + 1

  omega12_star = omega12_z1 = omega_star = z1 * omega12_star + (1-z1) * omega12_z1
  holdfrailsd12_1[z] = sd(omega12_star)
  holdfrailmean12_1[z] = mean(omega12_star)
  
  omega_cf = rnorm(n/2, omega_star, sd = frailtysd * sqrt(frailtysd - rhos^2))
  
  o12save1 = c(omega_star, omega_cf)
  
  holdomega12_1[,z] = omega12_z1[1 : 10]
  if(holdc) c12_1_star = 1
  beta12_1_star = rnorm(1, holdbeta12_1[z - 1], sd = proposalsd)
  shape12_1_star = rtruncnorm(n = 1, a = 0, b = 3, mean = holdshape12_1[z - 1], sd = proposalsd)
  beta12_1_star = 0

  if(holdshape) shape12_1_star = shape12_1
  
    scale12_1_star = rgamma(1, shape = 0.01 + sum(dat1$s12, na.rm = T),
    scale = (0.01 + 1/shape12_1_star *sum(dat1$y12 ^ shape12_1_star * 
    exp(c12_1_star * omega12_star), na.rm = T))^(-1))
    
  if(holdscale12) scale12_1_star = scale12_1

  c12prior_p = pnorm(c12_1_star, 0, sd = 3)
  beta12prior_p = pnorm(beta12_1_star, 0, sd = 3)
  shape12prior_p = pnorm(shape12_1_star, 0, sd = 3)
  shape12prior_p = ptruncnorm(shape12_1_star, a = 0, b = 3, mean = 0, sd = 3)
  
  c12prior_c = pnorm(holdc12_1[z - 1], 0, sd = 3)
  beta12prior_c = pnorm(holdbeta12_1[z - 1], 0, sd = 3)
  shape12prior_c = pnorm(holdshape12_1[z - 1], 0, sd = 3)
  shape12prior_c = ptruncnorm(holdshape12_1[z - 1], a = 0, b = 3, mean = 0, sd = 3)
  
    p = exp(like12_shape1(shape12_1 = shape12_1_star, scale12_1 = scale12_1_star,
            dat1 = dat1, c12_1 = c12_1_star, omega12_z1 = omega12_z1) + log(shape12prior_p) -
            (like12_shape1(shape12_1 = holdshape12_1[z - 1], scale12_1 = scale12_1_star,
            dat1 = dat1, c12_1 = c12_1_star, omega12_z1 = omega12_z1)+ log(shape12prior_c)))
    
  if(is.nan(p)) p = 0; if(is.na(p)) p = 0; if((p < 0)) p = 0;  if((p > 1)) p = 1
  if(scale12_1_star < 0 ) p = 0
   z1 = rbinom(1, 1, prob = p)

  if(z1 == 1) accept7 = accept7 + 1

  c12_1 = c12_1_star
  shape12_1 = z1 * shape12_1_star + (1-z1) *  holdshape12_1[z - 1]
  beta12_1 = z1 * beta12_1_star + (1-z1) *  holdbeta12_1[z - 1]
  
  holdbeta12_1[z] = beta12_1

  if(holdscale12)   z1 = 1
  scale12_1 = z1 * scale12_1_star + (1-z1) *  holdscale12_1[z - 1]

  holdscale12_1[z] = scale12_1
  holdc12_1[z] = c12_1
  holdshape12_1[z] = shape12_1
  
  z1 = p = rep(NA, (n/2))
  sap = sap2 = rep(NA, (n/2))
  
  for(i in 1:  (n/2)){
    sap[i] =   like13_omega1_i(omega13_star[i], i = i, scale13_1 = holdscale13_1[z - 1], dat1 = dat1, shape13_1 = holdshape13_1[z - 1],
                               c13_1 = holdc13_1[z - 1], shape23_1 = holdshape23_1[z - 1], scale23_1 = holdscale23_1[z - 1], theta23_1 = holdtheta23_1[z - 1],
                               c23_1 = holdc23_1[z - 1])[1]
    sap2[i] =   like13_omega1_i(omega13_z1[i], i = i, scale13_1 = holdscale13_1[z - 1], dat1 = dat1, shape13_1 = holdshape13_1[z - 1],
                                c13_1 = holdc13_1[z - 1], shape23_1 = holdshape23_1[z - 1], scale23_1 = holdscale23_1[z - 1], theta23_1 = holdtheta23_1[z - 1],
                                c23_1 = holdc23_1[z - 1])[1]
    
    p[i]  = exp(sap[i] - sap2[i] + sum(log(pnorm(omega13_star, 0, sd = frailtysd)))-
                   sum(log(pnorm(omega13_z1, 0, sd = frailtysd)) ))
    
    if(is.na(p[i])) p[i] = 0; if((p[i] < 0)) p[i] = 0;  if((p[i] > 1)) p[i] = 1
    z1[i] = rbinom(1, 1, prob = p[i])
    if(holdfrail13) z1[i] = 0
    
  }
  
  omega23_z1 = omega13_star = omega13_z1 = omega_star = z1 * omega13_star + (1-z1) * omega13_z1
  omega_cf = rnorm(n/2, omega_star, sd = frailtysd * sqrt(frailtysd - rhos^2))
  
  o13save1 = c(omega_star, omega_cf)

  if(sum(z1[1]) >= 1) accept8 = accept8 + 1
  
   holdfrailsd13_1[z] = sd(omega13_star)
   holdfrailmean13_1[z] = mean(omega13_star)
   
   z1 = p = rep(NA, (n/2))
   sap = sap2 = rep(NA, (n/2))
 
   o23save1 = o13save1
   
   holdomega23_0[,z] = omega23_z1[1 : 10]
   holdfrailsd23_1[z] = sd(omega23_star)
   holdfrailmean23_1[z] = mean(omega23_star)
   
   holdomega23_1[,z] = omega23_z1[1 : 10]
   
  holdomega13_1[,z] = omega13_z1[1 : 10]
  c13_1_star = rnorm(1, holdc13_1[z - 1], sd = proposalsd)

  if(holdc) c13_1_star = 1
  beta13_1_star = rnorm(1, holdbeta13_1[z - 1], sd = proposalsd)
  beta13_1_star = 0
  shape13_1_star = rtruncnorm(1, a = 0, b = 3, mean = holdshape13_1[z - 1], sd = proposalsd)
  if(holdshape) shape13_1_star = 1

  scale13_1_star = rgamma(1, shape = 0.01 + sum(dat1$s13, na.rm = T),
                          scale = (0.01 + 1/shape13_1_star *
                          sum(dat1$y13 ^ shape13_1_star * 
                          exp(c13_1_star * omega13_z1), na.rm = T))^(-1))
  
  if(holdscale13) scale13_1_star = scale13_1

  beta23_1_star = 0
  shape23_1_star = rnorm(1, holdshape23_1[z - 1], sd = proposalsd)
  
  shape23_1_star = rtruncnorm(1, a = 0, b = 3, mean = holdshape23_1[z - 1], sd = proposalsd)
  theta23_1_star = rnorm(1, -mod23_1$coefficients[2], sd = proposalsd)
  theta23_1_star = rnorm(1, holdtheta23_1[z - 1], sd = proposalsd)
  
  if(holdtheta){ theta23_1_star = theta23_1 }

  c23_1_star = rnorm(1, holdc23_1[z - 1], sd = proposalsd)
  c23_1_star = rtruncnorm(1, a=-3, b= 3, mean = holdc23_1[z - 1], sd = proposalsd)
  
  if(holdshape) shape23_1_star = 1
  if(holdc23) c23_1_star = 1
  
  scale23_1_star = rgamma(1, shape = 0.01 + sum(dat1$s23, na.rm = T),
                               scale = (0.01 + 1/shape23_1_star *
                               sum(dat1$y23 ^ shape23_1_star * 
                               exp(c23_1_star * omega23_z1 +
                               theta23_1_star * dat1$y12), na.rm = T))^(-1))

  if(holdscale23) scale23_1_star =  scale23_1

  c13prior_p = pnorm(c13_1_star, 0, sd = 3)
  beta13prior_p = pnorm(beta13_1_star, 0, sd = 3)
  shape13prior_p = pnorm(shape13_1_star, 0, sd = 3)
  shape13prior_p = ptruncnorm(shape13_1_star, a = 0, b = 3, mean = 0, sd = 3)
  
  c13prior_c = pnorm(holdc13_1[z - 1], 0, sd = 3)
  beta13prior_c = pnorm(holdbeta13_1[z - 1], 0, sd = 3)
  shape13prior_c = pnorm(holdshape13_1[z - 1], 0, sd = 3)
  shape13prior_c = ptruncnorm(holdshape13_1[z - 1], a = 0, b = 3, mean = 0, sd = 3)
  
  c23prior_p = pnorm(c23_1_star, 0, sd = 3)
  beta23prior_p = pnorm(beta23_1_star, 0, sd = 3)
  theta23prior_p = pnorm(theta23_1_star, 0, sd = 3)
  shape23prior_p = pnorm(shape23_1_star, 0, sd = 3)
  shape23prior_p = ptruncnorm(shape23_1_star, a = 0, b = 3, mean = 0, sd = 3)

  c23prior_c = pnorm(holdc23_1[z - 1], 0, sd = 3)
  beta23prior_c = pnorm(holdbeta23_1[z - 1], 0, sd = 3)
  shape23prior_c = pnorm(holdshape23_1[z - 1], 0, sd = 3)
  shape23prior_c = ptruncnorm(holdshape23_1[z - 1], a = 0, b = 3, mean = 0, sd = 3)
  theta23prior_c = pnorm(holdtheta23_1[z - 1], 0, sd = 3)
  
  theta23prop_p = pnorm(theta23_1_star, -mod23_1$coefficients[2], sd = proposalsd)
  theta23prop_c = pnorm(holdtheta23_1[z - 1], -mod23_1$coefficients[2], sd = proposalsd)
  
  theta23prop_p = pnorm(theta23_1_star, 0, sd = proposalsd)
  theta23prop_c = pnorm(holdtheta23_1[z - 1], 0, sd = proposalsd)

   p = exp(like13_shape1(shape13_1 = shape13_1_star, scale13_1 = holdscale13_1[z - 1], dat1 = dat1,
              omega13_z1 = omega13_z1, c13_1 = c13_1_star) + log(shape13prior_p)-
              (like13_shape1(shape13_1 = holdshape13_1[z - 1], scale13_1 = holdscale13_1[z - 1], dat1 = dat1,
              omega13_z1 = omega13_z1, c13_1 = c13_1_star)+ log(shape13prior_c)))

   p2 = exp( like23_theta1(theta23_1 = theta23_1_star, dat1 = dat1,c23_1 = holdc23_1[z - 1], omega13_z1 = omega13_z1,
                  shape23_1 = holdshape23_1[z - 1], scale23_1 = holdscale23_1[z - 1])+
                  log(theta23prior_p) + log(theta23prop_p) -
                  (like23_theta1(theta23_1 = holdtheta23_1[z - 1],dat1 = dat1, omega13_z1 = omega13_z1, c23_1 = holdc23_1[z - 1],
                  shape23_1 = holdshape23_1[z - 1], scale23_1 = holdscale23_1[z - 1])+ log(theta23prior_c) + log(theta23prop_c))+
                  like23_shape1(shape23_1 = shape23_1_star, dat1 = dat1, omega13_z1 = omega13_z1, c23_1 = holdc23_1[z - 1],
                  scale23_1 = holdscale23_1[z - 1], theta23_1 = holdtheta23_1[z - 1]) + log(shape23prior_p)-
                  (like23_shape1(shape23_1 = holdshape23_1[z - 1], dat1 = dat1, c23_1 = holdc23_1[z - 1], omega13_z1 = omega13_z1,
                  scale23_1 = holdscale23_1[z - 1], theta23_1 = holdtheta23_1[z - 1]) + log(shape23prior_c)))

   if(is.nan(p)) p = 0; if(is.na(p)) p = 0; if((p < 0)) p = 0;  if((p > 1)) p = 1
   if(is.nan(p2)) p2 = 0; if(is.na(p2)) p2 = 0; if((p2<0)) p2 = 0;  if((p2>1)) p2 = 1

   if(scale23_1_star < 0) {p2 = 0}
   if(scale13_1_star < 0) {p = 0}
   z1 = rbinom(1, 1, prob = p)
   z15 = rbinom(1, 1, prob = p2)
    
  if(z1 == 1) accept9 = accept9 + 1
  if(z15 == 1) accept10 = accept10 + 1
  
  scale13_1 = z1 * scale13_1_star + (1-z1) *  holdscale13_1[z - 1]
  scale13_1 = scale13_1_star
  
  c13_1 = c13_1_star
  beta13_1 = z1 * beta13_1_star + (1-z1) *  holdbeta13_1[z - 1]
  shape13_1 = z1 * shape13_1_star + (1-z1) *  holdshape13_1[z - 1]
  
  holdbeta13_1[z] = beta13_1 
  holdscale13_1[z] = scale13_1
  holdc13_1[z] = c13_1
  holdshape13_1[z] = shape13_1

  holdscale23_1[z] = scale23_1 = z15 * scale23_1_star + (1-z15) *  holdscale23_1[z - 1]
  holdscale23_1[z] = scale23_1 = scale23_1_star
  
  holdtheta23_1[z] = z15*theta23_1_star + (1-z15) *  holdtheta23_1[z - 1]
  holdc23_1[z] = c23_1 = z15 * c23_1_star + (1-z15) *  holdc23_1[z - 1]
  holdshape23_1[z] =  shape23_1 = z15 * shape23_1_star + (1-z15) *  holdshape23_1[z - 1]
  holdbeta23_1[z] = beta23_1 = z15 * beta23_1_star + (1-z15) *  holdbeta23_1[z - 1]
  
if(cep){ 
  
  o12save1flip = c(o12save1[(n/2 + 1) : n], o12save1[1 : (n/2)])
  o13save1flip = c(o13save1[(n/2 + 1) : n], o13save1[1 : (n/2)])

  o12save0flip = o12save0
  o13save0flip = o13save0
  
  Fw_0 = Fw = Fw_1 <- NULL
  j = 1
  xtrue = c(x_0, x_1)
  
  intfunction = function(j, i, t){
    exp(
      -cLambda13_frailty_lk(x = t, xdata1 = xtrue[i], omega2 = o13save0flip[i], scale13 = holdscale13_0[z], 
      shape13 = holdshape13_0[z],   c13 = holdc13_0[z], beta13_1 = holdbeta13_0[z]) - 
       cLambda12_frailty_lk(x = t, xdata1 = xtrue[i], 
      omega1 = o12save0flip[i],   scale12 = holdscale12_0[z], shape12 = holdshape12_0[z], c12 = holdc12_0[z], beta12_1 = holdbeta12_0[z])) *
      lambda12_frailty(t, xdata1 = xtrue[i], omega1 = o12save0[i], scale12 = holdscale12_0[z], shape12 = holdshape12_0[z],
      c12 = holdc12_0[z], beta12_1 = holdbeta12_0[z]) * 
      exp(-cLambda23_frailty_lk(x = tau_t - t, xdata1 = xtrue[i], omega2 = o13save0flip[i], scale23 = holdscale23_0[z], 
      shape23 = holdshape23_0[z],   c23 = holdc23_0[z], theta23 = holdtheta23_0[z], v_predict = t, beta23_1 = holdbeta23_0[z]
   ))
  }
  for(i in 1:  n){
    y =tryCatch({ integrate(f = intfunction, lower = 0, upper = tau_t, i = i, j = j)$val
    },
    error = function(error_condition) {
      NA}
    )  
    
    L12 = cLambda12_frailty_lk(x = tau_t, xdata1 = xtrue[i], shape12 = holdshape12_0[z],
                              omega1 = o12save0flip[i],   scale12 =holdscale12_0[z], c12 = holdc12_0[z], beta12_1 = holdbeta12_0[z])
    L13 = cLambda13_frailty_lk(x = tau_t, xdata1 = xtrue[i], omega2 = o13save0flip[i], scale13 = holdscale13_0[z], shape13 = holdshape13_0[z],
                              c13 = holdc13_0[z], beta13_1 = holdbeta13_0[z])
    
    Fw =  (exp(-L13 - L12) + y)

    Fw_0 = c(Fw_0, Fw)
  }
  
  intfunction = function(j, i, t){
    exp(
      -cLambda13_frailty_lk(x = t, xdata1 = xtrue[i], omega2 = o13save1flip[i], scale13 = holdscale13_1[z], 
        shape13 = holdshape13_1[z],   c13 = holdc13_1[z], beta13_1 = holdbeta13_1[z]) - 
        cLambda12_frailty_lk(x = t, xdata1 = xtrue[i], 
        omega1 = o12save1flip[i],   scale12 = holdscale12_1[z], shape12 = holdshape12_1[z], c12 = holdc12_1[z], beta12_1 = holdbeta12_1[z])) *
      lambda12_frailty(t, xdata1 = xtrue[i], omega1 = o12save1flip[i], scale12 = holdscale12_1[z], shape12 = holdshape12_1[z],
        c12 = holdc12_1[z], beta12_1 = holdbeta12_1[z])* 
      exp(-cLambda23_frailty_lk(x = tau_t - t, xdata1 = xtrue[i], omega2 = o13save1flip[i], scale23 = holdscale23_1[z], 
        shape23 = holdshape23_1[z], c23 = holdc23_1[z], theta23 = holdtheta23_1[z], v_predict = t, beta23_1 = holdbeta23_1[z]))
  }
  ## for z = 1
  for(i in 1:  n){
    y =tryCatch({ integrate(f = intfunction, lower = 0, upper = tau_t, i = i, j = j)$val
    },
    error = function(error_condition) {
      NA}
   )  
    L12 =cLambda12_frailty_lk(x = tau_t, xdata1 = xtrue[i], shape12 = holdshape12_1[z],
                             omega1 = o12save1flip[i],   scale12 =holdscale12_1[z], c12 = holdc12_1[z], beta12_1 = holdbeta12_1[z])
    L13 = cLambda13_frailty_lk(x = tau_t, xdata1 = xtrue[i], omega2 = o13save1flip[i], scale13 = holdscale13_1[z], 
                              shape13 = holdshape13_1[z],   c13 = holdc13_1[z], beta13_1 = holdbeta13_1[z])
    
    Fw =  (exp(-L13 - L12) + y)
    
    Fw_1 = c(Fw_1, Fw) 
  }
  
  s0cumulative = -pweibull(q = tau_s, scale=1/holdscale12_0[z], shape= holdshape12_0[z], lower = F, log = T) * exp(holdc12_0[z]*o12save0flip)
  s1cumulative = -pweibull(q = tau_s, scale=1/holdscale12_1[z], shape= holdshape12_1[z], lower = F, log = T) * exp(holdc12_1[z]*o12save1flip)
  
  dat = data.frame(cbind(s0cumulative/s1cumulative, Fw_1 - Fw_0))
  dat$X = log(s0cumulative/s1cumulative)
  dat$Y = c(Fw_1 - Fw_0) 

  reg<-lm(formula = Y ~ X,  data=dat)
  holdslope[z] = (summary(reg)$coef[2,1])
  holdint[z] = (summary(reg)$coef[1,1])
  
  saveCEPx[, z] = dat$X
  saveCEPy[, z] = dat$Y
}
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
                    theta23_0 = holdtheta23_0, theta23_1 = holdtheta23_1,
                    accept1 = accept1, accept2 = accept2, 
                    accept3 = accept3, accept4 = accept4, 
                    accept5 = accept5, accept6 = accept6, 
                    accept7 = accept7, accept8 = accept8, 
                    accept9 = accept9, accept10 = accept10, 
                    holdfrailsd12_0 = holdfrailsd12_0,
                    holdfrailsd13_0 = holdfrailsd13_0,
                    holdfrailsd23_0 = holdfrailsd23_0,
                    holdfrailsd12_1 = holdfrailsd12_1,
                    holdfrailsd13_1 = holdfrailsd13_1,
                    holdfrailsd23_1 = holdfrailsd23_1
                    
)

colMeans(params, na.rm = T)

result = list(params = params, saveCEPx = saveCEPx, saveCEPy = saveCEPy, 
              accept.rate1 = accept1/(SIM), 
              accept.rate2 = accept2/(SIM), 
              accept.rate3 = accept3/(SIM), 
              accept.rate4 = accept4/(SIM), 
              accept.rate5 = accept5/(SIM), 
              accept.rate6 = accept6/(SIM), 
              accept.rate7 = accept7/(SIM), 
              accept.rate8 = accept8/(SIM), 
              accept.rate9 = accept9/(SIM), 
              accept.rate10 = accept10/(SIM), 
              
              args = list(SIM = SIM, burnin = burnin, n = n))

return(result)

}