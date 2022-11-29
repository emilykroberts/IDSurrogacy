library(dplyr)
library(survival)
require(stats)
library(splines2)
library(pracma)  ## for numerical differentiation

## For estimating B-spline basis coefficients
obj.fn.1 <- function(eta.vec, spline.pred, h0.truth){ 
  ##
  spline.vals <- sapply(1:nrow(spline.pred), function(i) sum(spline.pred[i,]*(eta.vec)))
  sq.diffs <- (log(h0.truth)-spline.vals)^2
  return(sum(sq.diffs))
}

## Simulate data from illness-death, shared frailty
## Adapted from simID in SemiCompRisks package
## lt.type: draw left-truncation times from "unif" (uniform) or "norm" (normal)
##          - if "unif", lt=c(a, b) so that L ~ unif(n, a, b)
##          - if "norm", lt=c(mean, sd), so that L ~ Norm(mean, sd)
simID.LT <- function(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                     alpha1.true, alpha2.true, alpha3.true, frailtysd = frailtysd,
                     kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens, lt.type = "unif", lt)
{
  
  x1 = x2 = x3 = 0
  theta.true = 1
  n <- dim(x1)[1]
  p1 <- dim(x1)[2]
  p2 <- dim(x2)[2]
  p3 <- dim(x3)[2]
  alpha3.true = alpha2.true = alpha1.true = 1
  kappa3.true = kappa1.true = kappa2.true = 1
  if(theta.true >0)
  {
   # gamma.true <- rgamma(n, 1/theta.true, 1/theta.true)
    gamma.true <- exp(rnorm(n, mean = 0, sd = frailtysd))
  }
  if(theta.true == 0)
  {
    gamma.true <- rep(1, n)
  }
   # LP1	<- as.vector(beta1.true %*% t(x1))
#  LP2	<- as.vector(beta2.true %*% t(x2))
#  LP3	<- as.vector(beta3.true %*% t(x3))
  
  LP1= LP2 = LP3 = rep(0, n)
  
  Rind <- NULL
  R <- rweibull(n, shape = alpha1.true, scale = exp(-(log(kappa1.true) +
                                                        LP1 + log(gamma.true))/alpha1.true))
  D <- rweibull(n, shape = alpha2.true, scale = exp(-(log(kappa2.true) +
                                                        LP2 + log(gamma.true))/alpha2.true))
  
  yesR <- R < D
  
  D[yesR] <- R[yesR] + rweibull(sum(yesR), shape = alpha3.true,
                                scale = exp(-(log(kappa3.true) + LP3[yesR] + log(gamma.true[yesR]))/alpha3.true))
  delta1 <- rep(NA, n)
  delta2 <- rep(NA, n)
  y1 <- R
  y2 <- D
  #Cen <- runif(n, cens[1], cens[2])
  Cen <- rep(1, n)
  ind01 <- which(D < R & D < Cen)
  y1[ind01] <- D[ind01]
  delta1[ind01] <- 0
  delta2[ind01] <- 1
  ind10 <- which(R < D & R < Cen & D >= Cen)
  y2[ind10] <- Cen[ind10]
  delta1[ind10] <- 1
  delta2[ind10] <- 0
  ind00 <- which(R >= Cen & D >= Cen)
  y1[ind00] <- Cen[ind00]
  y2[ind00] <- Cen[ind00]
  delta1[ind00] <- 0
  delta2[ind00] <- 0
  ind11 <- which(R < Cen & D < Cen & R < D)
  delta1[ind11] <- 1
  delta2[ind11] <- 1
  
  lt.type = "norm"
  lt = rep(1, n)
  if (lt.type == "unif") L <- runif(n, lt[1], lt[2])
  if (lt.type == "norm") L <- rnorm(n, lt[1], lt[2])
  ret <- data.frame(cbind(y1, delta1, y2, delta2, L))
  return(ret)
  
  print(ret)
}


#################
## b-spline BH ##
#################

## Log-likelihood function for illness-death with shared frailty applied to left-truncated data
## bspline baseline hazard !!
## b.1: bspline basis function for the 1-transition (corresponds to y1)
## b.2: bspline basis function for the 2-transition (corresponds to y2)
## b.3.y2my1: bsspline basis function for the 3-transition (corresponds to y2-y1)
## log h0,1(t) = phi.1,0 * B_1,0(t) + ... + phi.1,k1 * B_1,k1(t), where B_1,l(t) are the bspline basis functions
## log h0,2(t) = phi.2,0 * B_2,0(t) + ... + phi.2,k1 * B_2,k1(t)
## log h0,3(t) = phi.3,0 * B_3,0(t) + ... + phi.3,k1 * B_3,k1(t)
## para: vector of true
##       - c(phi.1.vector, phi.2.vector, phi.3.vector, log(theta), beta1, beta2, beta3)
logLike.SCR.SM.LT.bSpline123.dropPrevCases <- function(para, y1, y2, delta1, delta2, l, Xmat1=NULL, Xmat2=NULL, Xmat3=NULL, frailty=TRUE,
                                                       b.1,  
                                                       b.2,  
                                                       b.3.y2my1)
{
  ##
  if (is.vector(Xmat1)==T) Xmat1 = matrix(Xmat1, nrow=1)
  if (is.vector(Xmat2)==T) Xmat2 = matrix(Xmat2, nrow=1)
  if (is.vector(Xmat3)==T) Xmat3 = matrix(Xmat3, nrow=1)
  num.Bspline.params.1 <- ncol(b.1)
  num.Bspline.params.2 <- ncol(b.2)
  num.Bspline.params.3 <- ncol(b.3.y2my1)
  phi1 <- para[1:(1+num.Bspline.params.1-1)]
  phi2 <- para[(1+num.Bspline.params.1):(1+num.Bspline.params.1 + num.Bspline.params.2 - 1)]
  phi3 <- para[(1+num.Bspline.params.1 + num.Bspline.params.2):(1+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1)]
  
  if(frailty == TRUE){
    theta    <- exp(para[(2+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1)])
    thetaInv <- 1 / theta
  }
  ##
  nP.0 <- ifelse(frailty, (2+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1), (1+num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 - 1))
  nP.1 <- ncol(Xmat1)
  nP.2 <- ncol(Xmat2)
  nP.3 <- ncol(Xmat3)
  ##
  eta.1 <- as.vector(Xmat1 %*% para[nP.0 + c(1:nP.1)])
  eta.2 <- as.vector(Xmat2 %*% para[nP.0 + nP.1 + c(1:nP.2)])
  eta.3 <- as.vector(Xmat3 %*% para[nP.0 + nP.1 + nP.2 + c(1:nP.3)])
  ##
  type1 <- as.numeric(delta1 == 1 & delta2 == 1 & l < y1)
  type2 <- as.numeric(delta1 == 0 & delta2 == 1 & l < y1)
  type3 <- as.numeric(delta1 == 1 & delta2 == 0 & l < y1)
  type4 <- as.numeric(delta1 == 0 & delta2 == 0 & l < y1)
  ##
  ## log.h1star.y1 <- log(alpha1) + log(kappa1) + (alpha1 - 1) * log(y1) + eta.1
  ## Uses h1 - need to bSpline
  ## log.h1,0(t) = phi0 * B0(t) + ... + phik * Bk(t)
  ##############################################################################
  B0.1.y1 <- as.vector(matrix(phi1, nrow=1) %*% t(as.matrix(b.1)))
  log.h1star.y1 <- B0.1.y1 + eta.1
  
  ## log.h2star.y1 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y1) + eta.2
  ##############################################################################
  B0.2.y1 <- as.vector(matrix(phi2, nrow=1) %*% t(as.matrix(predict(b.2, y1))))
  log.h2star.y1 <- B0.2.y1 + eta.2
  
  ## log.h2star.y2 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y2) + eta.2
  ##############################################################################
  B0.2.y2 <- as.vector(matrix(phi2, nrow=1) %*% t(as.matrix(b.2)))
  log.h2star.y2 <- B0.2.y2 + eta.2
  
  ## log.h3star.y2 <- log(alpha3) + log(kappa3) + (alpha3 - 1) * log(y2-y1) + eta.3
  ##################################################################################
  B0.3.y2my1 <- as.vector(matrix(phi3, nrow=1) %*% t(as.matrix(b.3.y2my1)))
  log.h3star.y2 <- B0.3.y2my1 + eta.3
  
  ## q.y1 <- kappa1*(y1)^alpha1 * exp(eta.1) + kappa2*(y1)^alpha2 * exp(eta.2)
  ## Uses h1 - needs bSpline
  ## Essentially uses H0.1.y1 and H0.2.y1
  ## exp(phi.3.truth %*% t(b.3.event))
  ############################################################################
  h0.1.y1 <- exp(B0.1.y1)
  h0.1.y1.interpolate.func <- approxfun(y1, h0.1.y1, rule=2)
  H0.1.y1 <- sapply(y1, function(x) integrate(h0.1.y1.interpolate.func, 0, x, stop.on.error = F)$value)
  
  h0.2.y1 <- exp(B0.2.y1)
  h0.2.y1.interpolate.func <- approxfun(y1, h0.2.y1, rule=2)
  H0.2.y1 <- sapply(y1, function(x) integrate(h0.2.y1.interpolate.func, 0, x, stop.on.error = F)$value)
  q.y1 <- H0.1.y1 * exp(eta.1) + H0.2.y1 * exp(eta.2)
  
  ## q.y2 <- kappa1*(y2)^alpha1 * exp(eta.1) + kappa2*(y2)^alpha2 * exp(eta.2)
  ## Essentially uses H0.1.y2 and H0.2.y2
  ############################################################################
  B0.1.y2 <- as.vector(matrix(phi1, nrow=1) %*% t(as.matrix(predict(b.1, y2))))
  h0.1.y2 <- exp(B0.1.y2)
  h0.1.y2.interpolate.func <- approxfun(y2, h0.1.y2, rule=2)
  H0.1.y2 <- sapply(y2, function(x) integrate(h0.1.y2.interpolate.func, 0, x, stop.on.error = F)$value)
  
  h0.2.y2 <- exp(B0.2.y2)
  h0.2.y2.interpolate.func <- approxfun(y2, h0.2.y2, rule=2)
  H0.2.y2 <- sapply(y2, function(x) integrate(h0.2.y2.interpolate.func, 0, x, stop.on.error = F)$value)
  q.y2 <- H0.1.y2 * exp(eta.1) + H0.2.y2 * exp(eta.2)
  
  ## q.l <- kappa1*(l)^alpha1 * exp(eta.1) + kappa2*(l)^alpha2 * exp(eta.2)
  ## Essentially uses H0.1.l and H0.2.l
  ############################################################################
  B0.1.l <- as.vector(matrix(phi1, nrow=1) %*% t(as.matrix(predict(b.1, l))))
  h0.1.l <- exp(B0.1.l)
  h0.1.l.interpolate.func <- approxfun(l, h0.1.l, rule=2)
  H0.1.l <- sapply(l, function(x) integrate(h0.1.l.interpolate.func, 0, x, stop.on.error = F)$value)
  
  B0.2.l <- as.vector(matrix(phi2, nrow=1) %*% t(as.matrix(predict(b.2, l))))
  h0.2.l <- exp(B0.2.l)
  h0.2.l.interpolate.func <- approxfun(l, h0.2.l, rule=2)
  H0.2.l <- sapply(l, function(x) integrate(h0.2.l.interpolate.func, 0, x, stop.on.error = F)$value)
  q.l <- H0.1.l * exp(eta.1) + H0.2.l * exp(eta.2)
  
  ## w.y1.y2 <- kappa3*(y2-y1)^alpha3 * exp(eta.3)
  ## This is H0.3(y2-y1)
  ############################################################################
  h0.3.y2my1 <- exp(B0.3.y2my1)
  h0.3.y2my1.interpolate.func <- approxfun(y2-y1, h0.3.y2my1, rule=2)
  H0.3.y2my1 <- sapply(y2-y1, function(x) integrate(h0.3.y2my1.interpolate.func, 0, x, stop.on.error = F)$value)
  w.y1.y2 <- H0.3.y2my1 * exp(eta.3)
  ##
  k1 <- w.y1.y2
  k2.y1 <- q.y1 - q.l
  k2.y2 <- q.y2 - q.l
  ##
  if(frailty == TRUE)
  {
    logLike1 <- log.h1star.y1 + log.h3star.y2 + log(1+theta) - ((thetaInv + 2) * log(1 + (theta * (k1 + k2.y1))))
    logLike2 <- log.h2star.y1 - ((thetaInv + 1) * log(1 + (theta * k2.y1)))  ## Making in terms of y1
    logLike3 <- log.h1star.y1 - ((thetaInv + 1) * log(1 + (theta * (k1 + k2.y1))))
    logLike4 <- - thetaInv * log(1 + (theta * k2.y1))  ## Making in terms of y1
  }
  if(frailty == FALSE)
  {
    logLike1 <- log.h1star.y1 + log.h3star.y2 - (k1 + k2.y1)
    logLike2 <- log.h2star.y1 - k2.y1  ## Making in terms of y1
    logLike3 <- log.h1star.y1 - (k1 + k2.y1)
    logLike4 <- - k2.y1  ## Making in terms of y1
  }
  ##
  loglh <- sum(logLike1[type1==1]) + sum(logLike2[type2==1]) + sum(logLike3[type3==1]) + sum(logLike4[type4==1]) 
  return(-loglh)
}


## Fitting illness-death, shared frailty, on left-truncated data
## bspline baseline hazard functions
## Note: startVals must be specified (see get.Bspline.startVals)
## Semi-Markov model
FreqID.LT.bSpline123.v3 <- function(Y, lin.pred, data, startVals, frailty=TRUE, 
                                    b.1, b.2, b.3.y2my1, 
                                    bdy.knots.b.1, 
                                    bdy.knots.b.2, 
                                    bdy.knots.b.3.y2my1,
                                    method = "optim")
{	
  
  ##
  y1     <- as.vector(Y[,1])
  delta1 <- as.vector(Y[,2])
  y2     <- as.vector(Y[,3])
  delta2 <- as.vector(Y[,4])
  l      <- as.vector(Y[,5])
  Xmat1  <- as.matrix(model.frame(lin.pred[[1]], data=data))  
  Xmat2  <- as.matrix(model.frame(lin.pred[[2]], data=data))  
  Xmat3  <- as.matrix(model.frame(lin.pred[[3]], data=data))
  ##
  cat("This is going to take quite awhile...~1 hour for 4,000 observations")
  if (method == "optim"){
    logLike <- function(p) logLike.SCR.SM.LT.bSpline123.dropPrevCases(p, y1=y1, y2=y2, delta1=delta1, delta2=delta2, l=l,
                                                                      Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3, frailty=frailty,
                                                                      b.1 = b.1,  
                                                                      b.2 = b.2, 
                                                                      b.3.y2my1 = b.3.y2my1)
    
    optim.control = list(REPORT = 50, maxit = 2000)  
    
    fit1 <- optim(startVals,   
                  logLike, hessian = TRUE, method="Nelder-Mead", control = optim.control)
    H = fit1$hessian
    est = fit1$par
    logLike=-fit1$value
    code = fit1$convergence
    
    if (method == "nlm"){    
      fit2 <- nlm(logLike.SCR.SM.LT.bSpline123.dropPrevCases, p=startVals, y1=y1, y2=y2,
                  delta1=delta1, delta2=delta2, l=l,
                  Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3, frailty=frailty,
                  b.1 = b.1,
                  b.2 = b.2,
                  b.3.y2my1 = b.3.y2my1, hessian=T)
      H = fit2$hessian
      logLike = -fit2$minimum
      est = fit2$estimate
      code = fit2$code
    }
  }
  ##
  myLabels <- c("phi1", "phi2", "phi3")
  if(frailty == TRUE) myLabels <- c(myLabels, "log(theta)")
  myLabels <- c(myLabels, colnames(Xmat1), colnames(Xmat2), colnames(Xmat3))
  nP <- c(ncol(Xmat1), ncol(Xmat2), ncol(Xmat3))
  
  value <- list(estimate=est, H=H, logLike=logLike, myLabels=myLabels, frailty=frailty, nP=nP, 
                code = code)
  
  class(value) <- c("Freq", "ID", "Ind", "WB", "semi-Markov")
  
  return(value)
  ##
  invisible()
}


## Getting bspline starting values for each baseline hazard; estimates bspline coefs and betas
## It's convoluted, but it works okay.
get.Bspline.startVals <- function(y, delta, lin.pred, data, b, knot.loc, Bspline.degree, method="nls"){
  ## 
  num.Bspline.params <- length(knot.loc) + Bspline.degree + 1
  
  fitCox <- coxph(as.formula(paste("Surv(y, delta)", paste(lin.pred, collapse = ""))), data=data)
  H0.info <- basehaz(fitCox, centered = F) ## Get cumulative hazard from Cox
  H0.est.linInterpol <- approxfun(H0.info$time, H0.info$hazard, rule = 2) ## Get lin interpolated function
  yvals <- seq(min(y), max(y), length=500) ## get plenty of y-values
  h0.est.func2 <- numdiff(H0.est.linInterpol, yvals) ## estimate h0 using numerical differentiation
  h0.est.smooth.xy <- loess.smooth(yvals, h0.est.func2, evaluation=500) ## Smooth out using loess
  
  ## Now using least squares to get estimated phi coefficients based on loess-smoothed h0
  h0.est.smooth <- h0.est.smooth.xy$y
  yvals.bspline <- data.frame(predict(b, yvals))
  
  obj.fn.1 <- function(eta.vec, spline.pred, h0.truth){ 
    ##
    spline.vals <- sapply(1:nrow(spline.pred), function(i) sum(spline.pred[i,]*(eta.vec)))
    sq.diffs <- (log(h0.truth)-spline.vals)^2
    return(sum(sq.diffs))
  }
  if (method == "nlm"){
    phi.h0.est.smooth=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params), 
                          spline.pred=yvals.bspline, 
                          h0.truth=h0.est.smooth)$estimate
    return(c(phi.h0.est.smooth, fitCox$coef))
  }
  if (method == "nls"){
    ## Alternatively -- this seems to give better results? 
    g <- function(bSpline.mat, phi.vec){
      phi.list <- as.list(phi.vec)
      mat.list <- as.list(data.frame(bSpline.mat))
      return(apply(do.call(rbind, Map('*', phi.list, mat.list)), 2, sum))
    }
    phi.h0.est.smooth.ls <- nls(h0.est.smooth ~ exp(g(yvals.bspline, phi.vec)), 
                                start = list(phi.vec = rep(-1, ncol(yvals.bspline))))
    phi.h0.est.smooth <- summary(phi.h0.est.smooth.ls)$parameters[,1]  ## These will be good for starting values!
    return(c(phi.h0.est.smooth, fitCox$coef))
  }
}


################
## Weibull BH ##
################
## Log-likelihood illness-death, shared frailty, left-truncated data
## Weibull baseline hazards
logLike.weibull.SCR.SM.LT <- function(para, y1, y2, delta1, delta2, l, Xmat1=NULL, Xmat2=NULL, Xmat3=NULL, frailty=TRUE)
{
  ##
  kappa1    <- exp(para[1])
  alpha1 <- exp(para[2])
  kappa2    <- exp(para[3])
  alpha2 <- exp(para[4])
  kappa3    <- exp(para[5])
  alpha3 <- exp(para[6])
  if(frailty == TRUE){
    theta    <- exp(para[7])
    thetaInv <- 1 / theta
  }
  ##
  nP.0 <- ifelse(frailty, 7, 6)
  nP.1 <- ncol(Xmat1)
  nP.2 <- ncol(Xmat2)
  nP.3 <- ncol(Xmat3)
  ##
  eta.1 <- as.vector(Xmat1 %*% para[nP.0 + c(1:nP.1)])
  eta.2 <- as.vector(Xmat2 %*% para[nP.0 + nP.1 + c(1:nP.2)])
  eta.3 <- as.vector(Xmat3 %*% para[nP.0 + nP.1 + nP.2 + c(1:nP.3)])
  ##
  type1 <- as.numeric(delta1 == 1 & delta2 == 1 & l < y1)
  type2 <- as.numeric(delta1 == 0 & delta2 == 1 & l < y1)
  type3 <- as.numeric(delta1 == 1 & delta2 == 0 & l < y1)
  type4 <- as.numeric(delta1 == 0 & delta2 == 0 & l < y1)
  type5 <- as.numeric(delta1 == 1 & delta2 == 1 & y1 <= l & l < y2)
  type6 <- as.numeric(delta1 == 1 & delta2 == 0 & y1 <= l & l < y2)
  ##
  log.h1star.y1 <- log(alpha1) + log(kappa1) + (alpha1 - 1) * log(y1) + eta.1
  log.h2star.y1 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y1) + eta.2
  log.h2star.y2 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y2) + eta.2
  log.h3star.y2 <- log(alpha3) + log(kappa3) + (alpha3 - 1) * log(y2-y1) + eta.3
  ##
  q.y1 <- kappa1*(y1)^alpha1 * exp(eta.1) + kappa2*(y1)^alpha2 * exp(eta.2)
  q.y2 <- kappa1*(y2)^alpha1 * exp(eta.1) + kappa2*(y2)^alpha2 * exp(eta.2)
  q.l <- kappa1*(l)^alpha1 * exp(eta.1) + kappa2*(l)^alpha2 * exp(eta.2)
  ##
  w.y1.y2 <- kappa3*(y2-y1)^alpha3 * exp(eta.3)
  w.y1.l  <- kappa3*(l-y1)^alpha3 * exp(eta.3)
  ##
  k1 <- w.y1.y2
  k2.y1 <- q.y1 - q.l
  k2.y2 <- q.y2 - q.l
  k3 <- w.y1.y2 - w.y1.l
  ##
  if(frailty == TRUE)
  {
    logLike1 <- log.h1star.y1 + log.h3star.y2 + log(1+theta) - ((thetaInv + 2) * log(1 + (theta * (k1 + k2.y1))))
    logLike2 <- log.h2star.y1 - ((thetaInv + 1) * log(1 + (theta * k2.y1)))  ## Making in terms of y1
    logLike3 <- log.h1star.y1 - ((thetaInv + 1) * log(1 + (theta * (k1 + k2.y1))))
    logLike4 <- - thetaInv * log(1 + (theta * k2.y1))  ## Making in terms of y1
    logLike5 <- log.h3star.y2 - ((thetaInv + 1) * log(1 + (theta * k3)))
    logLike6 <- - thetaInv * log(1 + (theta * k3))
  }
  if(frailty == FALSE)
  {
    logLike1 <- log.h1star.y1 + log.h3star.y2 - (k1 + k2.y1)
    logLike2 <- log.h2star.y1 - k2.y1  ## Making in terms of y1
    logLike3 <- log.h1star.y1 - (k1 + k2.y1)
    logLike4 <- - k2.y1  ## Making in terms of y1
    logLike5 <- log.h3star.y2 - k3
    logLike6 <- - k3
  }
  ##
  loglh <- sum(logLike1[type1==1]) + sum(logLike2[type2==1]) + sum(logLike3[type3==1]) + sum(logLike4[type4==1]) + sum(logLike5[type5==1]) + sum(logLike6[type6==1])
  ##
  return(-loglh)
}

## Fit model: illness-death, shared frailty, left-truncated data
## Weibull baseline hazards
FreqID.LT <- function(Y, lin.pred, data, model = "semi-Markov", startVals, frailty=TRUE, method="nlm")
{	
  
  
  Y = dat1
  
  ##
  y1 <- rowsums(cbind(dat1$y12, dat1$y13), na.rm = T)
  y1 <- dat1$y12
  delta1 <- dat1$s12
  y2     <- dat1$y13
  delta2 <- as.vector(dat1$s13)
  l      <- as.vector(Y[,5])
  lin.pred = 0
  Xmat1  <- as.matrix(rep(0, n/2))  
  Xmat2  <- as.matrix(rep(0, n/2))  
  Xmat3  <- as.matrix(rep(0, n/2))
  ##
  fit.survreg.1 <- survreg(as.formula(Surv(y1, delta1) ~
                                            1), dist="weibull", data=dat1)
  fit.survreg.2 <- survreg(as.formula(Surv(y2, delta2) ~1), dist="weibull", data=dat1)
  
  dat1$delta2= delta2
  
data.delta1_1 = dat1[delta1==1,]

plot(y1, y2)

  data.delta1_1$y2.m.y1 = y2[delta1==1] - y1[delta1==1]
  #fit.survreg.3 <- survreg(as.formula(Surv(y2.m.y1, delta2) ~ 1), dist="weibull", data=data.delta1_1)
  fit.survreg.3 <- survreg(as.formula(Surv(y23, s23) ~ 1), dist="weibull", data=dat1)
  
  
  alpha1      <- 1 / fit.survreg.1$scale
  alpha2      <- 1 / fit.survreg.2$scale
  alpha3     	<- 1 / fit.survreg.3$scale
  
  startVals =  NULL
  frailty = T
  
  if (is.null(startVals)==T){
    startVals     <- c(-alpha1*coef(fit.survreg.1)[1], log(alpha1),
                       -alpha2*coef(fit.survreg.2)[1], log(alpha2),
                       -alpha3*coef(fit.survreg.3)[1], log(alpha3))
    if(frailty == TRUE) startVals <- c(startVals, 0.5)
    startVals     <- c(startVals,
                       -coef(fit.survreg.1)[-1] * alpha1,
                       -coef(fit.survreg.2)[-1] * alpha2,
                       -coef(fit.survreg.3)[-1] * alpha3)
  }
  
  model = "semi-Markov"
  method = "nlm"
  if(model == "semi-Markov")
  {
    if (method == "optim"){
      cat("Fitting illness-death model with Weibull baseline hazards ... this should take < 1 min \n")
      logLike <- function(p) logLike.weibull.SCR.SM.LT(p, y1=y1, y2=y2, delta1=delta1, delta2=delta2, l=l,
                                                       Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3, frailty=frailty)
      optim.control = list(REPORT = 50)
      fit1 <- optim(startVals, #* runif(length(startVals), 0.9, 1.1), 
                    logLike, hessian = TRUE, method="Nelder-Mead", control = optim.control)
      value <- list(estimate=fit1$par, H=fit1$hessian, logLike=-fit1$value, code=fit1$convergence)#, Xmat=list(Xmat1, Xmat2, Xmat3))
    }
    if (method == "nlm"){
      cat("Fitting illness-death model with Weibull baseline hazards ... this should take < 1 min \n")
      fit1  <- suppressWarnings(nlm(logLike.weibull.SCR.SM.LT, p=startVals,
                                    y1=y1, delta1=delta1, y2=y2, delta2=delta2, Xmat1=as.matrix(Xmat1), Xmat2=as.matrix(Xmat2), Xmat3=as.matrix(Xmat3),
                                    l = l, frailty=frailty,
                                    iterlim=1000, hessian=TRUE))
      value <- list(estimate=fit1$est, H=fit1$hessian, logLike=-fit1$minimum, code=fit1$code)
    }
  }
  ##
  if(model == "semi-Markov")
  {
    class(value) <- c("Freq", "ID", "Ind", "WB", "semi-Markov")
  }
  
  return(value)
  ##
  invisible()
}



## For plotting WB baseline hazard function
WB.haz <- function(log.kappa, log.alpha, t){
  kappa = exp(log.kappa); alpha = exp(log.alpha)
  h = alpha * kappa * t^(alpha - 1)
  return(h)
}

## Generating simulated data used in simulation study
genWB.simData <- function(n = 5000, frailty = T, seed = 1){
  
  theta.true = ifelse(frailty == T, 0.25, 0)
  
  log.kappa1.true = -9.98  
  log.alpha1.true = 1.05 
  log.kappa2.true = -10.01
  log.alpha2.true = 1.15  
  log.kappa3.true = -5.92  
  log.alpha3.true = 0.92  
  ##
  kappa1.true = exp(log.kappa1.true)
  alpha1.true = exp(log.alpha1.true)
  kappa2.true = exp(log.kappa2.true)
  alpha2.true = exp(log.alpha2.true)
  kappa3.true = exp(log.kappa3.true)
  alpha3.true = exp(log.alpha3.true)
  
  beta1.true = c(-0.03021803)  
  beta2.true = c(-0.33064510)  
  beta3.true = c(-0.10652843)
  
  ## left-truncation params - drawing from uniform(Q1.l, Q3.l)
  Q1.l = 65  
  Q3.l = 78
  
  ## Male 
  prev.gender = 0.57
  
  cens = c(30, 30)
  
  ## B-spline baseline hazard specifications
  n.internalKnots.1 <- 1
  Bspline.degree.1 <- 1
  num.Bspline.params.1 <- (n.internalKnots.1) + (Bspline.degree.1 + 1)
  
  n.internalKnots.2 <- 1
  Bspline.degree.2 <- 1
  num.Bspline.params.2 <- (n.internalKnots.2) + (Bspline.degree.2 + 1)
  
  n.internalKnots.3 <- 1  ## Trying to add more flexibility - the previously estimated BH was off!
  Bspline.degree.3 <- 2
  num.Bspline.params.3 <- (n.internalKnots.3) + (Bspline.degree.3 + 1)
  
  lin.pred = list(as.formula(~gender), as.formula(~gender), as.formula(~gender))
  
  num.Bspline.params.1 = (n.internalKnots.1) + (Bspline.degree.1 + 1)
  num.Bspline.params.2 = (n.internalKnots.2) + (Bspline.degree.2 + 1)
  num.Bspline.params.3 = (n.internalKnots.3) + (Bspline.degree.3 + 1)
  
  n.params <- ifelse(frailty == T, num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 + length(beta1.true) + length(beta2.true) + length(beta3.true) + 1, 
                     num.Bspline.params.1 + num.Bspline.params.2 + num.Bspline.params.3 + length(beta1.true) + length(beta2.true) + length(beta3.true))
  
  BSparamNames <- c(paste0("phi1.", c(0:(num.Bspline.params.1-1))), paste0("phi2.", c(0:(num.Bspline.params.2-1))), paste0("phi3.", c(0:(num.Bspline.params.3-1))))
  WBparamNames <- c("log(kappa1)", "log(alpha1)", "log(kappa2)", "log(alpha2)", "log(kappa3)", "log(alpha3)")
  if (frailty == T){
    BSparamNames <- c(BSparamNames, "log(theta)")
    WBparamNames <- c(WBparamNames, "log(theta)") 
  }
  BSparamNames <- c(BSparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
  WBparamNames <- c(WBparamNames,  paste0("beta1.", 1:length(beta1.true)), paste0("beta2.", 1:length(beta2.true)), paste0("beta3.", 1:length(beta3.true)))
  
  ## DATA
  set.seed(seed)  
  gender = rbinom(n, 1, prev.gender)
  data0 = data.frame(gender)
  x1 = x2 = x3 = as.matrix(data0)
  
  lt = c(Q1.l, Q3.l) - 65  # Origin starts at age 65
  
  Y.tmp = data.frame(simID.LT(x1, x2, x3, beta1.true, beta2.true, beta3.true,
                              alpha1.true, alpha2.true, alpha3.true,
                              kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens, lt = lt), data0)
  ## Make sure left-truncated data
  Y = Y.tmp %>% filter(y1 > L)
  data = data.frame(Y)
  
  delta1 = Y$delta1; delta2 = Y$delta2
  y1 = Y[,1]; delta1 = Y[,2]; y2 = Y[,3]; delta2 = Y[,4]; l=Y[,5]
  
  ## Estimate B-spline params
  ## log h0(t) = phi0 * B0(t) + ... + phik * Bk(t)
  bdy.knots.b.1 = c(0, max(y1)) 
  bdy.knots.b.2 = c(0, max(y2))
  bdy.knots.b.3.y2my1 = c(0, max(y2-y1))
  
  ## 1-transition
  knot.loc.1 <- quantile(y1[delta1==1], ppoints(n.internalKnots.1))
  b.1.event <- bSpline(y1[delta1==1], knots = knot.loc.1, degree = Bspline.degree.1, intercept=TRUE, Boundary.knots=bdy.knots.b.1)
  b.1 <- predict(b.1.event, y1)
  h0.1.truth.event <- alpha1.true*kappa1.true*y1[delta1==1]^(alpha1.true-1)
  phi.1.truth=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.1), spline.pred=b.1.event, h0.truth=h0.1.truth.event)$estimate
  
  ## 2-transition
  knot.loc.2 <- quantile(y2[delta2==1], ppoints(n.internalKnots.2))
  b.2.event <- bSpline(y2[delta2==1], knots = knot.loc.2, degree = Bspline.degree.2, intercept=TRUE, Boundary.knots=bdy.knots.b.2)
  b.2 <- predict(b.2.event, y2)
  h0.2.truth.event <- alpha2.true*kappa2.true*y2[delta2==1]^(alpha2.true-1)
  phi.2.truth=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.2), spline.pred=b.2.event, h0.truth=h0.2.truth.event)$estimate
  
  ## 3-transition
  knot.loc.3 <- quantile((y2-y1)[delta1==1], ppoints(n.internalKnots.3))
  b.3.event <- bSpline((y2-y1)[delta1==1], knots = knot.loc.3, degree = Bspline.degree.3, intercept=TRUE, Boundary.knots=bdy.knots.b.3.y2my1)
  b.3.y2my1 <- predict(b.3.event, y2-y1)
  h0.3.truth.event <- alpha3.true*kappa3.true*(y2-y1)[delta1==1]^(alpha3.true-1)
  phi.3.truth=nlm(f=obj.fn.1, p=rep(-1,num.Bspline.params.3), spline.pred=b.3.event, h0.truth=h0.3.truth.event)$estimate
  
  ## Start values for b-spline simulations    
  startVals = c(phi.1.truth,
                phi.2.truth,
                phi.3.truth)
  if (frailty == T) startVals = c(startVals, log(theta.true))
  startVals = c(startVals, beta1.true, beta2.true, beta3.true)
  
  ## Params used to generate data 
  WBpara = c(log(kappa1.true), log(alpha1.true),
             log(kappa2.true), log(alpha2.true),
             log(kappa3.true), log(alpha3.true)) 
  if (frailty == T) WBpara = c(WBpara, log(theta.true))
  WBpara = c(WBpara, beta1.true, beta2.true, beta3.true)
  
  t1 <- t2 <- 0:30
  t3 <- 0:15
  
  
  return(list(Y = Y,
              lin.pred = lin.pred,
              data = data,
              startVals = startVals,
              frailty = frailty,frailty,
              b.1=b.1, b.2=b.2, b.3.y2my1=b.3.y2my1,
              bdy.knots.b.1=bdy.knots.b.1, 
              bdy.knots.b.2=bdy.knots.b.2, 
              bdy.knots.b.3.y2my1=bdy.knots.b.3.y2my1, 
              BSparamNames = BSparamNames,
              WBparamNames = WBparamNames,
              num.Bspline.params.1 = num.Bspline.params.1,
              num.Bspline.params.2 = num.Bspline.params.2,
              num.Bspline.params.3 = num.Bspline.params.3,
              WBpara = WBpara, 
              plot.WB1 = list(x = t1, y = WB.haz(log.kappa1.true, log.alpha1.true, t1)),
              plot.WB2 = list(x = t2, y = WB.haz(log.kappa2.true, log.alpha2.true, t2)),
              plot.WB3 = list(x = t3, y = WB.haz(log.kappa3.true, log.alpha3.true, t3))))  
}  




