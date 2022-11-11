cLambda13_frailty_lk<-function(x,xdata1, omega2, scale13, shape13, c13, beta13_1) {
  return(-pweibull(x, scale=1/scale13, shape = shape13, lower = F, log = T)*exp(beta13_1*xdata1 + c13*omega2))}
