cLambda23_frailty_lk<-function(x,xdata1,v_predict, omega2, shape23, scale23, omega1, theta23, c23, beta23_1) {
  return(-pweibull(x, scale=1/scale23, shape = shape23, lower = F, log = T)*exp(-(beta23_1*xdata1+c23*omega2+theta23*v_predict)))}
