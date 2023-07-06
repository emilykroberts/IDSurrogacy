#' create true cep curve
#'
#' @description create cep true based on true values
#'
#' @param dat0 data for placebo arm
#' @param dat1 data for treatment arm
#' @param write logical value if results should be written to a file
#' @param params_list true values of parameters
#' @param plotwrite logical value if plot should be saved to a file
#'
#' @return true values of CEP curve
#'
#' @examples 
#' example(true_cep(dat0 = dat0, dat1 = dat1, write = T, params_list = 
#' params_list, plotwrite = T))
true_cep = function(dat0, dat1, write, params_list, plotwrite){
  
  params_list = attach(params_list)
  Fsave1 = Fsave = Fw_0 = Fw = Fw_1 = NULL
  x = rep(0, n)
  
  intfunction = function(j, i, t){
    exp(
      - cLambda13_frailty_lk(x = t, xdata1 = x[i], omega2 = omega13true0[i], scale13 = scale13_0, 
                             shape13 = shape13_0, c13 = c13_0, beta13_1 = beta13_0) - 
        cLambda12_frailty_lk(x = t, xdata1 = x[i],  omega1 = omega12true0[i], scale12 = scale12_0, 
                             shape12 = shape12_0, c12 = c12_0, beta12_1 = beta12_0)) * 
      lambda12_frailty(t, xdata1 = x[i], omega1 = omega12true0[i], scale12 = scale12_0, shape12 = shape12_0, 
                       c12 = c12_0, beta12_1 = beta12_0) * 
      exp( - cLambda23_frailty_lk(x = tau_t - t, xdata1 = x[i], omega2 = omega23true0[i], scale23 = scale23_0, 
                                  shape23 = shape23_0, c23 = c23_0, theta23 = theta23_0, v_predict = t, beta23_1 = beta23_0))
  }
  
  
  for(i in 1:n){
    y = tryCatch({ integrate(f = intfunction, lower = 0, upper = tau_t, i = i, j = j)$val
    }, 
    error = function(error_condition) {
      NA}
    ) 
    
    L12 = cLambda12_frailty_lk(x = tau_t, xdata1 = x[i], shape12 = shape12_0, 
                               omega1 = omega12true0[i], scale12 = scale12_0, c12 = c12_0, beta12_1 = beta12_0)
    L13 = cLambda13_frailty_lk(x = tau_t, xdata1 = x[i], omega2 = omega13true0[i], scale13 = scale13_0, shape13 = shape13_0, 
                               c13 = c13_0, beta13_1 = beta13_0)
    
    Fw = (exp( - L13 - L12) + y)
    Fw_0 = c(Fw_0, Fw)
  }
  
  intfunction = function(j, i, t){
    exp(
      - cLambda13_frailty_lk(x = t, xdata1 = x[i], omega2 = omega13true1[i], scale13 = scale13_1, 
                             shape13 = shape13_1, c13 = c13_1, beta13_1 = beta13_1) - 
        cLambda12_frailty_lk(x = t, xdata1 = x[i], omega1 = omega12true1[i], scale12 = scale12_1, 
                             shape12 = shape12_1, c12 = c12_1, beta12_1 = beta12_1)) * 
      lambda12_frailty(t, xdata1 = x[i], omega1 = omega12true1[i], scale12 = scale12_1, shape12 = shape12_1, 
                       c12 = c12_1, beta12_1 = beta12_1) * 
      exp( - cLambda23_frailty_lk(x = tau_t - t, xdata1 = x[i], omega2 = omega23true1[i], scale23 = scale23_1, 
                                  shape23 = shape23_1, 
                                  c23 = c23_1, theta23 = theta23_1, v_predict = t, beta23_1 = beta23_1))
  }
  
  for(i in 1:n){
    y = tryCatch({ integrate(f = intfunction, lower = 0, upper = tau_t, i = i, j = j)$val
    }, 
    error = function(error_condition) { NA}
    ) 
    
    L12 = cLambda12_frailty_lk(x = tau_t, xdata1 = x[i], shape12 = shape12_1, 
                               omega1 = omega12true1[i], scale12 = scale12_1, c12 = c12_1, beta12_1 = beta12_1)
    L13 = cLambda13_frailty_lk(x = tau_t, xdata1 = x[i], omega2 = omega13true1[i], scale13 = scale13_1, shape13 = shape13_1, 
                               c13 = c13_1, beta13_1 = beta13_1)
    
    Fw = (exp( - L13 - L12) + y)
    if(Fw > 1) Fw = NA
    Fw_1 = c(Fw_1, Fw)
  }
  
  z = 1
  
  s0cumulative = - pweibull(q = tau_s, scale = 1 / scale12_0, shape = shape12_0, lower = F, log = T) * exp(beta12_0 * x + c12_0 * omega12true0)
  s1cumulative = - pweibull(q = tau_s, scale = 1 / scale12_1, shape = shape12_1, lower = F, log = T) * exp(beta12_1 * x + c12_1 * omega12true1)
  
  
  dat = data.frame(cbind(s0cumulative / s1cumulative, Fw_1 - Fw_0))
  dat$X = s0cumulative / s1cumulative
  dat$Y = c(Fw_1 - Fw_0)
  dat$x = x
  
  reg_true = reg = lm(formula = Y ~ log(X), data = dat) 
  
  fname = paste('cep', '.n', n, array_id, '.txt', sep = "")
  res = cbind(summary(reg)$coef[1, 1], summary(reg)$coef[2, 1])
  print(res)
  if(write) write.table(res, file = fname, sep = "\t", row.names = F, col.names = T)
  
  d2 = ggplot(dat, aes(log(X), Y, alpha = 0.01)) + geom_point(alpha = 0.5) + theme_classic() + 
    ggtitle("Illness - Death CEP Curve with True Values" ) + 
    theme(legend.position = "none") + 
    xlab("Delta S_i") + ylim( - 1, 1) + 
    ylab(TeX("$\\Delta T_i = P(T_i(1) > \\tau_T|x_i, \\omega_{.i}^1) - P(T_i(0) > \\tau_T |x_i, \\omega_{.i}^0)$$ ")) + 
    xlab(TeX("$\\Delta S_i = log \\frac{\\cLambda_{12}^{0} (\\tau_S| x_i, \\omega_{12i}^0)}{\\cLambda_{12}^{1} (\\tau_S| x_i, \\omega_{12i}^1)}$$ ")) + 
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
    geom_hline(yintercept = mean(dat$Y, na.rm = T), linetype = "dashed", col = 'red') + geom_vline(xintercept = mean(log(dat$X), na.rm = T), linetype = "dashed", col = 'red') + 
    theme(text = element_text(size = 16)) + geom_smooth(method = 'lm', level = 0.95) + coord_cartesian(xlim = c( - 5, 5)); 
  
  d2 = d2 + theme(axis.title.x = element_text(
    margin = margin(t = 5, r = 0, b = 50, l = 0, unit = "mm")))
  line_1 = expr(paste("Intercept ", !!c(round(summary(reg)$coef[1, 1], 3)), " Slope ", !!c(round(summary(reg)$coef[2, 1], 3)), " R\U00B2 = ", !!round(summary(reg_true)$r.squared, 3), 
                      sep = ""))
  
  line_2 = expr(paste("Generating Parameters ", alpha[12] ^ 0, " = ", !!shape12_0, ", ", 
                      alpha[13] ^ 0, " = ", !!shape13_0, ", ", 
                      alpha[23] ^ 0, " = ", !!shape23_0, ", ", 
                      alpha[12] ^ 1, " = ", !!shape12_1, ", ", 
                      alpha[13] ^ 1, " = ", !!shape13_1, ", ", 
                      alpha[23] ^ 1, " = ", !!shape23_1, 
                      sep = ""))
  line_3 = expr(paste(gamma[12] ^ 0, " = ", !!scale12_0, ", ", 
                      gamma[13] ^ 0, " = ", !!scale13_0, ", ", 
                      gamma[23] ^ 0, " = ", !!scale23_0, ", ", 
                      gamma[12] ^ 1, " = ", !!scale12_1, ", ", 
                      gamma[13] ^ 1, " = ", !!scale13_1, ", ", 
                      gamma[23] ^ 1, " = ", !!scale23_1, 
                      sep = ""))
  line_4 = expr(paste(kappa[12] ^ 0, " = ", !!c12_0, ", ", 
                      kappa[13] ^ 0, " = ", !!c13_0, ", ", 
                      kappa[23] ^ 0, " = ", !!c23_0, ", ", 
                      kappa[12] ^ 1, " = ", !!c12_1, ", ", 
                      kappa[13] ^ 1, " = ", !!c13_1, ", ", 
                      kappa[23] ^ 1, " = ", !!c23_1, 
                      sep = ""))
  
  line_5 = expr(paste( theta[23] ^ 0, " = ", !!theta23_0, ", ", 
                       theta[23] ^ 1, " = ", !!theta23_1, ", ", 
                       rho[s], " = ", !!rhos, ", ", 
                       rho[t], " = ", !!rhot, ", ", 
                       tau[s], " = ", !!tau_s, ", ", 
                       tau[t], " = ", !!tau_t, 
                       sep = ""))
  
  d3 = ggdraw(d2) + 
    draw_label(line_1, x = 0.55, y = 0.215) + # use relative coordinates for positioning
    draw_label(line_2, x = 0.55, y = 0.175) + 
    draw_label(line_3, x = 0.55, y = 0.135) + 
    draw_label(line_4, x = 0.55, y = 0.095) + 
    draw_label(line_5, x = 0.55, y = 0.055)
  print(d3)
  if(plotwrite) ggsave(paste0("CEPtrue", ", ", scenario, ", ", array_id, ".jpeg"), d3)
  
}
