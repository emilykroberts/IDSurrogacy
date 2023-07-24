#' create true cep curve
#'
#' @description create cep true based on true values
#'
#' @param dat0 data for placebo arm
#' @param dat1 data for treatment arm
#' @param write logical value if results should be written to a file
#' @param params_list true values of parameters
#' @param plotwrite logical value if plot should be saved to a file
#' @param fignum figure number for saving
#' @param tau_s time of evaluation for S
#' @param tau_t time of evaluation for T
#'
#' @return true values of CEP curve
#'
#' @examples 
#' example(true_cep(dat0 = dat0, dat1 = dat1, write = T, params_list = 
#' params_list, plotwrite = T))
true_cep = function(dat0, dat1, write, params_list, plotwrite, fignum, tau_s, tau_t){
  
  c12_0 = params_list$c12_0 
  c13_0 = params_list$c13_0 
  c23_0 = params_list$c23_0 
  c12_1 = params_list$c12_1 
  c13_1 = params_list$c13_1 
  c23_1 = params_list$c23_1
  shape12_1 = params_list$shape12_1 
  shape13_1 = params_list$shape13_1 
  shape23_1 = params_list$shape23_1 
  theta23_0 = params_list$theta23_0 
  theta23_1 = params_list$theta23_1 
  beta12_0 = params_list$beta12_0 
  beta13_0 = params_list$beta13_0 
  beta23_0 = params_list$beta23_0 
  beta12_1 = params_list$beta12_1 
  beta13_1 = params_list$beta13_1
  beta23_1 = params_list$beta23_1
  scale12_0 = params_list$scale12_0 
  scale13_0 = params_list$scale13_0 
  scale23_0 = params_list$scale23_0 
  scale12_1 = params_list$scale12_1 
  scale13_1 = params_list$scale13_1 
  scale23_1 = params_list$scale23_1 
  shape12_0 = params_list$shape12_0 
  shape13_0 = params_list$shape13_0 
  shape23_0 = params_list$shape23_0
  Fsave1 = Fsave = Fw_0 = Fw = Fw_1 = NULL
  x = rep(0, n)
  
  intfunction = function(j, i, t){
    exp(
      - cLambda13_frailty_lk(x = t, xdata1 = x[i], omega2 = omega13true0[i], scale13 = scale13_0, 
                             shape13 = shape13_0, c13 = c13_0, beta13_1 = beta13_0) - 
        cLambda12_frailty_lk(x = t, xdata1 = x[i], omega1 = omega12true0[i], scale12 = scale12_0, 
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
  
  s0cumulative = - pweibull(q = tau_s, scale = 1 / scale12_0, shape = shape12_0, lower = F, log = T) * exp(beta12_0 * x + c12_0 * omega12true0)
  s1cumulative = - pweibull(q = tau_s, scale = 1 / scale12_1, shape = shape12_1, lower = F, log = T) * exp(beta12_1 * x + c12_1 * omega12true1)
  
  dat = data.frame(cbind(s0cumulative / s1cumulative, Fw_1 - Fw_0))
  dat$X = log(s0cumulative / s1cumulative)
  dat$Y = c(Fw_1 - Fw_0)
  dat$x = x
  
  time23 = rowSums(cbind(dat0[, 1], dat0[, 3]), na.rm = T)
  time23_1 = rowSums(cbind(dat1[, 1], dat1[, 3]), na.rm = T)
  
  time23[dat0[, 2] == 0] = NA
  time23_1[dat1[, 2] == 0] = NA
  
  reg = lm(formula = Y ~ X, data = dat) 
  
  fname = paste('cep', '.n', n, array_id, '.txt', sep = "")
  res = cbind(summary(reg)$coef[1, 1], summary(reg)$coef[2, 1])
  print(res)
  if(write) write.table(res, file = fname, sep = "\t", row.names = F, col.names = T)
  
  d2 = ggplot(dat, aes((X), Y, alpha = 0.01)) + geom_point(alpha = 0.5) + theme_classic() + 
    ggtitle("Illness - Death CEP Curve with True Values" ) + 
    theme(legend.position = "none") + 
    xlab("Delta S_i") + ylim( - 1, 1) + 
    ylab(TeX("$\\Delta T_i = P(T_i(1) > \\tau_T|\\omega_{.i}^1) - P(T_i(0) > \\tau_T |\\omega_{.i}^0)$$ ")) + 
    xlab(TeX("$\\Delta S_i = log \\frac{\\Lambda_{12}^{0} (\\tau_S|\\omega_{12i}^0)}{\\Lambda_{12}^{1} (\\tau_S|\\omega_{12i}^1)}$$ ")) + 
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
    geom_hline(yintercept = mean(dat$Y, na.rm = T), linetype = "dashed", col = 'red') + geom_vline(xintercept = mean((dat$X), na.rm = T), linetype = "dashed", col = 'red') + 
    theme(text = element_text(size = 16)) + geom_smooth(method = 'lm', level = 0.95) + coord_cartesian(xlim = c( - 5, 5)); 
  
  d2 = d2 + theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 50, l = 0, unit = "mm")))
  line_1 = expr(paste("Intercept ", !!c(round(summary(reg)$coef[1, 1], 3)), " Slope ", !!c(round(summary(reg)$coef[2, 1], 3)), " R\U00B2 = ", !!round(summary(reg)$r.squared, 3), 
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
  
  line_6 = expr(paste( "by ", tau[s], " n12(0) ", !!sum(dat0[dat0[, 1] < tau_s, 2], na.rm = T), ", ",
                       "n13(0) ", !!sum(dat0[dat0[, 3] < tau_s, 4], na.rm = T), ", ",
                       "n23(0) ", !!sum(dat0[time23 < tau_s, 6], na.rm = T), ", ",
                       "n12(1) ", !!sum(dat1[dat1[, 1] < tau_s, 2], na.rm = T), ", ",
                       "n13(1) ", !!sum(dat1[dat1[, 3] < tau_s, 4], na.rm = T), ", ",
                       "n23(1) ", !!sum(dat1[time23_1 < tau_s, 6], na.rm = T),
                       sep = ""))
  
  line_7 = expr(paste( "by ", tau[t], " n12(0) ", !!sum(dat0[dat0[, 1] < tau_t, 2], na.rm = T), ", ",
                       "n13(0) ", !!sum(dat0[dat0[, 3] < tau_t, 4], na.rm = T), ", ",
                       "n23(0) ", !!sum(dat0[time23 < tau_t, 6], na.rm = T), ", ",
                       "n12(1) ", !!sum(dat1[dat1[, 1] < tau_t, 2], na.rm = T), ", ",
                       "n13(1) ", !!sum(dat1[dat1[, 3] < tau_t, 4], na.rm = T), ", ",
                       "n23(1) ", !!sum(dat1[time23_1 < tau_t, 6], na.rm = T),
                       sep = ""))
  
  d3 = ggdraw(d2) +
    draw_label(line_1, x = 0.55, y = 0.215) + # use relative coordinates for positioning
    draw_label(line_2, x = 0.55, y = 0.185) +
    draw_label(line_3, x = 0.55, y = 0.145) +
    draw_label(line_4, x = 0.55, y = 0.105) +
    draw_label(line_5, x = 0.55, y = 0.065) +
    draw_label(line_6, x = 0.55, y = 0.035) +
    draw_label(line_7, x = 0.55, y = 0.010)
  
  print(d3)
  if(plotwrite) ggsave(paste0("Figure", fignum, ".jpeg"), d3)
  return(d3)
  
}
#' create true cep curve
#'
#' @description create cep true based on true values
#'
#' @param dat0 data for placebo arm
#' @param dat1 data for treatment arm
#' @param write logical value if results should be written to a file
#' @param params_list true values of parameters
#' @param plotwrite logical value if plot should be saved to a file
#' @param fignum figure number for saving
#' @param tau_s time of evaluation for S
#' @param tau_t time of evaluation for T
#'
#' @return true values of CEP curve
#'
#' @examples 
#' example(true_cep(dat0 = dat0, dat1 = dat1, write = T, params_list = 
#' params_list, plotwrite = T))
true_cep = function(dat0, dat1, write, params_list, plotwrite, fignum, tau_s, tau_t){
  
  c12_0 = params_list$c12_0 
  c13_0 = params_list$c13_0 
  c23_0 = params_list$c23_0 
  c12_1 = params_list$c12_1 
  c13_1 = params_list$c13_1 
  c23_1 = params_list$c23_1
  shape12_1 = params_list$shape12_1 
  shape13_1 = params_list$shape13_1 
  shape23_1 = params_list$shape23_1 
  theta23_0 = params_list$theta23_0 
  theta23_1 = params_list$theta23_1 
  beta12_0 = params_list$beta12_0 
  beta13_0 = params_list$beta13_0 
  beta23_0 = params_list$beta23_0 
  beta12_1 = params_list$beta12_1 
  beta13_1 = params_list$beta13_1
  beta23_1 = params_list$beta23_1
  scale12_0 = params_list$scale12_0 
  scale13_0 = params_list$scale13_0 
  scale23_0 = params_list$scale23_0 
  scale12_1 = params_list$scale12_1 
  scale13_1 = params_list$scale13_1 
  scale23_1 = params_list$scale23_1 
  shape12_0 = params_list$shape12_0 
  shape13_0 = params_list$shape13_0 
  shape23_0 = params_list$shape23_0
  Fsave1 = Fsave = Fw_0 = Fw = Fw_1 = NULL
  x = rep(0, n)
  
  intfunction = function(j, i, t){
    exp(
      - cLambda13_frailty_lk(x = t, xdata1 = x[i], omega2 = omega13true0[i], scale13 = scale13_0, 
                             shape13 = shape13_0, c13 = c13_0, beta13_1 = beta13_0) - 
        cLambda12_frailty_lk(x = t, xdata1 = x[i], omega1 = omega12true0[i], scale12 = scale12_0, 
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
  
  s0cumulative = - pweibull(q = tau_s, scale = 1 / scale12_0, shape = shape12_0, lower = F, log = T) * exp(beta12_0 * x + c12_0 * omega12true0)
  s1cumulative = - pweibull(q = tau_s, scale = 1 / scale12_1, shape = shape12_1, lower = F, log = T) * exp(beta12_1 * x + c12_1 * omega12true1)
  
  dat = data.frame(cbind(s0cumulative / s1cumulative, Fw_1 - Fw_0))
  dat$X = log(s0cumulative / s1cumulative)
  dat$Y = c(Fw_1 - Fw_0)
  dat$x = x
  
  time23 = rowSums(cbind(dat0[, 1], dat0[, 3]), na.rm = T)
  time23_1 = rowSums(cbind(dat1[, 1], dat1[, 3]), na.rm = T)
  
  time23[dat0[, 2] == 0] = NA
  time23_1[dat1[, 2] == 0] = NA
  
  reg = lm(formula = Y ~ X, data = dat) 
  
  fname = paste('cep', '.n', n, array_id, '.txt', sep = "")
  res = cbind(summary(reg)$coef[1, 1], summary(reg)$coef[2, 1])
  print(res)
  if(write) write.table(res, file = fname, sep = "\t", row.names = F, col.names = T)
  
  d2 = ggplot(dat, aes(log(X), Y, alpha = 0.01)) + geom_point(alpha = 0.5) + theme_classic() + 
    ggtitle("Illness - Death CEP Curve with True Values" ) + 
    theme(legend.position = "none") + 
    xlab("Delta S_i") + ylim( - 1, 1) + 
    ylab(TeX("$\\Delta T_i = P(T_i(1) > \\tau_T|\\omega_{.i}^1) - P(T_i(0) > \\tau_T |\\omega_{.i}^0)$$ ")) + 
    xlab(TeX("$\\Delta S_i = log \\frac{\\Lambda_{12}^{0} (\\tau_S|\\omega_{12i}^0)}{\\Lambda_{12}^{1} (\\tau_S|\\omega_{12i}^1)}$$ ")) + 
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
    geom_hline(yintercept = mean(dat$Y, na.rm = T), linetype = "dashed", col = 'red') + geom_vline(xintercept = mean(log(dat$X), na.rm = T), linetype = "dashed", col = 'red') + 
    theme(text = element_text(size = 16)) + geom_smooth(method = 'lm', level = 0.95) + coord_cartesian(xlim = c( - 5, 5)); 
  
  d2 = d2 + theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 50, l = 0, unit = "mm")))
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
  
  line_6 = expr(paste( "by ", tau[s], " n12(0) ", !!sum(dat0[dat0[, 1] < tau_s, 2], na.rm = T), ", ",
                       "n13(0) ", !!sum(dat0[dat0[, 3] < tau_s, 4], na.rm = T), ", ",
                       "n23(0) ", !!sum(dat0[time23 < tau_s, 6], na.rm = T), ", ",
                       "n12(1) ", !!sum(dat1[dat1[, 1] < tau_s, 2], na.rm = T), ", ",
                       "n13(1) ", !!sum(dat1[dat1[, 3] < tau_s, 4], na.rm = T), ", ",
                       "n23(1) ", !!sum(dat1[time23_1 < tau_s, 6], na.rm = T),
                       sep = ""))
  
  line_7 = expr(paste( "by ", tau[t], " n12(0) ", !!sum(dat0[dat0[, 1] < tau_t, 2], na.rm = T), ", ",
                       "n13(0) ", !!sum(dat0[dat0[, 3] < tau_t, 4], na.rm = T), ", ",
                       "n23(0) ", !!sum(dat0[time23 < tau_t, 6], na.rm = T), ", ",
                       "n12(1) ", !!sum(dat1[dat1[, 1] < tau_t, 2], na.rm = T), ", ",
                       "n13(1) ", !!sum(dat1[dat1[, 3] < tau_t, 4], na.rm = T), ", ",
                       "n23(1) ", !!sum(dat1[time23_1 < tau_t, 6], na.rm = T),
                       sep = ""))
  
  d3 = ggdraw(d2) +
    draw_label(line_1, x = 0.55, y = 0.215) + # use relative coordinates for positioning
    draw_label(line_2, x = 0.55, y = 0.185) +
    draw_label(line_3, x = 0.55, y = 0.145) +
    draw_label(line_4, x = 0.55, y = 0.105) +
    draw_label(line_5, x = 0.55, y = 0.065) +
    draw_label(line_6, x = 0.55, y = 0.035) +
    draw_label(line_7, x = 0.55, y = 0.010)
  
  print(d3)
  if(plotwrite) ggsave(paste0("Figure", fignum, ".jpeg"), d3)
  return(d3)
  
}
