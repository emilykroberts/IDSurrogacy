library(truncnorm); library(latex2exp); library(cowplot); library(xtable); library(here); library(ggpubr)
library(ggplot2); library(survival); library(survminer); library(mvtnorm); library(frailtypack)
setwd("./R")

######## sensitivity analyses and supplemental files that do not require many simulations run in parallel
array_id = 10
#array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) # allows for running on cluster
#array_id = as.numeric(commandArgs(TRUE)) # allows for running on cluster
if(is.na(array_id)) array_id = 10
n = 600 # set various parameters to run simulations (see details in corresponding functions)
rhost = rhos = rhot = 0.5
specify = T
effecttheta = 0
equalfrail = T
independent = T
diffscale1323 = T
frailtysd = 0.4

tau_s = 1
tau_t = 2

ss = rep(1:8, 200)
scenario = ss[array_id]

write = T
plotwrite = T

# keep source files in a folder called R directory
for(i in 1:(length(list.files(pattern = "\\.R$")))){
  source(list.files(pattern = "\\.R$")[i]
  )
}

setwd("./results/supplement") # save results about true CEP curve and simulations in results/supplement folder
set.seed(10)

# simulate data, generate true CEP curve with generative values for Supplemental Figure 1a
dat = sim_data(n = n, array_id = array_id, scenario = 1, effecttheta = 0, 
               rhos = 0.5, rhot = 0.5, rhost = 0.5, frailtysd = frailtysd, 
               diffscale1323 = diffscale1323, specify = specify, equalfrail = equalfrail, independent = independent)

dat0 = dat$dat0; dat1 = dat$dat1; params_list = dat$params
omega12true0 = dat$o12save0; omega12true1 = dat$o12save1; omega13true0 = dat$o13save0
omega13true1 = dat$o13save1; omega23true0 = dat$o23save0; omega23true1 = dat$o23save1

true_cep(dat0 = dat0, dat1 = dat1, write = F, params_list = params_list, plotwrite = T, fignum = 'S1a',
         tau_s = tau_s, tau_t = tau_t)

# simulate data, generate true CEP curve with generative values for Supplemental Figure 1b
dat = sim_data(n = n, array_id = array_id, scenario = 2, effecttheta = 0, 
               rhos = 0.5, rhot = 0.5, rhost = 0.5, frailtysd = frailtysd, 
               diffscale1323 = diffscale1323, specify = specify, equalfrail = equalfrail, independent = independent)

dat0 = dat$dat0; dat1 = dat$dat1; params_list = dat$params
omega12true0 = dat$o12save0; omega12true1 = dat$o12save1; omega13true0 = dat$o13save0
omega13true1 = dat$o13save1; omega23true0 = dat$o23save0; omega23true1 = dat$o23save1

true_cep(dat0 = dat0, dat1 = dat1, write = F, params_list = params_list, plotwrite = T, fignum = 'S1b',
         tau_s = tau_s, tau_t = tau_t)

# simulate data, generate true CEP curve with generative values for Supplemental Figure 1c
dat = sim_data(n = n, array_id = array_id, scenario = 3, effecttheta = 0, 
               rhos = 0.5, rhot = 0.5, rhost = 0.5, frailtysd = frailtysd, 
               diffscale1323 = diffscale1323, specify = specify, equalfrail = equalfrail, independent = independent)

dat0 = dat$dat0; dat1 = dat$dat1; params_list = dat$params
omega12true0 = dat$o12save0; omega12true1 = dat$o12save1; omega13true0 = dat$o13save0
omega13true1 = dat$o13save1; omega23true0 = dat$o23save0; omega23true1 = dat$o23save1

true_cep(dat0 = dat0, dat1 = dat1, write = F, params_list = params_list, plotwrite = T, fignum = 'S1c',
         tau_s = tau_s, tau_t = tau_t)

# simulate data, generate true CEP curve with generative values for Supplemental Figure 1d
dat = sim_data(n = n, array_id = array_id, scenario = 4, effecttheta = 0, 
               rhos = 0.5, rhot = 0.5, rhost = 0.5, frailtysd = frailtysd, 
               diffscale1323 = diffscale1323, specify = specify, equalfrail = equalfrail, independent = independent)

dat0 = dat$dat0; dat1 = dat$dat1; params_list = dat$params
omega12true0 = dat$o12save0; omega12true1 = dat$o12save1; omega13true0 = dat$o13save0
omega13true1 = dat$o13save1; omega23true0 = dat$o23save0; omega23true1 = dat$o23save1

true_cep(dat0 = dat0, dat1 = dat1, write = F, params_list = params_list, plotwrite = T, fignum = 'S1d',
         tau_s = tau_s, tau_t = tau_t)

# simulate data, generate true CEP curve with generative values for Supplemental Figure 2a
dat = sim_data(n = n, array_id = array_id, scenario = 5, effecttheta = 0, 
               rhos = 0.5, rhot = 0.5, rhost = 0.5, frailtysd = frailtysd, 
               diffscale1323 = diffscale1323, specify = specify, equalfrail = equalfrail, independent = independent)

dat0 = dat$dat0; dat1 = dat$dat1; params_list = dat$params
omega12true0 = dat$o12save0; omega12true1 = dat$o12save1; omega13true0 = dat$o13save0
omega13true1 = dat$o13save1; omega23true0 = dat$o23save0; omega23true1 = dat$o23save1

true_cep(dat0 = dat0, dat1 = dat1, write = F, params_list = params_list, plotwrite = T, fignum = 'S2a',
         tau_s = tau_s, tau_t = tau_t)

# simulate data, generate true CEP curve with generative values for Supplemental Figure 2b
dat = sim_data(n = n, array_id = array_id, scenario = 6, effecttheta = 0, 
               rhos = 0.5, rhot = 0.5, rhost = 0.5, frailtysd = frailtysd, 
               diffscale1323 = diffscale1323, specify = specify, equalfrail = equalfrail, independent = independent)

dat0 = dat$dat0; dat1 = dat$dat1; params_list = dat$params
omega12true0 = dat$o12save0; omega12true1 = dat$o12save1; omega13true0 = dat$o13save0
omega13true1 = dat$o13save1; omega23true0 = dat$o23save0; omega23true1 = dat$o23save1

true_cep(dat0 = dat0, dat1 = dat1, write = F, params_list = params_list, plotwrite = T, fignum = 'S2b',
         tau_s = tau_s, tau_t = tau_t)

# simulate data, generate true CEP curve with generative values for Supplemental Figure 2c
dat = sim_data(n = n, array_id = array_id, scenario = 7, effecttheta = 0, 
               rhos = 0.5, rhot = 0.5, rhost = 0.5, frailtysd = frailtysd, 
               diffscale1323 = diffscale1323, specify = specify, equalfrail = equalfrail, independent = independent)

dat0 = dat$dat0; dat1 = dat$dat1; params_list = dat$params
omega12true0 = dat$o12save0; omega12true1 = dat$o12save1; omega13true0 = dat$o13save0
omega13true1 = dat$o13save1; omega23true0 = dat$o23save0; omega23true1 = dat$o23save1

true_cep(dat0 = dat0, dat1 = dat1, write = F, params_list = params_list, plotwrite = T, fignum = 'S2c',
         tau_s = tau_s, tau_t = tau_t)

# simulate data, generate true CEP curve with generative values for Supplemental Figure 2d
dat = sim_data(n = n, array_id = array_id, scenario = 8, effecttheta = 0, 
               rhos = 0.5, rhot = 0.5, rhost = 0.5, frailtysd = frailtysd, 
               diffscale1323 = diffscale1323, specify = specify, equalfrail = equalfrail, independent = independent)

dat0 = dat$dat0; dat1 = dat$dat1; params_list = dat$params
omega12true0 = dat$o12save0; omega12true1 = dat$o12save1; omega13true0 = dat$o13save0
omega13true1 = dat$o13save1; omega23true0 = dat$o23save0; omega23true1 = dat$o23save1

true_cep(dat0 = dat0, dat1 = dat1, write = F, params_list = params_list, plotwrite = T, fignum = 'S2d',
         tau_s = tau_s, tau_t = tau_t)

set.seed(10)

# simulate data, generate true CEP curve with generative values for Supplemental Figure 3a
dat = sim_data(n = n, array_id = array_id, scenario = 2, effecttheta = 0, 
               rhos = 0.95, rhot = 0.95, rhost = 0.5, frailtysd = frailtysd, 
               diffscale1323 = diffscale1323, specify = specify, equalfrail = equalfrail, independent = independent)

dat0 = dat$dat0; dat1 = dat$dat1; params_list = dat$params
omega12true0 = dat$o12save0; omega12true1 = dat$o12save1; omega13true0 = dat$o13save0
omega13true1 = dat$o13save1; omega23true0 = dat$o23save0; omega23true1 = dat$o23save1

true_cep(dat0 = dat0, dat1 = dat1, write = F, params_list = params_list, plotwrite = T, fignum = 'S3a',
         tau_s = tau_s, tau_t = tau_t)

# simulate data, generate true CEP curve with generative values for Supplemental Figure 3b
dat = sim_data(n = n, array_id = array_id, scenario = 2, effecttheta = 0, 
               rhos = 0.0, rhot = 0.0, rhost = 0.5, frailtysd = frailtysd, 
               diffscale1323 = diffscale1323, specify = specify, equalfrail = equalfrail, independent = independent)

dat0 = dat$dat0; dat1 = dat$dat1; params_list = dat$params
omega12true0 = dat$o12save0; omega12true1 = dat$o12save1; omega13true0 = dat$o13save0
omega13true1 = dat$o13save1; omega23true0 = dat$o23save0; omega23true1 = dat$o23save1

true_cep(dat0 = dat0, dat1 = dat1, write = F, params_list = params_list, plotwrite = T, fignum = 'S3b',
         tau_s = tau_s, tau_t = tau_t)

set.seed(10)
# simulate data, generate true CEP curve with generative values for Supplemental Figure 4a
dat = sim_data(n = n, array_id = array_id, scenario = 2, effecttheta = -1, 
               rhos = 0.5, rhot = 0.5, rhost = 0.5, frailtysd = frailtysd, 
               diffscale1323 = diffscale1323, specify = specify, equalfrail = equalfrail, independent = independent)

dat0 = dat$dat0; dat1 = dat$dat1; params_list = dat$params
omega12true0 = dat$o12save0; omega12true1 = dat$o12save1; omega13true0 = dat$o13save0
omega13true1 = dat$o13save1; omega23true0 = dat$o23save0; omega23true1 = dat$o23save1

true_cep(dat0 = dat0, dat1 = dat1, write = F, params_list = params_list, plotwrite = T, fignum = 'S4a',
         tau_s = tau_s, tau_t = tau_t)

# simulate data, generate true CEP curve with generative values for Supplemental Figure 4b
dat = sim_data(n = n, array_id = array_id, scenario = 2, effecttheta = 0, 
               rhos = 0.5, rhot = 0.5, rhost = 0.5, frailtysd = frailtysd, 
               diffscale1323 = diffscale1323, specify = specify, equalfrail = equalfrail, independent = independent)

dat0 = dat$dat0; dat1 = dat$dat1; params_list = dat$params
omega12true0 = dat$o12save0; omega12true1 = dat$o12save1; omega13true0 = dat$o13save0
omega13true1 = dat$o13save1; omega23true0 = dat$o23save0; omega23true1 = dat$o23save1

true_cep(dat0 = dat0, dat1 = dat1, write = F, params_list = params_list, plotwrite = T, fignum = 'S4b',
         tau_s = tau_s, tau_t = tau_t)

# simulate data, generate true CEP curve with generative values for Supplemental Figure 4c
dat = sim_data(n = n, array_id = array_id, scenario = 2, effecttheta = 1, 
               rhos = 0.5, rhot = 0.5, rhost = 0.5, frailtysd = frailtysd, 
               diffscale1323 = diffscale1323, specify = specify, equalfrail = equalfrail, independent = independent)

dat0 = dat$dat0; dat1 = dat$dat1; params_list = dat$params
omega12true0 = dat$o12save0; omega12true1 = dat$o12save1; omega13true0 = dat$o13save0
omega13true1 = dat$o13save1; omega23true0 = dat$o23save0; omega23true1 = dat$o23save1

true_cep(dat0 = dat0, dat1 = dat1, write = F, params_list = params_list, plotwrite = T, fignum = 'S4c',
         tau_s = tau_s, tau_t = tau_t)

