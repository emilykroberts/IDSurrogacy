library(truncnorm); library(latex2exp); library(cowplot); library(xtable); library(here); library(ggpubr)
library(ggplot2); library(survival); library(survminer); library(frailtyEM); library(frailtypack); library(mvtnorm)

setwd("./R")

array_id = 10; array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) # allows for running on cluster
if(is.na(array_id)) array_id = 10
n = 600 # set various parameters to run simulations (see details in corresponding functions)
rhost = rhos = rhot = 0.5
SIM = 3000
specify = T
effecttheta = 0
equalfrail = T
independent = T
diffscale1323 = T

proposalsd = 0.1
proposalsdfrail = 0.003
frailtysd = 0.4

tau_s = 1
tau_t = 2

ss = rep(1:8, 100)
scenario = ss[array_id]

write = T
plotwrite = T

# keep source files in a folder called R directory
for(i in 1:(length(list.files(pattern = "\\.R$")))){
  source(list.files(pattern = "\\.R$")[i]
  )
}

### step one: run main simulation study many times
set.seed(1 + array_id)

# simulate data
dat = sim_data(n = n, array_id = array_id, scenario = scenario, effecttheta = 0, 
               rhos = rhos, rhot = rhot, rhost = rhost, frailtysd = frailtysd, 
               diffscale1323 = diffscale1323, specify = specify, equalfrail = equalfrail, independent = independent)

dat0 = dat$dat0
dat1 = dat$dat1
omega12true0 = o12save0 = dat$o12save0
omega12true1 = o12save1 = dat$o12save1
omega13true0 = o13save0 = dat$o13save0
omega13true1 = o13save1 = dat$o13save1
omega23true0 = o23save0 = dat$o23save0
omega23true1 = o23save1 = dat$o23save1

true_params = dat$params

setwd("./results") # save results about true CEP curve and simulations in results folder

params_list = true_params

# run simulation
set.seed(1 + array_id)
params_res = run_sim(SIM = SIM, array_id = array_id, rhos = rhos, rhot = rhot, frailtysd = frailtysd, 
                     params_list = true_params, tau_s = tau_s, tau_t = tau_t, dat0 = dat0, dat1 = dat1,
                     proposalsdfrail = proposalsdfrail, proposalsd = proposalsd, MFS = F,
                     independent = independent, equalfrail = equalfrail, holdshape = F, holdtheta = F
)

final_results(params_matrix = params_res, write = write)
# saves results based on scenario and simulation ID (array_id) for later aggregation across many replications
# creates corresponding named Figure3.jpeg (Figure 3) in results for manuscript
plot_results(params_matrix = params_res, write = plotwrite, fignum = 3)

# after simulations are run, post-processing to aggregate results and recreate tables in manuscript
# for example, this individual file corresponds to results that can eventually contribute to Table 2 of aggregated simulation results
if(F){
  
  b = NULL
  
  fil = NULL
  for(i in 1:8){
    
    params <- do.call(rbind, sapply(list.files(path=fil, pattern=paste0("params,", i), full.names=TRUE), read.table, head = TRUE, simplify = F))
    
    a = ( rbind(rbind(round(colMeans(params[c(
      "int", "slope", "shape12_0", "shape13_0", "shape23_0",
      "shape12_1", "shape13_1", "shape23_1",
      "scale12_0", "scale13_0", "scale23_0",
      "scale12_1", "scale13_1", "scale23_1",
      "theta23_0", "theta23_1"
    )], na.rm = T), 3)),
    rbind(round(colMeans(params[c(
      "intse", "slopese", "shape12_0SE", "shape13_0SE", "shape23_0SE",
      "shape12_1SE", "shape13_1SE", "shape23_1SE",
      "scale12_0SE", "scale13_0SE", "scale23_0SE",
      "scale12_1SE", "scale13_1SE", "scale23_1SE",
      "theta23_0SE", "theta23_1SE"
    )], na.rm = T), 3)),
    rbind(round(apply(params[c(
      "int", "slope", "shape12_0", "shape13_0", "shape23_0",
      "shape12_1", "shape13_1", "shape23_1",
      "scale12_0", "scale13_0", "scale23_0",
      "scale12_1", "scale13_1", "scale23_1",
      "theta23_0", "theta23_1"
    )], 2, FUN = sd, na.rm = T), 3))))
    
    b = rbind(b, a)
    
  }
  
  # can be saved to produce tables such as table2results.csv
  write.csv(b, file = "table2results.csv")
}

######## sensitivity analyses and supplemental files that do not require many simulations run in parallel (running set to FALSE)
if(F){
  setwd("./supplement") # save results about true CEP curve and simulations in results/supplement folder
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
  
}
