library(truncnorm); library(latex2exp); library(cowplot); library(xtable); library(here); library(ggpubr)
library(ggplot2); library(survival); library(survminer); library(mvtnorm); library(frailtypack)
setwd("./R")

############################ Run Simulation Iteration ############################
array_id = 10
#array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) # allows for running on cluster
#array_id = as.numeric(commandArgs(TRUE)) # allows for running on cluster
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

ss = rep(1:8, 200)
scenario = ss[array_id]

write = T
plotwrite = T

# keep source files in a folder called R directory
for(i in 1:(length(list.files(pattern = "\\.R$")))){
  source(list.files(pattern = "\\.R$")[i]
  )
}
setwd("./results/Table2") # save results of simulations in results Table 2 folder
set.seed(1 + array_id)

dat = sim_data(n = n, array_id = array_id, scenario = scenario, effecttheta = 0, 
               rhos = rhos, rhot = rhot, rhost = rhost, frailtysd = frailtysd, 
               diffscale1323 = diffscale1323, specify = specify, equalfrail = equalfrail, 
               independent = independent) # simulate data using sim_data

true_params = dat$params
set.seed(1 + array_id)
params_res = run_sim(SIM = SIM, array_id = array_id, rhos = rhos, rhot = rhot, frailtysd = frailtysd, 
                     params_list = true_params, tau_s = tau_s, tau_t = tau_t, dat0 = dat$dat0, dat1 = dat$dat1,
                     proposalsdfrail = proposalsdfrail, proposalsd = proposalsd, MFS = F,
                     independent = independent, equalfrail = equalfrail, holdshape = F, holdtheta = F
) # run simulation using run_sim function

# write final results params,2,10.txt denoting scenario and simulation ID (array_id) for simulated dataset
final_results(params_matrix = params_res, write = write)
# create Figure3.jpeg (Figure 3) in manuscript
# note: figure in manuscript corresponds to array_id = 10
if(array_id == 10){
  plot_results(params_matrix = params_res, write = plotwrite, fignum = 3)
}
