library(truncnorm); library(latex2exp); library(cowplot); library(xtable); library(here); library(ggpubr)
library(ggplot2); library(survival); library(mvtnorm); library(frailtypack); library(survminer); 
setwd("./R")

############################ Run Simulation Iteration ############################
array_id = 402
#array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) # allows for running on cluster
#array_id = as.numeric(commandArgs(TRUE)) # allows for running on cluster
if(is.na(array_id)) array_id = 402
n = 600 # set various parameters to run simulations (simulation is run first in this file)
rhost = rhos = rhot = 0.5
SIM = 3000
specify = F
effecttheta = 0
equalfrail = T
independent = T
diffscale1323 = T

proposalsd = 0.1
proposalsdfrail = 0.003
frailtysd = 0.4

tau_s = 1
tau_t = 2

scenario = 2

write = T
plotwrite = F

# keep source files in a folder called R directory
for(i in 1:(length(list.files(pattern = "\\.R$")))){
  source(list.files(pattern = "\\.R$")[i]
  )
}

setwd("./results/supplement/TableS3") # save supplemental simulations results in TableS3 folder

# corresponds to Supplement Table 3 tableS3results.csv - misspecify distribution
# simulate data
set.seed(1 + array_id)
dat = sim_data(n = n, array_id = array_id, scenario = scenario, effecttheta = effecttheta, 
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

params_list = true_params = dat$params

# run simulation
set.seed(1 + array_id)
params_res = run_sim(SIM = SIM, array_id = array_id, rhos = rhos, rhot = rhot, frailtysd = frailtysd, 
                     params_list = true_params, tau_s = tau_s, tau_t = tau_t, dat0 = dat0, dat1 = dat1,
                     proposalsdfrail = proposalsdfrail, proposalsd = proposalsd, MFS = F,
                     independent = independent, equalfrail = equalfrail, holdshape = F, holdtheta = F
)

final_results(params_matrix = params_res, write = write)

