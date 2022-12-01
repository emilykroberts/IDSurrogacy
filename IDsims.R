library(truncnorm); library(latex2exp); library(cowplot); library(xtable)
library(ggplot2); library(survival); library(survminer); library(frailtyEM); library(frailtypack); library(mvtnorm)
library(here)

setwd("./R")

array_id = 1; array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) # allows for running on cluster
if(is.na(array_id)) array_id = 1
n = 600 # set various parameters to run simulations
rhost = rhos = rhot = .5
SIM = 3000

effecttheta = F
holdfrail12 = F 
holdfrail13 = F
holdtheta = F
holdscale12 =
holdscale23 = 
holdscale13 = F
holdc = T
holdc23 = T
holdshape = F

equalfrail = T
independent = T
diffscale1323 = T

proposalsd = 0.1
proposalsdfrail = 0.003
frailtysd = 0.4

cep = T
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

set.seed(1 + array_id)

# simulate data
dat = sim_data(n = n, array_id = array_id, scenario = scenario, effecttheta = effecttheta,
      rhos = rhos, rhot = rhot, rhost = rhost, frailtysd = frailtysd, diffscale1323 = diffscale1323)

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

true_cep(dat0 = dat0, dat1 = dat1, write = write, params_list = true_params, plotwrite = plotwrite)

params_list = true_params

# run simulation
set.seed(1  + array_id)
params_res = run_sim(SIM = SIM, rhos = rhos, rhot = rhot, frailtysd = frailtysd, params_list = true_params)

# view results
plot_traceplots(params_matrix = params_res, variable = "int")
plot_traceplots(params_matrix = params_res, variable = "slope")

final_results(params_matrix = params_res, write = write)

plot_results(params_matrix = params_res, write = plotwrite)

### after simulations are run, can read in results or analyze pseudo data example based on where results are stored
if(F){
# read in results from many simulations
params <- do.call(rbind, sapply(list.files(pattern="params,1", full.names=TRUE), read.table, head = TRUE, simplify = F)); round(colMeans(params)[c(T,F)], 3) # means

xtable( rbind(rbind(round(colMeans(params[c(
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

####
setwd('..')

# create plots for real data
setwd('..')
rtog = read.csv("simulated_data.csv")
ST = (cbind(rtog$metastatic_prostate_cancer_years, rtog$survival_years - rtog$metastatic_prostate_cancer_years,
     rtog$survival_years))
status = cbind(as.numeric(rtog$metastatic_prostate_cancer == 1), as.numeric(rtog$metastatic_prostate_cancer == 1 &
     rtog$survival == 1), as.numeric(rtog$survival == 1 & rtog$metastatic_prostate_cancer == 2))

ST = data.frame(cbind(ST, status, rtog$trt))
names(ST) = c("y12", "y23", "y13", "s12", "s23", "s13")

dat0 = ST[ST$trt == 0,] # put data in form for run_sim function for analysis
dat1 = ST[ST$trt == 1,]
  
# cumulative incidence of S
ci = mstate::Cuminc(time = rtog$metastatic_prostate_cancer_years, status = rtog$metastatic_prostate_cancer,
                    data = rtog, group = rtog$trt, failcodes = c(1,2))
ci1 = ci[ci$group == 0,]; ci2 = ci[ci$group == 1,]
plot(c(0, ci1$time, max(ci1$time + 0.1)), c(0, ci1$CI.2, max(ci1$CI.2)),type = 's', ylim = c(0,1), lty = 'dashed', 
       ylab = "Probability", xlab = "Time from Randomization",
       main = "Estimates based on the cumulative incidence functions\nBlack dashed line indicates z = 1, solid red line indicates z = 0",
       lwd = 2); lines(c(0, ci2$time, max(ci2$time) + 0.1), c(0, ci2$CI.2, max(ci2$CI.2)),type = 's', col = 'red')

# KM plots
rtog$metastatic_prostate_cancer_cause = as.numeric(rtog$metastatic_prostate_cancer == 1)
rtog$metastatic_prostate_cancer_years_comp = as.numeric(rtog$metastatic_prostate_cancer_years != 0)
rtogs = rtog[rtog$metastatic_prostate_cancer == 1,]

a = ggsurvplot(survfit(Surv(rtog$metastatic_prostate_cancer_years, rtog$metastatic_prostate_cancer_cause) ~ rtog$trt), 
               data = rtog, risk.table = F, ggtheme = theme_classic2(base_size=20), legend.title = "Treatment", 
               legend.labs=c("With antiandrogen\ntherapy","Without antiandrogen\ntherapy ")) + 
               ggtitle("KM Curve of Intermediate Outcome S\n(Individuals are Censored at T)") + ylab("Freedom from S") + xlab("Time")
b = ggsurvplot(survfit(Surv(rtog$survival_years, rtog$survival) ~ rtog$trt), data = rtog, risk.table = F, 
               ggtheme = theme_classic2(base_size=20), legend.title = "Treatment", 
               legend.labs=c("With antiandrogen\ntherapy","Without antiandrogen\ntherapy ")) + 
  ggtitle("KM Curve of True Outcome T") + ylab("Freedom from T") + xlab("Time")
c = ggsurvplot(survfit(Surv((rtogs$survival_years - rtogs$metastatic_prostate_cancer_years), rtogs$survival) ~ rtogs$trt), 
               data = rtogs, risk.table = F, ggtheme = theme_classic2(base_size=20), legend.title = "Treatment", 
               legend.labs=c("With antiandrogen\ntherapy","Without antiandrogen\ntherapy ")) + 
               ggtitle("KM Curve of Time between S to T\nFor Those who Experienced S") + ylab("Freedom from T Survival (Post S)") + xlab("Time")

#ggsave(file = "KM_S.jpeg", a$plot, width = 8, height = 8) # save KM curves in manuscript
#ggsave(file = "KM_T.jpeg", b$plot, width = 8, height = 8)
#ggsave(file = "KM_ST.jpeg", c$plot, width = 8, height = 8)
}



