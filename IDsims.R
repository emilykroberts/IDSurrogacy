library(MASS); library(mvtnorm); library(MCMCpack); library(LaplacesDemon)
library(ggplot2); library(survival); library(survminer); library(frailtyEM); library(frailtypack)
library(truncnorm); library(latex2exp); library(kableExtra); library(Rfast)
library(xtable); library(ggforce); library(wesanderson); library(cowplot); library(here)

array_id = 1; array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(is.na(array_id)) array_id = 1
{n = 600
rhost = rhos = .5; rhot = .5
SIM = 3000
holdtheta = F
cep = T
holdscale12 =
holdscale23 = 
holdscale13 = F

equalfrail = T
independent = T
diffscale1323 = T
holdshape = F

holdc = T
holdc23 = T

proposalsd = 0.1
proposalsdfrail = 0.003

effecttheta = F
holdfrail12 = F 
holdfrail13 = F
tau_s = 1
tau_t = 2
frailtysd = 0.4

ss = rep(1:8, 100)
scenario = ss[array_id]

write = T
plotwrite = as.numeric(array_id == 10)
}

setwd("/R") # keep source files in a folder called R directory
list.files()
for(i in 1:(length(list.files()))){
  source(  list.files()[i]
)
}

set.seed(1 + array_id)

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

true_cep(dat0 = dat0, dat1 = dat1, write = write, params_list = true_params, plotwrite = plotwrite)

params_list = true_params

set.seed(1  + array_id)
params_res = run_sim(SIM = SIM, rhos = rhos, rhot = rhot, frailtysd = frailtysd, params_list = true_params)

plot_traceplots(params_matrix = params_res, variable = "int")
plot_traceplots(params_matrix = params_res, variable = "slope")

final_results(params_matrix = params_res, write = write)

 
plot_results(params_matrix = params_res, write = plotwrite)


# read in results from many simulations
library(xtable)

params <- do.call(rbind, sapply(list.files(path=fil, pattern="params,1", full.names=TRUE), read.table, head = TRUE, simplify = F)); round(colMeans(params)[c(T,F)], 3) # means

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
# create plots for real data

rtog = read.csv("~/Dropbox (University of Michigan)/rtog9601/NCT00002874-D1-Dataset.csv")
 rtog = rtog[rtog$include_in_analysis == 1,]
 ST = (cbind(rtog$metastatic_prostate_cancer_years, rtog$survival_years - rtog$metastatic_prostate_cancer_years,
       rtog$survival_years))
 status = cbind(as.numeric(rtog$metastatic_prostate_cancer == 1), as.numeric(rtog$metastatic_prostate_cancer == 1 &
                                        rtog$survival == 1), as.numeric(rtog$survival == 1 & rtog$metastatic_prostate_cancer == 2))
 
 ST = data.frame(cbind(ST, status))
 names(ST) = c("y12", "y23", "y13", "s12", "s23", "s13")
 
 ST$trt = rtog$rx - 1
 dat0 = data.frame(ST[ST$"trt" == 0,])
 dat1 = data.frame(ST[ST$"trt" == 1,])
 
 dat0 = dat0[1:(n/2), ]
 dat1 = dat1[1:(n/2), ]
 
 # run prentice criteria
 time23 = rowSums(cbind(ST[,1], ST[,2]), na.rm = T)
 time23[status[,1] == 0] = NA
 Ttime = apply(cbind(ST[,3], time23), 1, min, na.rm = T) 
 Tstat = apply(cbind(status[,3], status[,1]), 1, max, na.rm = T) 
 
 if(F){
 # 0 = censored, 1 = outcome of interest, 2 = competing event
 a = cmprsk::cuminc(rtog$metastatic_prostate_cancer_years, fstatus = rtog$metastatic_prostate_cancer, 
           cencode = 0)
 plot(a)

 ci_fit = 
  cmprsk::cuminc(
   ftime = rtog$metastatic_prostate_cancer_years, 
   fstatus = rtog$metastatic_prostate_cancer, 
   group = rtog$rx - 1,
   cencode = 0
  )

 
 # cumulative incidence of T
 a = cmprsk::cuminc(rtog$survival_years, fstatus = rtog$survival, group = rtog$rx, cencode = 0)
 plot(a)
 
 # cumulative incidence of S
 mstatcumul = mstate::Cuminc(time = rtog$metastatic_prostate_cancer_years, status = rtog$metastatic_prostate_cancer , data = rtog)
 
 plot(mstatcumul, conf.int = .95, xlab = "Time (years)",use.ggplot = TRUE, main = "Cumulative Incidence for Surrogate Endpoint S", lty = 1:4)
 
 ci = mstate::Cuminc(time = rtog$metastatic_prostate_cancer_years, status = rtog$metastatic_prostate_cancer , data = rtog, group = rtog$rx , failcodes = c(1,2))
 
 ci1 = ci[ci$group == 1,]
 ci2 = ci[ci$group == 2,]
 
 plot(c(0, ci1$time, 15), c(0, ci1$CI.1, max(ci1$CI.1)),
    type = 's', xlim = c(0, 14), ylim = c(0,
                       1), lty = 'dashed', ylab = "Probability", xlab = "Years from Randomization",
    main = "Estimates based on the cumulative incidence functions\nBlack dashed line indicates z = 1, solid red line indicates z = 0",
    lwd = 2); lines(c(0, ci2$time, 15), c(0, ci2$CI.2, max(ci2$CI.2)),
            type = 's', col = 'red')
 }
 
 t_mod = coxph(Surv(Ttime, Tstat) ~ trt, dat = ST)
 ts_mod = coxph(Surv(Ttime, Tstat) ~ y12 + trt, dat = ST)
 s_mod = coxph(Surv(y12, s12) ~ trt, dat = ST)
 
 as.numeric(summary(s_mod)$coefficients[,"Pr(>|z|)"] < 0.05)
 summary(t_mod)
 summary(ts_mod)
 
 exp(-coef(ts_mod)["trt"])
 exp(-coef(t_mod)["trt"])
 ST$id = 1:(nrow(ST))

 ST$Ttime = Ttime
 ST$Tstat = Tstat
 
 new = tmerge(data1 = ST, data2 = ST, id = id, tstop = rtog$survival_years)
 new2 = tmerge(data1 = new, data2 = ST, id = id, ev = event(y12))
 
 ts_mod2 = coxph(Surv(Ttime, Tstat) ~ y12 + trt, dat = new2, cluster = id)
 exp(-coef(ts_mod2)["trt"])
 
 # KM plots
 
 rtog$metastatic_prostate_cancer_cause = as.numeric(rtog$metastatic_prostate_cancer == 1)
 rtog$metastatic_prostate_cancer_years_comp = as.numeric(rtog$metastatic_prostate_cancer_years != 0)

 theme_set( theme_classic(base_size = 16))

 a = ggsurvplot(survfit(Surv(rtog$metastatic_prostate_cancer_years, rtog$metastatic_prostate_cancer_cause) ~ rtog$rx), data = rtog, risk.table = F, ggtheme = theme_classic2(base_size=20), legend.title = "Treatment", legend.labs=c("With antiandrogen\ntherapy","Without antiandrogen\ntherapy ")) + ggtitle("KM Curve of Intermediate Outcome S\n(Individuals are Censored at T)") + ylab("Freedom from S") + xlab("Time (years)")
 a

 b = ggsurvplot(survfit(Surv(rtog$survival_years, rtog$survival) ~ rtog$rx), data = rtog, risk.table = F, ggtheme = theme_classic2(base_size=20), legend.title = "Treatment", legend.labs=c("With antiandrogen\ntherapy","Without antiandrogen\ntherapy ")) + ggtitle("KM Curve of True Outcome T") + ylab("Freedom from T") + xlab("Time (years)")
 b

 rtogs = rtog[rtog$metastatic_prostate_cancer == 1,]
 d = ggsurvplot(survfit(Surv((rtogs$survival_years - rtogs$metastatic_prostate_cancer_years), rtogs$survival) ~ rtogs$rx), data = rtogs, risk.table = F, ggtheme = theme_classic2(base_size=20), legend.title = "Treatment", legend.labs=c("With antiandrogen\ntherapy","Without antiandrogen\ntherapy ")) + ggtitle("KM Curve of Time between S to T\nFor Those who Experienced S") + ylab("Freedom from T Survival (Post S)") + xlab("Time (years)")
 d
 
 if(F){
  ggsave(file = "KM_S.jpeg", a$plot, width = 8, height = 8)
  ggsave(file = "KM_T.jpeg", b$plot, width = 8, height = 8)
  ggsave(file = "KM_ST.jpeg", d$plot, width = 8, height = 8)
 }
 



