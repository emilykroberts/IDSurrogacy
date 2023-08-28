library(truncnorm); library(latex2exp); library(cowplot); library(xtable); library(here); library(ggpubr)
library(ggplot2); library(survival); library(survminer); library(frailtypack); library(mvtnorm)

setwd("./R")

array_id = 1
# array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) # allows for running on cluster
if(is.na(array_id)) array_id = 1
# set various parameters to run data example
rhost = rhos = rhot = 0.5
SIM = 3000
specify = T
equalfrail = T
independent = T

proposalsd = 0.1
proposalsdfrail = 0.003
frailtysd = 0.4

tau_s = 5
tau_t = 8

write = T
plotwrite = T

# keep source files in a folder called R directory
for(i in 1:(length(list.files(pattern = "\\.R$")))){
  source(list.files(pattern = "\\.R$")[i]
  )
}

setwd('./data')
rtog = read.csv("simulated_data.csv")
setwd('./results')

set.seed(1 + array_id)

# form data for run_sim function for analysis
ST = (cbind(rtog$metastatic_prostate_cancer_years, rtog$survival_years - rtog$metastatic_prostate_cancer_years, 
            rtog$survival_years))
status = cbind(as.numeric(rtog$metastatic_prostate_cancer == 1), 
               as.numeric(rtog$metastatic_prostate_cancer == 1 & rtog$survival == 1), 
               as.numeric(rtog$survival == 1 & rtog$metastatic_prostate_cancer == 2))

ST = data.frame(cbind(ST, status, rtog$trt))
names(ST) = c("y12", "y23", "y13", "s12", "s23", "s13", "trt")
ST$y23[ST$y23 == 0] = NA

dat0 = ST[ST$trt == 0, ] 
dat1 = ST[ST$trt == 1, ]

# cumulative incidence of S (Figure S6)
ci = mstate::Cuminc(time = rtog$metastatic_prostate_cancer_years, status = rtog$metastatic_prostate_cancer, 
                    data = rtog, group = rtog$trt, failcodes = c(1, 2))
ci1 = ci[ci$group == 0, ]; ci2 = ci[ci$group == 1, ]

jpeg("FigureS6.jpeg", width = 650, height = 650) 
plot(c(0, ci1$time, max(ci1$time + 0.1)), c(0, ci1$CI.2, max(ci1$CI.2)), type = 's', ylim = c(0, 1), lty = 'dashed', 
     ylab = "Probability", xlab = "Time from Randomization", 
     main = "Estimates based on the cumulative incidence functions\nBlack dashed line indicates z = 1, solid red line indicates z = 0", 
     lwd = 2); lines(c(0, ci2$time, max(ci2$time) + 0.1), c(0, ci2$CI.2, max(ci2$CI.2)), type = 's', col = 'red')
dev.off()

# create and save KM curves in manuscript (Figure 5)
rtog$metastatic_prostate_cancer_cause = as.numeric(rtog$metastatic_prostate_cancer == 1)
rtog$metastatic_prostate_cancer_years_comp = as.numeric(rtog$metastatic_prostate_cancer_years != 0)
rtogs = rtog[rtog$metastatic_prostate_cancer == 1, ]

a = ggsurvplot(survfit(Surv(rtog$metastatic_prostate_cancer_years, rtog$metastatic_prostate_cancer_cause) ~ rtog$trt), 
               data = rtog, risk.table = F, ggtheme = theme_classic2(base_size = 20), legend.title = "Treatment", 
               legend.labs = c("z = 1", "z = 0")) + 
  ggtitle("KM Curve of Intermediate Outcome S\n(Individuals are Censored at T)") + ylab("Freedom from S") + xlab("Timee (years)")
b = ggsurvplot(survfit(Surv(rtog$survival_years, rtog$survival) ~ rtog$trt), data = rtog, risk.table = F, 
               ggtheme = theme_classic2(base_size = 20), legend.title = "Treatment", 
               legend.labs = c("z = 1", "z = 0")) + 
  ggtitle("KM Curve of True Outcome T") + ylab("Freedom from T") + xlab("Time (years)")
c = ggsurvplot(survfit(Surv((rtogs$survival_years - rtogs$metastatic_prostate_cancer_years), rtogs$survival) ~ rtogs$trt), 
               data = rtogs, risk.table = F, ggtheme = theme_classic2(base_size = 20), legend.title = "Treatment", 
               legend.labs = c("z = 1", "z = 0")) + 
  ggtitle("KM Curve of Time between S to T\nFor Those who Experienced S") + ylab("Freedom from T Survival (Post S)") + xlab("Timee (years)")

multi.page = ggarrange(a$plot, b$plot, c$plot,ncol = 2, labels = c("", ""), nrow = 2) 
ggexport(multi.page, filename = "Figure5.jpeg", width = 1100, height = 1100)

# run method on data
params_res = run_sim(SIM = SIM, array_id = array_id, rhos = rhos, rhot = rhot, frailtysd = frailtysd, 
                     params_list = true_params, tau_s = tau_s, tau_t = tau_t, dat0 = dat0, 
                     dat1 = dat1, proposalsdfrail = proposalsdfrail, proposalsd = proposalsd, MFS = F,
                     independent = independent, equalfrail = equalfrail, holdshape = F, holdtheta = F
)

final_results(params_matrix = params_res, write = write)

params = params_res$params
# create results table - produces nice Latex code for manuscript from single run of run_sim
xtable( rbind(rbind(round(colMeans(params[c(
  "int", "slope", "shape12_0", "shape13_0", "shape23_0", 
  "shape12_1", "shape13_1", "shape23_1", 
  "scale12_0", "scale13_0", "scale23_0", 
  "scale12_1", "scale13_1", "scale23_1", 
  "theta23_0", "theta23_1"
)], na.rm = T), 3)), 
rbind(round(apply(params[c(
  "int", "slope", "shape12_0", "shape13_0", "shape23_0", 
  "shape12_1", "shape13_1", "shape23_1", 
  "scale12_0", "scale13_0", "scale23_0", 
  "scale12_1", "scale13_1", "scale23_1", 
  "theta23_0", "theta23_1"
)], 2, FUN = sd, na.rm = T), 3))))

# corresponds to Table 3
tab3 = (rbind(rbind(round(colMeans(params[c(
  "int", "slope", "shape12_0", "shape13_0", "shape23_0", 
  "shape12_1", "shape13_1", "shape23_1", 
  "scale12_0", "scale13_0", "scale23_0", 
  "scale12_1", "scale13_1", "scale23_1", 
  "theta23_0", "theta23_1"
)], na.rm = T), 3)), 
rbind(round(apply(params[c(
  "int", "slope", "shape12_0", "shape13_0", "shape23_0", 
  "shape12_1", "shape13_1", "shape23_1", 
  "scale12_0", "scale13_0", "scale23_0", 
  "scale12_1", "scale13_1", "scale23_1", 
  "theta23_0", "theta23_1"
)], 2, FUN = sd, na.rm = T), 3))))
write.csv(tb3, file = "table3results.csv")

# save results - named Figure6.jpeg (Figure 6) in results for manuscript
plot_results(params_matrix = params_res, fignum = 6, write = plotwrite)

# set MFS to true to use MFS as the surrogate endpoint
params_res = run_sim(SIM = SIM, array_id = array_id, rhos = rhos, rhot = rhot, frailtysd = frailtysd, 
                     params_list = true_params, tau_s = tau_s, tau_t = tau_t, dat0 = dat0, 
                     dat1 = dat1, proposalsdfrail = proposalsdfrail, proposalsd = proposalsd, MFS = T,
                     independent = independent, equalfrail = equalfrail, holdshape = F, holdtheta = F
)

params = params_res$params
plot_results(params_matrix = params_res, fignum = 'S7', write = plotwrite) # file FigureS7.jpeg in the supplement

