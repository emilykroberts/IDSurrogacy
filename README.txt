# Illness-Death-Surrogacy

-----------------------------------------------------------------------------------------------------
Manuscript details
----------------------------------------------------------------------------------------------------- 

This R code provides functions to produce simulation studies and results in "Surrogacy Validation for Time-to-Event Outcomes with Illness-Death Frailty Models".
Manuscript Authors: Emily K. Roberts, Michael R. Elliott, Jeremy M. G. Taylor
Code Inquiries: Please contact Emily K. Roberts (emily-roberts-1@uiowa.edu) with any questions about this code.

-----------------------------------------------------------------------------------------------------
Functionality, structure, files, and functions
-----------------------------------------------------------------------------------------------------

### Core functionality

The main functions in this repository are written to simulate data and run analysis (either of simulated or real data) for surrogate endpoint validation of two time-to-event outcomes in a clinical trial.

### Purpose of functions found within the 'R' folder

Data in the illness death format for a randomized trial can be generated using the sim_data.R file. This will provide two datasets ("dat0" and "dat1") that follow the three transitions for each treatment arm of the illness death model.

true_cep.R creates a Causal Effect Predictiveness (CEP) plot based on the true values of all model parameters (no estimation is done).

The run_sim.R function runs the simulations/MCMC and performs estimation and surrogacy validation based on our proposed methods.

The following files contain likelihood components and calculations for parameters in the proposed models: cLambda12_frailty_lk.R (cumulative hazard for the baseline to surrogate transition), cLambda13_frailty_lk.R (cumulative hazard for the baseline to death transition), cLambda23_frailty_lk.R (cumulative hazard for the surrogate to death transition), lambda12_frailty.R (hazard for the baseline to surrogate transition), like12_omega_i.R (frailty for baseline to surrogate), like12_shape.R (Weibull shape for baseline to surrogate), like13_omega_i.R (hazard for the baseline to death transition), like13_shape.R (Weibull shape for baseline to death), like23_shape.R (Weibull shape for surrogate to death), like23_theta.R (Weibull time-varying covariate). These files are called within the run_sim function to run the Bayesian estimation scheme. 

Based on the estimation in run_sim, the draws of parameters from the MCMC and other values and saved, which can be displayed and written to a file using final_results.R.
      
plot_results.R displays the estimated CEP curve from a dataset, and plot_traceplots.R is a diagnostic tool to examine the traceplots of the parameter draws as desired by the user.

## Overview of scripts to Reproduce Results

IDsims.R is a file to generate the main simulation results; .txt files are created with the posterior results (one file per dataset). The simulations in the manuscript are run in parallel on a high performance computer cluster (a .slurm file is provided for completeness), though this is not necessary for checking of individual file results. A single simulation can be run by setting the array_id variable to a fixed number and potentially changing the number of MCMC iterations (SIM) to a smaller quantity. Results from individual files can be aggregated within IDaggregate.R to create a table suitable for LaTeX or for output to a .csv file. Since the output from many (pre-run) simulation replications are already provided for reproducibility, aggregating 200 simulation replications per scenario can be run immediately using this file. 

IDdataexample.R provides code for the results and plots for the data example, and a sample data set to mimic features of the data is included in simulated_data.csv within the 'data' folder.

Other settings within the supplement are found in scripts called IDsimsS1.R, IDsimsS2_1_3.R, IDsimsS2_4_6.R, IDsimsS3.R, and IDsimsS4.R. Further instructions for replication are found below.

All results can be found within the 'results' file structure for the reproducibility check.

-----------------------------------------------------------------------------------------------------
More detailed steps to reproduce all figures and tables
-----------------------------------------------------------------------------------------------------

To reproduce figures and tables corresponding to the main and supplemental results, please follow the steps for each R script below (the entire script can be run as-is). Note that any simulation replication requires two major steps (found in separate scripts) in order to run the simulation many times and aggregate the results (described below):

1. Use IDsims.R to run primary simulations

One dataset is simulated and analyzed via MCMC if this code is run once. plot_results.R can be run on a given dataset to produce a CEP curve (Figure3.jpeg; note that this figure is only generated for the manuscript when array_id = 10). By running this file many times over different array_id's (recommended to be done using a computing cluster), multiple .txt files are created with the results (one file per simulated dataset). These files are named in the following manner: 'params', the scenario (to denote which data generating scenario was used), and the array_id (to denote the simulation replication). For example, the current file creates params,2,10.txt, and 8 scenarios x 200 simulations have been pre-run. These results from this main simulation will be aggregated in Step 3 to produce Table 2 (found in table2results.csv). 
Note: due to large file size and space constraints, all pre-run (ie intermediate) files are available on figshare here: https://doi.org/10.6084/m9.figshare.23933013.v2
In the online version, empty folders are provided to designate where these files belong and populate to.

2. Use supplemental simulation files to run other additional and sensitivity analysis simulations located in the supplemental material. Note there are five supplemental material scripts in total that can each be run independently.

The supplemental figures that do not require MCMC replications are provided in the file IDsimsS1.R. These files include the components of Supplemental Figure 1 (in 4 parts as FigureS1a.jpeg, FigureS1b.jpeg, FigureS1c.jpeg, FigureS1d.jpeg), Supplemental Figure 2 (in 4 parts as FigureS2a.jpeg, FigureS2b.jpeg, FigureS2c.jpeg, FigureS2d.jpeg), Supplemental Figure 3 (in 2 parts as FigureS3a.jpeg, FigureS3b.jpeg), and Supplemental Figure 4 (in 3 parts as FigureS4a.jpeg, FigureS4b.jpeg, FigureS4c.jpeg).

Variations of the main simulation results (i.e. results that require running simulations many times) in the supplementary material are located in individual scripts corresponding to the supplement numbering system. The code to run these results are found in IDsimsS2_1_3.R (tableS2results.csv: Table S2 rows 1-3 and FigureS5.jpeg; note that this figure is only generated for the manuscript when array_id = 12), IDsimsS2_4_6.R (tableS2results.csv: Table S2 rows 4-6), IDsimsS3.R (tableS3results.csv: Table S3), and IDsimsS4.R (tableS4results.csv: Table S4), respectively (note Supplemental Table 2 is comprised of two separate sets of analyses that are divided into rows 1-3 and 4-6 in Table S2). 

3. Aggregate simulation results to produce simulation summary tables

The files produced by steps 1 and 2 above may be aggregated in IDaggregate.R over many simulation replications. These results from this main simulation in step 1 will appear in the manuscript main text in Table 2 (found in table2results.csv). Similarly, the file will aggregate all replications from step 2 to produce Table S2, Table S3, and Table S4.

4. Use functions and IDdataexample.R to run the analysis via MCMC on a 'real' dataset. Note a pseudo dataset is provided within the 'data' folder.

Results from the true data analysis appear in the manuscript main text in Table 3 (found in table3results.csv). plot_results.R is run on the data to produce a CEP curve (found in Figure6.jpeg). Code is also given to produce the Kaplan-Meier plots found in Figure 5 (Figure5.jpeg) and cumulative incidence plots code (FigureS6.jpeg). Variations in the supplementary material for the dataset are also calculated in this file to produce plots for Metastasis-Free Survival (FigureS7.jpeg). 

-----------------------------------------------------------------------------------------------------
Simulated dataset notes
-----------------------------------------------------------------------------------------------------

Since we are unable to share the real data due to confidentiality concerns, all code may be run with a simulated dataset that we have provided.

Variables:
y12 time S
s12 censoring indicator for S  
y13 time T
s13 censoring indicator for T  
y23 time T-S (difference between the events occurring)
s23 censoring indicator for T occurring after S
trt treatment group
metastatic_prostate_cancer_years intermediate (surrogate) time to event
survival_years intermediate (surrogate) censor
survival true clinical time to event
metastatic_prostate_cancer true clinical censor

Note that the last four variables are clinical data, and the first six variables are transformed into gap times for compatibility with the proposed method.

-----------------------------------------------------------------------------------------------------
Version notes
-----------------------------------------------------------------------------------------------------

## System Information after loading required packages on a local computer

sessionInfo()R version 4.1.2 (2021-11-01)Platform: x86_64-apple-darwin17.0 (64-bit)Running under: macOS Monterey 12.6Matrix products: defaultLAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dyliblocale:[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8attached base packages:[1] stats     graphics  grDevices utils     datasets  methods   base     other attached packages: [1] here_1.0.1        frailtypack_3.5.0 doBy_4.6.13       survC1_1.0-3      MASS_7.3-58.1     boot_1.3-28       [7] frailtyEM_1.0.1   survminer_0.4.9   ggpubr_0.4.0      survival_3.4-0    ggplot2_3.4.0     xtable_1.8-4     [13] cowplot_1.1.1     latex2exp_0.9.5   truncnorm_1.0-8  loaded via a namespace (and not attached): [1] Rcpp_1.0.9           mvtnorm_1.1-3        lattice_0.20-45      tidyr_1.2.1          zoo_1.8-11           [6] assertthat_0.2.1     rprojroot_2.0.3      utf8_1.2.2           R6_2.5.1             backports_1.4.1     [11] rootSolve_1.8.2.3    pillar_1.8.1         rlang_1.0.6          rstudioapi_0.14      data.table_1.14.2   [16] car_3.1-0            Matrix_1.5-1         textshaping_0.3.6    labeling_0.4.2       splines_4.1.2       [21] statmod_1.4.37       stringr_1.4.1        munsell_0.5.0        broom_1.0.1          compiler_4.1.2      [26] numDeriv_2016.8-1.1  Deriv_4.1.3          xfun_0.33            systemfonts_1.0.4    pkgconfig_2.0.3     [31] microbenchmark_1.4.9 mgcv_1.8-40          tidyselect_1.2.0     tibble_3.1.8         gridExtra_2.3       [36] km.ci_0.5-6          fansi_1.0.3          dplyr_1.0.10         withr_2.5.0          grid_4.1.2          [41] nlme_3.1-160         gtable_0.3.1         lifecycle_1.0.3      DBI_1.1.3            magrittr_2.0.3      [46] KMsurv_0.1-5         scales_1.2.1         cli_3.4.1            stringi_1.7.8        carData_3.0-5       [51] farver_2.1.1         ggsignif_0.6.4       ragg_1.2.3           survMisc_0.5.6       generics_0.1.3      [56] vctrs_0.5.1          expint_0.1-7         tools_4.1.2          glue_1.6.2           mstate_0.3.2        [61] purrr_0.3.5          abind_1.4-5          colorspace_2.0-3     rstatix_0.7.0        knitr_1.40  