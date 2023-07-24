# Illness-Death-Surrogacy

This R code provides functions to produce simulation studies as described in Surrogacy Validation for "Time-to-Event Outcomes with Illness-Death Frailty Models."

This package can be installed using the following R code:

devtools::install_github("emilykroberts/IDSurrogacy", build = TRUE, build_opts = c()) 

library(IDSurrogacy) 

### Resources

* [Ask a question/Open an issue: coming soon](https://github.com/emilykroberts) (GitHub issues for bug reports, feature requests)


-----------------------------------------------------------------------------------------------------
Functionality, structure, files, and functions
-----------------------------------------------------------------------------------------------------

### Core functionality

The main functions in this repository are written to simulate data and run analysis (either of simulated or real data) for surrogate endpoint validation of two time-to-event outcomes in a randomized clinical trial.

### Purpose of functions

First, data in the illness death format for a randomized trial can be generated using the sim_data.R file. This will provide two datasets ("dat0" and "dat1") that follow the three transitions for each treatment arm described in the manuscript.

The following files contain the likelihood components or calculations for different parameters in the proposed models based on the Bayesian framework:

cLambda12_frailty_lk.R, cLambda13_frailty_lk.R, cLambda23_frailty_lk.R,  lambda12_frailty.R, like12_omega_i.R, like12_shape.R, like13_omega_i.R, like13_shape.R, like23_shape.R, like23_theta.R

based on the parameter (shape of the cumulative hazard, frailty terms, Weibull shape parameters, Weibull scale parameters, and Weibull time-dependent coefficients, respectively as described in the manuscript). These files are needed within the run_sim.R function to run the Bayesian estimation scheme. 

true_cep.R creates a Causal Effect Predictiveness (CEP) plot based on the true values of all parameters and frailties (no estimation is done).

The run_sim.R function runs the simulations/MCMC and performs estimation and surrogacy validation based on our proposed methods.

Based on the estimation in run_sim.R, the draws of parameters from the MCMC and other values and saved, which can be shown and written to a file using final_results.R.
      
plot_results.R displays the estimated CEP curve from a dataset, and plot_traceplots.R  is a diagnostic tool to examine the traceplots of the parameter draws as desired by the user.

IDsims.R is a file to set up parameter values and generate the simulation results found in the main paper. Once this file is run, .txt files are created with the results (one file per dataset over "SIM" iterations of the MCMC). These results can be read in using code also within the file to create a table suitable for LaTeX or for output to a .csv file. The simulations in the manuscript are run in parallel on a high performance computer cluster. A slurm file is provided to run many iterations in parallel, though a single simulation can be run by setting the array_id variable to a fixed number and changing the number of MCMC iterations (SIM).

Within IDdataexample.R, code is provided for the results and plots for the data example, and a sample data set to mimic features of the data is included in simulated_data.csv within the data folder.

### Contributing 

If you are interested in contributing to the development please open an issue to request.

### References and other literature

Roberts, E.K., Elliott, M.R., Taylor, J.M.G. Surrogacy Validation for Time-to-Event Outcomes with Illness-Death Frailty Models. (submitted; preprint https://arxiv.org/abs/2211.15826).
