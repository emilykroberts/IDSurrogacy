# Illness-Death-Surrogacy

This R code provides functions to produce simulation studies as described in Surrogacy Validation for "Time-to-Event Outcomes with Illness-Death Frailty Models."

This package can be installed using the following R code:

devtools::install_github("emilykroberts/IDSurrogacy", build = TRUE, build_opts = c()) 

library(IDSurrogacy) 

### Resources

* [Ask a question/Open an issue: coming soon](https://github.com/emilykroberts) (GitHub issues for bug reports, feature requests)

### Core functionality

The main functions in this repository are written to simulate and run analysis for surrogacy validation of two time-to-event outcomes in a randomized clinical trial.

### Purpose of functions

Data in the illness death format for a randomized trial can be generated using the sim_data.R file. This will provide two datasets that follow the three transitions for each treatment arm described in the manuscript.

The following files contain the likelihood components or calculations for different parameters in the proposed models:

cLambda12_frailty_lk.R, cLambda13_frailty_lk.R, cLambda23_frailty_lk.R, lambda12_frailty.R, like12_omega_i .R, like12_omega1_i.R, like12_shape.R, like12_shape1.R, like13_omega_i.R, like13_omega1_i.R, like13_shape.R, like13_shape1.R, like23_shape.R, like23_shape1.R, like23_theta.R, like23_theta1.R

based on the parameter (shape of the Weibull, frailties, or cumulative hazard). 

true_cep.R creates a CEP plot based on the true values of all parameters and frailties (no estimation is done during this function).

The run_sim.R function runs the simulations and performs estimation and surrogacy validation based on our proposed methods.

Based on the estimation in run_sim.R, the draws of parameters from the MCMC and other values and saved, which can be shown in final_results.R       
      
plot_results.R displays the estimated CEP curve from a dataset, and plot_traceplots.R  is a diagnostic tool to examine the traceplots of the parameter draws.

Finally, IDsims.R is a file intended to set up parameter values and generate simulation results found in the main paper. Once this file is run, .txt files are created with the simulation results. These results can be read in using code also found in the file to create a table suitable for LaTeX. Code is also provided for the plots for the data example, and a sample data set to mimic features of the data is available.

### Contributing 

If you are interested in contributing to the development please open an issue to request.

### References and other literature

Roberts, E.K., Elliott, M.R., Taylor, J.M.G. Surrogacy Validation for Time-to-Event Outcomes with Illness-Death Frailty Models. (submitted; preprint https://arxiv.org/abs/2211.15826).
