# GLOFhazard

This repository contains the source codes to estimate GLOF peak discharges from the Himalayan lake-size distribution, 
and to estimate return periods of GLOF volumes and peak discharges in the entire Himalayas and seven subregions.
A mandatory prerequisite to run the scripts is to install R (https://www.r-project.org/) and 
RStudio (tps://rstudio.com) on your machine. The comments within the scripts provide further details on model dependencies
and usage of functions. 


## Scripts

### bayes_lm_piecewise_regression_stan.R

Bayesian piecewise regression model to learn peak discharge Qp from eta, the product of flood volume and the breach rate. Written in R.

### lm_piecewise_const.stan

This model calls a script to learn the posterior distributions of the model parameters, written in STAN.

### R_Script_Hazard_from_GLOFs_PNAS_supp.R

Script to estimate regional GLOF hazard from empirical GLOF rates and predicted GLOF volumes and discharges. Contains all commands to reproduce the figures from the associated manuscript.

### HDIofMCMC.R

Parts of our codes are based on the HDI function from Kruschke, J. K. (2014). Doing Bayesian Data Analysis: 
A Tutorial with R, JAGS, and Stan. 2nd Edition. Academic Press / Elsevier. We redistribute this function for ease of use, though the original function can be found here: https://sites.google.com/site/doingbayesiandataanalysis/software-installation


## Input data

... can be found here: add Zenodo


## References

Veh, G., Korup, O. & A. Walz: Hazard from Himalayan glacier lake outburst floods. In review.

## Contact

Georg Veh

Working group on natural hazards

University of Potsdam

georg.veh@uni-potsdam.de

https://www.uni-potsdam.de/en/umwelt/research/natural-hazards.html
