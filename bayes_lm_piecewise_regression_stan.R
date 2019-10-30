#####################################################################################
#
# Bayesian Robust Piecewise Linear Regression Model using STAN
#
# Code by Oliver Korup, 11 Nov 2018
#
# University of Potsdam, Germany (korup@uni-potsdam.de)
#
# Parts of this script and the associated STAN file are
# motivated by code by Max Joseph (mbjoseph/lm_segs.stan) at
# https://gist.github.com/mbjoseph/25d649e46602a419f9765638d5a2bfbc
#
# This R-script requires additional libraries to run the STAN
# probabilistic programming language. This script also requires
# three separate files to be located in the working directory:
#
# (1) A *.csv file named
# o_connor-and-beebee-2009.csv
#
# (2) A *.stan file named
# lm_piecewise_const.stan
#
# (3) A *.r file named
# HDIofMCMC.R
#
# 
# Disclaimer: 
# By using these scripts you acknowledge that you do so at your own risk.
# The author of this script does not accept any liability
# for any inconvenience, damage or loss arising from the use of these
# scripts.
#
####################################################################################

# PRELIMS---------------------------------------------------------------------------

# Set working directory:
# Change this to suit a working directory on your machine:
# Make sure all other necessary files listed below are within the same
# working directory
setwd("~/MAIN/RC/RSTAN")

# Load additional R libraries;
# use require() where needed
library(foreign)
library(rstan)
library(gridExtra)
library(bayesplot)
library(LaplacesDemon)

# Select multiple cores as suggested in the STAN manual
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# LOAD DATA--------------------------------------------------------------------------
# Load data on natural dam breaks taken from O'Connor and Beebee (2009):
# J. E. O’Connor, R. A. Beebee, “Floods from natural rock-material dams” 
# in Megaflooding on Earth and Mars, D. M. Burr, P. A. Carling, V. R. Baker, Eds.
# (Cambridge University Press, 2009), pp. 128–171.
dams <- read.csv2("o_connor-and-beebee-2009.csv", header = TRUE, dec = ",")

# Compute dimensionless discharge following O'Connor and Beebee (2009)
Qp_star <- dams$Qp / (9.81 ^ 0.5 * dams$Breach_d ^ 2.5)

# Compute eta following O'Connor and Beebee (2009)
eta <- (dams$V0 / dams$Breach_d ^ 3) * dams$Breach_Rate / 
  (9.81 ^ 0.5 * dams$Breach_d ^ 0.5)

# Create dataset of log10-transformed variables of Qp_star and eta
d <- as.data.frame(cbind(log10(eta), log10(Qp_star)))
colnames(d) <- c("eta", "Qp_star")

# Trap for missing data
d <- na.omit(d)


# Specify locations for prediction for STAN---------------------------------------------
# Set number of predicted data (= not included in the original data)
n_pred <- 60

# Set equal-spaced values of eta for prediction
x_pred <- seq(-3, 2, length.out = n_pred)


# Specify degrees of freedom for t-distributed noise
# (this parameter can also be learned by the model if necessary)
nu <- 10

# Standardise log10-transformed data to improve sampling in STAN
# and create a list of all necessary data and parameters for the
# regression model
stan_d <- list(y = c(scale(d$Qp_star)), 
               x = c(scale(d$eta)),
               n = nrow(d),
               nu = nu,
               n_tilde = n_pred,
               x_tilde = (x_pred - mean(d$eta)) / sd(d$eta)
               )


# Model for piecewise linear regression where the second interval is constant------
# Initial model for running a quick test
minit <- stan("lm_piecewise_const.stan",   # separate STAN file
              data = stan_d, 
              chains = 1, 
              iter = 2,
              control = list(adapt_delta = 0.99))

# CALL AND RUN MAIN MODEL
modfit <- stan(fit = minit , 
             iter = 2000,    # number of iterations
             data = stan_d, 
             chains = 5,     # number of chains
             cores = 4)

# Display tabulated posterior information
modfit

# Run (and check) diagnostics of sampling; re-run model if any of these diagnostics
# throws an error or a warning
check_divergences(modfit)
check_energy(modfit)
check_hmc_diagnostics(modfit)
check_treedepth(modfit)

# Extract posterior information
post <- rstan::extract(modfit)

saveRDS(post, "GLOF-dimless-posterior.rds")


# PLOT ON ORIGINAL SCALE ----------------------------------------------------------

# Re-transform standardised parameters to original scale:
# Change point
chpt <- post$chgpoint * sd(d$eta) + mean(d$eta)

# Intercept before change point
alpha_1 <- post$alpha * sd(d$Qp_star) + mean(d$Qp_star) + 
  post$beta * sd(d$Qp_star) / sd(d$eta) * 
  (-mean(d$eta))

# Intercept after change point
alpha_2 <- post$alpha * sd(d$Qp_star) + mean(d$Qp_star) + 
  post$beta * sd(d$Qp_star) / sd(d$eta) * 
  (chpt - mean(d$eta))

# Slope before change point
beta_1 <- post$beta * sd(d$Qp_star) / sd(d$eta) * 1

# Slope after change point (= 0)
beta_2 <- post$beta * sd(d$Qp_star) / sd(d$eta) * 0

# Plot original data (cf. Fig. 8.7 in O'Connor and Beebee, 2009)
par(mar = c(5, 5, 1, 1), mfcol = c(1, 1))
plot(d$eta, d$Qp_star, 
     xlab = expression(eta), 
     ylab = expression("Q"[p]^"*"),
     cex = 1.2, cex.lab = 1.2,
     las = 1,
     lab = c(5, 4, 7))
grid()

# Add credible piecewise regression lines, colour-coded by segment
for (i in seq_along(post$lp__)) {
  # Before change point
  curve(alpha_1[i] + beta_1[i] * x - 
          beta_1[i] * (x - chpt[i]) * ifelse(x < chpt[i], 0, 1),
        from = 1.2 * min(d$eta), to = chpt[i],
        add = TRUE,
        col = rgb(77/255, 106/255, 1, 0.05))
  # After change point
  curve(alpha_2[i] + beta_2[i] * x - 
          beta_2[i] * (x - chpt[i]) * ifelse(x < chpt[i], 0, 1),
        from = chpt[i], to = max(d$eta) * 1.2,
        add = TRUE,
        col = rgb(0, 213/255 , 1, 0.05))
}

# Add posterior estimate of change-point location
rug(chpt, col = "darkgrey")
polygon(density(chpt)$x, 0.5 * density(chpt)$y + min(d$Qp_star), 
        col = rgb(0.1, 0.1, 0.1, 0.25), border = NA)


# ADD POSTERIOR PREDICTIVE ---------------------------------------------------------

# Set new prediction points x_tilde in values of log10(eta)
x_tilde <- seq(1.1 * min(d$eta), 
               1.1 * max(d$eta), 
               length.out = 1000)

# Standardise x_tilde by empirical mean and s.d. of data in
# O'Connor and Beebee (2009)
x_tilde_stan <- (x_tilde - mean(d$eta)) / sd(d$eta)

# Set number of samples from predictive posterior
smp <- 5000

# Create output container for x_tilde, median(y_tilde), and 95% HDI of y_tilde
post_pred <- data.frame(cbind(x_tilde,
                              med_y_tilde = NA,
                              lwrHDI_y_tilde = NA,
                              uprHDI_y_tilde = NA))


# Sample from standard t-distribution, which models the data noise of each pointwise
# estimate
y_tilde_stan <- apply(t(x_tilde_stan), 2, function(x) rt(smp, df = nu) * post$sigma +
                        post$alpha + post$beta * x - post$beta *
                        (x - post$chgpoint) * ifelse(x < post$chgpoint, 0, 1))


# Re-transform standardised variables to original scale:
# Each column of yy contains the posterior predictive of log10(Qp_star)
yy <- y_tilde_stan * sd(d$Qp_star) + mean(d$Qp_star)

# Load HDI function
# Code based on the HDI function from
# Kruschke, J. K. (2014). Doing Bayesian Data Analysis: 
# A Tutorial with R, JAGS, and Stan. 2nd Edition. Academic Press / Elsevier.
source("HDIofMCMC.R")

# Collect point estimates of posterior of y_tilde measured in values of log10(Qp_star)
post_pred$med_y_tilde <- apply(yy, 2, median)    # median predictive posterior
HDI <- apply(yy, 2, function(x) HDIofMCMC(x, credMass = 0.95))  # 95% HDI
post_pred$lwrHDI_y_tilde <- HDI[1, ]
post_pred$uprHDI_y_tilde <- HDI[2, ]

# Add pointwise estimates to plot
lines(x_tilde, post_pred$med_y_tilde)
lines(x_tilde, post_pred$lwrHDI_y_tilde, lty = 2, col = "darkgrey")
lines(x_tilde, post_pred$uprHDI_y_tilde, lty = 2, col = "darkgrey")

# Replot original data points
points(d$eta, d$Qp_star, pch = 21, cex = 1.2, col = "white", bg = "black")

# FIN----------------------------------------------------------------------------------