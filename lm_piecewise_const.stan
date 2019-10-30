//
// Bayesian Piecewise Linear Regression
// Code by Oliver Korup, 11 November 2018
// University of Potsdam, Germany (korup@uni-potsdam.de)
//
// This STAN code learns a piecewise regression model with one change point separating a linear from a constant segment. The model is to be called from within the accompanying file named bayes_lm_piecewise_regression_stan.R
//
// Make sure the last line of this file is blank.
//
// Disclaimer:
// By using these scripts you acknowledge that you do so at your own risk. The author of this script does not accept any liability for any inconvenience, damage or loss arising from the use of these scripts.
//
// Parts of this script are motivated by code by Max Joseph (mbjoseph/lm_segs.stan) at https://gist.github.com/mbjoseph/25d649e46602a419f9765638d5a2bfbc
//
data {
  int<lower=0> n;
  real x[n];
  real y[n];
  real<lower=0> nu;  // degrees of freedom for robust noise
  int<lower=0> n_tilde;  // number of data points to predict
  vector[n_tilde] x_tilde;  // new x_tilde to predict y_tilde
}
parameters {
  real alpha;  // regression intercept
  real beta;   // regression slope
  real<lower=0> sigma;  // data noise
  real chgpoint;   // bounds of chgpoint location need explicit truncation in hyperprior below
}
transformed parameters {
  vector[n] x2;   // indicator variable that tracks whether x_i is above or below the change point

  for (i in 1:n) {   // determine value of indicator function
    if (x[i] < chgpoint) {
      x2[i] = 0;
    } else {
      x2[i] = 1;
    }
  }
}
model {
  real mu[n];
  alpha ~ normal(0, 1);  // prior on intercept alpha
  beta ~ normal(0, 1);   // prior on slope beta
  sigma ~ cauchy(0, 2.5);  // prior on sigma (half-Cauchy due to lower=0 above)
  chgpoint ~ normal(-0.1, 0.1) T[-1, 0.5]; // truncated prior on change-point location informed by physical considerations outlined in O'Connor and Beebee (2009)

  for(i in 1:n){
    mu[i] = alpha + beta * x[i] - beta * (x[i] - chgpoint) * x2[i];
  }
  // Gaussian noise option
  // y ~ normal(mu, sigma);
  // Robust noise option
  y ~ student_t(nu, mu, sigma);
}
generated quantities {
  vector[n_tilde] y_tilde;   // generate predictions from the posterior
  for (j in 1:n_tilde)
    y_tilde[j] = student_t_rng(nu, alpha + beta * x_tilde[j] - beta * (x_tilde[j] - chgpoint) * x2[j], sigma);
}
