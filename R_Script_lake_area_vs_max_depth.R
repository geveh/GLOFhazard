####################################################################################################
#######             Estimate the maximum depth of Himalyan glacier lakes from          #############
#######                               glacier lake area                                #############
#######                                  by Georg Veh                                  #############
#######                                 30 July, 2019                                  #############
#######                   final version for PNAS 25 October, 2019                      #############
####################################################################################################

### IMPORTANT NOTE I:  MAKE SURE TO DROP ALL FILES TO RUN THE CODE IN ONE SINGLE FOLDER. USE     ###    
### setwd() to link R's working directory to exactly that folder. Adapt the command below to     ###
### link to the folder where you dropped the data                                                ###

setwd("D:/nrc_user/veh/LW_F/Hazard_from_GLOFs_PNAS")

### IMPORTANT NOTE II:  We here use rjags to generate posterior samples of lake depth,           ###    
### predicted from lake area. This script builds on code from John K. Kruschke, which makes use  ###
### of a Gibbs Sampler called JAGS. YOU MUST INSTALL THE STANDALONE PROGRAM 'JAGS' in order to   ###
### call the R package 'rjags'. Carefully follow the instructions how to install JAGS and rjags  ###
### given here: https://sites.google.com/site/doingbayesiandataanalysis/software-installation    ###

# If you have done so, call the package 'rjags'

require(rjags)

# Read the data.frame with empirically measured lake depths and volumes. Data and sources
# are given in SI Appendix, Table S1.

area.depth <- readRDS("area_vs_depth.rds")

# Setup for the robust regression model of lake area vs. lake depth. y (lake depth) is 
# normally distributed. See supplementary figure S1a for a graphical illustration of the
# individual model parameters. 

# We need to write a model string to disk, which is then called from JAGS

modelstring = "
model {
for( i in 1 : Ndata ) {
y[i] ~ dnorm( mu[i] , tau )
mu[i] <- beta0 + beta1 * x[i]
}
beta0 ~ dnorm( 0 , 1.0E-12 )
beta1 ~ dnorm( 0 , 1.0E-12 )
tau ~ dgamma( 0.001 , 0.001 )
}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

# We build the model in the log-log space, so that we need to log-transfrom
# x (lake area) and y (depth)

x <- log10(area.depth$area)
y <- log10(area.depth$depth)

# Number of data points

Ndata <- length(y)

# Specify data, as a list, needed for the JAGS sampler

dataList <- list(
  x = x,
  y = y,
  Ndata = Ndata)

# Initialize the Markov Chain Monte Carlo (MCMC):
# Use R's built-in least-squares regression to get plausible initial values:

lmInfo <- lm( dataList$y ~ dataList$x ) 
b0Init <- lmInfo$coef[1]
bInit  <- lmInfo$coef[2]
tauInit <- length(dataList$y) / sum(lmInfo$res^2)

initsList <- list(
  beta0 = b0Init,
  beta1 = bInit,
  tau = tauInit)


# Run the chains

parameters <- c("beta0" , "beta1" , "tau")  # The parameter(s) to be monitored.
adaptSteps <- 1000                          # Number of steps to "tune" the samplers.
burnInSteps <- 2000                         # Number of steps to "burn-in" the samplers.
nChains <- 3                                # Number of chains to run in parallel.
numSavedSteps <- 100000                     # Total number of steps in chains to save.
thinSteps <- 1                              # Number of steps to "thin" (1=keep every step).
nPerChain <- ceiling(( numSavedSteps * thinSteps ) / nChains) # Steps per chain.

# Create, initialize, and adapt the model:

jagsModel <- jags.model( "model.txt" , data = dataList , inits = initsList , 
                        n.chains = nChains , n.adapt = adaptSteps )

# Burn-in:

cat("Burning in the MCMC chain...\n")
update(jagsModel, n.iter=burnInSteps )

# The saved MCMC chain:

cat("Sampling final MCMC chain...\n")
codaSamples <- coda.samples(jagsModel, variable.names = parameters, 
                            n.iter = nPerChain, thin = thinSteps )

# Convert coda-object codaSamples to a matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]

mcmcChain <- as.matrix(codaSamples)
chainLength <- NROW(mcmcChain)

# For convenience later, append a column with tau converted to sigma:

sigma <- 1 / sqrt(mcmcChain[, "tau" ]) # Convert precision to SD
mcmcChain <- cbind(mcmcChain , sigma)

# We integrate out these parameters to simulate draws from the predictive 
# posterior distribution for unobserved inputs. We thus use these chains to
# predict lake depth for any lake in the Himalayas. Save this object to disk.

saveRDS(mcmcChain, "mcmcChain.rds")

##### DONE! ####

# This is the code for SI Appendix, Fig. S2

# Posterior prediction: Specify x values for which predicted y's are wanted:

extrapolationExtent <- 0.5*(range(x)[2]-range(x)[1])
lowX <- range(x)[1] - extrapolationExtent
highX <- range(x)[2] + extrapolationExtent
xPostPred <- seq(lowX, highX, length = 500)

# Define matrix for recording posterior predicted y values at each x value.
# One row per x value, with each row holding random predicted y values.

yPostPred <- matrix(0, nrow = length(xPostPred), ncol = chainLength)

# Define matrix for recording HDI limits of posterior predicted y values.

yHDIlim <- matrix(0, nrow = length(xPostPred), ncol = 2)

# Generate posterior predicted y values.
# This gets only one y value, at each x, for each step in the chain.

for (chainIdx in 1:chainLength) {
  yPostPred[,chainIdx] <- rnorm(length(xPostPred) ,
                               mean = (mcmcChain[chainIdx,"beta0"] + mcmcChain[chainIdx,"beta1"] * xPostPred) ,
                               sd = rep(sigma[chainIdx] , length(xPostPred)))
}

# We here use an umodified program from John K. Kruschke to caluculate the Highest density interval  
# (HDI) for predicted y values at each x location. The original code is distributed
# here: https://sites.google.com/site/doingbayesiandataanalysis/software-installation

source("HDIofMCMC.R")

for (xIdx in 1:length(xPostPred)) {
  yHDIlim[xIdx,] <- HDIofMCMC( yPostPred[xIdx,] )
}

# Display data with believable regression lines and posterior predictions.

pdf(file = "area_vs_maxdepth.pdf", width = 5.5, height = 3)

par(mar = c(4,4,0.5,0.5) , mgp = c(2.3,0.85,0), las = 1 )

# Plot the original data: Specify range of values (i.e. the figure margins) to be plotted.

xRang = max(x) - min(x)
yRang = max(y) - min(y)
limMult = 0.25
xLim= c( min(x) - limMult*xRang , max(x)+limMult*xRang )
yLim= c( min(y) - limMult*yRang , max(y)+limMult*yRang )

# We first create an empty plot without showing the data, so that we can add
# the original data on top of the believable regression lines.

plot(x, y, 
     cex  = 1.5,  lwd = 1.5,   col = "black", 
     xlim = xLim, ylim = yLim,
     lend = 1, cex.lab = 1, 
     xlab = expression(log[10]*"(lake area) [m²]") ,
     ylab = expression(log[10]*"(max. lake depth) [m]"), 
     type = "n")

# Impose a grid for better readability.

grid()

# Plot a smattering of 500 sampled, believable regression lines.

xComb <- seq(xLim[1] - 1, xLim[2] + 1, length = 201)

for ( i in seq(from=1, to= chainLength, by = 500)) {
  
  lines(xComb , 
        mcmcChain[i,"beta0"] + mcmcChain[i,"beta1"] * xComb , 
        col = rgb(0, 0, 0, 0.05))

  }


# Superimpose posterior predicted 95% HDIs:

lines(xPostPred, yHDIlim[,1], lwd = 2, col = "skyblue", lend = 1)
lines(xPostPred, yHDIlim[,2], lwd  =2, col = "skyblue", lend = 1)

# Now add the original data values

points(x , y , cex=1.25 , lwd=1.25 , col="black", pch = 21, bg = "skyblue")

# Close plot - done!

dev.off()