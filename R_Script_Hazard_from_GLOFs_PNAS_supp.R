####################################################################################################
#######             Hazard from Himalayan moraine-dammed glacial lakes                 #############
#######                                                                                #############
#######                                  by Georg Veh                                  #############
#######                                 10 July, 2019                                  #############
#######                     more comments added 09 October, 2019                       #############
#######                   final version for PNAS 30 October, 2019                      #############
####################################################################################################

### IMPORTANT NOTE I: PARTS OF THIS CODE CONSUME RAM WELL-BEYOND 100 GB (depending on the        ###    
### number of cores you use). DATA ARE WRITTEN TO DISK AND MAY CONSUME SEVERAL TENS OF GBs.      ###
### MAKE SURE YOU HAVE THE NECESSARY SYSTEM REQUIREMENTS.                                        ###

### IMPORTANT NOTE II: By using these scripts you acknowledge that you do so at your own risk.   ###                     
### The author of this script does not accept any liability for any inconvenience, damage or     ###
### loss arising from the use of these scripts.                                                  ###


### IMPORTANT NOTE III: MAKE SURE THAT YOU DROP ALL FILES TO RUN THE CODE IN ONE SINGLE FOLDER.  ###    
### Use setwd() to link R's working directory to exactly that folder. Adapt the command below    ###
### to link to the folder where you dropped the data                                             ###

setwd("D:/nrc_user/veh/LW_F/Hazard_from_GLOFs_PNAS")


### Load packages to your session. If you have not installed some of these packages yet, write
### install.packages("packagename").

require(sp)
require(fitdistrplus)
require(rgeos)
require(pbapply)
require(stringr)
require(parallel)
require(doParallel)
require(abind)
require(extRemes)
require(sfsmisc)

### Load data

# We use a data-optimized Oblique mercator projection, which minimizes high distortions
# given the oblique stretch of the Himalayas. All shapefiles that we use in the following 
# come in this projection

# We need the historic GLOF count to obtain regional GLOF rates. In the past three decades 38
# GLOFs have been previously reported (glofs.3dec.omerc) or newly detected (landsat.glofs),
# as reported in G. Veh, O. Korup, S. von Specht, S. Roessner, A. Walz, Unchanged frequency of  
# moraine-dammed glacial lake outburst floods in the Himalaya. Nat. Clim. Change 9, 379-383 (2019).

# Load previously unreported GLOFs, detected from Landsat images.

landsat.glofs.omerc <- readRDS("landsat_omerc_glofs.rds")

# Load reported GLOFs.

glofs.3dec.omerc <- readRDS("glofs_3dec_omerc.rds")

# Load the lake inventory of Maharjan et al. (2018). Here we use a subset of moraine-dammed lakes,
# reprojected to the oblique mercator projection. Further details on this inventory are given here:
# Maharjan, S.B., Mool, P.K., Lizong, W., Xiao, G., Shrestha, F., Shrestha, R.B., Khanal, N.R.,
# Bajracharya, S.R., Joshi, S., Shai, S., Baral, P. (2018). The status of glacial lakes in the Hindu
# Kush Himalaya. ICIMOD Research Report 2018/1. Kathmandu: ICIMOD
# http://lib.icimod.org/record/33736/files/icimodRR1-018.pdf

present.lakes.omerc <- readRDS("present_lakes_omerc.rds")

# Load the study region, including the seven subregions Hindu Kush, Karakoram, Western Himalaya,
# Central Himalaya, Eastern Himalaya, Nyainqentanglha, and Hengduan Shan

aoi.buffer.int.omerc <- readRDS("aoi_buffer_int_omerc.rds")

# In the following, we mostly define the tails of our distributions by the Highest density interval (HDI).
# We here use an umodified program from John K. Kruschke to caluculate HDI. The original code is distributed
# here: https://sites.google.com/site/doingbayesiandataanalysis/software-installation

source("HDIofMCMC.R")

####################################################################################################
###############         Estimate flood volumes and peak discharge per lake        ##################
####################################################################################################

# Need to import objects: mcmcChain, post, d

# mcmcChain is the posterior predictive of lake depths vs. area, see the script
# "R_Script_lake_area_vs_max_depth.R"

mcmcChain <- readRDS("mcmcChain.rds")

# post is a dataframe containing the posterior predictive of the piece-wise
# regression model of Qp vs. eta.

post <- readRDS("GLOF-dimless-posterior.rds")

# d is the original data from O'Connor and Beebee (2009) where we based this model on.

d <- readRDS("beebee-oconnor-2009-data.rds")


# Fit a log-normal distribution to empirical breach rates from the table on natural
# dam failures in O'Connor and Beebee (2009), pp.148-162.
# These are all reported natural dams breaks.

published.k <- c(3.7 * 10^-4, 1.4  * 10^-3, 1.9 * 10^-3, 1.11 * 10^-2,
                 1.5 * 10^-3, 2.54 * 10^-2, 6   * 10^-4, 2.5  * 10^-2,
                 3.4 * 10^-3, 1.13 * 10^-2, 1.7 * 10^-3, 2.2  * 10^-3,
                 2.5 * 10^-3, 3.2  * 10^-3, 1.1 * 10^-1, 3.8  * 10^-3,
                 4.3 * 10^-3, 3.17 * 10^-2, 3.7 * 10^-4)


# Fit lognormal distribution to k values.

k.lnorm.fit <- fitdist(published.k, distr = "lnorm", method = "mle")

# Generate 100,000 realizations from this distribution, which we use for 
# calculating eta

all.k <- rlnorm(n = 100000, 
                meanlog = k.lnorm.fit$estimate[1], 
                sdlog = k.lnorm.fit$estimate[2])

# We allow for twice the range of documented breach rates to account for the effects of 
# failure mechanisms that may have remained unobserved. 

all.k <- all.k[all.k < max(published.k)*2]

# Define gravitational acceleration g

g <- 9.80665

# Calculate lake area from OMERC-projected lakes

present.lakes.omerc$Area <- gArea(present.lakes.omerc, byid = T)

# Delete lakes < 0.01

lakes.sub <- present.lakes.omerc[present.lakes.omerc$Area >= 10000, ]

# Generate container for posterior predictives of maximum lake depths.
# This matrix has a size of m lakes and n and posterior predicted values.

### We iterate over all lakes and estimate their depth, flood volume and peak discharge
### We sample 100 maximum lake depths from the posterior predictives of lake area vs. maximum depth.
### We discretised lake depth into 1 to 100 % breaches and converted these to flood volumes, assuming
### a circular lake area. For each flood volume (per lake level drop, per 100 estimated lake depth, for all lakes)
### we assumed 100 breach rates. We predict Qp from the posterior predictive of 100 models of Qp versus eta.  
### We do this sampling scheme 10 times for computational efficiency.
### See SI Appendix, Fig. S1 b) and c) for a graphical notation of our workflow.


for (run in 1:10) {
  
  yPostPred = matrix(0 , nrow = nrow(lakes.sub), ncol = 100)
  
  # Define a matrix for recording HDI limits of posterior predicted lake depths.
  
  yHDIlim = matrix(0 , nrow = nrow(lakes.sub) , ncol = 2)
  
  # Generate posterior predicted lake depths.
  # This gets one depth value, at each lake area, for each step in the chain.
  # We sample 1,000 (i.e. 10 x 100) value pairs of beta0, beta1 and sigma from the chains.
  # We assumed the mean of the robust regession to be normally distributed.
  # We learned the model in log-log-space, so that lake areas must be log-transformed
  # when making predictions of lake depth. The latter are thus also in log-space.
  
  chainSamp <- sample(1: nrow(mcmcChain), size = 100, replace = F)
  
  for ( chainIdx in 1:length(chainSamp)) {
    
    yPostPred[ , chainIdx] <- rnorm(length(lakes.sub$Area) ,
                                    mean = (mcmcChain[chainSamp[chainIdx], "beta0"] + 
                                               mcmcChain[chainSamp[chainIdx], "beta1"] * 
                                               log10(lakes.sub$Area)) ,
                                    sd = rep(mcmcChain[chainSamp[chainIdx], "sigma"], length(lakes.sub$Area)) )
  }
  
  # For each lake, generate the 95 % HDI of credible depth values.
    
  for ( xIdx in 1:nrow(lakes.sub)) {
    
    yHDIlim[xIdx,] <- HDIofMCMC(yPostPred[xIdx,] )
  
    }
  
  # Delete samples outside the 95 % HDI.
    
  depthPostSamp <- t(sapply(1:nrow(yPostPred), function (i) {
    
    return(yPostPred[i, yPostPred[i, ] > yHDIlim[i, 1] & yPostPred[i, ] < yHDIlim[i, 2]])
    
  }))
  
  # Depth values are predicted in log10. Convert to original space and generate a list
  # of depth values for each lake.
    
  depth.total <- 10^depthPostSamp
  depth.total.list <- lapply(seq_len(nrow(depth.total)), function(i) depth.total[i, ])
  
  # Assume that all breach depths from 0% to 100% are possible. 
    
  perc.breach <- seq(0.01, 1, by = 0.01)
  
  # Generate intervals of lake depth according to 0-100% breach. This corresponds
  # to the height between lake bed and the breach of the dam.
    
  lake.depth <- lapply(depth.total.list, function(x) outer(x, perc.breach))  
  
  # Breach depth is the inverse of lake depth, measured from the lake surface to 
  # the sill.
  
  breach.depth <- lapply(lake.depth, function (x) {
    
    # the entry at 100 is the brim-full lake
    apply(x, 2, function(m) x[, 100] - m)
    
  })
  
  # Assume that each lake has a circular area. Calculate the radius as r = sqrt(a / pi)
  
  r <- sqrt(lakes.sub$Area/pi)
  r.list <- lapply(r, function (x) x)
  
  # Estimate the volume of each lake, assuming a bowl-shaped lake bathymetry. 
  # This corresponds to a half ellipsoid with radius r and depth c.
  
  vol.tot <- ((4/3) * pi * depth.total * r^2)   / 2
  vol.tot.list <- lapply(seq_len(nrow(vol.tot )), function(i) vol.tot [i, ])
  
  # Calculate the lake volume, measured from the bottom until the breach depth.
  # This is the cap of the ellipsoid, i.e. the water that remains after the breach.
  
  vol.lake <- mapply(function (r, h, c) {(pi * r^2) * (((2*h)/3) - c + ((c^3)/ (3 * h^2)))}, 
                     r.list, depth.total.list, breach.depth, SIMPLIFY = F)
  
  # Subtract the total volume from the remaining lake volume after the outburst
  # to obtain the flood volume V0.
  
  released.vol  <- mapply(function (x, y) { 
    
    out <- apply(x, MARGIN = 2, function (h) {y - h})
    
    return(out[, -100])
    
  }, vol.lake, vol.tot.list, SIMPLIFY = F)
  
  # Generate an ID per run, which we use to save the simulations to a specified folder.
  
  run.id <- str_pad(run, width = 2, side = "left", pad = 0)
  
  if (run == 1) {dir.create("Vol_present")}
  
  # Create a directory where to save the predicted flood volumes.
  
  vol.dir <- paste0("Vol_present/Run_", run.id)
  dir.create(vol.dir)
  
  # Write the predicted flood volumes to disk.
  
  saveRDS(released.vol, paste0(vol.dir, "/Run_", run.id, ".rds"))
  
  # Setup the cluster to predict peak discharges. We use here 50 cores. 
  # Use detectCores() to obtain the number of cores that you can use on your machine, 
  # but leave some cores free for other processes running on your machine at the same time.
  
  cl <- makeCluster(50)
  registerDoParallel(cl, cores = 50)
  
  # Export the necessary data to the cluster.
  
  clusterExport(cl, c("all.k", "breach.depth", "released.vol", 
                      "d", "g", "run", "post" ))
                         
  # We iterate over all lakes (here the object 'breach.depth' is a list, which each entry
  # representing one lake) and predict peak discharge with samples of k and V0.
  # We use pbapply package, which creates a nice progress bar for the console.
  # Create a directory, where predicted peak discharges will be stored
  
  dir.create("Qp_present")
  
  pblapply(1:length(breach.depth), cl = cl, function (xxx) {
  
    # Sample 100 breach rates
     
    k <- sample(all.k, 100)
    
    h <- breach.depth[[xxx]][, -100]
    
    # Calculate k_star with the formula given in O'Connor and Beebee (2009), p.141
    
    k_star  <- outer(X = h, Y = k , FUN = function (x, y) y / ((g^0.5) * (x^0.5)))
    
    # We extract the corresponding flood volume for the breach depth.
    
    vol <- released.vol[[xxx]]
    
    # Calculate V_star with the formula given in O'Connor and Beebee (2009), p.141
    
    V_star <- vol/ h^3
    
    # Set new prediction points x_tilde, measured as log10(eta)
    # Note 'eta' is log10(eta) from the original data. Multiply k_star with V_star
    # to obtain eta.
    
    eta.list <- list()
    
    for(i in 1:(dim(k_star)[3])) {
      
      eta.list[[i]] <- log10(k_star[, , i] * V_star)
      
    }
    
    eta <- do.call(abind::abind, list(eta.list, rev.along = 0, use.dnns = TRUE))
    
    # For convenience, we name eta x_tilde and standardise it.
   
    x_tilde <- eta
    x_tilde_stan <- (x_tilde - mean(d$eta)) / sd(d$eta) 
    
    # Set number of samples from predictive posterior
    
    smp <- 100
    smp.idx <- sample(1:5000, 100, replace = F)
    
    # Sample from a t-distribution with posterior location and scale.
    # post.igma is the variance in the t-distribution.
    # post.alpha is the intercept in the piecewise regression model.
    # post.beta is the slope in the piecewise regression model.
    # And post.cutpoint is the learned location of the breakpoint in the model,
    # linking the two segments
    
    post.sigma <- post$sigma[smp.idx]
    post.alpha <- post$alpha[smp.idx]
    post.beta  <- post$beta[smp.idx]
    post.cutpoint <- post$cutpoint[smp.idx]
    
    Qp <- sapply(1:99, function (i) {
      
      eta.over.depth <- x_tilde_stan[, i, ]
      
      # Set fixed degrees of freedom.
      
      nu <- 10
      
      # Notation of the piecewise model. The ifelse-statement is the indicator function
      # to specify, whether the predicted point is beyond the breakpoint.
      
      y_tilde_stan <- apply(t(as.numeric(eta.over.depth)), 2, function(m) rt(smp, df = nu) * post.sigma +
                              post.alpha + post.beta * m - post.beta *
                              (m - post.cutpoint) * ifelse(m < post.cutpoint, 0, 1))
      
      # Re-transform to original scale:
      # Each column of yy contains the posterior predictive of log10(Qp_star)
      # Note that 'Qp_star' is log10(Qp_star) from the original data.
      
      yy <- (y_tilde_stan * sd(d$Qp_star)) + mean(d$Qp_star)
      
      th <- h[, i]
      
      # Re-transform Qp_star to the original scale 
      
      Qp <- apply(yy, 1, function (z) (10^ z * g^0.5 * rep(th, 100)^2.5))
      
    })
    
    # We sample 1,000 values from the predicted posterior. We could instead report all predicted
    # posterior values, which would come with the price of larger memory consumption.
    
    o <- apply(Qp, 2, sample, 1000, replace = F)
    
    lake.id <- stringr::str_pad(xxx, width = 4, side = "left", pad = 0)
    
    # As before, create a folder write the predictions to disk.
    
    lake.dir <- paste0("Qp_present/Lake_", lake.id)
    dir.create(lake.dir)
    saveRDS(o, paste0(lake.dir, "/", "Lake_", lake.id, "_run", run, ".rds"))
    
    # Clear workspace.
    
    gc()
  
  })
  
  # Stop the clusters.
  
  stopCluster(cl)
  
  # Report the progress on the console.
  
  message(paste0("run ", run, " completed at ", Sys.time()))

} # end run


####################################################################################################
###############         Estimate return periods of FLOOD VOLUMES for              ##################
###############          the entire Himalayas and seven subregions                ##################
####################################################################################################

# Read the predicted flood volumes back to memory. List all directories where we had dropped
# the files before.

Vol.dirs <- list.files(path = "Vol_present", full.names = T, recursive = T)

all.Vol <- lapply(Vol.dirs, function (xxx) {
  
  o <- readRDS(xxx)
  o <- lapply(o, as.numeric)
  o.b <- do.call(rbind, o)
})

all.Vol2 <- do.call(cbind, all.Vol)

# Truncate the distribution of predicted flood volumes to the 95% HDI for each lake.
# Set up the cluster and export the flood volumes to the cores.

cl <- makeCluster(10)
registerDoParallel(cl, cores = 10)
clusterExport(cl, list("all.Vol2", "HDIofMCMC" ))

all.Vol3 <- pbapply(all.Vol2, cl = cl, MARGIN = 1,  function (x) {
  
  # The function HDIofMCMC() returns the lower and upper bound of the HDI,
  # which we use to clip the vector of predicted flood volumes.
  
  hdi <- HDIofMCMC(x)
  lake.hdi <- x[x >= hdi[1] & x <= hdi[2]]
  
})

stopCluster(cl)

# transpose the matrix to obtain the original data structure.

all.Vol <- t(all.Vol3)

# Delete the intermediate steps and clear cache.

rm(all.Vol2, all.Vol3)
gc()



####### Return periods of FLOOD VOLUMES for the ENTIRE HIMALAYAS ########


# Calculate the posterior GLOF rate in the entire Himalayas: 
# The formula is (#GLOFs + 1) / (length of the study period + 1)

# Number of years

sample_number <- 30

# Calculate the posterior GLOF rate in the entire Himalayas as
# (Number of GLOFs + 1) / (Length of study period + 1)

lambda_mu_post = (nrow(glofs.3dec.omerc) + nrow(landsat.glofs.omerc) + 1) / (sample_number + 1) 


# Generate a list with 1,000 time series, each 10,000 years long.
# The number of samples per year is assumed to follow a Poisson distribution with
# the rate being the posterior mean annual GLOF rate ('lambda_mu_post').

# Specify number of years.

years <- 10000

Vol.present.container <- list()

for (i in 1: 1000) {
  
  Vol.present.container[[i]] <- sample(x = all.Vol,
                                       size = sum(rpois(years, lambda_mu_post)), 
                                       replace = F)
}

# Setup the cluster - remember the instructions from above.

cl <- makeCluster(50)
registerDoParallel(cl, cores = 50)

clusterExport(cl, c("Vol.present.container", "lambda_mu_post"))

# We iterate over each time series, fit extreme value distributions, and extract
# exceedance probabilities. These are the return levels for a given return period.
# Here we extract the 2- to 250-year event, but could be adapted to any desired return level.  
# Like the input, the output is a list with each entry showing the estimated return levels
# for the defined return periods, including their confidence intervals.

Vol.levels.present.HKKHN <- parLapply(Vol.present.container, cl = cl, function(x) {
  
  smpnum <- 10000
  
  # Obtain mean rate of events per year
  
  t_u <- lambda_mu_post
  
  # Define the above-threshold for the flood peaks, here the 80th percentile.
  
  th <- 0.8
  
  # Fit the Bayesian Poisson Generalized Pareto model.
  
  fit_pp <- extRemes::fevd(x,
                           threshold = quantile(x, th),
                           type = "PP",
                           method = "Bayes",
                           time.units = paste0(t_u, "/year"),
                           units = "cumecs",
                           priorParams = list(v = c(1, 10, 10)),
                           use.phi = TRUE,
                           verbose = TRUE)
  
  # Extract the return levels.
  
  reg_levels <- extRemes::return.level(fit_pp, do.ci = TRUE, return.period = c(2:250))
  
  return(reg_levels)
  
})

stopCluster(cl)

# You can save the estimated return periods to disk.

save(list = "Vol.levels.present.HKKHN", file =  "Vol_levels_present_HKKHN.RData")


# Plot the estimated return levels for a given return period ('i.e. hockey stick curves')
# The following part includes a lot of graphical adjustment, and expects a basic knowledge
# on generating and adjusting figures in R. 

options(scipen = 999)

# Note that figure closes approx. 250 lines later with command dev.off().

png(filename = "Return_level_plots_all_volume.png", units = "in", res = 500,
    width = 7.2, height = 4)

# The alignment of the figure panels corresponds to the layout SI Appendix, Fig. S6

l <- layout(matrix(c(3,2,1,
                    4,9,10,
                    5,11,12,
                    6,7,8), 3, 4))

par(omi = c(0.3,0.1,0.1,0.1),
    mar = c(0.3,2.5,1.5,0.3), 
    lend = 1, las = 1)

for (i in 1: length(Vol.levels.present.HKKHN)) {
  
  # Open an empty plot without showing any data(type = "n")
  # Set graphical limits, log x-axis, and so on.
  # Return levels are in million cumecs.
  
  if (i == 1) {plot(2:250, 
                    Vol.levels.present.HKKHN[[i]][, 2]/ 10^6, 
                    xlim = c(2, 208),
                    log  = "x", type = "n", 
                    ylim = c(0, 80), las = 1,
                    cex.axis = 0.6, cex.lab = 0.6, 
                    mgp = c(3.15, 0.5, 0),
                    yaxt = "n", tck = -0.04, 
                    ylab = "",  xlab = "")
    
    axis(side = 2, at = seq(0, 80, by = 20),  
         cex.axis=0.6, cex.lab = 0.6, las = 1,
         mgp = c(3.15, 0.5, 0), tck = -0.04,
         # the following command adds a comma as a thousands separator
         labels = formatC(seq(0, 80, by = 20), format = "d", big.mark = ','))
    
  }
  
  # Add the  estimated return periods from the 1,000 iterations.
  # The second column contains the mean estimates of the return levels.
  
  lines(2:250, Vol.levels.present.HKKHN[[i]][, 2]/ 10^6,
        col = rgb(red = 0, green = 0, blue = 100, alpha = 1, maxColorValue = 250) )
  
}

# Add a thick darkblue line representing the mean return level from all 1,000 model runs.

o <- sapply(Vol.levels.present.HKKHN, function (x) x[, 2])
Vol.mean <- apply(o, 1, mean)/ 10^6
lines(2:250, Vol.mean , lwd = 1.5, col = "darkblue")

# Calculate the variance for the estimated return levels from the 1,000 model runs.

HDIs <- apply(o, 1, HDIofMCMC)

box()
grid()

# Add return periods of empirical flood volumes as small ticks ('rug') to the plot.

rug(x = sapply(as.numeric(glofs.3dec.omerc$RelVol[!is.na(glofs.3dec.omerc$RelVol)]),
               function (x) {which.min(abs(x- Vol.mean))+1 }),
    side = 1, lwd = 2, col = "brown", ticksize = 0.05)

# Extract the 100-year flood volume and its 95% HDI.

Vol.100 <- round(Vol.mean["100-year return level"], digits = 2)

# Calculate the lower and the upper bound of the 100-year GLOF volume.

Vol.100.lower <- HDIs[1, ]["100-year return level"] / 10^6
Vol.100.upper <- HDIs[2, ]["100-year return level"] / 10^6

minus.bound <- round(Vol.100-Vol.100.lower, digits = 1)
plus.bound <- round(Vol.100.upper-Vol.100, digits = 1)

# Add the regional code to the plot.

text(2, 70,  pos = 4, "All regions",
     font = 2, cex = 0.9, col = "black")

# Add the 100-year flood volume and its 95% HDI to the plot.

text(2, 50,  pos = 4, 
     font = 3, cex = 0.9, col = "darkblue",
     substitute(paste(a[b]^c), list(a = Vol.100, 
                                    b = -1*minus.bound,
                                    c = paste0("+",plus.bound))))

# Add the posterior GLOF rate from the entire Himalayas to the plot.

text(2, 30, pos = 4, 
     labels = round(lambda_mu_post, digits = 2), font = 1, cex = 0.9, col = "brown")



####### Return periods of FLOOD VOLUMES for ALL SEVEN SUBREGIONS. ########

# To estimate regional GLOF return periods, we again generate a list, in which each entry
# holds 1,000 time series of 10,000 years length. The number of samples per year is 
# informed by the historical GLOF count in the past 30 years.

# We now iterate over each subregion, and generate a list with seven entries, one for each 
# subregion.

Vol.present.per.region.container <- list()

for (i in 1: nrow(aoi.buffer.int.omerc)) {
  
  # Generate an index to identify the lakes in a given subregion.
  
  int.pres <- gIntersects(aoi.buffer.int.omerc[i, ], lakes.sub, byid = T)
  
  # Extract the corresponding flood volumes from the matrix of predicted flood volumes.
  
  Vol.region.present <- all.Vol[, int.pres[, 1]]
  Vol.region.present <- as.numeric(Vol.region.present)
  
  # Extract all previously reported and newly detected GLOFs that happened in or 
  # 'intersect with' this region, respectively.
  
  landsat.int <- gIntersects(aoi.buffer.int.omerc[i, ], landsat.glofs.omerc, byid = T)
  glof.hist.int <- gIntersects(aoi.buffer.int.omerc[i, ], glofs.3dec.omerc, byid = T)
  
  # Calculate total sum of previously reported and newly detected GLOFs.
  
  counts_sum <- sum(c(landsat.int[,1], glof.hist.int[,1]))
  
  # Calculate posterior GLOF rate.
  
  lambda_mu_post = (counts_sum  + 1) / (sample_number + 1)
  
  # We now sample 1,000 time series of 10,000 years length, informed by the regional posterior 
  # GLOF rate. This part is identical to the approach for the entire Himalayas above.
  
  Vol.present.per.region <- list()
  
  for (idx in 1:1000) { 
    
    Vol.present.per.region[[idx]] <- sample(x = Vol.region.present,
                                            size = sum(rpois(years, lambda_mu_post)),
                                            replace = F)
  }
  
  Vol.present.per.region.container[[i]] <-  Vol.present.per.region
  
}

# Add names to the entries in the list to avoid confusion.

names(Vol.present.per.region.container) <- aoi.buffer.int.omerc$Region


# Generate a dataframe that holds the posterior GLOF rates for the seven
# subregions of the Himalaya.

rates.mu.posterior <- data.frame(matrix(NA, nrow = 1, ncol = 7))
colnames(rates.mu.posterior ) <- aoi.buffer.int.omerc$Region

# Define the length of study period, here 30 years.

sample_number <- 30

# Iterate over the study regions, extract the number of historical and newly detected 
# GLOFs in each region, sum them up, and add them to the dataframe.

for (i in 1:nrow(aoi.buffer.int.omerc)) {
  
  landsat.int <- gIntersects(aoi.buffer.int.omerc[i, ], landsat.glofs.omerc, byid = T)
  glof.hist.int <- gIntersects(aoi.buffer.int.omerc[i, ], glofs.3dec.omerc, byid = T)
  
  counts_sum <- sum(c(landsat.int, glof.hist.int)) 
  
  lambda_mu_post <- (counts_sum  +1) / (sample_number + 1) 
  
  rates.mu.posterior[1, i] <- lambda_mu_post
}


# As above, estimate return levels for the subregions, We here reorder the loop
# to sort the subregions from East to West. The final result is a list with seven entries, 
# each holding the return levels (2- to 250-year event) per region.

Vol.levels.present.per.region <- list()

for (p in c(3,4,5,6,7,1,2)) {
  
  # Extract the 1,000 time series per subregion.
  
  pres.sub <- Vol.present.per.region.container[[p]]
  
  # We iterate over the 1,000 time series in parallel and fit the extreme-value-models.
  
  cl <- makeCluster(50)
  registerDoParallel(cl, cores = 50)
  
  clusterExport(cl, c("pres.sub", "p", "rates.mu.posterior"))
  
  Vol.level.pres <- parLapply(pres.sub, cl = cl, function(x) {
    
    # Define number of years per time series.
    
    smpnum <- 10000
    
    # Obtain mean rate of events per year

    t_u <- rates.mu.posterior[1, p]
    
    # Define the threshold for the POT method, here the 80th percentile.
    
    th <- 0.8
    
    # Fit model the Poisson Generalized Pareto model.
    
    fit_pp <- extRemes::fevd(x,
                             threshold = quantile(x, th),
                             type = "PP",
                             method = "Bayes",
                             time.units = paste0(t_u, "/year"),
                             units = "cumecs",
                             priorParams = list(v = c(1, 10, 10)),
                             use.phi = TRUE,
                             verbose = TRUE)
    
    # Extract the regional return levels.
    
    reg_levels <- extRemes::return.level(fit_pp, do.ci = TRUE, return.period = c(2:250))
    
    return(reg_levels)
    
  })
  
  # Close the cluster; clear cache.
  
  stopCluster(cl)
  gc()
  
  # Add the estimated return levels to the list.
  
  Vol.levels.present.per.region[[p]] <- Vol.level.pres
  
  # Same code as above to plot return levels. See comments above.
  
  for (m in 1: length(Vol.level.pres)) {
    
    if (m == 1) {plot(2:250, 
                      Vol.level.pres[[m]][, 2]/ 10^6, 
                      xlim = c(2, 208),
                      log = "x", type = "n", 
                      ylim = c(0, 80), las = 1,
                      cex.axis = 0.6, cex.lab = 0.6, 
                      mgp = c(3.15, 0.5, 0),
                      yaxt = "n", tck = -0.04, 
                      ylab = "", xlab = "")
      axis(side = 2, at = seq(0, 80, by = 20),  cex.axis = 0.6, cex.lab = 0.6, las = 1,
           mgp = c(3.15, 0.5, 0), tck = -0.04,
           labels = formatC(seq(0, 80, by = 20), format = "d", big.mark = ','))
      
    }
    
    
    lines(2:250, Vol.level.pres[[m]][, 2]/ 10^6,
          col = rgb(red = 0, green = 0, blue = 100, alpha = 1, maxColorValue = 250))
    
  }
  
  o <- sapply(Vol.level.pres, function (x) x[, 2])
  Vol.mean <- apply(o, 1, mean)
  lines(2:250,   Vol.mean/ 10^6 , lwd = 1.5, col = "darkblue")
  
  HDIs <- apply(o, 1, HDIofMCMC)

  box()
  grid()
  
  Vol.100 <- round(Vol.mean["100-year level"] / 10^6, digits = 2)
  
  Vol.100.lower <- HDIs[1, ]["100-year level"] / 10^6
  Vol.100.upper <- HDIs[2, ]["100-year level"] / 10^6
  
  minus.bound <- round(Vol.100-Vol.100.lower, digits = 1)
  plus.bound <- round(Vol.100.upper-Vol.100, digits = 1)
  
  # Add the name of the subregion to the plot.
  
  text(2, 70, pos = 4, aoi.buffer.int.omerc$Region[p],
       font = 2, cex = 0.9, col = "black")
  
  
  text(2, 50, pos = 4, 
       font = 3, cex = 0.9, col = "darkblue",
       substitute(paste(a[b]^c), list(a= Vol.100, 
                                      b= -1*minus.bound ,
                                      c = paste0("+", plus.bound))))
  
  # Add the regional posterior GLOF rate to the plot.
  
  text(2, 30,  pos = 4, round(rates.mu.posterior[1,p], digits = 2), 
       font = 1, cex = 0.9, col = "brown")
  
}


# You can save the estimated return periods to disk.

save(list = "Vol.levels.present.per.region", file =  "Vol_levels_present_per_region.RData")

dev.off()


####################################################################################################
###############         Estimate return periods of PEAK DISCHARGE for             ##################
###############          the entire Himalayas and seven subregions                ##################
####################################################################################################

# List all directories with the predicted peak discharges per lake. 
# Note that the folder structure here is slightly different to the approach above: 
# here, each lake has its own folder of predicted peak discharges.

Qp.dirs <- list.dirs(path = "Qp_present", full.names = T, recursive = F)

# We iterate over the lake folders in parallel and read the predicted peak discharges
# to the memory. 

cl <- makeCluster(50)
registerDoParallel(cl, cores = 50)

clusterExport(cl, list("Qp.dirs", "HDIofMCMC" ))

all.Qp <- pbsapply(X = Qp.dirs, cl = cl, function (xxx) {
  
  runs <- list.files(xxx, full.names = T, recursive = F, include.dirs = T)
  
  o <- lapply(runs, readRDS)
  
  # We convert the predicted peak discharges to log-space, concatenate them
  # to one long vector, and truncate them to the 95% HDI.
  
  o2 <- as.numeric(log10(unlist(o)))
  hdi <- HDIofMCMC(o2)
  lake.hdi <- o2[o2 >= hdi[1] & o2 <= hdi[2]]
  return(lake.hdi)
  
  })

# Stop cluster and clear cache.

stopCluster(cl)
gc()


####### Return periods of PEAK DISCHARGE for the ENTIRE HIMALAYAS ########

# Define length of the study period.

sample_number <- 30

# Posterior GLOF rate for samples: Sum up the number of historical and newly detected 
# GLOFs and divide by the length of the study period.

lambda_mu_post = (nrow(glofs.3dec.omerc) + nrow(landsat.glofs.omerc) + 1) / (sample_number + 1) 

years <- 10000

# Generate a list with 1,000 entries, each holding a time series of 10,000 years.
# The number of GLOF peak discharges per year follows a Poisson distribution informed
# by the regional GLOF rate.

Qp.present.container <- list()

for (i in 1: 1000) {

  # Note that we log10-transformed Qp above, so that we need to re-transform it to the
  # original space.
  
  Qp.present.container[[i]] <- 10^sample(x = all.Qp,
                  size = sum(rpois(years, lambda_mu_post)), 
                  replace = F)
}

# The approach for estimating return periods of GLOF peak discharges is, in essence,
# equivalent for the GLOF flood volumes above: We iterate over the 1,000 time series,
# define an above-threshold, and fit the Poisson Generalized Pareto distribution.

cl <- makeCluster(50)
registerDoParallel(cl, cores = 50)

clusterExport(cl, c("Qp.present.container", "lambda_mu_post"))


Qp.levels.present.HKKHN <- parLapply(Qp.present.container, cl = cl, function(x) {
  
  # Define the length of the time series.
  
  smpnum <- 10000
  
  # Obtain mean rate of GLOFs per year.
  
  t_u <- lambda_mu_post
  
  # Define the above-threshold for the flood peaks, here the 80th percentile. 
  
  th <- 0.8
  
  # Fit the Poisson Generalized Pareto model.
  
  fit_pp <- extRemes::fevd(x,
                           threshold = quantile(x, th),
                           type = "PP",
                           method = "Bayes",
                           time.units = paste0(t_u, "/year"),
                           units = "cumecs",
                           priorParams = list(v = c(1, 10, 10)),
                           use.phi = TRUE,
                           verbose = TRUE)
  
  # Extract the return level for a (continous) range of return periods.
  
  Qp_levels <- extRemes::return.level(fit_pp, do.ci = TRUE, return.period = c(2:250))
  
  return(Qp_levels)
  
})

# Stop the cluster; clear cache.

stopCluster(cl)
gc()

# Save the estimated return periods for the entire Himalayas to disk.

save(list = "Qp.levels.present.HKKHN", file = "Qp_levels_present_HKKHN.RData")


# Now plot the estimated return period of GLOF peak discharges for the entire Himalayas.
# The code here is identical to the code for the flood volumes, so that it is left here
# largely uncommmented.

options(scipen = 999)

png(filename = "Peak_discharge_return_level_plots_all.png", units = "in", res = 500,
    width = 7.2, height = 4)

l <- layout(matrix(c(3, 2, 1,
                     4, 9, 10,
                     5, 11,12,
                     6, 7, 8), 3, 4))

par(omi = c(0.3, 0.1, 0.1, 0.1),
    mar = c(0.3, 2.5, 1.5, 0.3), 
    lend = 1, las = 1)

for (i in 1: length(Qp.levels.present.HKKHN)) {
  
  if (i == 1) {plot(2:250, 
                    Qp.levels.present.HKKHN[[i]][, 2], 
                    xlim = c(2, 208),
                    log  = "x", type = "n", 
                    ylim = c(0, 40000), las = 1,
                    cex.axis = 0.6, cex.lab = 0.6, 
                    mgp = c(3.15, 0.35, 0),
                    yaxt = "n", tck = -0.04, 
                    ylab = "", xlab = "")
    
    axis(side = 2, at = seq(0, 40000, by = 10000),  cex.axis = 0.6, cex.lab = 0.6, las = 1,
         mgp = c(3.15, 0.5, 0), tck = -0.04,
         labels = formatC(seq(0, 40000, by = 10000), format = "d", big.mark = ','))
    }
  
  
 lines(2:250, Qp.levels.present.HKKHN[[i]][, 2],
       col = rgb(red = 0, green = 0, blue = 100, alpha = 1, maxColorValue = 250) )
  
}

o <- sapply(Qp.levels.present.HKKHN, function (x) x[, 2])
qp.mean <- apply(o, 1, mean)
lines(2:250, qp.mean , lwd = 1.5, col = "darkblue")

HDIs <- apply(o, 1, HDIofMCMC)

box()
grid()

# Add small ticks of estimated return periods from historical peak discharges.

rug(x = sapply(as.numeric(glofs.3dec.omerc$QP[!is.na(glofs.3dec.omerc$QP)]), 
               function (x) {which.min(abs(x- qp.mean))+1 }),
    side = 1, lwd = 2, col = "brown", ticksize = 0.05)

# Extract mean 100-year GLOF peak discharge from the 1,000 simulations...

qp.100 <- round(qp.mean["100-year return level"])

# ... and its 95% HDI.

qp.100.lower <- HDIs[1, ]["100-year return level"]
qp.100.upper <- HDIs[2, ]["100-year return level"]

minus.bound <- round(qp.100-qp.100.lower)
plus.bound <- round(qp.100.upper-qp.100)

# Add the regional ID to the panel.

text(2,35000,  pos = 4, "All regions",
     font = 2, cex = 0.9, col = "black")

# Add the 100-year GLOF and its 95% HDI to the plot.

text(2,25000,  pos = 4, 
     font = 3, cex = 0.9, col = "darkblue",
     substitute(paste(a[b]^c), list(a = formatC(qp.100, format = "d", big.mark = ","), 
                                    b = formatC(-1*minus.bound, format = "d", big.mark = ",") ,
                                    c = paste0("+", formatC(plus.bound, format = "d", big.mark = ",")))))

# Add the posterior GLOF rate from the entire Himalayas to the panel.

text(2,15000,  pos = 4, round(lambda_mu_post, digits = 2), font = 1, cex = 0.9, col = "brown")


####### Return periods of PEAK DISCHARGE for the SEVEN SUBREGIONS ########

# As for GLOF volumes, we generate a list with seven entries for each subregion.
# For each subregion, sample 1,000 time series of 10,000 years length, which are
# informed by the regional GLOF rate (number of GLOFs divided by the length of the
# study period, i.e. 30 years).

qp.present.per.region.container <- list()

for (i in 1:nrow(aoi.buffer.int.omerc)) {
  
  int.pres <- gIntersects(aoi.buffer.int.omerc[i, ], lakes.sub, byid = T)
  
  Qp.region.present <- all.Qp[, int.pres[, 1]]
  Qp.region.present <- as.numeric(Qp.region.present)
  
  landsat.int <- gIntersects(aoi.buffer.int.omerc[i, ], landsat.glofs.omerc, byid = T)
  glof.hist.int <- gIntersects(aoi.buffer.int.omerc[i, ], glofs.3dec.omerc, byid = T)
  
  counts_sum <- sum(c(landsat.int[,1], glof.hist.int[,1])) 
  
  lambda_mu_post = (counts_sum + 1) / (sample_number + 1)
  
  qp.present.per.region <- list()
  
  for (idx in 1:1000) { 
    
    # Note that QP samples were converted into log space, so that we need to re-transform
    # to the original scale.
    
    qp.present.per.region[[idx]] <- 10^sample(x = Qp.region.present,
                                              size = sum(rpois(years, lambda_mu_post)),
                                              replace = F)
  }
  
  qp.present.per.region.container[[i]] <- qp.present.per.region
  
}

# Add names to list entries.

names(qp.present.per.region.container) <- aoi.buffer.int.omerc$Region
 
# Estimate regional return levels of GLOF peak discharge.
# The output is a list with seven entries, with each element holding
# 1,000 estimated GLOF return periods.

Qp.levels.present.per.region <- list()

for (p in c(3, 4, 5, 6, 7, 1, 2)) {
  
  pres.sub <- qp.present.per.region.container[[p]]

  cl <- makeCluster(50)
  registerDoParallel(cl, cores = 50)

  clusterExport(cl, c("pres.sub", "p", "rates.mu.posterior"))

  Qp.level.pres <- parLapply(pres.sub,cl = cl, function(x) {

    # Set number of years.
    
    smpnum <- 10000
    
    # Obtain mean regional rate of GLOFs per year

    t_u <- rates.mu.posterior[1, p]

    # The the above-threshold for the time series, here the 80th percentile.
    
    th <- 0.8

    # Fit Poisson Generalized Pareto model.
    
    fit_pp <- extRemes::fevd(x,
                             threshold = quantile(x, th),
                             type = "PP",
                             method = "Bayes",
                             time.units = paste0(t_u, "/year"),
                             units = "cumecs",
                             priorParams = list(v = c(1, 10, 10)),
                             use.phi = TRUE,
                             verbose = TRUE)
    
    # Extract the return levels for a given return period, here the 2- to 250-year event.

    Qp_levels <- extRemes::return.level(fit_pp, do.ci = TRUE, return.period = c(2:250))

    return(Qp_levels)

  })

  stopCluster(cl)
  
  # Add all estimated return levels from the 1,000 simulations to the list.

  Qp.levels.present.per.region[[p]] <- Qp.level.pres}

# Save the estimated regional return periods of GLOF peak discharges to disk.

save(list = "Qp.levels.present.per.region",file =  "Qp_levels_present_per_region.RData")
  

# Now plot the return periods of peak dischages in the same style as for the flood volumes.
# The following lines of code use exactly the same notation as above, except for axis limits and
# annotation, so that we skip any comments for the sake of brevity.

for(p in 1: length(Qp.levels.present.per.region)) {

  Qp.level.pres <-  Qp.levels.present.per.region[[p]]
  
  for (m in 1: length(Qp.level.pres)) {
    
    if (m == 1) {plot(2:250, 
                      Qp.level.pres[[m]][, 2], 
                      xlim = c(2, 208), cex.axis = 0.6, 
                      cex.lab = 0.6, yaxt = "n",
                      mgp = c(3.15, 0.35, 0),
                      tck = -0.04, log = "x", 
                      type = "n", ylim = c(0, 40000), las = 1,
                      ylab = "", xlab = "")
    
    axis(side = 2, at = seq(0, 40000, by = 10000),  
         cex.axis = 0.6, cex.lab = 0.6, las = 1,
         mgp = c(3.15, 0.5, 0), tck = -0.04,
         labels = formatC(seq(0, 40000, by = 10000), format = "d", big.mark = ','))
    }
    
    lines(2:250, Qp.level.pres[[m]][, 2],
          col = rgb(red = 0, green = 0, blue = 100, alpha = 1, maxColorValue = 250)  )
  
  }
  
  o <- sapply(Qp.level.pres, function (x) x[, 2])
  qp.mean <- apply(o, 1, mean)
  lines(2:250,   qp.mean , lwd = 2, col = "darkblue")
  
  HDIs <- apply(o, 1, HDIofMCMC)
 
  box()
  grid()
  
  
  qp.100 <- round(qp.mean["100-year level"])
  
  qp.100.lower <- HDIs[1,]["100-year level"]
  qp.100.upper <- HDIs[2,]["100-year level"]
  
  minus.bound <- round(qp.100-qp.100.lower)
  plus.bound <- round(qp.100.upper-qp.100)
  
  
  text(2,35000,  pos = 4, aoi.buffer.int.omerc$Region[p],
       font = 2, cex = 0.9, col = "black")
  
  
  text(2,25000,  pos = 4, 
       font = 3, cex = 0.9, col = "darkblue",
       substitute(paste(a[b]^c), list(a = formatC(qp.100, format = "d", big.mark = ","), 
                                      b = formatC(-1*minus.bound, format = "d", big.mark = ",") ,
                                      c = paste0("+", formatC(plus.bound, format = "d", big.mark = ",")))))
  
  text(2,15000,  pos = 4, round(rates.mu.posterior[1,p], digits = 2), 
       font = 1, cex = 0.9, col = "brown")
  
}

dev.off()



####################################################################################################
##############         Plot frequency densities of flood volumes and  peak      ####################
##############    discharges from all moraine-dammed lakes in the  Himalayas    ####################
####################################################################################################


# Plot Fig. 2 with frequency-density curves of flood volumes and peak discharges estimated from 
# the present lake-size distribution in the Himalayas.
# Write a PDF with eight subpanels: a) for the complete Himalaya, b-h) for the seven RGI regions.

# Convert matrices with samples of Vol and Qp to vectors, and sort the shapefile with the subregions
# from West to East.

all.Qp.num <- as.numeric(all.Qp)
all.Vol.num <- as.numeric(all.Vol)
aoi.buf.sort <- aoi.buffer.int.omerc[c(3, 4, 5, 6, 7, 1, 2), ] 

# Open the PDF with some layout specifications.

pdf(file = "Vol_Qp_for_all_subregions_freq_density.pdf", 
    width = 7.2, height = 4.4)

# We use narrow figure margins and small distances between the panels.

par(mfrow = c(2, 4), 
    oma =   c(3, 3, 2, 2),
    mai =   c(0.01, 0.01, 0.01, 0.01), 
    mar =   c(0.3, 0.3, 0.3, 0.3), 
    lend =  1)


# Define a function that creates equally spaced bins on log scale. We use this function
# to calculate frequency densities, i.e. the number of GLOFs of given bin size, normalized
# by the size of the bin.

lseq <- function(from = 1, to = 100000, length.out = 6) {
  exp(seq(log(from), log(to), length.out = length.out))}

# Define the range, for which histogram bins shall be calculated. To avoid empty bins,
# we adapt the range to the distribution of the data.

min.all.vol <- 10^(log10(min(all.Vol.num))-0.1)
max.all.vol <- 10^(log10(max(all.Vol.num))+0.1)

brea <- lseq(from = min.all.vol, 
             to = max.all.vol, 
             length.out = 45) # we obtain 45 bins


# Plot the first panel (a) for the entire Himalayas. hist() is a useful comment: it splits
# the data into bins that we provide with the lseq() function, and returns the mids of 
# the bins that we need for plotting.

### GLOF VOLUMES FOR THE ENTIRE HIMALAYAS

h1 <- hist(all.Vol.num, breaks = brea, plot = FALSE)

x = h1$mids

# Normalize frequency per bin by the bin width, and normalize all bins, so that the
# the density adds to 1.

y =  h1$counts/ diff(brea)
y = y / sum(y)

# Define the range of the x- and y-axis.

XLIM = c(5*10^-2,  0.5*10^10)
YLIM = c(5*10^-12, 0.5*10^0)

# Plot the frequency density as points in log-log space.

plot(x, y, 
     log= "xy", pch = 21, 
     type = "p", cex = 1,
     xlim = XLIM, ylim = YLIM,
     col = "black", bg =  "navy",
     lwd = 0.25, axes = F,
     xlab = "", ylab = "", 
     las = 1, xaxt = "n", yaxt = "n")

# Add some nice logarithmic axes on all sides. 'eaxis()' is a flexible extension for axis
# labels in the 'sfsmisc' package.

eaxis(side = 1, cex.axis = 0.65, cex.lab = 0.65, mgp = c(3, 0.45, 0),
      tck = -0.02, labels = F, at.small = FALSE)

eaxis(side = 2, cex.axis = 0.65, cex.lab = 0.65, mgp = c(1, 0.45, 0), 
      tck = -0.03, at.small=FALSE)

eaxis(side = 3, cex.axis = 0.65, cex.lab = 0.65, mgp = c(3, 0.3, 0),
      tck = -0.03, at.small=FALSE)

eaxis(side = 4, cex.axis = 0.65, cex.lab = 0.65, mgp = c(3, 0.45, 0),
      tck = -0.02, labels = F, at.small = FALSE)

# Add a bounding box around the plot...

box()

# ... and some additional lines for orientation.

abline(v = 10^seq(-2, 10, by = 2), lty = 2, col = "grey50", lwd = 0.4)
abline(h = 10^seq(-12, 0, by = 2), lty = 2, col = "grey50", lwd = 0.4)

# Add axes labels for the entire plot to the outer figure margins

mtext(side = 1, line = 2, 
      expression('Peak discharge Q'['p']*' [m'^3*' s'^-1*'] / Flood volume V'[0]*' [m'^3*']'), 
      cex = 0.9, outer = T)

mtext(side = 2, line = 1.5, 
      expression('Frequency density'), 
      cex = 0.9, outer = T)


### PEAK DISCHARGES FOR THE ENTIRE HIMALAYAS

# Estimating the frequency densities for Qp is literally the same approach as above for
# GLOF volumes. We take 100 mio samples of Qp to prevent the hist() function to crash. 

qp.samp <- 10^sample(all.Qp.num, 10^8, replace = F)


min.all.vol <- 10^(log10(min(qp.samp))-0.1)
max.all.vol <- 10^(log10(max(qp.samp))+0.1)

brea <- lseq(from = min.all.vol, 
             to = max.all.vol, 
             length.out = 45)

h1 <- hist(qp.samp, breaks = brea, plot = FALSE)

x = h1$mids

# Normalize the frequency per bin by the bin width.

y =  h1$counts/ diff(brea)
y = y/ sum(y)

# Add the frequency densities of Qp to the existing plot by using the
# function 'points()' .

points(x, y,  pch = 21, cex = 1, col = "black",
       bg = "darkorange", lwd = 0.25)


# Label the plot with the letter "a" in the topright corner.

text(x = 0.5*max(XLIM), y =  0.5*max(YLIM), 
  labels = "a", font = 2, cex = 1.5, adj = c(0.5, NA)) 

# Put the region in the lowerleft corner

text(x = 0.9*min(XLIM), y = 10^-11, 
     labels = "All regions", 
     font = 2, cex = 0.9, pos = 4) 

# Add the number of glacier lakes in this region.

text(x = 0.9*min(XLIM), 
     y = 5*10^-8, 
     labels = length(lakes.sub), 
     col = "gray40",
     font = 3, cex = 0.9, pos = 4) 

# Add the area of glacier lakes [in square kilometres] in this region.

text(x = 0.9*min(XLIM), 
     y =  10^-9, 
     labels = round(sum(lakes.sub$Area)/10^6, digits = 1), 
     col = "gray40",
     font = 3, cex = 0.9, pos = 4) 


# THE SEVEN SUBREGIONS

for (i in seq_along(aoi.buf.sort)) {
  
  # Select a subregion. 
  
  int <- gIntersects(aoi.buf.sort[i, ], lakes.sub, byid = T)
  
  ### FLOOD VOLUMES.
  
  region.vol <- all.Vol[int[,1 ], ]
  
  # Do the same calculations for the entire Himalayas above:define the bin sizes.
  
  min.region.vol <- 10^(log10(min(region.vol))-0.1)
  max.region.vol <- 10^(log10(max(region.vol))+0.1)
  
  brea <- lseq(from = min.region.vol, 
               to = max.region.vol, 
               length.out = 45)
  
  
  # Use hist() to calculate frequencies per bin.
  
  h1 <- hist(as.numeric(region.vol), breaks = brea, plot = FALSE)
  x = h1$mids
  
  # Normalize frequency per bin by the bin width.
  
  y =  h1$counts/ diff(brea)
  y = y/ sum(y)
  
  # Plot the frequency densities as points on log-log scale.
  
  plot(x, y, 
       log ="xy", pch=21, type = "p", 
       xlim = XLIM, ylim = YLIM,
       cex = 1, col = "black",
       bg=  "navy", lwd = 0.25,
       axes = F, las = 1,
       xlab = "", ylab = "", 
       xaxt = "n", yaxt = "n")
  
  
  # Now there comes a lot of axis labeling:
  # Panels should only have labels at the outer margins on all sides,
  # but not at the interior. There could be smarter solutions than the one below.
  
  if (i  %in% c(1:3)) {
    eaxis(1, cex.axis = 0.65, cex.lab = 0.65, mgp = c(1, 0.3, 0), 
          tck = -0.02, labels = F, at.small = FALSE) 
    
    eaxis(3, cex.axis = 0.65, cex.lab = 0.65, mgp = c(1, 0.3, 0), 
          tck=-0.03, at.small = FALSE) 
    
  } else { eaxis(1, cex.axis = 0.65, cex.lab = 0.65, mgp = c(3, 0.45, 0),
                 tck = -0.03, at.small = FALSE)
    
    eaxis(3, cex.axis = 0.65, cex.lab = 0.65, mgp = c(3, 0.45, 0),
                   tck = -0.02, at.small = FALSE, labels = F)
  }
  
  if (i  %in% 4) { eaxis(2, cex.axis = 0.65, cex.lab = 0.65, mgp = c(3, 0.45, 0),
                         tck = -0.03, at.small = FALSE) }
  
  if (i  %in% c(3,7)) { eaxis(4, cex.axis = 0.65, cex.lab = 0.65, mgp = c(1, 0.3, 0), 
                              tck = -0.03, at.small = FALSE) }
  
  if (i  %in% c(1,2,4,5,6)) { eaxis(4, cex.axis = 0.65, cex.lab = 0.65, 
                                    mgp = c(1, 0.3, 0),  tck = -0.02, labels = F,
                                    at.small = FALSE) }
  
  if (i  %in% c(1,2,3,5,6,7)) { eaxis(3, cex.axis = 0.65, cex.lab = 0.65, 
                                      mgp = c(1, 0.3, 0), tck = -0.02, labels = F,
                                      at.small = FALSE) }
  
  if (i  %in% c(1,2,3,5,6,7)) { eaxis(2, cex.axis = 0.65, cex.lab = 0.65, 
                                      mgp = c(1, 0.45, 0), tck = -0.02, 
                                      labels = F, at.small = FALSE)}
 
  # Add a bounding box and some additional lines for orientation.
   
  box()
  abline(v = 10^seq(-2, 10, by = 2), lty = 2, col = "grey50", lwd = 0.4)
  abline(h = 10^seq(-12, 0, by = 2), lty = 2, col = "grey50", lwd = 0.4)
  
  
  ### PEAK DISCHARGE
  
  # Now we do the same calculations for Qp.
  
  region.qp <- 10^all.Qp[,int[,1 ] ]
  
  if(length(region.qp) > 10^8) {
    region.qp <- sample(region.qp, 10^8)
    
  }

  min.region.qp <- 10^(log10(min(region.qp))-0.1)
  max.region.qp <- 10^(log10(max(region.qp))+0.1)
  
  brea <- lseq(from = min.region.qp, 
               to = max.region.qp, 
               length.out = 50)
  
  h1 <- hist(as.numeric(region.qp), breaks = brea, plot = FALSE)
  
  x = h1$mids
  
  # Normalize frequency per bin by the bin width.
  
  y =  h1$counts/ diff(brea)
  y = y/ sum(y)
  
  # Plot the points.
  
  points(x, y,  pch=21, cex = 1,col = "black",
         bg =  "darkorange",
         lwd = 0.25)

  # Label the subregions sequentially from b-h in the topright corner.
  
  text(x = 0.5*max(XLIM), y = 0.55*max(YLIM), 
       labels = letters[i+1], font = 2, cex = 1.5, adj = c(0.5, NA)) 
  
  # Add the name of the subregion to the topright corner.
  
  text(x = 0.9*min(XLIM), 
       y = 10^-11, 
       labels = aoi.buf.sort$Region[i], 
       font = 2, cex = 0.9, pos = 4) 
  
  # Add the number of glacier lakes in the subregion.
  
  text(x = 0.9*min(XLIM), 
       y = 5*10^-8, 
       labels = sum(int[, 1]), 
       col = "gray40",
       font = 3, cex = 0.9, pos = 4) 
  
  # Add  the area of glacier lakes in the subregion.
  
  text(x = 0.9*min(XLIM), 
       y =  10^-9, 
       labels = round(sum(lakes.sub$Area[int[, 1]])/10^6, digits = 1), 
       col = "gray40",
       font = 3, cex = 0.9, pos = 4) 
  }

dev.off()

# Plot finished!

####################################################################################################
##############         Explore influence of lake size distribution by           ####################
##############         assuming a fixed GLOF rate of 1 GLOF per year            ####################
####################################################################################################

# The model setup is identical to the workflow for estimating the regional GLOF hazard above. 
# The only difference is that we exchange the regional GLOF rate (derived from historical GLOFs
# in a given region) by a fixed GLOF rate of 1 GLOF per year.

### First: RELATIVE GLOF VOLUMES

# We take the usual model setup: 1,000 time series of 10,000 years length. This
# time, we fix the GLOF rate in the Poisson process at 1 in the entire Himalayas 
# and all subregions.

years <- 10000

Vol.present.HKKHN.lambda1 <- list()

for (i in 1: 1000) {
  
  Vol.present.HKKHN.lambda1[[i]] <- sample(x = all.Vol,
                                       size = sum(rpois(years, 1)), # lambda = 1
                                       replace = F)
}


cl <- parallel::makeCluster(50)
registerDoParallel(cl, cores = 50)
clusterExport(cl, c("Vol.present.HKKHN.lambda1"))

Vol.levels.HKKHN.lambda1 <- parSapply(Vol.present.HKKHN.lambda1, cl = cl, function(x) {
  
  # Define length of time series.
  
  smpnum <- 10000
  
  # Fix the mean rate of GLOFs per year at 1.
  
  t_u <- 1
  
  # Set above-threshold.
  
  th <- 0.8
  
  # Fit the Poisson Generalized Pareto model.
  
  fit_pp <- extRemes::fevd(x,
                           threshold = quantile(x, th),
                           type = "PP",
                           method = "MLE",
                           time.units = paste0(t_u, "/year"),
                           units = "cumecs",
                           priorParams = list(v = c(1, 10, 10)),
                           use.phi = TRUE,
                           verbose = TRUE)
  
  # We only extract the 100-year return level for simplicity.
  
  vol_levels <- as.numeric(extRemes::return.level(fit_pp,  return.period = c(100)))
  
  return(vol_levels)
  
})

stopCluster(cl)

# Save the list with the 100-year flood volumes, assuming a rate of 1 GLOF per year, to disk.

save(list = "Vol.levels.HKKHN.lambda1",
     file =  "Vol_levels_HKKHN_lambda1 .RData")


vol100.mean.lambda1  <- mean(Vol.levels.HKKHN.lambda1)/ 10^6

# Generate a list with seven entries, for each subregion.
# Sample 1,000 time series of GLOF volumes for each region using a rate
# of 1 GLOF per year, and add the samples to the list.

Vol.per.region.lambda1 <- list()

for (i in 1:nrow(aoi.buffer.int.omerc)) {
  
  int.pres <- gIntersects(aoi.buffer.int.omerc[i, ], lakes.sub, byid = T)
  
  Vol.region.present <- all.Vol[ int.pres[, 1], ]
  Vol.region.present <- as.numeric(Vol.region.present)
  
  Vol.present.per.region <- list()
  
  for (idx in 1:1000) {
    
    Vol.present.per.region[[idx]] <- sample(x = Vol.region.present,
                                            size = sum(rpois(years, 1)), # lambda = 1
                                            replace = F)
  }
  
  Vol.per.region.lambda1[[i]] <-  Vol.present.per.region
  
}


names(Vol.per.region.lambda1) <- aoi.buffer.int.omerc$Region

# Estimate the 100-year GLOF volume for each subregion using a fixed rate of 1 GLOF
# per year.

Vol100.region.lambda1 <- list()

for (p in 1:length(Vol.per.region.lambda1)) {
  
  pres.sub <- Vol.per.region.lambda1[[p]]
  
  cl <- makeCluster(50)
  registerDoParallel(cl, cores = 50)
  
  clusterExport(cl, c("pres.sub"))
  
  Vol.level.pres <- parLapply(pres.sub,cl = cl, function(x) {
    
    smpnum <- 10000
    
    # Keep the GLOF rate fixed at 1 GLOF per year.
    
    t_u <- 1
    th <- 0.8
    
    # Fit the Poisson GP model.
    
    fit_pp <- extRemes::fevd(x,
                             threshold = quantile(x, th),
                             type = "PP",
                             method = "MLE",
                             time.units = paste0(t_u, "/year"),
                             units = "cumecs",
                             priorParams = list(v = c(1, 10, 10)),
                             use.phi = TRUE,
                             verbose = TRUE)
    
    vol_levels <- as.numeric(extRemes::return.level(fit_pp,  return.period = c(100)))
    
    return(vol_levels)
    })
  
  stopCluster(cl)
  
  Vol100.region.lambda1[[p]] <- Vol.level.pres
  
}

# Concatenate the seven list elements to a matrix with seven rows,
# where the columns are the 1,000 estimated 100-GLOF volumes. 

vol100.l1 <- t(sapply(Vol100.region.lambda1, unlist ))

rownames(vol100.l1) <- aoi.buffer.int.omerc$Region

# Sort matrix positions from West to East, convert the GLOF volumes to 
# million cubicmetres.

vol100.l1 <- vol100.l1[c(3, 4, 5, 6, 7, 1, 2), ]
vol100.l1 <- vol100.l1 / 10^6

# Divide all estimated 100-year flood volumes by the mean estimate from the entire
# Himalayas to obtain an estimate of the relative regional GLOF hazard.

vol100.l1.rel <- vol100.l1/ vol100.mean.lambda1 

# Convert this matrix into a 3-column dataframe, 
# with a column 'value' containing the 100-year GLOFs relative to the entire Himalayas,
# a column 'region' denoting to which region the 100-year GLOFs belong,
# and a column 'type' saying that the column contains GLOF VOLUMES.

vol100.rel.df <- data.frame(
  
  value = as.numeric(t(vol100.l1.rel))*100,
  region = rep(rownames(vol100.l1.rel), each = ncol(vol100.l1.rel)),
  type = rep("vol", length(vol100.l1.rel))
  
)

### Second: RELATIVE GLOF PEAK DISCHARGES

# Comments only where necessary.

Qp.HKKHN.container.lambda1 <- list()

for (i in 1: 1000) {
  
  Qp.HKKHN.container.lambda1[[i]] <- 10^sample(x = all.Qp,
                                         size = sum(rpois(years, 1)), # lambda = 1
                                         replace = F)
}


cl <- makeCluster(50)
registerDoParallel(cl, cores = 50)
clusterExport(cl, c("Qp.HKKHN.container.lambda1"))

qp100.HKKHN.lambda1 <- parLapply(Qp.HKKHN.container.lambda1, cl = cl, function(x) {
  
  # Number of years in the time series
  
  smpnum <- 10000

  # Fix mean annual rate of GLOFs at 1.
  
  t_u <- 1
  th <- 0.8
  
  # Fit the Poisson GP model.
  
  fit_pp <- extRemes::fevd(x,
                           threshold = quantile(x, th),
                           type = "PP",
                           method = "MLE",
                           time.units = paste0(t_u, "/year"),
                           units = "cumecs",
                           priorParams = list(v = c(1, 10, 10)),
                           use.phi = TRUE,
                           verbose = TRUE)
  
  # Extract only the 100-year GLOF peak discharge.
  
  qp_levels <- as.numeric(extRemes::return.level(fit_pp,  return.period = c(100)))
  
  return(qp_levels)
  
})

stopCluster(cl)

save(list = "qp100.HKKHN.lambda1", file =  "qp100_HKKHN_lambda1.RData")

qp100.mean.lambda1 <- mean(unlist(qp100.HKKHN.lambda1))


# Estimate the 100-year GLOF volume for each subregion using a fixed rate of 1 GLOF
# per year.

qp.present.per.region.lambda1 <- list()

for (i in 1:nrow(aoi.buffer.int.omerc)) {
  
  int.pres <- gIntersects(aoi.buffer.int.omerc[i, ], lakes.sub, byid = T)
  
  Qp.region.present <- all.Qp[ , int.pres[, 1]]
  Qp.region.present <- as.numeric(Qp.region.present)
  
  qp.present.per.region <- list()
  
  for (idx in 1:1000) { 
    
    qp.present.per.region[[idx]] <- 10^sample(x = Qp.region.present,
                                              size = sum(rpois(years, 1)), # lambda = 1
                                              replace = F)
  }
  
  qp.present.per.region.lambda1[[i]] <-  qp.present.per.region
  
}


names(qp.present.per.region.lambda1) <- aoi.buffer.int.omerc$Region

# Estimate the 100-year GLOF peak discharge for each subregion using a fixed rate of 1 GLOF
# per year.

qp100.per.region.lambda1 <- list()

for (p in 1:length(qp.present.per.region.lambda1)) {
  
  pres.sub <- qp.present.per.region.lambda1[[p]]
  
  cl <- makeCluster(50)
  registerDoParallel(cl, cores = 50)
  
  clusterExport(cl, c("pres.sub"))
  
  qp.level.pres <- parSapply(pres.sub,cl = cl, function(x) {
    
    smpnum <- 10000
    
    # Keep the GLOF rate fixed at 1 GLOF per year
    
    t_u <- 1
    th <- 0.8
    
    # Fit Poisson GP model.
    
    fit_pp <- extRemes::fevd(x,
                             threshold = quantile(x, th),
                             type = "PP",
                             method = "MLE",
                             time.units = paste0(t_u, "/year"),
                             units = "cumecs",
                             priorParams = list(v = c(1, 10, 10)),
                             use.phi = TRUE,
                             verbose = TRUE)
    
    qp_levels <- as.numeric(extRemes::return.level(fit_pp,  return.period = c(100)))
    
    return(qp_levels)
    
  })
  
  stopCluster(cl)
  
  qp100.per.region.lambda1[[p]] <- qp.level.pres
 
}


qp100.l1 <- do.call(rbind, qp100.per.region.lambda1)

# Sort matrix positions from West to East.

qp100.l1 <- qp100.l1[c(3,4,5,6,7,1,2), ]

# Divide all estimated 100-year peak discharges by the mean estimate from the entire
# Himalayas to obtain an estimate of the relative regional GLOF hazard.

qp100.l1.rel <- qp100.l1/ qp100.mean.lambda1 
row.names(qp100.l1.rel) <- aoi.buffer.int.omerc$Region[c(3,4,5,6,7,1,2)]

# Convert this matrix into a 3-column dataframe, 
# with a column 'value' containing the 100-year GLOFs relative to the entire Himalayas,
# a column 'region' denoting to which region the 100-year GLOFs belong,
# and a column 'type' saying that the column contains GLOF PEAK DISCHARGES.

qp100.rel.df <- data.frame(
  
  value = as.numeric(t(qp100.l1.rel))*100,
  region = rep(rownames(qp100.l1.rel), each = ncol(qp100.l1.rel)),
  type = rep("Qp", length(qp100.l1.rel)),
  stringsAsFactors = F
  
)

# Concatenate the dataframes of regional relative GLOF volumes and peak discharges.

rel.vol.qp <- rbind(vol100.rel.df, qp100.rel.df)

# Convert the colum 'region' to a factor.

rel.vol.qp$region <- factor(rel.vol.qp$region, 
                            levels =  c("Hindu Kush", 
                                        "Karakoram", 
                                        "Western Himalaya", 
                                        "Central Himalaya", 
                                        "Eastern Himalaya",  
                                        "Nyainqentanglha" ,
                                        "Hengduan Shan"))

# Plot Fig. 4: Open a new PDF.

pdf("Fractions_of_hazard.pdf",
    width = 7,
    height = 2.94)

# Set some graphical parameters.

par(mar = c(0.5, 5.5, 0.5, 0.5),
    lend = 1,  mgp = c(2.5, 0.6, 0))

# We now generate a boxplot that is grouped I) by both the 'type' or metric of hazard, i.e.
# flood volume or peak discharge, and II) by the seven subregions.
# We here also modify the basic arguments in the boxplot() function, making
# the median line white, thicker, etc.

boxplot(value ~ type + region, 
        data = rel.vol.qp,
        col = c("navy", "darkorange"), # Volumes are blue, peak discharges are orange
        outline = F, pars = list(boxwex = 1),
        ylab = "Fraction of regional\n vs. Himalaya-wide hazard [%]",
        las = 1, whisklty = 1,
        medcol = "white", medlwd = 1,
        at = (1:20)[-seq(3, 18, by = 3)],
        xlab = "", cex.lab = 0.9,
        cex.axis = 0.8, xaxt = 'n',
        ylim = c(0, 220))

# Add a dotted line to separte the regions.

abline(v = seq(3, 18, by = 3), lty = 3, lwd = 0.5, col = "gray60")

# Write the names of the regions below the boxplot

text(x = seq(1.5, 19.5, by = 3),
     y = -5,
     adj = c(0.5, 0),
     cex = 0.75,
     labels =  c("Hindu Kush", 
                 "Karakoram" , 
                 "Western\n Himalaya", 
                 "Central\n Himalaya", 
                 "Eastern\n Himalaya",  
                 "Nyainqen-\ntanglha" ,
                 "Hengduan\n Shan"))

# Draw a grey horizontal line at 100% and label it accordingly.

abline(h = 100, lty = 1, lwd = 2, col = "gray30")
text(x = 0, y = 110, 
     pos = 4, col = "gray30", cex = 0.8,
     labels =  expression('V'[100]* ' / Q'['p100']*' (entire study area)' ))

# Finally add a legend to the top left corner.

legend(x = 0, y = 220, 
       x.intersp = 0.5,
       cex = 0.8,
       bty = "n",
       fill = c("navy", "darkorange"), 
       legend = c( expression('V'['100']),
                   expression('Q'['p100'])))
  
dev.off()

# Done!
