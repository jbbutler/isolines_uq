# Script to simulate coverage rates for the confidence region construction process in Mammen and Polonik (2013),
# assuming we are an oracle which has access to the true survival function.
#
# Jimmy Butler, December 2024

library(mvtnorm)
library(mnormt)
library(dplyr)
library(foreach)
library(doSNOW)
library(parallel)

# loading up helper functions
source('~/isolines_uq/scripts/R/auxiliary_scripts/utils.R')
source('~/isolines_uq/scripts/R/auxiliary_scripts/distributionIsolines.R')

### simulation specs ###

# specifications for the survival function grid
lbs <- c(-2,-2)
ubs <- c(10, 10)
n_boundary_coords <- 100
n_interior_coords <- 5

# specifications for the confidence regions
ps <- c(0.1)
alpha <- 0.05
ns <- c(1000, 10000, 50000)
isolines <- list()

for (i in 1:length(ps)) {
    isolines[[as.character(ps[i])]] <- drawBivtIsoline(numCoords=400, gridUbs=c(10,10), gridLbs=c(-2,-2), prob=ps[i])
}

# list of betas to choose from (ranging from small to large)
# highest beta chosen as 1+p because this beta covers all of R2, any higher beta will still cover this
n_betas <- 64
betas <- 1/exp(seq(3, 20, length.out=n_betas))

# parallelization specifications
n_cores <- 64
# number of simulations for each beta (to determine coverage rate)
n_sims <- 5000
# number of monte carlo simulations for determining bhat (for a particular beta)
n_monte_carlo <- 10000

### helper functions **PUT THIS IN A SEPARATE SCRIPT EVENTUALLY** ###

# function that draws a grid of points over which we wish to evaluate the sup difference between
# the true survival function and the empirical estimate
drawSupRegion <- function(n_boundary_coords, n_interior_coords, gridLbs, gridUbs, p, beta) {
    
    # function that we will use our root-finding algorithm for
    pmvtFixedCoord <- function(radius, angle, sigma, df, prob, xCenter, yCenter) {
        xCoord <- xCenter + radius*cos(angle)
        yCoord <- yCenter + radius*sin(angle)
        exceedanceProb <- pmvt(lower=c(xCoord, yCoord), upper=c(Inf, Inf), df=df, sigma=sigma)
        return(exceedanceProb - prob)
    }
    
    radii <- rep(NA, n_boundary_coords*n_interior_coords)
    full_angles <- rep(NA, n_boundary_coords*n_interior_coords)
    angles <- seq(0, pi/2, length.out=n_boundary_coords)
    maxRad <- sqrt(sum((gridUbs-gridLbs)**2))

    # for each angle in first quadrant, find radius that gives a point with desired exceedance probability
    # by finding roots of pmvtFixedCoord given the angle
    # find both the lower and upper radius, and then store a mesh of them
    for (i in 1:n_boundary_coords) {
        angle <- angles[i]
        lower_radius <- uniroot(pmvtFixedCoord, interval=c(0, maxRad), angle=angle,
                             sigma=matrix(c(1, 0.7, 0.7, 1), nrow=2), df=4, prob=p+beta, xCenter=gridLbs[1], yCenter=gridLbs[2])$root
        upper_radius <- uniroot(pmvtFixedCoord, interval=c(0, maxRad), angle=angle,
                             sigma=matrix(c(1, 0.7, 0.7, 1), nrow=2), df=4, prob=p-beta, xCenter=gridLbs[1], yCenter=gridLbs[2])$root
        
        radii_per_angle <- seq(lower_radius, upper_radius, length.out=n_interior_coords)
        angles_per_angle <- rep(angle, length.out=n_interior_coords)
        
        radii[((i-1)*n_interior_coords + 1):(i*n_interior_coords)] <- radii_per_angle
        full_angles[((i-1)*n_interior_coords + 1):(i*n_interior_coords)] <- angles_per_angle
        
        
    }

    # convert back to cartesian coordinates
    xs <- radii*cos(full_angles) + gridLbs[1]
    ys <- radii*sin(full_angles) + gridLbs[2]
    
    return(cbind(xs, ys))
    
}

empSurv <- function(pt, dat) {
    point <- as.numeric(pt)
    return(mean((point[1] < dat[,1]) & (point[2] < dat[,2])))
}

### the actual simulating... ###

# big function to simulate coverage rates for a particular beta (we will parallelize the loop over beta)
parallelizedCode <- function(ind) {
    
    beta <- betas[ind]
    
    region <- drawSupRegion(n_boundary_coords, n_interior_coords, lbs, ubs, p, beta)
    true_survfunc <- pmt(x=-as.matrix(region), mean=rep(0,2), S=matrix(c(1,0.7,0.7,1), 2,2), df=4)
    
    ### compute bhat ###
    Z_betas <- rep(NA, n_monte_carlo)
    for (i in 1:n_monte_carlo) {
        sample_dat <- data.frame(rmvt(n, sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2), df = 4))
        sample_survfunc <- apply(region, 1, empSurv, dat=sample_dat)
        Z_betas[i] <- max(abs(sample_survfunc - true_survfunc))
    }
    bhat <- as.numeric(quantile(Z_betas, probs = 1-alpha))
    
    ### compute coverage rate with this bhat ###
    is_covereds <- rep(NA, n_sims)
    
    for (i in 1:n_sims) {
        dat <- data.frame(rmvt(n, sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2), df = 4))
    
        findEmpSurv <- function(row) {
            return(mean((dat[,1] > row[[1]]) & (dat[,2] > row[[2]])))
        }
        survs <- apply(isolines[[as.character(p)]], 1, findEmpSurv)
        is_covereds[i] <- all((survs <= p + bhat) & (survs >= p - bhat))
    } 
    cov_rate <- mean(is_covereds)
    return(c(beta, bhat, cov_rate))
}

### parallelized loop over betas ###

for (i in 1:length(ns)) {
    
    for (j in 1:length(ps)) {
    
        n <- ns[i]
        p <- ps[j]

        clust <- makeSOCKcluster(n_cores)
        registerDoSNOW(clust)
        pb <- txtProgressBar(min = 1, max = n_betas, style = 3)
        progress <- function(n) setTxtProgressBar(pb ,n)
        opts <- list(progress = progress)

        results <- foreach(k = 1:n_betas, 
                      .options.snow = opts, 
                      .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt'), 
                      .combine='rbind') %dopar% parallelizedCode(k)

        close(pb)
        stopCluster(clust)

        saveRDS(results, file=paste0('/pscratch/sd/j/jbbutler/oracle_analysis/bivt/p', p, '_alpha', alpha, '_n', n, '.RData'))
        
    }
    
}
