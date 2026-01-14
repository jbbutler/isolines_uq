# Script to simulate draws from the Z(beta) random variable from Mammen and Polonik, assuming
# an oracle situation where the true survival function and true sampling distribution is known.
#
# Jimmy Butler, January 2025

library(mvtnorm)
library(mnormt)
library(dplyr)
library(foreach)
library(doSNOW)
library(parallel)

# loading up helper functions
source('~/isolines_uq/scripts/R/auxiliary_scripts/utils.R')

# save path
save_path <- '/pscratch/sd/j/jbbutler/oracle_analysis/bivt/zbetas/'

### simulation specs ###

# specifications for the survival function grid
lbs <- c(-2,-2)
ubs <- c(10, 10)
n_boundary_coords <- 100
n_interior_coords <- 5

# specifications for the confidence regions
ps <- c(0.01)
ns <- c(50000)
isolines <- list()

# list of betas to choose from (ranging from small to large)
n_betas <- 64
betas <- 1/exp(seq(5, 20, length.out=n_betas))

# parallelization specifications
n_cores <- 64
# number of monte carlo simulations for determining bhat (for a particular beta)
n_monte_carlo <- 10000

# function that draws a grid of points over which we wish to evaluate the sup difference between
# the true survival function and the empirical estimate
drawSupRegion <- function(n_boundary_coords, n_interior_coords, gridLbs, gridUbs, p, beta) {
    
    # function that we will use our root-finding algorithm for
    pmvtFixedCoord <- function(radius, angle, sigma, prob, xCenter, yCenter) {
        xCoord <- xCenter + radius*cos(angle)
        yCoord <- yCenter + radius*sin(angle)
        exceedanceProb <- pmvt(lower=c(xCoord, yCoord), upper=c(Inf, Inf), sigma=sigma, df=4)
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
                             sigma=matrix(c(1, 0.7, 0.7, 1), nrow=2), prob=p+beta, xCenter=gridLbs[1], yCenter=gridLbs[2])$root
        upper_radius <- uniroot(pmvtFixedCoord, interval=c(0, maxRad), angle=angle,
                             sigma=matrix(c(1, 0.7, 0.7, 1), nrow=2), prob=p-beta, xCenter=gridLbs[1], yCenter=gridLbs[2])$root
        
        radii_per_angle <- seq(lower_radius, upper_radius, length.out=n_interior_coords)
        angles_per_angle <- rep(angle, length.out=n_interior_coords)
        
        radii[((i-1)*n_interior_coords + 1):(i*n_interior_coords)] <- radii_per_angle
        full_angles[((i-1)*n_interior_coords + 1):(i*n_interior_coords)] <- angles_per_angle
        
        
    }

    # convert back to cartesian coordinates
    X1 <- radii*cos(full_angles) + gridLbs[1]
    X2 <- radii*sin(full_angles) + gridLbs[2]
    
    return(data.frame(cbind(X1, X2)))
    
}

# function to to compute the empirical survival function (given observed data) on the sup region
computeSupregEmpsurv <- function(sup_reg, dat) {
    
    sup_reg_xs_unique <- sort(unique(sup_reg$X1))
    sup_reg_ys_unique <- sort(unique(sup_reg$X2))
    full_grid <- expand.grid(X1=sup_reg_xs_unique, X2=sup_reg_ys_unique)
    surv <- fastEmpSurv(full_grid, dat)

    full_surv_res <- data.table(X1=full_grid$X1, X2=full_grid$X2, surv=surv)
    sup_reg <- data.table(sup_reg)

    res <- full_surv_res[sup_reg, on=c('X1', 'X2')]
    
    return(res$surv)
}

# function to simulate monte carlo draws of Z(beta), given a particular beta
parallelizedCode <- function(ind) {
    
    beta <- betas[ind]
    
    region <- drawSupRegion(n_boundary_coords, n_interior_coords, lbs, ubs, p, beta)
    true_survfunc <- pmt(x=-as.matrix(region), mean=rep(0,2), S=matrix(c(1,0.7,0.7,1), 2,2), df=4)
    
    ### compute bhat ###
    Z_betas <- rep(NA, n_monte_carlo)
    for (i in 1:n_monte_carlo) {
        sample_dat <- data.frame(rmvt(n, sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2), df=4))
        sample_survfunc <- computeSupregEmpsurv(region, sample_dat)
        Z_betas[i] <- max(abs(sample_survfunc - true_survfunc))
    }
    
    beta_res <- list()
    beta_res$beta <- beta
    beta_res$zbetas <- Z_betas
    
    saveRDS(beta_res, file=paste0(ppath, 'beta_', ind, '.RData'))
}

### parallelized loop over betas ###
for (i in 1:length(ns)) {
    n <- ns[i]
    npath <- paste0(save_path, 'n', n, '/')
    dir.create(npath)
    
    for (j in 1:length(ps)) {
        p <- ps[j]
        ppath <- paste0(npath, 'p', p, '/')
        dir.create(ppath)

        clust <- makeSOCKcluster(n_cores)
        registerDoSNOW(clust)
        pb <- txtProgressBar(min = 1, max = n_betas, style = 3)
        progress <- function(n) setTxtProgressBar(pb ,n)
        opts <- list(progress = progress)

        foreach(k = 1:n_betas, 
                      .options.snow = opts, 
                      .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt')) %dopar% parallelizedCode(k)

        close(pb)
        stopCluster(clust)
        
    }
    
}
