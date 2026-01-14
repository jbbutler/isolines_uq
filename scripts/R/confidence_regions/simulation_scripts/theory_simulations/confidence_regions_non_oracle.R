# Script to simulate coverage rates for the confidence region construction process in Mammen and Polonik (2013),
# assuming no oracle (i.e., we don't know the underlying distribution so we must use bootstrap techniques to get bhat)
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
source('~/isolines_uq/scripts/R/auxiliary_scripts/distributionIsolines.R')

### simulation specs ###

# specifications for the survival function grid
lbs <- c(-2,-2)
ubs <- c(10, 10)
n_boundary_coords <- 70
n_interior_coords <- 4

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

# parallelization specifications (parallelizing over the different draws of data for each beta)
n_cores <- 64

# number of simulations for each beta (to determine coverage rate)
n_sims <- 256
# number of bootstrap replicates for each beta
B <- 500

### helper functions **PUT THIS IN A SEPARATE SCRIPT EVENTUALLY** ###

# function that draws a grid of points over which we wish to evaluate the sup difference between
# the true survival function (in bootstrap world) and the empirical estimate (bootstrap survival function in bootstrap world)
drawSupRegion <- function(dat, n_boundary_coords, n_interior_coords, gridLbs, gridUbs, p, beta) {
    
    # function that we will use our root-finding algorithm for
    empSurvFixedCoord <- function(radius, angle, prob, xCenter, yCenter) {
        xCoord <- xCenter + radius*cos(angle)
        yCoord <- yCenter + radius*sin(angle)
        exceedanceProb <- mean((dat[,1] > xCoord) & (dat[,2] > yCoord))
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
        lower_radius <- uniroot(empSurvFixedCoord, interval=c(0, maxRad), angle=angle,
                            prob=p+beta, xCenter=gridLbs[1], yCenter=gridLbs[2])$root
        upper_radius <- uniroot(empSurvFixedCoord, interval=c(0, maxRad), angle=angle,
                             prob=p-beta, xCenter=gridLbs[1], yCenter=gridLbs[2])$root
        
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

### the actual simulating... ###

# big function to simulate bhat for a particular random draw of data (we will parallelize the loop over the draws of data)
# bhat computed using bootstrap (i.e. no knowledge of underlying distribution)
parallelizedCode <- function(ind) {
    
    sample_dat <- data.frame(rmvt(n, sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2), df = 4))
    region <- drawSupRegion(sample_dat, n_boundary_coords, n_interior_coords, lbs, ubs, p, beta)

    emp_survfunc <- computeSupregEmpsurv(region, sample_dat)
    
    ### compute bhat ###
    Z_betas <- rep(NA, B)
    for (i in 1:B) {
        boot_dat <- sample_dat %>% sample_frac(1, replace = TRUE)
        boot_survfunc <- computeSupregEmpsurv(region, boot_dat)
        Z_betas[i] <- max(abs(emp_survfunc - boot_survfunc))
    }
    bhat <- as.numeric(quantile(Z_betas, probs = 1-alpha))
    return(bhat)
}

### parallelized loop over betas ###
for (i in 1:length(ns)) {
    n <- ns[i]
    
    for (j in 1:length(ps)) {
        p <- ps[j]
        bhats_res <- vector(mode='list', length=n_betas)
        
        for (k in 1:n_betas) {
            beta <- betas[k]

            clust <- makeSOCKcluster(n_cores)
            registerDoSNOW(clust)
            pb <- txtProgressBar(min = 1, max = n_sims, style = 3)
            progress <- function(n) setTxtProgressBar(pb ,n)
            opts <- list(progress = progress)

            results <- foreach(l = 1:n_sims, 
                          .options.snow = opts, 
                          .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt'), 
                          .combine='c') %dopar% parallelizedCode(l)

            close(pb)
            stopCluster(clust)
            
            bhats_res[[k]] <- results

        }
        
        res_lst <- list()
        res_lst$betas <- betas
        res_lst$bhats <- bhats_res
        
        saveRDS(res_lst, file=paste0('/pscratch/sd/j/jbbutler/non_oracle_analysis/bivt/p', p, '_alpha', alpha, '_n', n, '.RData'))
    }
    
}
