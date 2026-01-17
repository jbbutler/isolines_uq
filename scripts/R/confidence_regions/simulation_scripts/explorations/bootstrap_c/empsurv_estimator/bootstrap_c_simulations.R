# It appears that beta=0 is the right choice actually: if we assume we are able to get control on the distribution
# of worst case differences specifically over the true isoline, we can construct a confidence tube for it.
# More specifically, if we can estimate the 1-alpha quantile of the distribution of worst case differences over
# the true isoline and do so consistently, then we will have asymptotically exact coverage.
# 
# This script is all about simulating bhats via the bootstrap and computing confidence tubes, and evaluating coverage
#
# Jimmy Butler
# February 2025

library(dplyr)
library(mvtnorm)
library(mnormt)
library(foreach)
library(doSNOW)
library(parallel)

source('~/isolines_uq/scripts/R/confidence_regions_procedure/auxiliary_scripts/utils.R')
source('~/isolines_uq/scripts/R/confidence_regions_procedure/auxiliary_scripts/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions_procedure/auxiliary_scripts/confidenceRegions.R')

save_path <- '/pscratch/sd/j/jbbutler/bootstrap_c_tubes/'

# simulation parameters
ns <- c(1000, 5000, 10000, 50000, 100000)
alphas <- c(0.05, 0.1, 0.01)
ps <- c(0.1, 0.01, 0.001, 0.0001)
dists <- c('bivt', 'bivgauss')

n_cores <- 64

# number of simulations, i.e. number of bootstrap estimates of c_{1-alpha}
n_sims <- 500
# number of bootstrap replicates
B <- 1000

# part of loop that will be parallelized: the computation of c for each simulation
parallelizedCode <- function(ind) {

    dat <- sampling_func(n)
    regions <- drawBaseRegionsBetaZero(dat, alphas, p, B, lbs)

    return(regions)
}

# loop over the distributions

for (i in 1:length(dists)) {
    dist <- dists[i]
    sampling_func <- loadSamplingFunction(dist)
    
    distpath <- paste0(save_path, 'dist/')
    dir.create(distpath)
    distpath <- paste0(distpath, dist, '/')
    dir.create(distpath)
    print('Starting new distribution')
    
    for (j in 1:length(ps)) {
        p <- ps[j]
        ppath <- paste0(distpath, 'p', p, '/')
        dir.create(ppath)

        print('Starting new p')
        
        for (k in 1:length(ns)) {
            n <- ns[k]

            print('Starting new n')

            clust <- makeSOCKcluster(n_cores)
            registerDoSNOW(clust)
            pb <- txtProgressBar(min = 1, max = n_sims, style = 3)
            progress <- function(n) setTxtProgressBar(pb ,n)
            opts <- list(progress = progress)

            samp_regions <- foreach(l = 1:n_sims, 
                    .options.snow = opts, 
                    .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt')) %dopar% parallelizedCode(l)

            close(pb)
            stopCluster(clust)

            saveRDS(samp_regions, file=paste0(ppath, 'n', n, '_tubes.RData'))
        }
    }
}