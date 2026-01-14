# Script to compute Monte Carlo realizations of the empirical process, the 1-alpha quantile of which defines c_{1-alpha}
# Do so for various probability distributions with different isoline exceedance probabilities p, using Monte Carlo
#
# Jimmy Butler
# February 7, 2025

library(dplyr)
library(mvtnorm)
library(mnormt)
library(foreach)
library(doSNOW)
library(parallel)

source('~/isolines_uq/scripts/R/auxiliary_scripts/utils.R')
source('~/isolines_uq/scripts/R/auxiliary_scripts/distributionIsolines.R')

# simulation parameters
ns <- c(1000, 5000, 10000, 50000, 100000, 1000000)
ps <- c(0.1, 0.01, 0.05)
dists <- c('bivt', 'bivgauss')
n_monte_carlo <- 10000

# parallelization specs
n_cores <- 64

save_path <- '/pscratch/sd/j/jbbutler/true_cs/'

parallelizedCode <- function(ind) {

    dat <- samplingFunc(n)
    est_survfunc <- computeEmpSurvIrregular(isoline, dat)
    monte_carlo <- max(abs(est_survfunc - p))

    return(monte_carlo)
}

for (i in 1:length(dists)) {
    dist <- dists[i]
    samplingFunc <- loadSamplingFunction(dist=dist)
    
    distpath <- paste0(save_path, dist, '/')
    dir.create(distpath)

    for (j in 1:length(ps)) {
        p <- ps[j]
        isoline <- drawIsoline(dist=dist, numCoords=200, gridUbs=c(20,20), gridLbs=c(-2,-2), prob=p)   
        
        ppath <- paste0(distpath, 'p', p, '/')
        dir.create(ppath)
        
        for (k in 1:length(ns)) {
            n <- ns[k]
            savepath <- paste0(ppath, 'n', n, '.RData')
            
            clust <- makeSOCKcluster(n_cores)
            registerDoSNOW(clust)
            pb <- txtProgressBar(min = 1, max = n_monte_carlo, style = 3)
            progress <- function(n) setTxtProgressBar(pb ,n)
            opts <- list(progress = progress)

            bootstrapped_cs <- foreach(l = 1:n_monte_carlo, 
                    .options.snow = opts, 
                    .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt'), 
                    .combine='c') %dopar% parallelizedCode(l)

            saveRDS(bootstrapped_cs, file=savepath)

            close(pb)
            stopCluster(clust)
        }

    }
}