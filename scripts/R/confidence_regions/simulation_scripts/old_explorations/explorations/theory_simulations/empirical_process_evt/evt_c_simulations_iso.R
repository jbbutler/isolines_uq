# Script to compute Monte Carlo realizations of the empirical process, the 1-alpha quantile of which defines c_{1-alpha}
# Do so for various probability distributions with different isoline exceedance probabilities p, using Monte Carlo
# This is specifically for the extreme-value estimator of the survival function, using regular variation.
#
# Jimmy Butler
# April 14, 2025

library(dplyr)
library(mvtnorm)
library(mnormt)
library(foreach)
library(doSNOW)
library(parallel)

source('~/isolines_uq/scripts/R/confidence_regions_procedure/auxiliary_scripts/utils.R')
source('~/isolines_uq/scripts/R/confidence_regions_procedure/auxiliary_scripts/distributionIsolines.R')

# simulation parameters
## sample sizes to try
ns <- c(1000, 5000, 10000, 50000)
## rate at which we look out into the tail (q_n = n^-gamma)
gamma <- 1/2
## extreme value index of probability distribution
xi <- 1/4
## sequence of extreme probabilities
pn <- function(n) { 5/n }
plab <- '5_div_n'
dists <- c('bivt')
n_sims <- 500
B <- 1000

# parallelization specs
n_cores <- 64

save_path <- '/pscratch/sd/j/jbbutler/extreme_trueiso_c_tubes/'

# function to take care of the parallelization
parallelizedCode <- function(ind) {

    dat <- samplingFunc(n)
    boot_process <- rep(NA, B)

    for (l in 1:B) {
        boot_dat <- dat %>% sample_frac(1, replace = TRUE)
        est_survfunc <- apply(isoline, 1, blendedSurvivalFunc, dat=boot_dat, gamma=gamma, xi=xi)
        boot_process[l] <- max(abs(est_survfunc - p))
    }
    
    return(boot_process)
}

for (i in 1:length(dists)) {
    dist <- dists[i]
    samplingFunc <- loadSamplingFunction(dist=dist)
    
    distpath <- paste0(save_path, dist, '/')
    dir.create(distpath)

    ps <- pn(ns)

    for (j in 1:length(ps)) {
        p <- ps[j]
        n <- ns[j]

        isoline <- drawIsoline(dist, 200, c(50,50), c(0,0), p)
        savepath <- paste0(distpath, 'n', n, '_', 'p', plab, '.RData')
            
        clust <- makeSOCKcluster(n_cores)
        registerDoSNOW(clust)
        pb <- txtProgressBar(min = 1, max = n_sims, style = 3)
        progress <- function(n) setTxtProgressBar(pb ,n)
        opts <- list(progress = progress)

        empirical_process_draws <- foreach(k = 1:n_sims, 
                    .options.snow = opts, 
                    .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt')) %dopar% parallelizedCode(k)

        saveRDS(empirical_process_draws, file=savepath)

        close(pb)
        stopCluster(clust)
    }
}