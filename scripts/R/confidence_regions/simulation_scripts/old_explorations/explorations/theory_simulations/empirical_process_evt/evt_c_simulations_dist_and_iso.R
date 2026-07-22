# Script to compute Monte Carlo realizations of the empirical process, the 1-alpha quantile of which defines c_{1-alpha}
# Do so for various probability distributions with different isoline exceedance probabilities p, using Monte Carlo
# This is specifically for the extreme-value estimator of the survival function, using regular variation.
#
# Jimmy Butler
# March 31, 2025

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
xi <- 1
## sequence of extreme probabilities
pn <- function(n) { 5/n }
plab <- '5_div_n'
dists <- c('bivt_copula_pareto_margins')
n_monte_carlo <- 10000

# parallelization specs
n_cores <- 64

save_path <- '/pscratch/sd/j/jbbutler/true_cs/'

rdist <- function(n) {
    dat <- data.frame(rmvt(n=n, sigma=matrix(c(1, 0.7, 0.7, 1), nrow=2), df=1))
    margY1 <- (1-pt(dat$X1, df=1))^(-1)
    margY2 <- (1-pt(dat$X2, df=1))^(-1)
    Y <- data.frame(Y1=margY1, Y2=margY2) - 1

    return(Y)
}

# function to take care of the parallelization
parallelizedCode <- function(ind) {

    dat <- rdist(n)
    est_survfunc <- apply(isoline, 1, blendedSurvivalFunc, dat=dat, gamma=gamma, xi=xi)
    monte_carlo <- max(abs(est_survfunc - p))

    return(monte_carlo)
}

for (i in 1:length(dists)) {
    dist <- dists[i]
    
    distpath <- paste0(save_path, dist, '/')
    dir.create(distpath)

    ps <- pn(ns)

    for (j in 1:length(ps)) {
        p <- ps[j]
        n <- ns[j]
        isoline <- drawParetoBivtIsoline(numCoords=200, gridUbs=c(10000, 10000), gridLbs=c(0, 0), prob=p, df=1)
    
        savepath <- paste0(distpath, 'n', n, '_', 'p', plab, '.RData')
            
        clust <- makeSOCKcluster(n_cores)
        registerDoSNOW(clust)
        pb <- txtProgressBar(min = 1, max = n_monte_carlo, style = 3)
        progress <- function(n) setTxtProgressBar(pb ,n)
        opts <- list(progress = progress)

        empirical_process_draws <- foreach(l = 1:n_monte_carlo, 
                    .options.snow = opts, 
                    .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt'), 
                    .combine='c') %dopar% parallelizedCode(l)

        saveRDS(empirical_process_draws, file=savepath)

        close(pb)
        stopCluster(clust)
    }
}