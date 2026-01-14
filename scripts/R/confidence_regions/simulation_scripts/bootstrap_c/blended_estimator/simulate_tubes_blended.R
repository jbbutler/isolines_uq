# Script to simulate confidence tubes constructed using the extreme value estimator, assuming regular variation with a known index of regular variation.
#
# Jimmy Butler
# March 2025

set.seed(3434)

library(foreach)
library(doSNOW)
library(parallel)
library(dplyr)
library(mvtnorm)

source('~/isolines_uq/scripts/R/confidence_regions_procedure/auxiliary_scripts/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions_procedure/auxiliary_scripts/utils.R')
source('~/isolines_uq/scripts/R/confidence_regions_procedure/auxiliary_scripts/confidenceRegions.R')

save_path <- '/pscratch/sd/j/jbbutler/extreme_tubes/marginals_unknown_edf_only/'

rdist <- function(n) {
    dat <- data.frame(rmvt(n=n, sigma=matrix(c(1, 0.7, 0.7, 1), nrow=2), df=4))
    return(dat)
}

ns <- c(1000, 5000, 10000, 50000)
dist <- 'copula_tdf4_marginals_tdf4'
alphas <- c(0.05, 0.1, 0.01)
gamma <- 1/2
xi <- 1
pn <- function(n) { 5/n }

# number of cores for parallelization
n_cores <- 64
# number of simulations, i.e. number of bootstrap estimates of c_{1-alpha}
n_sims <- 500
# number of bootstrap replicates
B <- 1000

## Parameters for drawing isolines to check coverage
# number of coords for isoline drawing
numCoords <- 200
ub <- max(qt(1-(5/ns), df=4) + 100)
lvlset_ubs <- c(ub, ub)
lvlset_lbs <- c(0, 0)

# part of loop that will be parallelized: the computation of c for each simulation and evaluation of coverage
parallelizedCode <- function(ind) {

    dat <- rdist(n)
    # convert marginals to be draws from a lomax distribution (pareto starting at 0)
    # shape and scale both 1

    # marginal transformation with estimated marginal cdfs (blended empirical and GPD)
    transform <- function(pts, dat) {
        transformed_pts <- 1/(1-est_cdf(pts, dat, gamma)) - 1
        return(transformed_pts)
    }
    #
    # marginal transformation with known cdfs
    #transform <- function(pts, dat) {
    #    transformed_pts <- 1/(1-pt(pts, df=1)) - 1
    #    return(transformed_pts)
    #}
    # marginal transformation with estimated cdfs (empirical only)
    #transform <- function(pts, dat) {
    #    edf <- ecdf(dat)
    #    transformed_pts <- 1/(1-edf(pts)) - 1
    #    return(transformed_pts)
    #}

    transformed_dat_X1 <- transform(dat[,1], dat[,1])
    transformed_dat_X2 <- transform(dat[,2], dat[,2])
    transformed_dat <- data.frame(X1=transformed_dat_X1, X2=transformed_dat_X2)
    
    regions <- drawExtremeRegion(transformed_dat, alphas, p, B, gamma, xi)

    # evaluate coverage
    transformed_isoline <- data.frame(X1=transform(isoline[,1], dat[,1]), X2=transform(isoline[,2], dat[,2]))
    est_survfunc <- apply(transformed_isoline, 1, blendedSurvivalFunc, dat=transformed_dat, gamma=gamma, xi=xi)
    is_covereds <- list()
    cs <- regions$c_estimates

    for (i in 1:length(alphas)) {
        alpha <- alphas[i]
        c <- cs[[as.character(alpha)]]
        is_covereds[as.character(alpha)] <- all((est_survfunc <= p + c) & (est_survfunc >= p - c))
    }
    regions$is_covereds <- is_covereds
    return(regions)
}


#sampling_func <- loadSamplingFunction(dist)
distpath <- paste0(save_path, dist, '/')
dir.create(distpath)

ppath <- paste0(distpath, 'p5_div_n/')
dir.create(ppath)

for (k in 1:length(ns)) {
    n <- ns[k]
    p <- pn(n)

    isoline <- drawBivtIsoline(numCoords=numCoords, gridUbs=lvlset_ubs, gridLbs=lvlset_lbs, prob=p, df=4)

    print('Starting new n')

    clust <- makeSOCKcluster(n_cores)
    registerDoSNOW(clust)
    pb <- txtProgressBar(min = 1, max = n_sims, style = 3)
    progress <- function(n) setTxtProgressBar(pb ,n)
    opts <- list(progress = progress)

    samp_regions <- foreach(l = 1:n_sims, 
            .options.snow = opts, 
            .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt'),
            .errorhandling='stop') %dopar% parallelizedCode(l)

    close(pb)
    stopCluster(clust)

    saveRDS(samp_regions, file=paste0(ppath, 'n', n, '_tubes.RData'))
}
