# I have a hypothesis that the reason the marginal transformation is giving bad coverage is because
# when we transform the marginals with the estimate, we are generating confidence tubes which are
# correctly covering the isoline of the estimated transformed distribution, but not the original
# distribution. In this script, I'm testing this out by also evaluating coverage against the isoline
# of the distribution whose marginals are the estimated marginals and whose copula is bivariate t.
# This will hopefully let us know the effect of marginal estimation, and how we can remedy it.
#
# Jimmy Butler
# August 2025

set.seed(3434)

library(foreach)
library(doSNOW)
library(parallel)
library(dplyr)
library(mvtnorm)

source('~/isolines_uq/scripts/R/confidence_regions_procedure/auxiliary_scripts/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions_procedure/auxiliary_scripts/utils.R')
source('~/isolines_uq/scripts/R/confidence_regions_procedure/auxiliary_scripts/confidenceRegions.R')

save_path <- '/pscratch/sd/j/jbbutler/extreme_tubes/diagnosing_undercoverage/wrong_isoline/'

rdist <- function(n) {
    dat <- data.frame(rmvt(n=n, sigma=matrix(c(1, 0.7, 0.7, 1), nrow=2), df=1))
    return(dat)
}

ns <- c(10000)
dist <- 'copula_tdf1_marginals_tdf1'
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
ub <- max(qt(1-(5/ns), df=1) + 100)
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

    transformed_dat_X1 <- transform(dat[,1], dat[,1])
    transformed_dat_X2 <- transform(dat[,2], dat[,2])
    transformed_dat <- data.frame(X1=transformed_dat_X1, X2=transformed_dat_X2)
    
    regions <- drawExtremeRegion(transformed_dat, alphas, p, B, gamma, xi)

    # construct the wrong isoline that I hypothesize this procedure is covering with the right probability
    wrong_iso_X1 <- est_inv_cdf(pt(isoline[,1], df=1), dat[,1], gamma)
    wrong_iso_X2 <- est_inv_cdf(pt(isoline[,2], df=1), dat[,2], gamma)
    wrong_isoline <- data.frame(X1=wrong_iso_X1, X2=wrong_iso_X2)
        
    # evaluate coverage of both the true isoline and the wrong isoline
    transformed_isoline <- data.frame(X1=transform(isoline[,1], dat[,1]), X2=transform(isoline[,2], dat[,2]))
    transformed_wrong_iso <- data.frame(X1=transform(wrong_isoline[,1], dat[,1]), X2=transform(wrong_isoline[,2], dat[,2]))
    
    est_survfunc_true <- apply(transformed_isoline, 1, blendedSurvivalFunc, dat=transformed_dat, gamma=gamma, xi=xi)
    est_survfunc_wrong <- apply(transformed_wrong_iso, 1, blendedSurvivalFunc, dat=transformed_dat, gamma=gamma, xi=xi)    
    is_covereds_true <- list()
    is_covereds_wrong <- list()
    cs <- regions$c_estimates

    for (i in 1:length(alphas)) {
        alpha <- alphas[i]
        c <- cs[[as.character(alpha)]]
        is_covereds_true[as.character(alpha)] <- all((est_survfunc_true <= p + c) & (est_survfunc_true >= p - c))
        is_covereds_wrong[as.character(alpha)] <- all((est_survfunc_wrong <= p + c) & (est_survfunc_wrong >= p - c))
    }
    regions$is_covereds_true <- is_covereds_true
    regions$is_covereds_wrong <- is_covereds_wrong
    regions$orig_dat <- dat
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

    isoline <- drawBivtIsoline(numCoords=numCoords, gridUbs=lvlset_ubs, gridLbs=lvlset_lbs, prob=p, df=1)

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

