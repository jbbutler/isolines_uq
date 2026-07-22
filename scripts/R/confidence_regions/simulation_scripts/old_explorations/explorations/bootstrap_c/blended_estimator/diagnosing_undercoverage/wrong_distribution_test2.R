# An R script to do coverage simulations with transformations of the data using true and unknown marginals
# for the same simulated datasets. Hopefully side-by-side results will make it easier to debug what's going on..

set.seed(12345)

library(foreach)
library(doSNOW)
library(parallel)
library(dplyr)
library(mvtnorm)

source('~/isolines_uq/scripts/R/confidence_regions_procedure/auxiliary_scripts/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions_procedure/auxiliary_scripts/utils.R')
source('~/isolines_uq/scripts/R/confidence_regions_procedure/auxiliary_scripts/confidenceRegions.R')

save_path <- '/pscratch/sd/j/jbbutler/extreme_tubes/diagnosing_undercoverage/wrong_isoline2/'

rdist <- function(n) {
    dat <- data.frame(rmvt(n=n, sigma=matrix(c(1, 0.7, 0.7, 1), nrow=2), df=1))
    return(dat)
}

n <- 10000
n_datasets <- 5
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
ub <- qt(1-(5/n), df=1) + 100
lvlset_ubs <- c(ub, ub)
lvlset_lbs <- c(0, 0)

# part of loop that will be parallelized: the computation of c for each simulation and evaluation of coverage
parallelizedCode <- function(ind) {

    bivt_dat <- rdist(n)
    trans_dat_X1 <- est_inv_cdf(pt(bivt_dat[,1], df=1), orig_dat[,1], gamma)
    trans_dat_X2 <- est_inv_cdf(pt(bivt_dat[,2], df=1), orig_dat[,2], gamma)
    dat <- data.frame(X1=trans_dat_X1, X2=trans_dat_X2)
    # convert marginals to be draws from a lomax distribution (pareto starting at 0)
    # shape and scale both 1

    # marginal transformation with estimated marginal cdfs (blended empirical and GPD)
    transform_ecdf <- function(pts, dat) {
        transformed_pts <- 1/(1-est_cdf(pts, dat, gamma)) - 1
        return(transformed_pts)
    }

    transformed_dat_X1_ecdf <- transform_ecdf(dat[,1], dat[,1])
    transformed_dat_X2_ecdf <- transform_ecdf(dat[,2], dat[,2])
    transformed_dat_ecdf <- data.frame(X1=transformed_dat_X1_ecdf, X2=transformed_dat_X2_ecdf)

    regions <- drawExtremeRegion(transformed_dat_ecdf, alphas, p, B, gamma, xi)

    # evaluate coverage
    transformed_isoline_ecdf <- data.frame(X1=transform_ecdf(wrong_isoline[,1], dat[,1]), X2=transform_ecdf(wrong_isoline[,2], dat[,2]))
    est_survfunc_ecdf <- apply(transformed_isoline_ecdf, 1, blendedSurvivalFunc, dat=transformed_dat_ecdf, gamma=gamma, xi=xi)
    is_covereds <- list()
    cs <- regions$c_estimates

    for (i in 1:length(alphas)) {
        alpha <- alphas[i]
        c <- cs[[as.character(alpha)]]
        is_covereds[as.character(alpha)] <- all((est_survfunc_ecdf <= p + c) & (est_survfunc_ecdf >= p - c))
    }
    
    regions$is_covereds <- is_covereds
    regions$transform <- transform_ecdf
    regions$trans_dat <- dat
    regions$orig_dat <- orig_dat
    
    saveRDS(regions, file=paste0(orig_sim_path, 'tubes_', ind,'.RData'))
}


#sampling_func <- loadSamplingFunction(dist)
distpath <- paste0(save_path, dist, '/')
dir.create(distpath)

ppath <- paste0(distpath, 'p5_div_n/')
dir.create(ppath)

for (k in 1:n_datasets) {

    orig_sim_path <- paste0(ppath, 'n', n, '_dataset', k, '/')
    dir.create(orig_sim_path)

    print('Starting another dataset.')
    
    p <- pn(n)
    orig_dat <- rdist(n)
    
    isoline <- drawBivtIsoline(numCoords=numCoords, gridUbs=lvlset_ubs, gridLbs=lvlset_lbs, prob=p, df=1)
    wrong_iso_X1 <- est_inv_cdf(pt(isoline[,1], df=1), orig_dat[,1], gamma)
    wrong_iso_X2 <- est_inv_cdf(pt(isoline[,2], df=1), orig_dat[,2], gamma)
    # the "wrong" isoline we are trying to capture
    wrong_isoline <- data.frame(X1=wrong_iso_X1, X2=wrong_iso_X2)

    clust <- makeSOCKcluster(n_cores)
    registerDoSNOW(clust)
    pb <- txtProgressBar(min = 1, max = n_sims, style = 3)
    progress <- function(n) setTxtProgressBar(pb ,n)
    opts <- list(progress = progress)

    foreach(l = 1:n_sims, 
            .options.snow = opts, 
            .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt'),
            .errorhandling='stop') %dopar% parallelizedCode(l)

    close(pb)
    stopCluster(clust)

}