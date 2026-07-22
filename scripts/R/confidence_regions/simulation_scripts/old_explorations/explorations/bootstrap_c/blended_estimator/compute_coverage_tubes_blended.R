# Script to compute coverage rates for confidence tube construction procedure using extreme value theory
# Takes in tubes whose bounds have already been computed, and determines if true isolines are covered
#
# Jimmy Butler
# March 2025

library(dplyr)
library(mvtnorm)
library(mnormt)
library(foreach)
library(doSNOW)
library(parallel)

source('~/isolines_uq/scripts/R/confidence_regions_procedure/auxiliary_scripts/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions_procedure/auxiliary_scripts/utils.R')

save_path <- '/pscratch/sd/j/jbbutler/extreme_tubes/rv_index_unknown/'
dists <- c('copula_bivtdf1_marginals_tdf1')
p_labs <- c('p5_div_n')

numCoords <- 200
lvlset_ubs <- c(10000, 10000)
lvlset_lbs <- c(0, 0)

gamma <- 1/2
xi <- 1

n_cores <- 64

parallelizedCode <- function(k) { 

    tubes_sample <- n_res[[k]]
    dat <- tubes_sample$dat
    alphas <- tubes_sample$alphas
    
    est_survfunc <- apply(true_iso, 1, blendedSurvivalFunc, dat=dat, gamma=gamma, xi=xi)
    is_covereds <- rep(NA, length(alphas))
    c_alphas <- rep(NA, length(alphas))

    for (l in 1:length(alphas)) {

        c_alpha <- tubes_sample$c_estimates[[as.character(alphas[l])]]
        is_covereds[l] <- all((est_survfunc <= p + c_alpha) & (est_survfunc >= p - c_alpha))
    	c_alphas[l] <- c_alpha

    }
    return(cbind(n, p, alphas, c_alphas, is_covereds))
}

p_coverages <- vector(mode='list', length=length(p_labs))

for (k in 1:length(dists)) {
    dist <- dists[k]
    dist_path <- paste0(save_path, dist, '/')

    for (i in 1:length(p_labs)) {

        p_path <- paste0(dist_path, p_labs[i], '/')
        n_paths <- list.files(p_path)
        n_coverages <- vector(mode='list', length=length(n_paths))
    
        for (j in 1:length(n_paths)) {

            n_path <- paste0(p_path, n_paths[j])
            n_res <- readRDS(n_path)
            p <- n_res[[1]]$p
            n <- nrow(n_res[[1]]$dat)

            true_iso <- drawBivtIsoline(numCoords=numCoords, gridUbs=lvlset_ubs, gridLbs=lvlset_lbs, prob=p, df=1)

            clust <- makeSOCKcluster(n_cores)
            registerDoSNOW(clust)
            pb <- txtProgressBar(min = 1, max = length(n_res), style = 3)
            progress <- function(n) setTxtProgressBar(pb ,n)
            opts <- list(progress = progress)

            n_coverages[[j]] <- foreach(k = 1:length(n_res), 
                .options.snow = opts, 
                .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt'), 
                .combine='rbind') %dopar% parallelizedCode(k)

            close(pb)
            stopCluster(clust)
        }
        p_coverages[[i]] <- do.call(rbind.data.frame, n_coverages)
    }
}
full_coverages <- do.call(rbind.data.frame, p_coverages)
saveRDS(full_coverages, paste0(dist_path, 'full_coverage_results.RData'))