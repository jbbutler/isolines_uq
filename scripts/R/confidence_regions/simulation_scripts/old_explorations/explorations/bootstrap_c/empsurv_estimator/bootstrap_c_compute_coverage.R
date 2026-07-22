library(dplyr)
library(mvtnorm)
library(mnormt)
library(foreach)
library(doSNOW)
library(parallel)

source('~/isolines_uq/scripts/R/confidence_regions_procedure/auxiliary_scripts/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions_procedure/auxiliary_scripts/utils.R')

ps <- c(0.0001, 0.01, 0.001, 0.1)
save_path <- '/pscratch/sd/j/jbbutler/bootstrap_c_tubes/dist/'
dists <- c('bivt', 'bivgauss')
numCoords <- 200
lvlset_ubs <- c(20, 20)
lvlset_lbs <- c(-2, -2)

n_cores <- 64

parallelizedCode <- function(k) { 

    tubes_sample <- n_res[[k]]
    dat <- tubes_sample$dat
    alphas <- tubes_sample$alphas
    
    true_iso <- isolines[[as.character(p)]]
    est_survfunc <- computeEmpSurvIrregular(true_iso, dat)
    is_covereds <- rep(NA, length(alphas))
    c_alphas <- rep(NA, length(alphas))

    for (l in 1:length(alphas)) {

        c_alpha <- tubes_sample$c_estimates[[as.character(alphas[l])]]
        is_covereds[l] <- all((est_survfunc <= p + c_alpha) & (est_survfunc >= p - c_alpha))
    	c_alphas[l] <- c_alpha

    }
    return(cbind(n, p, alphas, c_alphas, is_covereds))
}

p_coverages <- vector(mode='list', length=length(ps))

for (k in 1:length(dists)) {
    dist <- dists[k]
    dist_path <- paste0(save_path, dist, '/')

    isolines <- list()
    for (i in 1:length(ps)) {
        isolines[[as.character(ps[i])]] <- drawIsoline(dist=dist, numCoords=numCoords, gridUbs=lvlset_ubs, gridLbs=lvlset_lbs, prob=ps[i]) 
    }

    for (i in 1:length(ps)) {

        p <- ps[i]
        p_path <- paste0(dist_path, 'p', p, '/')
        n_paths <- list.files(p_path)

        n_coverages <- vector(mode='list', length=length(n_paths))
    
        for (j in 1:length(n_paths)) {

            n_path <- paste0(p_path, n_paths[j])
            n_res <- readRDS(n_path)
            n <- nrow(n_res[[1]]$dat)

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
    full_coverages <- do.call(rbind.data.frame, p_coverages)
    saveRDS(full_coverages, paste0(dist_path, 'full_coverage_results.RData'))
}




