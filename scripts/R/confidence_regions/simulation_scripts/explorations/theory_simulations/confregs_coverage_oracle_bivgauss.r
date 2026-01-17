# Script to simulate coverage rates for the confidence region construction process in Mammen and Polonik (2013),
# assuming we are an oracle which has access to the true survival function.
#
# Jimmy Butler, December 2024

library(mvtnorm)
library(mnormt)
library(dplyr)
library(foreach)
library(doSNOW)
library(parallel)

# loading up helper functions
source('~/isolines_uq/scripts/R/auxiliary_scripts/utils.R')
source('~/isolines_uq/scripts/R/auxiliary_scripts/distributionIsolines.R')

path <- '/pscratch/sd/j/jbbutler/oracle_analysis/bivgauss/'

### simulation specs ###

# specifications for the confidence regions
ps <- c(0.1, 0.01)
alphas <- c(0.01, 0.05, 0.1)
ns <- c(1000, 10000, 50000)
isolines <- list()

for (i in 1:length(ps)) {
    isolines[[as.character(ps[i])]] <- drawBivGaussIsoline(numCoords=200, gridUbs=c(10,10), gridLbs=c(-2,-2), prob=ps[i])
}

# parallelization specifications
n_cores <- 64
# number of simulations for each beta (to determine coverage rate)
n_sims <- 5000

### the actual simulating... ###

# big function to simulate coverage rates for a particular beta (we will parallelize the loop over beta)
parallelizedCode <- function(ind) {
    
    zbetas_file <- readRDS(paste0(path, 'zbetas/n', n, '/p', p, '/beta_', ind, '.RData'))
    beta <- zbetas_file$beta
    zbetas <- zbetas_file$zbetas
    
    bhat <- as.numeric(quantile(zbetas, probs = 1-alpha))
    
    ### compute coverage rate with this bhat ###
    is_covereds <- rep(NA, n_sims)
    
    for (i in 1:n_sims) {
        dat <- data.frame(rmvnorm(n, sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2)))
    
        findEmpSurv <- function(row) {
            return(mean((dat[,1] > row[[1]]) & (dat[,2] > row[[2]])))
        }
        survs <- apply(isolines[[as.character(p)]], 1, findEmpSurv)
        is_covereds[i] <- all((survs <= p + bhat) & (survs >= p - bhat))
    } 
    cov_rate <- mean(is_covereds)
    return(c(beta, bhat, cov_rate))
}

### parallelized loop over betas ###

for (i in 1:length(ns)) {
    
    for (j in 1:length(ps)) {

        for (k in 1:length(alphas)) {
    
            n <- ns[i]
            p <- ps[j]
            alpha <- alphas[k]

            n_betas <- length(list.files(paste0(path, 'zbetas/n', n, '/p', p, '/')))

            clust <- makeSOCKcluster(n_cores)
            registerDoSNOW(clust)
            pb <- txtProgressBar(min = 1, max = n_betas, style = 3)
            progress <- function(n) setTxtProgressBar(pb ,n)
            opts <- list(progress = progress)

            results <- foreach(l = 1:n_betas, 
                          .options.snow = opts, 
                          .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt'), 
                          .combine='rbind') %dopar% parallelizedCode(l)

            results <- data.frame(results)
            colnames(results) <- c('beta', 'bhat', 'covrate')

            close(pb)
            stopCluster(clust)

            saveRDS(results, file=paste0(path, 'covrates/p', p, '_alpha', alpha, '_n', n, '.RData'))

        }
    }
    
}