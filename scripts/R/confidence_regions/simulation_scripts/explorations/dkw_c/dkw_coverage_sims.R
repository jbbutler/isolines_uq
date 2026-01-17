# Script to load up base confidence regions and compute their coverages
library(foreach)
library(doSNOW)
library(parallel)
library(mvtnorm)
library(data.table)
source('~/isolines_uq/scripts/R/auxiliary_scripts/distributionIsolines.R')

save_path <- '/pscratch/sd/j/jbbutler/dkw_analysis/'

n_cores <- 64
n_iter <- 1000
ns <- c(1000, 10000, 100000, 500000, 1000000)
distribution <- 'karachi'
numCoords <- 200
ps <- c(0.1, 0.05, 0.01)
alphas <- c(0.05, 0.1, 0.01)

# bounds of grid you used to draw the regions (will be used for drawing true isolines)
lvlset_lbs <- c(40, 0)
lvlset_ubs <- c(150, 100)

if (distribution=='bivt'){
    isolineFunc <- drawBivtIsoline
} else if (distribution=='bivgauss') {
    isolineFunc <- drawBivGaussIsoline
} else if (distribution=='karachi'){
    isolineFunc <- drawKarachiIsoline
}

# create list of isoline basded on desired ps
isolines <- list()
for (i in 1:length(ps)) {
    isolines[[as.character(ps[i])]] <- isolineFunc(numCoords=numCoords, gridUbs=lvlset_ubs, gridLbs=lvlset_lbs, prob=ps[i]) 
}

parallelizedCode <- function(ind) {

    # draw the data, depending on what distribution you're interested in!
    if (distribution == 'bivgauss') {
        dat <- data.frame(rmvnorm(n, mean = rep(0, 2), sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2)))
    } else if (distribution == 'bivt') {
        dat <- data.frame(rmvt(n, sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2), df = 4))
    } else if (distribution == 'karachi') {
        dat <- rKarachiBetaKDE(n=n, b=0.00073, ubs=lvlset_ubs, lbs=lvlset_lbs)
    }

    combos <- expand.grid(alpha=alphas, p=ps)
    res_lst <- vector(mode='list', length=nrow(combos))

    for (i in 1:nrow(combos)) {
        # compute the DKW-derived tube bound
        alpha <- combos$alpha[i]
        p <- combos$p[i]
        isoline <- isolines[[as.character(p)]]
        c <- sqrt((1/(2*n))*log((2*(n+1))/(1-alpha)))

        # evaluate coverage
        findEmpSurv <- function(row) {
        	return(mean((dat[,1] > row[[1]]) & (dat[,2] > row[[2]])))
        }

        # evaluate if each point on the isoline satisfies conditions of being in confidence set
        survs <- apply(isoline, 1, findEmpSurv)
        covered <- all((survs <= p + c) & (survs >= p - c))

        sim_lst <- list()
        sim_lst$n <- n
        sim_lst$alpha <- alpha
        sim_lst$p <- p
        sim_lst$c <- c
        sim_lst$covered <- covered

        res_lst[[i]] <- sim_lst
    }
    
    results_sim <- do.call(rbind.data.frame, res_lst)

    return(results_sim)
}

# simulate across n
total_results <- vector(mode='list', length=length(ns))
for (k in 1:length(ns)) {
    
    n <- ns[k]

    clust <- makeSOCKcluster(n_cores)
    registerDoSNOW(clust)
    pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
    progress <- function(n) setTxtProgressBar(pb ,n)
    opts <- list(progress = progress)
    nresult <- foreach(ind = 1:n_iter,
                   .packages = c('mvtnorm', 'dplyr'),
                   .options.snow = opts,
		   .combine = 'rbind') %dopar% {

        parallelizedCode(ind)
    } 
    total_results[[k]] <- nresult

    close(pb)
    stopCluster(clust)
}


total_results <- do.call(rbind.data.frame, total_results)
saveRDS(total_results, paste0(save_path, distribution, '_dkw_results.RData'))
