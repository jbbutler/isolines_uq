# Script to create survival functions, using either the integral estimate of the Gaussian KDE
# or the empirical survival function

library(foreach)
library(doSNOW)
library(parallel)
library(mvtnorm)
library(dplyr)

source('/global/homes/j/jbbutler/isolines_uq/scripts/R/confidenceRegions.R')
source('/global/homes/j/jbbutler/isolines_uq/scripts/R/utils.R')

path <- '/global/cscratch1/sd/jbbutler/sims/regions/empirical_survfuncs/'
set.seed(12345)
n_cores <- 64

ns <- c(1000)
B <- 500
asympIndep <- FALSE
kernelSurv <- FALSE
n_iter <- 250

lb <- 0
ub <- 5
gticks <- 250

grid <- expand.grid(X1 = seq(lb, ub, length.out = gticks),
                    X2 = seq(lb, ub, length.out = gticks))

parallelizedCode <- function(ind) {


    prefix <- paste0(n, 'n_', gticks, 'gticks_', lb, 'lb_', ub, 'ub_')

    if (asympIndep) {
        dat <- data.frame(rmvnorm(n, mean = rep(0, 2), sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2)))
        suffix <- 'bivgauss_'
    } else {
        dat <- data.frame(rmvt(n, sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2), df = 4))
        suffix <- 'bivt_'
    }

    if (kernelSurv) {

        survtag <- 'gausskernel'
	
        fname <- paste0(prefix, suffix, survtag)
        dir.create(path = paste0(path, fname), showWarnings = FALSE)
	dir.create(path = paste0(path, fname, '/simulation_', ind), showWarnings = FALSE)

    	bw <- c(bandwidth.nrd(dat[,1]), bandwidth.nrd(dat[,2]))
    	est_survfunc <- apply(grid, 1, kernSurv, dat=dat, bw=bw)

	saveRDS(est_survfunc, file = paste0(path, fname, '/simulation_', ind, '/simulation_', ind, '_estimate.RData'))
	saveRDS(dat, file = paste0(path, fname, '/simulation_', ind, '/simulation_', ind, '_data.RData'))
	dir.create(path = paste0(path, fname, '/simulation_', ind, '/bootstrap_estimates'), showWarnings = FALSE)

#       boot_survfuncs <- vector(mode='list', length=B)
    
        for (i in 1:B) {
            boot_dat <- dat %>% sample_frac(size=1, replace=TRUE)
            boot_bw <- c(bandwidth.nrd(boot_dat[,1]), bandwidth.nrd(boot_dat[,2]))
            #boot_survfuncs[[i]] <- apply(grid, 1, kernSurv, dat=boot_dat, bw=boot_bw)
	    boot_survfunc <- apply(grid, 1, kernSurv, dat=boot_dat, bw=boot_bw)

	    saveRDS(boot_survfunc, file = paste0(path,
                                   fname, '/simulation_', ind, '/bootstrap_estimates/', i, '_simulation_', ind, '_boot.RData'))
            
        }

    } else {

        survtag <- 'empirical'

        fname <- paste0(prefix, suffix, survtag)
        dir.create(path = paste0(path, fname), showWarnings = FALSE)
        dir.create(path = paste0(path, fname, '/simulation_', ind), showWarnings = FALSE)

        est_survfunc <- fastEmpSurv(grid, dat)

        saveRDS(est_survfunc, file = paste0(path, fname, '/simulation_', ind, '/simulation_', ind, '_estimate.RData'))
	saveRDS(dat, file = paste0(path, fname, '/simulation_', ind, '/simulation_', ind, '_data.RData'))

	boot_survfuncs <- vector(mode='list', length=B)

	for (i in 1:B) {
            boot_dat <- dat %>% sample_frac(size=1, replace=TRUE)
            #boot_survfuncs[[i]] <- apply(grid, 1, empSurv, dat=boot_dat)
	    boot_survfuncs[[i]] <- fastEmpSurv(grid, dat)
        }
    }

    saveRDS(boot_survfuncs, file = paste0(path,
                                   fname, '/simulation_', ind, '/simulation_', ind, '_boots.RData'))
}

for (i in 1:length(ns)) {

    n <- ns[i]

    clust <- makeSOCKcluster(n_cores)
    registerDoSNOW(clust)
    pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
    progress <- function(n) setTxtProgressBar(pb ,n)
    opts <- list(progress = progress)

    foreach(i = 1:n_iter,
                   .packages = c('data.table', 'mvtnorm', 'dplyr', 'MASS', 'purrr'),
                   .options.snow = opts) %dopar% {

        parallelizedCode(i)

    }
    close(pb)
    stopCluster(clust)
}

