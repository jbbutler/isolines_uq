library(foreach)
library(doSNOW)
library(parallel)
library(mvtnorm)
library(dplyr)

source('/global/homes/j/jbbutler/isolines_uq/scripts/R/confidenceRegions.R')
source('/global/homes/j/jbbutler/isolines_uq/scripts/R/utils.R')


n_cores <- 64

ns <- c(1000)
# additional number of bootstrap replicates
addB <- 100
n_iter <- 500

lb <- 0
ub <- 3
gticks <- 70

grid <- expand.grid(X1 = seq(lb, ub, length.out = gticks),
                    X2 = seq(lb, ub, length.out = gticks))

path <- '/global/cscratch1/sd/jbbutler/sims/regions/true_kde_survfuncs/'

parallelizedCode <- function(ind) {

    surv_lst <- readRDS(paste0(path, fname, '/', ind, '_', prefix, suffix, '.RData'))	
    est_survfunc <- surv_lst$est_survfunc
    dat <- surv_lst$data
    boot_survfuncs <- vector(mode='list', length=addB)

    for (i in 1:addB) {
	boot_dat <- dat %>% sample_frac(size=1, replace=TRUE)
        boot_bw <- c(bandwidth.nrd(boot_dat[,1]), bandwidth.nrd(boot_dat[,2]))
        boot_survfuncs[[i]] <- apply(grid, 1, kernSurv, dat=boot_dat, bw=boot_bw)
    }

    augmented_boot_survfuncs <- append(surv_lst$boot_survfuncs, boot_survfuncs)
    surv_lst$boot_survfuncs <- augmented_boot_survfuncs

    saveRDS(surv_lst, file = paste0(path,
                                   fname, '/', ind, '_', prefix, suffix, '.RData'))

}

for (i in 1:length(ns)) {

    n <- ns[i]

    prefix <- paste0(n, 'n_', gticks, 'gticks_', lb, 'lb_', ub, 'ub_')
    suffix <- 'bivt'
    fname <- paste0(prefix, suffix)

    clust <- makeSOCKcluster(n_cores)
    registerDoSNOW(clust)
    pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
    progress <- function(n) setTxtProgressBar(pb ,n)
    opts <- list(progress = progress)

    foreach(i = 1:n_iter,
                   .packages = c('mvtnorm', 'dplyr', 'MASS', 'purrr'),
                   .options.snow = opts) %dopar% {

        parallelizedCode(i)

    }
    close(pb)
    stopCluster(clust)
}

