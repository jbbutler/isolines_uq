library(foreach)
library(doSNOW)
library(parallel)
library(mvtnorm)

source('/global/homes/j/jbbutler/isolines_uq/scripts/R/confidence_regions_procedure/confidenceRegions.R')
source('/global/homes/j/jbbutler/isolines_uq/scripts/R/auxiliary_scripts/utils.R')
source('/global/homes/j/jbbutler/isolines_uq/scripts/R/auxiliary_scripts/karachiTools.R')

path <- '/pscratch/sd/j/jbbutler/sims/regions/base_confregs_empsurv_smoothboot/'

n_cores <- 64

ns <- c(1000, 10000, 100000, 500000, 1000000)
B <- 500
beta_funcs_dict <- list()
beta_funcs_dict[[as.character(1/2)]] <- function(n) {return((1/n)^(1/2))}
beta_funcs_dict[['sqrt(log(n)/n)']] <- function(n) {return((log(n)/n)^(1/2))}
distribution <- 'bivgauss'
n_iter <- 500

lbs <- c(-2, -2)
ubs <- c(5, 5)
ps <- c(0.05, 0.01, 0.1)
alphas <- c(0.01, 0.05, 0.1)

gticks <- 400
grid <- expand.grid(X1 = seq(lbs[1], ubs[1], length.out = gticks),
                    X2 = seq(lbs[2], ubs[2], length.out = gticks))
grid_obj <- list()
grid_obj$grid <- grid
grid_obj$ubs <- ubs
grid_obj$lbs <- lbs

# store simulation results for particular grid settings
subpath <- paste0(path, distribution, '/', gticks, 'x', gticks, '_on_[', lbs[1], ',', ubs[1], ']x[',
		  lbs[2], ',', ubs[2], ']')
# create the directory
dir.create(subpath)

parallelizedCode <- function(ind) {

    if (distribution == 'bivgauss') {
        dat <- data.frame(rmvnorm(n, mean = rep(0, 2), sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2)))
    } else if (distribution == 'bivt') {
        dat <- data.frame(rmvt(n, sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2), df = 4))
    } else if (distribution == 'karachi') {
        dat <- rKarachiBetaKDE(n=n, b=0.00073, ubs=ubs, lbs=lbs)
    }

    base_out <- drawBaseRegions_smoothboot(dat=dat, grid=grid_obj, beta_funcs_dict=beta_funcs_dict, alphas=alphas, ps=ps, B=B)
    prefix <- paste0(n, 'n_', B, 'B_', gticks, 'x', gticks, '_on_[', lbs[1], ',', ubs[1], ']x[',
                  lbs[2], ',', ubs[2], ']_')
    fname <- paste0(prefix, distribution, '_empirical')
    dir.create(path = paste0(subpath, '/', fname), showWarnings = FALSE)

    saveRDS(base_out, file = paste0(subpath, '/', 
				   fname, '/', 'simulation_', ind, '.RData'))

}

for (i in 1:length(ns)) {

    n <- ns[i]

    clust <- makeSOCKcluster(n_cores)
    registerDoSNOW(clust)
    pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
    progress <- function(n) setTxtProgressBar(pb ,n)
    opts <- list(progress = progress)

    foreach(i = 1:n_iter,
		   .packages = c('mvtnorm', 'purrr', 'data.table', 'dplyr', 'ismev', 'ks'),
		   .options.snow = opts) %dopar% {

        parallelizedCode(i)

    }
    close(pb)
    stopCluster(clust)
}
