simulation_36.RData# File to compuate confidence regions for simuluations of a pre-specified set of sample sizes from particular distributions
# Does so by loading up premade simulated datasets, along with their empirical and bootstrap survival functions

library(foreach)
library(doSNOW)
library(parallel)
library(mvtnorm)
library(sf)
library(sfheaders)

source('/global/homes/j/jbbutler/isolines_uq/scripts/R/confidence_regions_procedure/confidenceRegions.R')
source('/global/homes/j/jbbutler/isolines_uq/scripts/R/utils.R')

n_cores <- 64

# choose sample sizes whose survival functions you wish to load up
ns <- c(10000, 15000, 20000, 25000, 30000)
# path where you are looking for survival functions to load (all data from same distribution)
suffix <- 'bivt'
load_path <- paste0('/global/cscratch1/sd/jbbutler/sims/regions/empirical_survfuncs/', suffix)
save_path <- paste0('/global/cscratch1/sd/jbbutler/sims/regions/base_confregs_load/', suffix, '/')
# get a list of the possible survfunc simulations to load, all differing by sample size
simulated_survfuncs <- list.files(load_path)
# get the file names that match the desired ns above
match_inds <- match(ns, as.integer(do.call(rbind.data.frame, strsplit(simulated_survfuncs, 'n_'))[,1]))
desired_sims <- simulated_survfuncs[match_inds]

# just make sure this is the same as the number of bootstrap replicates of the empirical survfuncs you're loading
# (TO DO) this is not necessasry when loading up simulations, so find a way to not have to specify
B <- 500
# specify all the alphas, ps, and q's you wish to try out
alphas <- c(0.01, 0.05, 0.1)
ps <- c(0.01, 0.05) 
beta_funcs_dict <- list()
beta_funcs_dict[['1/sqrt(n)']] <- function(n) {return(1/sqrt(n))}
n_iter <- 500

# setup the grids, make sure it matches the grid specifications of the survival functions you're plugging in
lb <- 0
ub <- 5
gticks <- 250
grid <- expand.grid(X1 = seq(lb, ub, length.out = gticks), 
                    X2 = seq(lb, ub, length.out = gticks))

parallelizedCode <- function(ind) {

    sim_path <- paste0(sims_path, '/', 'simulation_', ind)
    boots <- readRDS(paste0(sim_path, '/simulation_', ind, '_boots.RData'))
    dat <- readRDS(paste0(sim_path, '/simulation_', ind, '_data.RData'))
    estimate <- readRDS(paste0(sim_path, '/simulation_', ind, '_estimate.RData'))
	
    base_out <- drawBaseRegions(dat=dat, grid=grid, beta_funcs_dict=beta_funcs_dict,
			     	alphas=alphas, ps=ps, B=B, emp_surv=estimate, boot_survs=boots)

    dir.create(path = paste0(save_path, sim), showWarnings = FALSE)

    saveRDS(base_out, file = paste0(save_path, 
				   sim, '/', 'simulation_', ind, '.RData'))

}
for (sim in desired_sims) {

    sims_path <- paste0(load_path, '/', sim)

    clust <- makeSOCKcluster(n_cores)
    registerDoSNOW(clust)
    pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
    progress <- function(n) setTxtProgressBar(pb ,n)
    opts <- list(progress = progress)

    foreach(i = 1:n_iter,
		   .packages = c('mvtnorm', 'purrr', 'data.table', 'dplyr', 'ismev', 'ks', 'sf', 'sfheaders', 'concaveman'),
		   .options.snow = opts) %dopar% {

        parallelizedCode(i)

    }
    close(pb)
    stopCluster(clust)
}
