# Script to load up base confidence regions and compute their coverages
library(foreach)
library(doSNOW)
library(parallel)
library(mvtnorm)
library(data.table)
source('~/isolines_uq/scripts/R/auxiliary_scripts/distributionIsolines.R')

###### USER-SPECIFIED PARAMETERS ######

ns <- c(1000, 10000, 50000, 100000)
distribution <- 'bivt'
numCoords <- 500
load_create_tag <- 'create_tubes'
base_ps <- c(0.001, 0.005, 0.0001)
alphas <- c(0.05, 0.1, 0.01)
beta_func_labs <- c('0.5')
proj_ps_lst <- list()
proj_ps_lst[[as.character(0.005)]]=c(0.001, 0.0001)
proj_ps_lst[[as.character(0.001)]]=c(0.0001, 0.00001)
proj_ps_lst[[as.character(0.0001)]]=c(0.00001, 0.000001)
projs_ps <- c(0.001, 0.0001, 0.00001, 0.000001)

# bounds of grid you used to draw the regions (will be used for drawing true isolines)
lbs <- c(0,0)
ubs <- c(15,15)
# grid resolution used to draw regions (not used in computations, solely for bookkeeping)
gticks <- 400

# extents of grid on which the isolines will be drawn
iso_lbs <- c(0,0)
iso_ubs <- c(35,35)

n_cores <- 15
#######################################

grid_tag <- paste0(gticks, 'x', gticks, '_on_[', lbs[1], ',', ubs[1], ']x[',
                  lbs[2], ',', ubs[2], ']')

if (distribution=='bivt'){
    isolineFunc <- drawBivtIsoline
} else if (distribution=='bivgauss') {
    isolineFunc <- drawBivGaussIsoline
} else if (distribution=='karachi'){
    isolineFunc <- drawKarachiIsoline
}

# create projected isoline
#isoline <- isolineFunc(numCoords=numCoords, gridUbs=iso_ubs, gridLbs=iso_lbs, prob=proj_p)

proj_isolines <- list()
for (p in projs_ps) {
    proj_isolines[[as.character(p)]] <- isolineFunc(numCoords=numCoords, gridUbs=iso_ubs, gridLbs=iso_lbs, prob=p)
}

# loading and saving paths
load_path <- paste0('/pscratch/sd/j/jbbutler/sims/regions/bivt_projexp/proj_tubes/',
                    distribution, '/', grid_tag, '/')
save_path <- paste0('/pscratch/sd/j/jbbutler/sims/regions/bivt_projexp/proj_tubes_results/')
dir.create(save_path)

# get a list of the possible survfunc simulations to load, all differing by sample size
simulated_survfuncs <- list.files(load_path)
# get the file names that match the desired ns above
match_inds <- match(ns, as.integer(do.call(rbind.data.frame, strsplit(simulated_survfuncs, 'n_'))[,1]))
desired_ns <- simulated_survfuncs[match_inds]

# get regular expressions correspdonding to above choices of n, alpha, etc.
# to be used to search for the right files when loading up confregs
base_ps <- paste(paste0('p', base_ps), collapse='|')
alphas <- paste(paste0('alpha', alphas), collapse='|')
beta_func_labs <- paste(paste0('beta', beta_func_labs), collapse='|')

parallelizedCode <-  function(ind) {

    # grab desired simulations
    sim_confregs <- readRDS(paste0(sims_path, '/simulation_', ind, '.RData'))
    confreg_names <- names(sim_confregs)
    mask <- grepl(base_ps, confreg_names) &
            grepl(alphas, confreg_names) &
            grepl(beta_func_labs, confreg_names)
    sim_confregs <- sim_confregs[mask]
    res <- vector(mode='list', length=length(sim_confregs)*length(proj_ps_lst[[1]]))
    # determine coverage for each desired confreg, for this simulation
    ct <- 1
    for (i in 1:length(sim_confregs)) {
        confreg <- sim_confregs[[i]]
        sim_lst <- list()
        sim_lst$n <- n
        sim_lst$alpha <- confreg$alpha
        sim_lst$base_p <- confreg$base_p
        sim_lst$beta_func <- confreg$beta_func
        sim_lst$sim_num <- ind
        
	proj_ps <- proj_ps_lst[[as.character(confreg$base_p)]]
        for (proj_p in proj_ps) {
	    isoline <- proj_isolines[[as.character(proj_p)]]
	    sim_lst$proj_p <- proj_p
	    proj_region <- confreg$proj_regions[[paste0('proj_p', proj_p)]]
	    tube_top <- proj_region$top
	    tube_bottom <- proj_region$bottom

            comp2xs <- CJ(iso_x=isoline[,1], bott_x=tube_bottom[,1], sorted=FALSE)
            comp2ys <- CJ(iso_y=isoline[,2], bott_y=tube_bottom[,2], sorted=FALSE)
            lower_cov <- !any((comp2xs[[1]] <= comp2xs[[2]]) & (comp2ys[[1]] <= comp2ys[[2]]))
            covered <- lower_cov

            if (!is.null(tube_top)) {
	        comp1xs <- CJ(top_x=tube_top[,1], iso_x=isoline[,1], sorted=FALSE)
                comp1ys <- CJ(top_y=tube_top[,2], iso_y=isoline[,2], sorted=FALSE)
	        upper_cov <- !any((comp1xs[[1]] <= comp1xs[[2]]) & (comp1ys[[1]] <= comp1ys[[2]]))
	        covered <- upper_cov & lower_cov
            }

            sim_lst$covered <- covered
	    res[[ct]] <- sim_lst
	    ct <- ct + 1
        }
    }

    results_sim <- do.call(rbind.data.frame, res)

    return(results_sim)
}

total_results <- vector(mode='list', length=length(desired_ns))
for (j in 1:length(desired_ns)) {

    nsim <- desired_ns[j]
    sims_path <- paste0(load_path, '/', nsim)
    n <- as.numeric(strsplit(nsim, 'n_')[[1]][1])

    n_iter <- length(list.files(sims_path))

    clust <- makeSOCKcluster(n_cores)
    registerDoSNOW(clust)
    pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
    progress <- function(n) setTxtProgressBar(pb ,n)
    opts <- list(progress = progress)

    nresult <- foreach(i = 1:n_iter,
                   .packages = c('mvtnorm', 'purrr', 'data.table', 'dplyr', 'ismev', 'ks'),
                   .options.snow = opts,
                   .combine = 'rbind') %dopar% {

        parallelizedCode(i)

    }
    total_results[[j]] <- nresult

    close(pb)
    stopCluster(clust)

}

save_name <- paste0('isolinenumCoords_', numCoords)

total_results <- do.call(rbind.data.frame, total_results)
saveRDS(total_results, paste0(save_path, save_name, '.RData'))
