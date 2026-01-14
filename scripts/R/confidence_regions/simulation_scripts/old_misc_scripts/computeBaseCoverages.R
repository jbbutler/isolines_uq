# Script to load up base confidence regions and compute their coverages
library(foreach)
library(doSNOW)
library(parallel)
library(mvtnorm)
library(data.table)
source('~/isolines_uq/scripts/R/auxiliary_scripts/distributionIsolines.R')

###### USER-SPECIFIED PARAMETERS ######

ns <- c(1000, 10000, 100000, 500000, 1000000)
distribution <- 'bivgauss'
numCoords <- 500
load_create_tag <- 'empsurv_smoothboot'
ps <- c(0.1, 0.05, 0.01)
alphas <- c(0.05, 0.1, 0.01)
beta_func_labs <- c('0.5', 'sqrt\\(log\\(n\\)/n\\)')

# bounds of grid you used to draw the regions (will be used for drawing true isolines)
lvlset_lbs <- c(-2, -2)
lvlset_ubs <- c(5, 5)

# bounds of grid on which region bhats were calculated
ubs <- c(5, 5)
lbs <- c(-2, -2)

# grid resolution used to draw regions (not used in computations, solely for bookkeeping)
gticks <- 400

n_cores <- 10
#######################################

# tag to keep track of which grid you used to construct your confidence regions
grid_tag <- paste0(gticks, 'x', gticks, '_on_[', lbs[1], ',', ubs[1], ']x[',
                  lbs[2], ',', ubs[2], ']')

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

# loading and saving paths
load_path <- paste0('/pscratch/sd/j/jbbutler/sims/regions/base_confregs_', 
		    load_create_tag, '/', distribution, '/', grid_tag, '/')
save_path <- paste0('/pscratch/sd/j/jbbutler/sims/regions/base_coverage_', 
                    load_create_tag, '/', distribution, '/', grid_tag, '/')
dir.create(save_path)

# get a list of the possible survfunc simulations to load, all differing by sample size
simulated_survfuncs <- list.files(load_path)
# get the file names that match the desired ns above
str_ns <- as.character(ns)
match_inds <- match(ns, as.integer(do.call(rbind.data.frame, strsplit(simulated_survfuncs, 'n_'))[,1]))
desired_ns <- simulated_survfuncs[match_inds]

# get regular expressions correspdonding to above choices of n, alpha, etc.
# to be used to search for the right files when loading up confregs
ps <- paste(paste0('p', ps), collapse='|')
alphas <- paste(paste0('alpha', alphas), collapse='|')
beta_func_labs <- paste(paste0('beta', beta_func_labs), collapse='|')

parallelizedCode <-  function(ind) {

    # grab desired simulations
    sim_confregs <- readRDS(paste0(sims_path, '/simulation_', ind, '.RData'))
    dat <- sim_confregs$data
    # get rid of the last name, since that's the data and not confidence region
    confreg_names <- names(sim_confregs)
    mask <- grepl(ps, confreg_names) & 
	    grepl(alphas, confreg_names) & 
	    grepl(beta_func_labs, confreg_names)
    sim_confregs <- sim_confregs[mask]
    res <- vector(mode='list', length=length(sim_confregs))
    # determine coverage for each desired confreg, for this simulation
    for (i in 1:length(res)) {

	confreg <- sim_confregs[[i]]
	sim_lst <- list()
        sim_lst$n <- n
        sim_lst$alpha <- confreg$alpha
        sim_lst$p <- confreg$p
        sim_lst$beta_func <- confreg$beta_func
	sim_lst$sim_num <- ind
        
        p_tag <- as.character(sim_lst$p)

        #tube_top <- confreg$tube_top
        #tube_bottom <- confreg$tube_bottom

        isoline <- isolines[[p_tag]]

        # old way of evaluating coverage
	#comp1xs <- CJ(top_x=tube_top[,1], iso_x=isoline[,1], sorted=FALSE)
        #comp1ys <- CJ(top_y=tube_top[,2], iso_y=isoline[,2], sorted=FALSE)
        #comp2xs <- CJ(iso_x=isoline[,1], bott_x=tube_bottom[,1], sorted=FALSE)
        #comp2ys <- CJ(iso_y=isoline[,2], bott_y=tube_bottom[,2], sorted=FALSE)

        #upper_cov <- !any((comp1xs[[1]] <= comp1xs[[2]]) & (comp1ys[[1]] <= comp1ys[[2]]))
        #lower_cov <- !any((comp2xs[[1]] <= comp2xs[[2]]) & (comp2ys[[1]] <= comp2ys[[2]]))
        #covered <- upper_cov & lower_cov

        # new way of determining coverage
	bhat <- confreg$bhat
	p <- confreg$p

	findEmpSurv <- function(row) {
    		return(mean((dat[,1] > row[[1]]) & (dat[,2] > row[[2]])))
	}
	# evaluate if each point on the isoline satisfies conditions of being in confidence set
	survs <- apply(isoline, 1, findEmpSurv)
	covered <- all((survs <= p + bhat) & (survs >= p - bhat))

	sim_lst$covered <- covered

	res[[i]] <- sim_lst

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

#ntag <- paste0('n', paste(ns, collapse='|'))
#ptag <- paste0('p', paste(ps, collapse='|'))
#atag <- paste0('alpha', paste(alphas, collapse='|'))
#btag <- paste0('beta_func', paste(beta_func_labs, collapse='|'))
#save_name <- paste(c(ntag, ptag, atag, btag), collapse='_')
n_tag <- paste(ns, collapse='_')
save_name <- paste0('n', n_tag, 'isolinenumCoords_', numCoords)

total_results <- do.call(rbind.data.frame, total_results)
saveRDS(total_results, paste0(save_path, save_name, '.RData'))

