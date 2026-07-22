# File that takes in base confidence regions and projects them into distribution tails

library(foreach)
library(doSNOW)
library(parallel)
library(mvtnorm)
source('~/isolines_uq/scripts/R/confidence_regions_procedure/confidenceRegions.R')

##### USER SPECIFIED PARAMETERS #####

distribution <- 'bivt'
asympIndep <- FALSE
lbs <- c(0,0)
ubs <- c(15,15)
gticks <- 400
# ps you are projecting from
base_ps <- c(0.005, 0.001, 0.0001)
proj_ps_lst <- list()
proj_ps_lst[[as.character(0.005)]]=c(0.001, 0.0001)
proj_ps_lst[[as.character(0.001)]]=c(0.0001, 0.00001)
proj_ps_lst[[as.character(0.0001)]]=c(0.00001, 0.000001)

n_cores <- 15

####################################

# set up paths to access base confidence regions to project
gtag <- paste0(gticks, 'x', gticks, '_on_[', lbs[1], ',', ubs[1],']x[', lbs[2], ',', ubs[2],']')
base_path <- '/pscratch/sd/j/jbbutler/sims/regions/bivt_projexp/base_tubes/'
base_path <- paste0(base_path, distribution, '/', gtag, '/')
nsize_collections <- list.files(base_path)

# set up directories to save projected uncertainty regions
proj_path <- '/pscratch/sd/j/jbbutler/sims/regions/bivt_projexp/proj_tubes/'
proj_path <- paste0(proj_path, distribution, '/')
dir.create(proj_path)
proj_path <- paste0(proj_path, gtag, '/')
dir.create(proj_path)


# define function to be run in parallel
# for a particular sample size, will load and project associated confidence regions of each sim
# all in parallel
parallelizedCode <- function(ind, nsize_collection, proj_path_n) {

    sim_load_path <- paste0(base_path, nsize_collection, '/simulation_', ind, '.RData')
    base_confreg_list <- readRDS(sim_load_path)
    # grab desired base confregs to project to
    proj_confreg_list <- list()
    for (base_p in base_ps) {

    	base_confreg_names <- grep(paste0('p', base_p), names(base_confreg_list), value=TRUE)
    	proj_ps <- proj_ps_lst[[as.character(base_p)]]
    	for (base_name in base_confreg_names) {

	    base_confreg <- base_confreg_list[[base_name]]
            #proj <- projectBaseRegion(base_out=base_confreg, proj_p=proj_p, asympIndep=asympIndep)

            proj_regions <- list()
        
            for (proj_p in proj_ps) {
	        projTube <- list()
	        projTube$top <- (base_confreg$tube_top)*((base_confreg$p/proj_p)**(1/4))
	        projTube$bottom <- (base_confreg$tube_bottom)*((base_confreg$p/proj_p)**(1/4))
                proj_regions[[paste0('proj_p', proj_p)]] <- projTube
	    }

	    proj_out <- list()
	    proj_out$proj_regions <- proj_regions
	    proj_out$n <- nrow(base_confreg$data)
	    proj_out$alpha <- base_confreg$alpha
	    proj_out$base_p <- base_confreg$p
	    proj_out$proj_p <- proj_ps
	    proj_out$beta_func <- base_confreg$beta_func
	    proj_out$sim_num <- ind

	    proj_confreg_list[[base_name]] <- proj_out

        }
    }

    sim_save_path <- paste0(proj_path_n, 'simulation_', ind, '.RData')
    saveRDS(proj_confreg_list, file=sim_save_path) 

}

for (nsize_collection in nsize_collections) {

    num_sims <- length(list.files(paste0(base_path, nsize_collection, '/')))
    proj_path_n <- paste0(proj_path, nsize_collection, '/')
    dir.create(proj_path_n)

    clust <- makeSOCKcluster(n_cores)
    registerDoSNOW(clust)
    pb <- txtProgressBar(min = 1, max = num_sims, style = 3)
    progress <- function(n) setTxtProgressBar(pb ,n)
    opts <- list(progress = progress)

    # do parallelized loop
    foreach(i = 1:num_sims,
                   .packages = c('mvtnorm', 'dplyr', 'ismev', 'ks'),
                   .options.snow = opts) %dopar% {

            parallelizedCode(i, nsize_collection, proj_path_n)

    }
    
    close(pb)
    stopCluster(clust)

}
