library(foreach)
library(doSNOW)
library(parallel)
library(mvtnorm)

source('/global/u1/j/jbbutler/isolines_uq/scripts/R/confidenceRegions.R')
source('/global/u1/j/jbbutler/isolines_uq/scripts/R/utils.R')

proj_p <- 0.01
tol <- 0.0004

asympIndep <- TRUE

setwd('/global/cscratch1/sd/jbbutler/sims/regions/')

all_folders <- list.files(path = 'bases')
specifiers <- c('0.25q', 'bivgauss', '0.05pbase')

n_cores <- 64

parallelizedCode <- function(ind, folder, path, save_folder) {

	base_out <- readRDS(paste0(path, ind, '_', folder, '.RData'))
        proj_out <- projectBaseRegion(base_out$orig_out, proj_p, asympIndep=asympIndep)

        if ('bivgauss' %in% specifiers) {
            isoline <- drawBivGaussIsoline(p = proj_p, tol = tol, grid = proj_out$proj_grid)
	}
	if ('bivt' %in% specifiers) {
	    isoline <- drawBivtIsoline(p = proj_p, tol = tol, grid = proj_out$proj_grid)
	}

	is_covered <- nrow(setdiff(isoline, proj_out$proj_region)) == 0


        cov_lst <- list()
        cov_lst$proj_out <- proj_out
        cov_lst$isoline <- isoline
        cov_lst$is_covered <- is_covered

        saveRDS(cov_lst, file = paste0('projections/', save_folder, '/', ind, '_', save_folder, '.RData'))

}

for (folder in all_folders) {

    folder_specified <- all(as.logical(lapply(specifiers, function(x) grepl(x, folder))))

    if (folder_specified) {

	path <- paste0('bases/', folder, '/')

        identifiers <- strsplit(folder, '_')[[1]]
        fname <- paste(append(identifiers, paste0(proj_p, 'pproj'), after = 3), collapse = '_')
        dir.create(paste0('projections/', fname), showWarnings = FALSE)

        n_iter <- length(list.files(path))
	clust <- makeSOCKcluster(n_cores)
	registerDoSNOW(clust)
        pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
        progress <- function(n) setTxtProgressBar(pb ,n)
        opts <- list(progress = progress)

        foreach(i = 1:n_iter,
                   .packages = c('mvtnorm', 'dplyr', 'ismev', 'ks'),
                   .options.snow = opts,
		   .errorhandling = 'pass') %dopar% {

            parallelizedCode(i, folder, path, fname)

    }

    close(pb)
    stopCluster(clust)

    }
}
