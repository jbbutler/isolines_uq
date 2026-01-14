library(MASS)
library(ismev)
library(evd)
library(ggplot2)
library(dplyr)
library(mvtnorm)
library(parallel)
library(doSNOW)
library(foreach)

ns <- c(5000)
n_boot_samps <- 1000
n_samps <- 1000
n_cores <- 64

source('/global/homes/j/jbbutler/isolines_uq/scripts/R/orig_isolines.R')

parallelizedCode <- function(n, n_boot, j) {

    boot_out_n_samp <- list()
    gauss_dat <- data.frame(rmvnorm(n, mean = rep(0, 2), sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2)))
    boot_out_n_samp$data <- gauss_dat

    boot_samps <- vector(mode = 'list', length = n_boot)

    for (k in 1:n_boot) {
        boot_out <- list()
        boot_dat <- gauss_dat %>% sample_frac(1, replace = TRUE)
        boot_out_full <- suppressWarnings(xContours(dat = boot_dat, faster = TRUE, asympIndep = TRUE, projContourLevels = c(0.001)))
        boot_out$boot_out_isolines <- c(list(boot_out_full$contourOrig), boot_out_full$projContours)
        boot_out$boot_out_levels <- c(boot_out_full$setup$baseContourLevel, boot_out_full$setup$projContourLevels)

        boot_samps[[k]] <- boot_out
    }
    boot_out_n_samp$bootstrap_info <- boot_samps

    path <- '/global/cscratch1/sd/jbbutler/sims/bootstrap_experiments/bivariate_gauss/'

    saveRDS(boot_out_n_samp, file = paste0(path, n, 'sims/', n_boot, 'boot/', j , '_', n, 'bigger_bivgauss_bootstrap.RData'))

}

for (i in 1:length(ns)) {

    print(paste('Starting n:', i))
    start <- proc.time()

    clust <- makeSOCKcluster(n_cores)
    registerDoSNOW(clust)
    pb <- txtProgressBar(min = 1, max = n_samps, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    boot_out_n_samps <- foreach (j = 1:n_samps, .packages = c('mvtnorm', 'dplyr', 'MASS', 'ismev', 'evd', 'ks'), .options.snow = opts, .errorhandling = 'pass') %dopar% {

	parallelizedCode(ns[i], n_boot_samps, j)

    }	
    close(pb)
    stopCluster(clust)

    print(paste('Finished n:', i))
    print(proc.time() - start)
 
}
