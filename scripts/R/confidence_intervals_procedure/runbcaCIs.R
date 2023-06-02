library(foreach)
library(doSNOW)
library(parallel)

source('/global/u1/j/jbbutler/isolines_uq/scripts/R/confidenceIntervals.R')

n_cores <- 64
big_dat <- readRDS('/global/u1/j/jbbutler/isolines_uq/sims/
		   bootstrap_experiments/bivariate_t/
		   1000bigger_bivt_bootstrap.RData')
n <- 1000
alpha <- 0.05
isoline_num <- 3
fixed_coord <- 5.63
ax <- 2
asympIndep <- FALSE

n_iter <- length(big_dat)

parallelizedCode <- function(ind) {

    boots <- big_dat[[ind]][[2]]
    orig <- big_dat[[ind]][[1]]

    res <- bcaMethod(boots, orig, alpha, isoline_num, fixed_coord, ax, asympIndep)
    return(c(res$bca_cI, res$per_cI))
}

start <- proc.time()

clust <- makeSOCKcluster(n_cores)
registerDoSNOW(clust)
pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

confints <- foreach(i = 1:10, 
		    .packages = c('mvtnorm', 'dplyr', 'MASS', 'ismev', 'evd', 'ks'), 
		    .options.snow = opts, .errorhandling = 'pass') %dopar% {
    parallelizedCode(i)
}

close(pb)
stopCluster(clust)

print(proc.time() - start)

confint_df <- as.data.frame(do.call(rbind, confints))
colnames(confint_df) <- c('bca_low', 'bca_high', 'per_low', 'per_high')

saveRDS(confint_df, file = paste0('/global/homes/j/jbbutler/isolines_uq/sims/bootstrap_experiments/bivariate_t/', n, '_bivt_confints.RData'))
