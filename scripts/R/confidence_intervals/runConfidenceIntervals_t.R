library(foreach)
library(doSNOW)
library(parallel)
library(bootstrap)

source('/global/u1/j/jbbutler/isolines_uq/scripts/R/confidenceIntervals.R')

n_cores <- 64

runs <- vector(mode = 'list', length = 1)

#runs[[1]] <- c(2500, 1000, 0.05, 2, 4.50, 1)
runs[[1]] <- c(2500, 1000, 0.05, 2, 0, 2)
#runs[[3]] <- c(2000, 1000, 0.05, 2, 4.50, 1)


for (k in 1:length(runs)) {

n <- runs[[k]][1]
n_boot <- runs[[k]][2]
alpha <- runs[[k]][3]
isoline_num <- runs[[k]][4]
fixed_coord <- runs[[k]][5]
ax <- runs[[k]][6]

asympIndep <- FALSE
n_iter <- 1000

parallelizedCode <- function(ind, n, n_boot) {

    path <- '/global/cscratch1/sd/jbbutler/sims/bootstrap_experiments/bivariate_t/'

    dat <- readRDS(paste0(path, n, 'sims/', n_boot, 'boot/', ind, '_', n,'bigger_bivt_bootstrap.RData'))

    boots <- dat[[2]]
    orig <- dat[[1]]

    res <- bcaMethod(boots, orig, alpha, isoline_num, fixed_coord, ax, asympIndep)
    return(res)
}

start <- proc.time()

print('Making cluster')
clust <- makeSOCKcluster(n_cores)
registerDoSNOW(clust)
pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
print('Finished cluster')
confints <- foreach(i = 1:n_iter,
                    .packages = c('mvtnorm', 'dplyr', 'MASS', 'ismev', 'evd', 'ks', 'bootstrap'),
                    .options.snow = opts, .errorhandling = 'pass') %dopar% {
    parallelizedCode(i, n, n_boot)
}

close(pb)
stopCluster(clust)

print('hi')

print(proc.time() - start)

prefix <- paste0(n, 'n_', n_boot, 'boot_', fixed_coord*100, '_', ax, '_', isoline_num)
saveRDS(confints, file = paste0('/global/homes/j/jbbutler/isolines_uq/results/', prefix, '_bivt_confints.RData'))
}

