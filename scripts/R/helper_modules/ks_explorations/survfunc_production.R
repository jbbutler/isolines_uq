# A script to produce datasets and their true survival functions, to be used
# to compare against binned survival functions for different binning sizes and get
# the sampling distributions for the error statistic

library(ks)
library(dplyr)
library(mvtnorm)
library(doSNOW)
library(foreach)
library(parallel)
library(MASS)

path <- '/global/cscratch1/sd/jbbutler/surv_funcs/' 

# kernel survival function from original analysis (gives true survival function)
kernSurv <- function(loc, dat, bw)
{
        p1 <- 1 - pnorm(loc[1], mean = dat[,1], sd = bw[1]/4)
        p2 <- 1 - pnorm(loc[2], mean = dat[,2], sd = bw[2]/4)
        return(mean(p1*p2))
}

# setting up parameters of the grid we are evaluating on
lb <- -5
ub <- 5
gticks <- 250

grid <- expand.grid(X1 = seq(lb, ub, length.out = gticks),
                    X2 = seq(lb, ub, length.out = gticks))

# the sample sizes we will be using
ns <- c(1000, 3000, 5000, 10000, 15000, 20000)
# the number of datasets to generate for each sample size
M <- 100

# number of cores on which to distribute parallel job
n_cores <- 15

# defining a function to be executed across multiple cores in parallel fashion
parallelizedCode <- function(ind) {

    samp_n <- list()
    sample <- data.frame(rmvt(n, sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2), df = 4))
    samp_n$data <- sample
    samp_n$grid <- grid
    bw <- c(bandwidth.nrd(sample[,1]), bandwidth.nrd(sample[,2]))
    surv <- apply(grid, 1, kernSurv, dat=sample, bw=bw)
    samp_n$surv <- surv

    fname <- paste0(n, 'n_', gticks, 'gticks_', 'bivt')
    dir.create(path = paste0(path, fname), showWarnings = FALSE)
    saveRDS(samp_n, file = paste0(path, fname, '/', ind, '_', fname, '.RData'))

}

# running everything

for (i in 1:length(ns)) {
    n <- ns[i]
    clust <- makeSOCKcluster(n_cores)
    registerDoSNOW(clust)
    pb <- txtProgressBar(min = 1, max = M, style = 3)
    progress <- function(n) setTxtProgressBar(pb ,n)
    opts <- list(progress = progress)

    foreach(j = 1:M,
                   .packages = c('mvtnorm', 'MASS'),
                   .options.snow = opts) %dopar% {

        parallelizedCode(j)

    }
    close(pb)
    stopCluster(clust)
}
