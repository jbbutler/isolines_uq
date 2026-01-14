library(mvtnorm)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(foreach)
library(doSNOW)
library(parallel)

source('/global/homes/j/jbbutler/isolines_uq/scripts/R/confidence_regions_procedure/auxiliary_scripts/utils.R')

dist <- 'bivt_copula_df1_marginals_df1'
# make the blended survival function estimator, given a dataset, a desired point, extreme value index, and a rate of exponential decay
rdist <- function(n) {
    dat <- data.frame(rmvt(n=n, sigma=matrix(c(1, 0.7, 0.7, 1), nrow=2), df=1))
    return(dat)
}

pdist <- function(point) {
    prob <- pmvt(lower=point, upper=Inf, df=1, sigma=matrix(c(1, 0.7, 0.7, 1), nrow = 2))
    return(prob)
}
# 1/n^gamma rate of going out into tail for empirical survival function estimation
gamma <- 1/2
# extreme value index (or 1/index of regular variation), assuming it's known
xi <- 1
# different rays through the origin we will extend out from
thetas <- c(pi/4)
thetas_labs <- c('pi/4')
res <- vector(mode='list')

# sample sizes
ns <- round(exp(seq(7, 14)))
n_sims <- 10000
# parallelization specs
n_cores <- 64

# function to run in parallel
parallelizedCode <- function(ind) {
    dat <- rdist(n)

    transform <- function(pts, dat) {
            transformed_pts <- 1/(1-est_cdf(pts, dat, gamma)) - 1
            return(transformed_pts)
    }

    transformed_dat_X1 <- transform(dat[,1], dat[,1])
    transformed_dat_X2 <- transform(dat[,2], dat[,2])
    transformed_dat <- data.frame(X1=transformed_dat_X1, X2=transformed_dat_X2)

    pt <- r*c(cos(theta), sin(theta))
    emp_surv <- mean((dat[,1] > pt[1]) & (dat[,2] > pt[2]))
    
    transformed_pt <- c(transform(pt[1], dat[,1]), transform(pt[2], dat[,2]))
    ext_surv <- blendedSurvivalFunc(transformed_pt, transformed_dat, gamma, xi)
    return(c(emp_surv, ext_surv))
}

for (i in 1:length(thetas)) {
    theta <- thetas[i]

    ext_mses <- rep(NA, length(ns))
    ext_biases <- rep(NA, length(ns))
    emp_mses <- rep(NA, length(ns))
    emp_biases <- rep(NA, length(ns))
    radii <- rep(NA, length(ns))

    for (j in 1:length(ns)) {
        n <- ns[j]
        p <- 5/n

        dist_lvlset <- function(r) {
            diff <- pdist(c(r*sin(theta), r*cos(theta)))
            return(diff-p)
        }

        r <- uniroot(dist_lvlset, interval=c(0, 100000))$root
        
        clust <- makeSOCKcluster(n_cores)
        registerDoSNOW(clust)
        pb <- txtProgressBar(min = 1, max = n_sims, style = 3)
        progress <- function(n) setTxtProgressBar(pb ,n)
        opts <- list(progress = progress)

        surv_estimates <- foreach(l = 1:n_sims, 
            .options.snow = opts, 
            .packages = c('mvtnorm'), 
            .combine='rbind') %dopar% parallelizedCode(l)

        close(pb)
        stopCluster(clust)

        ext_mses[j] <- mean((surv_estimates[,2] - p)**2)
        ext_biases[j] <- mean(surv_estimates[,2]) - p
        emp_mses[j] <- mean((surv_estimates[,1] - p)**2)
        emp_biases[j] <- mean(surv_estimates[,1]) - p
        radii[j] <- r

    }
    
    plt_df_wider <- data.frame(ext_mses, ext_biases, emp_mses, emp_biases, radii, ns)
    plt_df_wider$adjusted_ext_bias <- plt_df_wider$ext_biases*(plt_df_wider$ns)
    plt_df_wider$adjusted_emp_bias <- plt_df_wider$emp_biases*(plt_df_wider$ns)
    plt_df_wider$adjusted_ext_mse <- plt_df_wider$ext_mses*((plt_df_wider$ns)**2)
    plt_df_wider$adjusted_emp_mse <- plt_df_wider$emp_mses*((plt_df_wider$ns)**2)

    res[[thetas_labs[i]]] <- plt_df_wider

}

dir.create(paste0('~/isolines_uq/data/theoretical_data/regular_variation_probability_tests/rv_index_unknown/gamma_', gamma))

saveRDS(res, file=paste0('~/isolines_uq/data/theoretical_data/regular_variation_probability_tests/rv_index_unknown/gamma_', gamma, '/', dist, '.RData'))
