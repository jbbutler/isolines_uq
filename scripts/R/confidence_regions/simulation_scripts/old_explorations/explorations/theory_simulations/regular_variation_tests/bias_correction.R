library(mvtnorm)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(foreach)
library(doSNOW)
library(parallel)

estimateIsoR <- function(dat, p, gamma, xi, theta) {

    q <- nrow(dat)^(-gamma)

    empSurvFixedCoord <- function(radius, angle, prob) {
        xCoord <- radius*cos(angle)
        yCoord <- radius*sin(angle)
        exceedanceProb <- mean((dat[,1] > xCoord) & (dat[,2] > yCoord))
        return(exceedanceProb - prob)
    }

    maxRad <- max(10 + sqrt((dat[,1]**2) + (dat[,2]**2)))

    radius_q <- uniroot(empSurvFixedCoord, interval=c(0, maxRad), angle=theta,
                             prob=q)$root
    
    radius_p <- radius_q*((q/p)^(xi))

    return(radius_p)
}

B <- 1000
n <- 10000
p <- 5/n
gamma <- 1/2
xi <- 1/4
theta <- pi/4

n_cores <- 64

n_iter <- 500
est_rs <- rep(NA, n_iter)
bias_rs <- rep(NA, n_iter)

parallelizedCode <- function(i) {

    dat <- data.frame(rmvt(n, sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2), df = 4))
    est_r <- estimateIsoR(dat, p, gamma, xi, theta)
    boot_rs <- rep(NA, B)

    for (j in 1:B) {
        boot_dat <- dat %>% sample_frac(1, replace=TRUE)
        boot_rs[j] <- estimateIsoR(boot_dat, p, gamma, xi, theta)
    }
    
    bias_r <- est_r - mean(boot_rs)

    return(c(est_r, bias_r))

}

clust <- makeSOCKcluster(n_cores)
registerDoSNOW(clust)
pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
progress <- function(n) setTxtProgressBar(pb ,n)
opts <- list(progress = progress)

res <- foreach(l = 1:n_iter, 
            .options.snow = opts, 
            .packages = c('mvtnorm', 'dplyr'), 
            .combine='rbind') %dopar% parallelizedCode(l)

close(pb)
stopCluster(clust)

saveRDS(res, file=paste0('~/isolines_uq/data/theoretical_data/regular_variation_bias_correction/bivariate_t_4_dof.RData'))


