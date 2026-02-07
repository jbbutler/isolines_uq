# Script to explore bias and variance in using regular variation projection to estimate extreme probabilities.
#
# Jimmy Butler
# February 2025

library(mvtnorm)
library(dplyr)
library(ggplot2)
library(dplyr)
library(foreach)
library(doSNOW)
library(parallel)

set.seed(12345)

source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')

dist <- 'bivt_copula_df4_marginals_df4'

# defining sampling, survival functions to make our lives easier later..
rdist <- function(n) {
    dat <- data.frame(rmvt(n=n, sigma=matrix(c(1, 0.7, 0.7, 1), nrow=2), df=4))
    return(dat)
}
pdist <- function(point) {
    prob <- pmvt(lower=point, upper=Inf, df=4, sigma=matrix(c(1, 0.7, 0.7, 1), nrow = 2))
    return(prob)
}

# parameters
xi <- 1 # we're doing marginal transformation here
thetas <- c(pi/4, 0, pi/2)
thetas_labs <- c('pi/4', '0', 'pi/2')
gammas <- c(1/4, 1/3, 1/2, 2/3, 3/4)
gammas_labs <- c('1/4', '1/3', '1/2', '2/3', '3/4')


ns <- round(exp(seq(7, 13)))
n_sims <- 10000
n_cores <- 64

# initialize a list to store the data frames from each iteration
results_list <- list()
counter <- 1

# looping through theta, gamma, and n combos
for (i in 1:length(thetas)) {
    theta <- thetas[i]
    theta_lab <- thetas_labs[i]

    for (j in 1:length(gammas)) {
        gamma <- gammas[j]
        gamma_lab <- gammas_labs[j]

        for (k in 1:length(ns)) {
            n <- ns[k]
            p <- 5/n

            # calculate true radius r for this specific n and theta
            # estimate exceedance probability past this point
            dist_lvlset <- function(r) {
                diff <- pdist(c(r*sin(theta), r*cos(theta)))
                return(diff-p)
            }
            r <- uniroot(dist_lvlset, interval=c(0, 100000))$root
            
            # parallel code: monte carlo loop is parallelized
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

            clust <- makeSOCKcluster(n_cores)
            registerDoSNOW(clust)
            #pb <- txtProgressBar(min = 1, max = n_sims, style = 3)
            #progress <- function(n) setTxtProgressBar(pb ,n)
            #opts <- list(progress = progress)
            
            surv_estimates <- foreach(l = 1:n_sims, 
                .packages = c('mvtnorm'),
                .combine='rbind') %dopar% parallelizedCode(l)
            
            #close(pb)
            stopCluster(clust)

            # calculate raw metrics for this specific n, theta, gamma combination
            ext_mse_val <- mean((surv_estimates[,2] - p)**2)
            ext_bias_val <- mean(surv_estimates[,2]) - p
            emp_mse_val <- mean((surv_estimates[,1] - p)**2)
            emp_bias_val <- mean(surv_estimates[,1]) - p

            # store in a temporary 1-row dataframe
            temp_df <- data.frame(
                n = n,
                theta = theta,
                theta_lab = theta_lab,
                gamma = gamma,
                gamma_lab = gamma_lab,
                radius = r,
                target_p = p,
                ext_mse = ext_mse_val,
                ext_bias = ext_bias_val,
                emp_mse = emp_mse_val,
                emp_bias = emp_bias_val
            )
            
            # add to list and increment
            results_list[[counter]] <- temp_df
            counter <- counter + 1
            
            print(paste("Finished: Theta =", theta_lab, "| Gamma =", gamma_lab, "| n =", n))
        }
    }
}

# combine all list elements into one master dataframe
final_results <- bind_rows(results_list)

# Compute the adjusted columns
final_results <- final_results %>%
  mutate(
    adjusted_ext_bias = ext_bias * n,
    adjusted_emp_bias = emp_bias * n,
    adjusted_ext_mse = ext_mse * (n^2),
    adjusted_emp_mse = emp_mse * (n^2)
  )

# saving
out_dir <- paste0('~/isolines_uq/data/theoretical_data/regular_variation_probability_tests/rv_index_unknown')
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

write.csv(final_results, 
          file = paste0(out_dir, '/', dist, '_all_results.csv'), 
          row.names = FALSE)