# Script to compute coverage for confidence tubes in the style of Mammen and Polonik (2013)
# ADAPTED FOR A RAW-SCALE HIDDEN REGULAR VARIATION (HRV) SANITY CHECK
# Bypasses marginal transformation by exploiting symmetric heavy-tailed t4 margins.
# Estimates the interior tail index (xi_int = eta * xi) directly via pmin(X, Y).
# Strictly limits evaluation to the interior of the positive orthant.

set.seed(45678)

library(dplyr)
library(mvtnorm)
library(mnormt)
library(foreach)
library(doSNOW)
library(parallel)
library(argparse)
library(matrixStats)

source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/coverageEvaluation.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')

# parsing command line specified parameters
parser <- ArgumentParser(description = "Args for Raw Asymmetric Fractional Log HRV bias bounding simulations.")
parser$add_argument("--save_full_path", type = "character", default = NULL, required = FALSE,
                    help = "Optional file path to save full details about simulation (data sample, etc.).")
parser$add_argument("--n_cores", type = "integer", default = 64, 
                    help = "Number of cores for parallel processing.")
parser$add_argument("--save_df_path", type = "character", default = NULL, 
                    help = "Optional file path to save just the csv of coverage results for each simulation.")

args <- parser$parse_args()

if (is.null(args$save_full_path) && is.null(args$save_df_path)) {
    stop("Error: Both --save_full_path and --save_df_path are NULL. \n You must provide at least one output destination so results are saved.")
}

# simulation parameters
ns <- c(1000, 5000, 10000, 50000, 100000)
alphas <- c(0.05, 0.1, 0.01)
C <- 5

gamma1 <- 0.5
gamma2 <- 0.67
B <- 5000
n_sims <- 500
n_coords <- 50

# interior pad for angles
angular_pad <- 0.1

pn <- function(n){ C/n } 
dists <- c('gauss_copula_t4')
lbs <- c(0,0)
ubs <- c(200, 200)

n_cores <- args$n_cores
save_full_path <- args$save_full_path

# sampling function for a gaussian copula, t marginals (so, asymptotic independence with same marginal tails)
sample_gauss_copula_t4 <- function(n) {
    Sigma <- matrix(c(1, 0.7, 0.7, 1), 2, 2)
    Z <- mvtnorm::rmvnorm(n, mean = c(0,0), sigma = Sigma)
    U <- pnorm(Z)
    X <- qt(U, df = 4)
    return(data.frame(X))
}

# ----------------------------------------------------------------------
# Confidence Tube Generator adapted for Raw Scale HRV
# ----------------------------------------------------------------------
computeExtremeRegionBC_HRV_Raw <- function(dat, alphas, p, B, gamma1=1/2, gamma2=2/3, ncoords=50, lbs=c(0,0), progress_bar=FALSE, verbose=FALSE) {
    
    n_dat <- nrow(dat)
    q1 <- n_dat^(-gamma1)
    q2 <- n_dat^(-gamma2)

    if (verbose) print('Estimating interior tail index (xi_int) via Hill on pairwise minima.')
    # Raw minima isolates the interior joint tail decay rate (xi_int = eta * xi)
    base_minima <- pmin(dat[,1], dat[,2])
    xi_int_hat <- estimate_xi_hill(base_minima, gamma1)

    if (verbose) print('Evaluating Empirical Gap for raw fractional bias calibration.')
    
    # 1. Strictly Interior Angles
    thetas <- seq(angular_pad, (pi/2) - angular_pad, length.out = ncoords)
    
    # Project raw data strictly along the interior rays
    Z_static <- matrix(NA, nrow = n_dat, ncol = length(thetas))
    for(j in 1:length(thetas)) {
        Z_static[, j] <- pmin(dat[,1]/cos(thetas[j]), dat[,2]/sin(thetas[j]))
    }
    
    # Base empirical radius along each ray at anchor level q1
    rq1 <- colQuantiles(Z_static, probs=1-q1, type=1, drop=TRUE)
    
    # Extrapolate using the interior tail index
    rp_radii <- rq1 * (q1 / p)^xi_int_hat
    rq2_radii <- rq1 * (q1 / q2)^xi_int_hat
    
    isoline_p <- cbind(rp_radii * cos(thetas), rp_radii * sin(thetas))
    isoline_q2 <- cbind(rq2_radii * cos(thetas), rq2_radii * sin(thetas))
    
    # 2. Evaluate Empirical Survival Function exactly on the q2 HRV isoline
    hit_q2 <- (outer(dat[,1], isoline_q2[,1], '>') & outer(dat[,2], isoline_q2[,2], '>')) * 1.0
    emp_probs_q2 <- colMeans(hit_q2)
    
    # 3. Calculate Fractional Structural Bias
    delta_P <- emp_probs_q2 - q2
    C_theta <- delta_P / (q2 * log(q1 / q2))
    beta_frac_theta <- C_theta * log(q1 / p)
    
    # 4. Extract ASYMMETRIC Global Supremums
    beta_frac_pos <- max(0, max(beta_frac_theta))
    beta_frac_neg <- max(0, -min(beta_frac_theta))

    if (verbose) print('Estimating pure log-variance parameters via Centered Bootstrap.')

    hit_static <- (outer(dat[,1], isoline_p[,1], '>') & 
                   outer(dat[,2], isoline_p[,2], '>')) * 1.0

    qn <- n_dat^(-gamma1)
    boot_probs_mat <- matrix(NA, nrow=B, ncol=ncoords)

    if (progress_bar) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (k in 1:B) {
        idx <- sample.int(n_dat, n_dat, replace=TRUE)
        
        # A. RE-ESTIMATE INTERIOR TAIL INDEX FOR THIS DRAW
        boot_base_minima <- base_minima[idx]
        xi_int_boot <- estimate_xi_hill(boot_base_minima, gamma1)
        
        # Tabulate trick for memory performance
        w <- tabulate(idx, nbins = n_dat)
        emp_probs <- as.numeric(w %*% hit_static) / n_dat
        
        Z_boot <- Z_static[idx, , drop=FALSE]
        rq <- colQuantiles(Z_boot, probs=1-qn, type=1, drop=TRUE)
        rm(Z_boot) 
        
        # B. HRV Scale Extrapolation
        tail_probs <- ((rq / rp_radii)^(1/xi_int_boot)) * qn
        final_probs <- ifelse(emp_probs >= qn, emp_probs, tail_probs)
        
        boot_probs_mat[k, ] <- log(final_probs)
        if (progress_bar) utils::setTxtProgressBar(pb, k)
    }
    if (progress_bar) close(pb)
    
    rm(Z_static, hit_static) 

    # 5. Bootstrap Centering in Log-Space
    mean_log_boot_probs <- colMeans(boot_probs_mat)
    boot_log_devs <- sweep(boot_probs_mat, MARGIN = 2, STATS = mean_log_boot_probs, FUN = "-")
    
    boot_draws_Wplus <- apply(boot_log_devs, 1, max)
    boot_draws_Wminus <- apply(boot_log_devs, 1, min)

    c_plus_estimates <- list()
    c_minus_estimates <- list()
    raw_c_plus_estimates <- list()
    raw_c_minus_estimates <- list()
    
    for (i in 1:length(alphas)) {
        alpha <- alphas[i]        
        c_plus_log <- as.numeric(quantile(boot_draws_Wplus, probs = 1-(alpha/2)))
        c_minus_log <- as.numeric(quantile(boot_draws_Wminus, probs = (alpha/2)))
        
        raw_c_plus_estimates[as.character(alpha)] <- p * exp(c_plus_log) - p
        raw_c_minus_estimates[as.character(alpha)] <- p * exp(c_minus_log) - p
        
        c_plus_estimates[as.character(alpha)] <- p * exp(c_plus_log + beta_frac_neg) - p
        c_minus_estimates[as.character(alpha)] <- p * exp(c_minus_log - beta_frac_pos) - p
    }

    survFunc <- function(pts) {
        survProb <- vectorizedBlendedSurv(pts, dat, gamma1, xi_int_hat)
        return(survProb)
    }

    res_lst <- list(
        dat = dat,
        p = p,
        gamma1 = gamma1,
        gamma2 = gamma2,
        xi_hat = xi_int_hat,
        thetas = thetas,
        r_naive = rp_radii,
        beta_frac_pos = beta_frac_pos,
        beta_frac_neg = beta_frac_neg,
        c_plus_estimates = c_plus_estimates,
        c_minus_estimates = c_minus_estimates,
        raw_c_plus_estimates = raw_c_plus_estimates,
        raw_c_minus_estimates = raw_c_minus_estimates,
        survFunc = survFunc,
        transformed = FALSE
    )

    return(res_lst)
}

# loop over the distributions
for (i in 1:length(dists)) {
    dist <- dists[i]
    
    sampling_func <- sample_gauss_copula_t4
    
    if (!is.null(save_full_path)) {
        distpath <- paste0(save_full_path, dist, '/')
        if (!dir.exists(distpath)) dir.create(distpath, recursive = TRUE)
    }
    print(paste0('Starting new distribution: ', dist))
    
    for (k in 1:length(ns)) {
        n <- ns[k]
        p <- pn(n)

        gauss_isoline <- drawIsoline('bivgauss', 50, c(200,200), c(0,0), p, angular_extents=c(angular_pad, (pi/2)-angular_pad))
        copula_isoline <- data.frame(X1=pnorm(gauss_isoline[,1]), X2=pnorm(gauss_isoline[,2]))
        isoline <- data.frame(X1=qt(copula_isoline[,1], df=4), X2=qt(copula_isoline[,2], df=4)) 

        print(paste0('Starting new n: ', n))
        
        n_dirpath <- NULL
        if (!is.null(save_full_path)) {
            n_dirpath <- paste0(distpath, 'n', n, '/')
            if (!dir.exists(n_dirpath)) dir.create(n_dirpath, recursive = TRUE)
        }

        # set up parallel computations
        clust <- makeSOCKcluster(n_cores)
        registerDoSNOW(clust)
        clusterEvalQ(clust, {
                source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')
        })
            
        pb <- txtProgressBar(min = 1, max = n_sims, style = 3)
        progress <- function(iter) setTxtProgressBar(pb, iter)
        opts <- list(progress = progress)

        parallelizedCode <- function(ind) {
            dat <- sampling_func(n)
            
            regions <- computeExtremeRegionBC_HRV_Raw(
                dat=dat, alphas=alphas, p=p, B=B, 
                gamma1=gamma1, gamma2=gamma2, 
                ncoords=n_coords, lbs=lbs
            )

            is_covereds <- evaluateCoverage(regions, isoline)
            regions$is_covereds <- is_covereds
            
            if (!is.null(n_dirpath)) {
                saveRDS(regions, file=paste0(n_dirpath, 'simulation', ind, '.RData'))
            }
            
            res_df <- data.frame(
                dist = dist,
                p = p,
                n = n,
                alpha = names(regions$is_covereds),
                is_covered = unlist(regions$is_covereds, use.names=FALSE),
                c_plus_adj = unlist(regions$c_plus_estimates[names(regions$is_covereds)]),
                c_minus_adj = unlist(regions$c_minus_estimates[names(regions$is_covereds)]),
                c_plus_raw = unlist(regions$raw_c_plus_estimates[names(regions$is_covereds)]),
                c_minus_raw = unlist(regions$raw_c_minus_estimates[names(regions$is_covereds)]),
                beta_frac_pos = regions$beta_frac_pos,
                beta_frac_neg = regions$beta_frac_neg,
                xi_hat = regions$xi_hat  
            )
            
            return(res_df)
        }

        samp_df_list <- foreach(l = 1:n_sims, 
                .options.snow = opts,
                .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt', 'matrixStats')) %dopar% parallelizedCode(l)

        close(pb)
        stopCluster(clust)

        if (!is.null(args$save_df_path)) {
            coverage_df <- do.call(rbind, samp_df_list)

            save_fname <- paste0(args$save_df_path, 'experiment2.2.3.csv')
            file_exists <- file.exists(save_fname)
                
            write.table(coverage_df, 
                        file = save_fname, 
                        sep = ",", 
                        row.names = FALSE, 
                        col.names = !file_exists, 
                        append = file_exists)     
        }
    }
}