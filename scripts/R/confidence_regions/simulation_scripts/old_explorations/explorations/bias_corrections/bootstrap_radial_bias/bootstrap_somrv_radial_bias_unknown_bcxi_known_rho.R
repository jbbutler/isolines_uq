# Version of bootstrap_somrv_radial_bias_unknown_xi_known_rho.R that swaps the
# plain Hill estimator for the Caeiro-Gomes-Pestana Corrected Hill (CH)
# estimator. rho is still assumed known (oracle); only xi estimation changes.
# beta is estimated externally at a high level (k1 ~ n^0.999) and treated as
# fixed for the rest of the computation.
#
# Goal: check whether the SOMRV bias estimator better matches the true bias
# of the first-order MRV procedure when the Hill-bias source is removed.

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

parser <- ArgumentParser(description = "SOMRV bias with CH-estimator for xi.")
parser$add_argument("--n_cores", type = "integer", default = 64, 
                    help = "Number of cores for parallel processing.")
parser$add_argument("--save_df_path", type = "character", default = NULL, 
                    help = "File path to save the csv comparing true vs SOMRV biases.")
args <- parser$parse_args()


ns <- c(1000, 5000, 10000, 50000)
C <- 5
rho <- -1/2          # bivt with nu = 4 gives rho ~ -1/2 (assumed known)
gamma1 <- 1/2
gamma2 <- 2/3
gamma_rho <- 0.001   # high-level for beta external estimation

n_coords <- 50
n_monte_carlo <- 2000
B_boot <- 2000
n_sims <- 500

pn <- function(n){ C/n } 
dists <- c('bivt')
lbs <- c(0,0)
ubs <- c(200, 200)

n_cores <- args$n_cores


estimateSomrvBootstrapBias_CH <- function(dat, p, q1, q2, B, n_coords, 
                                          gamma1, gamma2, rho, gamma_rho, lbs) {
    
    base_radii <- sqrt(dat[,1]^2 + dat[,2]^2)
    
    # Estimate rho via FAGH at high level for the CH estimator only.
    # The oracle `rho` argument is preserved for use in M (the SOMRV multiplier).
    rho_hat_for_CH <- estimate_rho_FAGH(base_radii, gamma_rho = gamma_rho)
    
    # Estimate beta ONCE on the original sample at high level using rho_hat;
    # treat as fixed for the rest of the computation.
    base_beta_hat <- estimate_beta_GM(base_radii, rho_hat = rho_hat_for_CH, 
                                       gamma_rho = gamma_rho)
    
    # CH estimate of xi for the base sample, using the FAGH rho_hat
    base_ch <- estimate_xi_CH(base_radii, gamma = gamma1, 
                              rho_hat = rho_hat_for_CH,
                              beta_hat = base_beta_hat)
    base_xi_hat <- base_ch$xi_ch
    
    # Generate empirical q2 isoline
    isoline_q2_emp <- drawEmpiricalIsoline(dat=dat, n_coords=n_coords, 
                                           grid_lbs=lbs, p=q2)
    rq2_emp <- sqrt(isoline_q2_emp[,1]^2 + isoline_q2_emp[,2]^2)
    thetas <- atan2(isoline_q2_emp[,2], isoline_q2_emp[,1])
    
    n_dat <- nrow(dat)
    boot_rq2_proj_mat <- matrix(NA, nrow=B, ncol=n_coords)
    
    # Bootstrap loop: project from q1 to q2 using the base xi_hat (CH).
    # This mirrors the original script's behavior of using base_xi_hat in
    # the bootstrap, so the comparison to the Hill version is apples-to-apples.
    for (k in 1:B) {
        idx <- sample.int(n_dat, n_dat, replace=TRUE)
        boot_dat <- dat[idx, , drop=FALSE]
        
        boot_isoline_q2 <- drawExtremeIsoline(dat=boot_dat, p=q2, 
                                              n_coords=n_coords, 
                                              grid_lbs=lbs, gamma=gamma1, 
                                              xi=base_xi_hat)
        boot_rq2_proj_mat[k, ] <- sqrt(boot_isoline_q2[,1]^2 + boot_isoline_q2[,2]^2)
    }
    
    expected_rq2_boot <- colMeans(boot_rq2_proj_mat)
    in_sample_bias_est <- expected_rq2_boot - rq2_emp
    
    # M multiplier uses the ORACLE rho (passed in), not rho_hat_for_CH.
    # This isolates whether CH (with estimated rho) is enough to fix the
    # xi-bias source, while keeping M as clean as possible.
    M <- ((q2/p)^base_xi_hat) * (((q1/p)^rho - 1) / ((q1/q2)^rho - 1))
    
    somrv_bias_est <- M * in_sample_bias_est
    
    return(list(thetas = thetas, 
                rq2_emp = rq2_emp, 
                expected_rq2_boot = expected_rq2_boot,
                in_sample_bias_est = in_sample_bias_est,
                somrv_bias_est = somrv_bias_est,
                base_xi_hat = base_xi_hat,
                base_beta_hat = base_beta_hat,
                base_rho_hat = rho_hat_for_CH,
                base_xi_hill = base_ch$xi_hill))
}

for (i in 1:length(dists)) {
    dist <- dists[i]
    sampling_func <- loadSamplingFunction(dist)
    
    print(paste0('Starting new distribution: ', dist))
    
    for (k in 1:length(ns)) {
        n <- ns[k]
        p <- pn(n)
        q1 <- n^(-gamma1)
        q2 <- n^(-gamma2)

        # ---- True bias of the FIRST-ORDER projection with CH xi-plug-in ----
        # We re-derive the true bias under the CH estimator because the bias
        # of the first-order projection depends on which xi estimator is used.
        isoline <- drawIsoline(dist=dist, numCoords=n_coords, gridUbs=ubs, 
                               gridLbs=lbs, prob=p)
        thetas_true <- atan2(isoline[,2], isoline[,1])
        r_thetas_true <- sqrt(isoline[,1]^2 + isoline[,2]^2)
        r_hat_theta_mat <- matrix(NA, nrow=n_monte_carlo, ncol=length(thetas_true))

        print(paste0('Starting n: ', n, ' | Computing True MC Naive Radial Bias with CH...'))
        for (m in 1:n_monte_carlo) {
            dat_mc <- sampling_func(n)
            mc_radii <- sqrt(dat_mc[,1]^2 + dat_mc[,2]^2)
            
            # CH with FAGH-estimated rho, matching what the SOMRV estimator does
            mc_rho_hat <- estimate_rho_FAGH(mc_radii, gamma_rho = gamma_rho)
            mc_ch <- estimate_xi_CH(mc_radii, gamma = gamma1, 
                                     rho_hat = mc_rho_hat,
                                     gamma_rho = gamma_rho)
            mc_xi_hat <- mc_ch$xi_ch
            
            isoline_ext_estimate <- drawExtremeIsoline(dat=dat_mc, p=p, 
                                                       n_coords=n_coords, 
                                                       grid_lbs=lbs, 
                                                       gamma=gamma1, 
                                                       xi=mc_xi_hat)
            r_hat_theta_mat[m,] <- sqrt(isoline_ext_estimate[,1]^2 + 
                                         isoline_ext_estimate[,2]^2)
        }
        
        true_bias <- colMeans(r_hat_theta_mat) - r_thetas_true

        # ---- Parallel SOMRV bootstrap evaluation -------------------------
        clust <- makeSOCKcluster(n_cores)
        registerDoSNOW(clust)
        clusterEvalQ(clust, {
                source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')
                source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')
        })
            
        pb <- txtProgressBar(min = 1, max = n_sims, style = 3)
        progress <- function(iter) setTxtProgressBar(pb, iter)
        opts <- list(progress = progress)

        print(paste0('Evaluating CH-SOMRV radial bias over ', n_sims, ' simulations...'))
        
        parallelizedCode <- function(ind, true_bias, q1, q2) {
            dat <- sampling_func(n)
            
            somrv_res <- estimateSomrvBootstrapBias_CH(
                dat=dat, p=p, q1=q1, q2=q2, B=B_boot, 
                n_coords=n_coords, gamma1=gamma1, gamma2=gamma2, 
                rho=rho, gamma_rho=gamma_rho, lbs=lbs
            )
            
            res_df <- data.frame(
                dist = dist,
                p = p,
                n = n,
                sim_id = ind,
                theta = somrv_res$thetas,
                true_bias = true_bias,
                somrv_bias_est = somrv_res$somrv_bias_est,
                in_sample_bias_est = somrv_res$in_sample_bias_est,
                xi_hat_ch = somrv_res$base_xi_hat,
                xi_hat_hill = somrv_res$base_xi_hill,
                beta_hat = somrv_res$base_beta_hat,
                rho_hat = somrv_res$base_rho_hat
            )
            
            return(res_df)
        }

        samp_df_list <- foreach(l = 1:n_sims, 
                .options.snow = opts,
                .packages = c('mvtnorm', 'data.table', 'dplyr', 
                              'mnormt', 'matrixStats')) %dopar% 
            parallelizedCode(l, true_bias, q1, q2)

        close(pb)
        stopCluster(clust)

        coverage_df <- do.call(rbind, samp_df_list)
        save_fname <- paste0(args$save_df_path, 'somrv_simulation_CH_xi_truerho.csv')
        file_exists <- file.exists(save_fname)
            
        write.table(coverage_df, 
                    file = save_fname, 
                    sep = ",", 
                    row.names = FALSE, 
                    col.names = !file_exists, 
                    append = file_exists)     
    }
}