# Script to compute coverage rates of confidence tube procedure using the 
# Conservative Envelope Bounds (Limit as rho -> 0) to account for structural bias.
# Draws data from the bivariate t, evaluates xi dynamically using the Hill estimator,
# and interpolates the bootstrap over the structural bounds.
#
# Jimmy Butler
# March 2026

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
parser <- ArgumentParser(description = "Args for Conservative Envelope Bound survival function simulations.")
parser$add_argument("--save_full_path", type = "character", default = NULL, required = FALSE,
                    help = "Optional file path to save full details about simulation (data sample, etc.).")
parser$add_argument("--n_cores", type = "integer", default = 64, 
                    help = "Number of cores for parallel processing.")
parser$add_argument("--save_df_path", type = "character", default = NULL, 
                    help = "Optional file path to save just the csv of coverage results for each simulation.")

args <- parser$parse_args()

if (is.null(args$save_full_path) && is.null(args$save_df_path)) {
    stop("Error: Both --save_full_path and --save_df_path are NULL. 
          \n You must provide at least one output destination so results are saved.")
}

# simulation parameters
ns <- c(1000, 5000, 10000, 50000)
alphas <- c(0.05, 0.1, 0.01)
C <- 5

gamma1 <- 0.4 
gamma2 <- 0.5
B <- 2000
n_sims <- 500
n_lines <- 15
n_coords <- 50

pn <- function(n){ C/n } 
dists <- c('bivt')
lbs <- c(0,0)
ubs <- c(200, 200)

n_cores <- args$n_cores
save_full_path <- args$save_full_path

# function that constructs confidence tube(s) given a particular dataset 
# using the Conservative Envelope Bound (Limit as rho -> 0) to define the bias region.
computeExtremeRegionBC_Conservative <- function(dat, alphas, p, B, gamma1=1/2, gamma2=2/3, ncoords=50, n_lines=20, lbs=c(0,0), progress_bar=FALSE, verbose=FALSE) {
    
    n_dat <- nrow(dat)
    q1 <- n_dat^(-gamma1)
    q2 <- n_dat^(-gamma2)

    if (verbose) print('Estimating xi via Hill estimator.')
    base_radii <- sqrt(dat[,1]^2 + dat[,2]^2)
    xi_hat <- estimate_xi_hill(base_radii, gamma1)

    if (verbose) print('Constructing Conservative Envelope region.')
    
    # 1. Calculate the First-Order Naive Bound at p
    isoline_p_naive <- drawExtremeIsoline(dat=dat, p=p, n_coords=ncoords, grid_lbs=lbs, gamma=gamma1, xi=xi_hat)
    Rp_first_order <- sqrt(isoline_p_naive[,1]^2 + isoline_p_naive[,2]^2)
    thetas <- atan2(isoline_p_naive[,2], isoline_p_naive[,1])

    # 2. Empirical Anchor at q2
    isoline_q2_emp <- drawEmpiricalIsoline(dat=dat, n_coords=ncoords, grid_lbs=lbs, p=q2)
    Rq2_emp <- sqrt(isoline_q2_emp[,1]^2 + isoline_q2_emp[,2]^2)
    
    # 3. Naive Projection at q2
    isoline_q2_proj <- drawExtremeIsoline(dat=dat, p=q2, n_coords=ncoords, grid_lbs=lbs, gamma=gamma1, xi=xi_hat)
    Rq2_proj <- sqrt(isoline_q2_proj[,1]^2 + isoline_q2_proj[,2]^2)
    
    # 4. In-Sample Gap at q2
    in_sample_gap <- Rq2_proj - Rq2_emp
    
    # 5. Worst-Case Multiplier (rho -> 0 Limit)
    log_ratio <- log(q1 / p) / log(q1 / q2)
    M_limit <- ((q2 / p)^xi_hat) * log_ratio
    
    # 6. Maximal Second-Order Bound at p
    Rp_worst_case <- Rp_first_order - (M_limit * in_sample_gap)
    
    # 7. Construct Envelope Boundaries
    r_lower <- pmin(Rp_first_order, Rp_worst_case)
    r_upper <- pmax(Rp_first_order, Rp_worst_case)

    # 8. Stack points to form the full spatial region to evaluate in the bootstrap
    region_pts_list <- lapply(seq(0, 1, length.out = n_lines), function(weight) {
        r_interp <- r_lower + weight * (r_upper - r_lower)
        cbind(X1 = r_interp * cos(thetas), X2 = r_interp * sin(thetas))
    })
    ext_region_pts <- do.call(rbind, region_pts_list)

    if (verbose) print('Estimating tube size parameters via bootstrap over the Conservative Envelope.')

    # Set up hit matrix for the bootstrap loop over the constructed region
    region_thetas <- atan2(ext_region_pts[,2], ext_region_pts[,1])
    rp <- sqrt(rowSums(ext_region_pts^2))
    
    # SAFELY construct Z_static avoiding NaN generation on the axes
    Z_static <- matrix(NA, nrow = n_dat, ncol = length(region_thetas))
    pos_orthant_mask <- (dat[,1] > 0) & (dat[,2] > 0)
    
    ind_0 <- which(region_thetas == 0)
    if (length(ind_0) > 0) Z_static[, ind_0] <- dat[,1] * pos_orthant_mask
    
    ind_90 <- which(region_thetas == pi/2)
    if (length(ind_90) > 0) Z_static[, ind_90] <- dat[,2] * pos_orthant_mask
    
    inds_no_ax <- which(!(region_thetas == 0 | region_thetas == pi/2))
    if (length(inds_no_ax) > 0) {
        Z_static[, inds_no_ax] <- pmin(
            dat[,1] %o% (1/cos(region_thetas[inds_no_ax])), 
            dat[,2] %o% (1/sin(region_thetas[inds_no_ax]))
        )
    }

    hit_static <- (outer(dat[,1], ext_region_pts[,1], '>') & 
                   outer(dat[,2], ext_region_pts[,2], '>')) * 1.0

    qn <- n_dat^(-gamma1)
    
    boot_draws_Wplus <- rep(NA, B)
    boot_draws_Wminus <- rep(NA, B)

    if (progress_bar) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (k in 1:B) {
        idx <- sample.int(n_dat, n_dat, replace=TRUE)
        
        # Tabulate trick for massive memory savings
        w <- tabulate(idx, nbins = n_dat)
        emp_probs <- as.numeric(w %*% hit_static) / n_dat
        
        Z_boot <- Z_static[idx, , drop=FALSE]
        rq <- colQuantiles(Z_boot, probs=1-qn, type=1, drop=TRUE)

        rm(Z_boot) # explicit garbage collection
        
        tail_probs <- ((rq/rp)^(1/xi_hat))*qn
        final_probs <- ifelse(emp_probs >= qn, emp_probs, tail_probs)
        
        # Taking the max/min deviation over the entire conservative envelope region
        boot_draws_Wplus[k] <- max(final_probs - p)
        boot_draws_Wminus[k] <- min(final_probs - p)
        
        if (progress_bar) utils::setTxtProgressBar(pb, k)
    }
    if (progress_bar) close(pb)
    
    rm(Z_static, hit_static) # Final cleanup

    c_plus_estimates <- list()
    c_minus_estimates <- list()
    for (i in 1:length(alphas)) {
        alpha <- alphas[i]        
        c_plus_estimate <- as.numeric(quantile(boot_draws_Wplus, probs = 1-(alpha/2)))
        c_minus_estimate <- as.numeric(quantile(boot_draws_Wminus, probs = (alpha/2)))
        c_plus_estimates[as.character(alpha)] <- c_plus_estimate
        c_minus_estimates[as.character(alpha)] <- c_minus_estimate
    }

    survFunc <- function(pts) {
        survProb <- vectorizedBlendedSurv(pts, dat, gamma1, xi_hat)
        return(survProb)
    }

    # Compile the final result list
    res_lst <- list(
        dat = dat,
        p = p,
        gamma1 = gamma1,
        gamma2 = gamma2,
        xi_hat = xi_hat,
        thetas = thetas,
        r_naive = Rp_first_order,
        r_lower_env = r_lower,
        r_upper_env = r_upper,
        c_plus_estimates = c_plus_estimates,
        c_minus_estimates = c_minus_estimates,
        survFunc = survFunc,
        transformed = FALSE
    )

    return(res_lst)
}

# loop over the distributions
for (i in 1:length(dists)) {
    dist <- dists[i]
    sampling_func <- loadSamplingFunction(dist)
    
    if (!is.null(save_full_path)) {
        distpath <- paste0(save_full_path, dist, '/')
        if (!dir.exists(distpath)) dir.create(distpath, recursive = TRUE)
    }
    print(paste0('Starting new distribution: ', dist))
    
    for (k in 1:length(ns)) {
        n <- ns[k]
        p <- pn(n)

        # compute true p-isoline
        isoline <- drawIsoline(dist=dist, numCoords=n_coords, gridUbs=ubs, gridLbs=lbs, prob=p)

        thetas <- atan2(isoline[,2], isoline[,1])
        print(paste0('Starting new n: ', n))
        
        # create the n-specific directory BEFORE the parallel loop
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

        # define function to be executed in parallel
        parallelizedCode <- function(ind) {
            dat <- sampling_func(n)
            
            # use the new Conservative Envelope bounding function
            regions <- computeExtremeRegionBC_Conservative(
                dat=dat, alphas=alphas, p=p, B=B, 
                gamma1=gamma1, gamma2=gamma2, 
                ncoords=n_coords, n_lines=n_lines, lbs=lbs
            )

            is_covereds <- evaluateCoverage(regions, isoline)
            regions$is_covereds <- is_covereds
            
            # save directly to disk from the worker node
            if (!is.null(n_dirpath)) {
                saveRDS(regions, file=paste0(n_dirpath, 'simulation', ind, '.RData'))
            }
            
            # return only the concise row needed for the CSV
            res_df <- data.frame(
                dist = dist,
                p = p,
                n = n,
                alpha = names(regions$is_covereds),
                is_covered = unlist(regions$is_covereds, use.names=FALSE),
                c_plus_estimate = unlist(regions$c_plus_estimates[names(regions$is_covereds)]),
                c_minus_estimate = unlist(regions$c_minus_estimates[names(regions$is_covereds)]),
                xi_hat = regions$xi_hat  # Saved to track shape estimator variance
            )
            
            return(res_df)
        }

        # execute parallel loop, returning a list of simple data frames.
        samp_df_list <- foreach(l = 1:n_sims, 
                .options.snow = opts,
                .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt', 'matrixStats')) %dopar% parallelizedCode(l)

        close(pb)
        stopCluster(clust)

        # save CSV if save_df_path argument is provided
        if (!is.null(args$save_df_path)) {
            # combine concise data frame rows returned from parallel workers
            coverage_df <- do.call(rbind, samp_df_list)

            save_fname <- paste0(args$save_df_path, 'extremal_coverage_conservative_envelope_gamma1_0.4_gamma2_0.5.csv')
            file_exists <- file.exists(save_fname)
                
            # write to CSV using write.table for append support
            write.table(coverage_df, 
                        file = save_fname, 
                        sep = ",", 
                        row.names = FALSE, 
                        col.names = !file_exists, 
                        append = file_exists)     
        }
    }
}