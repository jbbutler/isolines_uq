# Script to compute coverage for confidence tubes in the style of Mammen and Polonik (2013)
# ADAPTED FOR A RAW-SCALE HIDDEN REGULAR VARIATION (HRV) SANITY CHECK,
# WITH A COOLEY ET AL. (2019)-STYLE SMOOTH TRANSITION from interior HRV scaling
# to first-order (marginal) MRV scaling near the axes.
#
# ORACLE-INDEX VARIANT: BOTH tail indices are FIXED at their true values everywhere --
# no Hill estimation on pmin/pmax, and no re-estimation in the bootstrap.
#   Interior/hidden index xi_int = eta * xi = 0.85 * 0.25 = 0.2125
#     (eta = (1 + rho)/2 = 0.85 for a Gaussian copula with correlation 0.7)
#   Marginal index        xi_marg = xi = 0.25  (t_4 margins)
#
# Purpose: with estimation removed, isolate the role of the second-order bias
# correction from any Hill-bias compensation (the asymptotic-independence analogue
# of the bivariate-t / log-gamma oracle runs). Since the uncorrected HRV tubes
# already slightly undercover, the corrected-vs-uncorrected contrast should be clean.
#
# MODIFIED: Includes computation of tube sizes (absolute and relative radial width).

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
parser <- ArgumentParser(description = "Args for Raw Asymmetric HRV ORACLE-index sims with Cooley smoothing.")
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
angular_pad <- 0    # full quadrant [0, pi/2]; Cooley smoothing handles the axis transition

# Cooley (2019) smoothing parameter (applied on the STANDARDIZED-scale weight)
beta_cooley <- 200

# --- ORACLE indices (fixed everywhere; no Hill) ---
# Gaussian copula with correlation rho_corr -> coefficient of tail dependence eta = (1 + rho_corr)/2
# t_4 margins -> marginal EVI xi = 1/4. Interior/hidden EVI = eta * xi.
rho_corr     <- 0.7
eta_true     <- (1 + rho_corr) / 2          # = 0.85
xi_marg_true <- 1/4                          # = 0.25 (marginal, t_4)
xi_int_true  <- eta_true * xi_marg_true      # = 0.2125 (interior / hidden)

pn <- function(n){ C/n }
dists <- c('gauss_copula_t4')
lbs <- c(0,0)
ubs <- c(200, 200)

n_cores <- args$n_cores
save_full_path <- args$save_full_path

# sampling function for a gaussian copula, t marginals (asymptotic independence, same marginal tails)
sample_gauss_copula_t4 <- function(n) {
    Sigma <- matrix(c(1, rho_corr, rho_corr, 1), 2, 2)
    Z <- mvtnorm::rmvnorm(n, mean = c(0,0), sigma = Sigma)
    U <- pnorm(Z)
    X <- qt(U, df = 4)
    return(data.frame(X))
}

# ----------------------------------------------------------------------
# Tube Size Evaluator
# ----------------------------------------------------------------------
computeTubeSize <- function(r_inner, r_outer, r_ref) {

  stopifnot(length(r_outer) == length(r_inner),
            length(r_ref)   == length(r_inner))

  w_abs <- r_outer - r_inner                # data-space distance
  w_rel <- w_abs / r_ref                    # dimensionless

  data.frame(
    width_abs_median = median(w_abs),
    width_abs_max    = max(w_abs),
    width_rel_median = median(w_rel),
    width_rel_max    = max(w_rel)
  )
}

# ----------------------------------------------------------------------
# Cooley-style smoothed effective exponent along a ray.
# ----------------------------------------------------------------------
cooley_xi_eff <- function(thetas, xi_int, xi_marg, beta_cooley) {
    a   <- 1 / xi_marg
    cs  <- cos(thetas)^a
    sn  <- sin(thetas)^a
    wmx <- pmax(cs, sn) / (cs + sn)
    m   <- 1 - wmx^beta_cooley
    m * xi_int + (1 - m) * xi_marg
}

# ----------------------------------------------------------------------
# Cooley-aware blended survival function
# ----------------------------------------------------------------------
vectorizedBlendedSurvCooley <- function(pts, dat, gamma, xi_int, xi_marg, beta_cooley) {
    pts <- as.matrix(pts)
    if (any(pts < 0)) {
        stop('Error: some points outside the nonnegative orthant.')
    }
    n_dat <- nrow(dat)
    qn    <- n_dat^(-gamma)

    match_x <- outer(dat[,1], pts[,1], ">")
    match_y <- outer(dat[,2], pts[,2], ">")
    empsurv_probs <- colMeans(match_x & match_y)

    theta  <- atan2(pts[,2], pts[,1])
    rp     <- sqrt(rowSums(pts^2))
    xi_eff <- cooley_xi_eff(theta, xi_int, xi_marg, beta_cooley)

    pos_orthant_mask <- (dat[,1] > 0) & (dat[,2] > 0)
    rq <- rep(NA_real_, length(theta))

    if (any(theta == 0)) {
        ind_0 <- which(theta == 0)
        M_0   <- dat[,1] * pos_orthant_mask
        rq[ind_0] <- quantile(M_0, probs = 1 - qn, names = FALSE, type = 1)
    }
    if (any(theta == pi/2)) {
        ind_90 <- which(theta == pi/2)
        M_90   <- dat[,2] * pos_orthant_mask
        rq[ind_90] <- quantile(M_90, probs = 1 - qn, names = FALSE, type = 1)
    }
    inds_no_ax <- which(!(theta == 0 | theta == pi/2))
    if (length(inds_no_ax) > 0) {
        angles_no_ax <- theta[inds_no_ax]
        inv_cos <- 1 / cos(angles_no_ax)
        inv_sin <- 1 / sin(angles_no_ax)
        proj_x  <- dat[,1] %o% inv_cos
        proj_y  <- dat[,2] %o% inv_sin
        M_int   <- pmin(proj_x, proj_y)
        rq[inds_no_ax] <- colQuantiles(M_int, probs = 1 - qn, type = 1, drop = TRUE)
    }

    tail_probs <- ((rq / rp)^(1 / xi_eff)) * qn
    ifelse(empsurv_probs >= qn, empsurv_probs, tail_probs)
}


# ----------------------------------------------------------------------
# Confidence Tube Generator: Raw Scale HRV + Cooley smoothing, ORACLE indices
# ----------------------------------------------------------------------
computeExtremeRegionBC_HRV_Raw <- function(dat, alphas, p, B, gamma1=1/2, gamma2=2/3,
                                           ncoords=50, lbs=c(0,0), beta_cooley=200,
                                           xi_int_fixed=0.2125, xi_marg_fixed=0.25,
                                           progress_bar=FALSE, verbose=FALSE) {

    n_dat <- nrow(dat)
    q1 <- n_dat^(-gamma1)
    q2 <- n_dat^(-gamma2)

    # ORACLE: fixed true indices, NO Hill estimation
    xi_int_hat  <- xi_int_fixed
    xi_marg_hat <- xi_marg_fixed

    # 1. Strictly Interior Angles
    thetas <- seq(angular_pad, (pi/2) - angular_pad, length.out = ncoords)
    xi_eff <- cooley_xi_eff(thetas, xi_int_hat, xi_marg_hat, beta_cooley)

    # standardized-scale transition weight (fixed); with oracle indices the bootstrap
    # effective exponent equals xi_eff exactly, so m_theta is only needed for parity/reference.
    a_marg <- 1 / xi_marg_hat
    cs <- cos(thetas)^a_marg; sn <- sin(thetas)^a_marg
    m_theta <- 1 - (pmax(cs, sn) / (cs + sn))^beta_cooley

    Z_static <- matrix(NA_real_, nrow = n_dat, ncol = length(thetas))
    pos_orthant_mask <- (dat[,1] > 0) & (dat[,2] > 0)
    for (j in 1:length(thetas)) {
        th <- thetas[j]
        if (th == 0) {
            Z_static[, j] <- dat[,1] * pos_orthant_mask
        } else if (th == pi/2) {
            Z_static[, j] <- dat[,2] * pos_orthant_mask
        } else {
            Z_static[, j] <- pmin(dat[,1]/cos(th), dat[,2]/sin(th))
        }
    }

    rq1 <- colQuantiles(Z_static, probs=1-q1, type=1, drop=TRUE)

    rp_radii  <- rq1 * (q1 / p)^xi_eff
    rq2_radii <- rq1 * (q1 / q2)^xi_eff

    isoline_p  <- cbind(rp_radii  * cos(thetas), rp_radii  * sin(thetas))
    isoline_q2 <- cbind(rq2_radii * cos(thetas), rq2_radii * sin(thetas))

    hit_q2 <- (outer(dat[,1], isoline_q2[,1], '>') & outer(dat[,2], isoline_q2[,2], '>')) * 1.0
    emp_probs_q2 <- colMeans(hit_q2)

    delta_P <- emp_probs_q2 - q2
    C_theta <- delta_P / (q2 * log(q1 / q2))
    beta_frac_theta <- C_theta * log(q1 / p)

    beta_frac_pos <- max(0, max(beta_frac_theta))
    beta_frac_neg <- max(0, -min(beta_frac_theta))

    if (verbose) print('Estimating pure log-variance via Centered Bootstrap (oracle indices fixed).')

    hit_static <- (outer(dat[,1], isoline_p[,1], '>') &
                   outer(dat[,2], isoline_p[,2], '>')) * 1.0

    qn <- n_dat^(-gamma1)
    boot_probs_mat <- matrix(NA, nrow=B, ncol=ncoords)

    if (progress_bar) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (k in 1:B) {
        idx <- sample.int(n_dat, n_dat, replace=TRUE)

        # ORACLE: indices fixed in the bootstrap too (no re-estimation) => xi_eff_boot == xi_eff
        w <- tabulate(idx, nbins = n_dat)
        emp_probs <- as.numeric(w %*% hit_static) / n_dat

        Z_boot <- Z_static[idx, , drop=FALSE]
        rq <- colQuantiles(Z_boot, probs=1-qn, type=1, drop=TRUE)
        rm(Z_boot)

        tail_probs <- ((rq / rp_radii)^(1/xi_eff)) * qn
        final_probs <- ifelse(emp_probs >= qn, emp_probs, tail_probs)

        boot_probs_mat[k, ] <- log(final_probs)
        if (progress_bar) utils::setTxtProgressBar(pb, k)
    }
    if (progress_bar) close(pb)

    rm(Z_static, hit_static)

    mean_log_boot_probs <- colMeans(boot_probs_mat)
    boot_log_devs <- sweep(boot_probs_mat, MARGIN = 2, STATS = mean_log_boot_probs, FUN = "-")

    boot_draws_Wplus  <- apply(boot_log_devs, 1, max)
    boot_draws_Wminus <- apply(boot_log_devs, 1, min)

    c_plus_estimates <- list(); c_minus_estimates <- list()
    raw_c_plus_estimates <- list(); raw_c_minus_estimates <- list()
    size_estimates <- list()

    for (i in 1:length(alphas)) {
        alpha <- alphas[i]
        akey <- as.character(alpha)

        c_plus_log  <- as.numeric(quantile(boot_draws_Wplus,  probs = 1-(alpha/2)))
        c_minus_log <- as.numeric(quantile(boot_draws_Wminus, probs = (alpha/2)))

        raw_c_plus_estimates[[akey]]  <- p * exp(c_plus_log)  - p
        raw_c_minus_estimates[[akey]] <- p * exp(c_minus_log) - p

        cplus  <- p * exp(c_plus_log  + beta_frac_neg) - p
        cminus <- p * exp(c_minus_log - beta_frac_pos) - p

        c_plus_estimates[[akey]]  <- cplus
        c_minus_estimates[[akey]] <- cminus

        s_lo <- p + cminus
        if (s_lo <= 0) {
            size_estimates[[akey]] <- data.frame(
                is_unbounded     = TRUE,
                width_abs_median = NA_real_,
                width_abs_max    = NA_real_,
                width_rel_median = NA_real_,
                width_rel_max    = NA_real_
            )
        } else {
            s_hi <- min(p + cplus, 0.999)

            # calculate extreme boundaries algebraically along the rays using power-law rules
            r_inner <- rp_radii * (p / s_hi)^xi_eff
            r_outer <- rp_radii * (p / s_lo)^xi_eff
            r_ref   <- rp_radii

            sz <- computeTubeSize(r_inner, r_outer, r_ref)
            size_estimates[[akey]] <- cbind(data.frame(is_unbounded = FALSE), sz)
        }
    }

    survFunc <- function(pts) {
        vectorizedBlendedSurvCooley(pts, dat, gamma1, xi_int_hat, xi_marg_hat, beta_cooley)
    }

    res_lst <- list(
        dat = dat,
        p = p,
        gamma1 = gamma1,
        gamma2 = gamma2,
        xi_hat = xi_int_hat,
        xi_marg_hat = xi_marg_hat,
        beta_cooley = beta_cooley,
        thetas = thetas,
        xi_eff = xi_eff,
        r_naive = rp_radii,
        beta_frac_pos = beta_frac_pos,
        beta_frac_neg = beta_frac_neg,
        c_plus_estimates = c_plus_estimates,
        c_minus_estimates = c_minus_estimates,
        raw_c_plus_estimates = raw_c_plus_estimates,
        raw_c_minus_estimates = raw_c_minus_estimates,
        size_estimates = size_estimates,
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
                ncoords=n_coords, lbs=lbs, beta_cooley=beta_cooley,
                xi_int_fixed=xi_int_true, xi_marg_fixed=xi_marg_true
            )

            is_covereds <- evaluateCoverage(regions, isoline)
            regions$is_covereds <- is_covereds

            if (!is.null(n_dirpath)) {
                saveRDS(regions, file=paste0(n_dirpath, 'simulation', ind, '.RData'))
            }

            alpha_names <- names(regions$is_covereds)
            size_df <- do.call(rbind, regions$size_estimates[alpha_names])

            res_df <- data.frame(
                dist = dist,
                p = p,
                n = n,
                alpha = alpha_names,
                is_covered = unlist(regions$is_covereds, use.names=FALSE),
                c_plus_adj = unlist(regions$c_plus_estimates[alpha_names]),
                c_minus_adj = unlist(regions$c_minus_estimates[alpha_names]),
                c_plus_raw = unlist(regions$raw_c_plus_estimates[alpha_names]),
                c_minus_raw = unlist(regions$raw_c_minus_estimates[alpha_names]),
                beta_frac_pos = regions$beta_frac_pos,
                beta_frac_neg = regions$beta_frac_neg,
                xi_hat = regions$xi_hat,
                xi_marg_hat = regions$xi_marg_hat,
                is_unbounded     = size_df$is_unbounded,
                width_abs_median = size_df$width_abs_median,
                width_abs_max    = size_df$width_abs_max,
                width_rel_median = size_df$width_rel_median,
                width_rel_max    = size_df$width_rel_max
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

            save_fname <- paste0(args$save_df_path, '1.0.1.2_gausscopula_t4margs_gamma1_0.5_gamma2_0.67_oraclexi_cooleysmooth.csv')
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