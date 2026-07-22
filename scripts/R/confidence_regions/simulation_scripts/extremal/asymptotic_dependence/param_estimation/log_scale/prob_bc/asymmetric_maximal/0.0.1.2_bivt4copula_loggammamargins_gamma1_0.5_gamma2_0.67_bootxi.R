# Coverage for confidence tubes (Mammen-Polonik style), extreme regime.
# Asymmetric Fractional Log-Probability bias correction + tube size.
#
# DISTRIBUTION: t_4-copula with LOG-GAMMA margins (xi = 1/4, rho = 0).
#   Slowly-varying (log x)^{a-1} tail factor => rho = 0 => slow convergence to
#   first-order MRV => the regime where bias correction should matter.
#
# Jimmy Butler, June 2026

set.seed(45678)

library(dplyr); library(mvtnorm); library(mnormt)
library(foreach); library(doSNOW); library(parallel)
library(argparse); library(matrixStats)

source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/coverageEvaluation.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')

parser <- ArgumentParser(description = "Asymmetric Fractional Log-Prob BC sims: t4-copula, log-gamma margins.")
parser$add_argument("--save_full_path", type = "character", default = NULL, required = FALSE)
parser$add_argument("--n_cores", type = "integer", default = 64)
parser$add_argument("--save_df_path", type = "character", default = NULL)
args <- parser$parse_args()

if (is.null(args$save_full_path) && is.null(args$save_df_path)) {
    stop("Error: Both --save_full_path and --save_df_path are NULL.")
}

# --- simulation parameters ---
ns <- c(1000, 5000, 10000, 50000, 100000)
alphas <- c(0.05, 0.1, 0.01)
C <- 5
gamma1 <- 0.5
gamma2 <- 0.67
B <- 5000
n_sims <- 500
n_coords <- 50

pn <- function(n){ C/n }

# --- log-gamma / copula parameters ---
# X = exp(G), G ~ Gamma(shape=a, rate=b):  Fbar(x) ~ C x^{-b} (log x)^{a-1}
#   tail index 1/xi = b  => b = 4 gives xi = 1/4
#   a != 1 gives a genuine slowly-varying log factor => rho = 0   (a = 1 => exact Pareto)
lg_shape <- 2      # a
lg_rate  <- 4      # b = 1/xi
copula_df <- 4     # t-copula degrees of freedom (nu); rho_cop = -2/nu = -1/2

dists <- c('t4copula_loggamma')
lbs <- c(0, 0)
ubs <- c(200, 200)

n_cores <- args$n_cores
save_full_path <- args$save_full_path

# ------------------------------------------------------------------
# Sampler: t_nu-copula + log-gamma margins, restricted to the positive
# orthant support the method assumes (margins are positive since X = exp(G) > 0).
# ------------------------------------------------------------------
sample_t4copula_loggamma <- function(n) {
    Sigma <- matrix(c(1, 0.7, 0.7, 1), 2, 2)
    # t-copula: multivariate t -> uniforms via its own marginal t CDF
    Tm <- mvtnorm::rmvt(n, sigma = Sigma, df = copula_df)      # nu must match the marginal t CDF below
    U  <- pt(Tm, df = copula_df)                               # t_nu-copula uniforms
    # log-gamma margins: X = exp( qgamma(U) )
    X1 <- exp(qgamma(U[,1], shape = lg_shape, rate = lg_rate))
    X2 <- exp(qgamma(U[,2], shape = lg_shape, rate = lg_rate))
    data.frame(X1 = X1, X2 = X2)
}

# marginal survival + quantile for the log-gamma (for the TRUE isoline)
qloggamma <- function(u) exp(qgamma(u, shape = lg_shape, rate = lg_rate))

computeExtremeRegionBC_FractionalLogAsymmetric <- function(dat, alphas, p, B, gamma1=1/2, gamma2=2/3,
                                                           ncoords=50, lbs=c(0,0),
                                                           progress_bar=FALSE, verbose=FALSE) {
    n_dat <- nrow(dat)
    q1 <- n_dat^(-gamma1); q2 <- n_dat^(-gamma2)

    base_radii <- sqrt(dat[,1]^2 + dat[,2]^2)
    xi_hat <- estimate_xi_hill(base_radii, gamma1)

    isoline_p  <- drawExtremeIsoline(dat=dat, p=p,  n_coords=ncoords, grid_lbs=lbs, gamma=gamma1, xi=xi_hat)
    isoline_q2 <- drawExtremeIsoline(dat=dat, p=q2, n_coords=ncoords, grid_lbs=lbs, gamma=gamma1, xi=xi_hat)

    hit_q2 <- (outer(dat[,1], isoline_q2[,1], '>') & outer(dat[,2], isoline_q2[,2], '>')) * 1.0
    emp_probs_q2 <- colMeans(hit_q2)

    delta_P <- emp_probs_q2 - q2
    C_theta <- delta_P / (q2 * log(q1 / q2))
    beta_frac_theta <- C_theta * log(q1 / p)

    beta_frac_pos <- max(0, max(beta_frac_theta))
    beta_frac_neg <- max(0, -min(beta_frac_theta))

    thetas <- atan2(isoline_p[,2], isoline_p[,1])
    rp <- sqrt(rowSums(isoline_p^2))

    Z_static <- matrix(NA, nrow = n_dat, ncol = length(thetas))
    pos_orthant_mask <- (dat[,1] > 0) & (dat[,2] > 0)
    ind_0  <- which(thetas == 0);    if (length(ind_0)  > 0) Z_static[, ind_0]  <- dat[,1] * pos_orthant_mask
    ind_90 <- which(thetas == pi/2); if (length(ind_90) > 0) Z_static[, ind_90] <- dat[,2] * pos_orthant_mask
    inds_no_ax <- which(!(thetas == 0 | thetas == pi/2))
    if (length(inds_no_ax) > 0) {
        Z_static[, inds_no_ax] <- pmin(
            dat[,1] %o% (1/cos(thetas[inds_no_ax])),
            dat[,2] %o% (1/sin(thetas[inds_no_ax]))
        )
    }

    hit_static <- (outer(dat[,1], isoline_p[,1], '>') & outer(dat[,2], isoline_p[,2], '>')) * 1.0
    qn <- n_dat^(-gamma1)
    boot_probs_mat <- matrix(NA, nrow=B, ncol=ncoords)

    if (progress_bar) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (k in 1:B) {
        idx <- sample.int(n_dat, n_dat, replace=TRUE)
        xi_boot <- estimate_xi_hill(base_radii[idx], gamma1)
        w <- tabulate(idx, nbins = n_dat)
        emp_probs <- as.numeric(w %*% hit_static) / n_dat
        Z_boot <- Z_static[idx, , drop=FALSE]
        rq <- colQuantiles(Z_boot, probs=1-qn, type=1, drop=TRUE)
        rm(Z_boot)
        tail_probs <- ((rq/rp)^(1/xi_boot))*qn
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
        c_plus_log  <- as.numeric(quantile(boot_draws_Wplus,  probs = 1-(alpha/2)))
        c_minus_log <- as.numeric(quantile(boot_draws_Wminus, probs = (alpha/2)))

        raw_c_plus_estimates[as.character(alpha)]  <- p * exp(c_plus_log)  - p
        raw_c_minus_estimates[as.character(alpha)] <- p * exp(c_minus_log) - p

        c_plus_estimates[as.character(alpha)]  <- p * exp(c_plus_log  + beta_frac_neg) - p
        c_minus_estimates[as.character(alpha)] <- p * exp(c_minus_log - beta_frac_pos) - p

        # --- tube size (bias-corrected walls) ------------------------------
        # walls -> survival levels -> radii via the same radial projection the
        # isoline used: r(s) = rp * (p / s)^{xi_hat}. Log-scale => s_lo > 0 always.
        s_hi <- p + as.numeric(c_plus_estimates[[as.character(alpha)]])   # inner (higher surv, smaller r)
        s_lo <- p + as.numeric(c_minus_estimates[[as.character(alpha)]])  # outer (lower surv, larger r)
        if (s_lo <= 0) {
            size_estimates[[as.character(alpha)]] <- data.frame(
                is_unbounded = TRUE, width_abs_median = NA_real_, width_abs_max = NA_real_,
                width_rel_median = NA_real_, width_rel_max = NA_real_)
        } else {
            r_inner <- rp * (p / s_hi)^xi_hat
            r_outer <- rp * (p / s_lo)^xi_hat
            size_estimates[[as.character(alpha)]] <- cbind(
                data.frame(is_unbounded = FALSE),
                computeTubeSize(r_inner, r_outer, rp))
        }
    }

    survFunc <- function(pts) vectorizedBlendedSurv(pts, dat, gamma1, xi_hat)

    list(
        dat = dat, p = p, gamma1 = gamma1, gamma2 = gamma2, xi_hat = xi_hat,
        thetas = thetas, r_naive = rp,
        beta_frac_pos = beta_frac_pos, beta_frac_neg = beta_frac_neg,
        c_plus_estimates = c_plus_estimates, c_minus_estimates = c_minus_estimates,
        raw_c_plus_estimates = raw_c_plus_estimates, raw_c_minus_estimates = raw_c_minus_estimates,
        size_estimates = size_estimates,
        survFunc = survFunc, transformed = FALSE
    )
}

for (i in 1:length(dists)) {
    dist <- dists[i]
    sampling_func <- sample_t4copula_loggamma

    if (!is.null(save_full_path)) {
        distpath <- paste0(save_full_path, dist, '/')
        if (!dir.exists(distpath)) dir.create(distpath, recursive = TRUE)
    }
    print(paste0('Starting new distribution: ', dist))

    for (k in 1:length(ns)) {
        n <- ns[k]
        p <- pn(n)

        # TRUE p-isoline: t-copula survival isoline -> log-gamma margins (same recipe as your gauss-copula scripts)
        tcop_isoline  <- drawIsoline(dist='bivt', numCoords=n_coords, gridUbs=ubs, gridLbs=lbs, prob=p)
        copula_isoline <- data.frame(X1 = pt(tcop_isoline[,1], df=copula_df),
                                     X2 = pt(tcop_isoline[,2], df=copula_df))
        isoline <- data.frame(X1 = qloggamma(copula_isoline[,1]),
                              X2 = qloggamma(copula_isoline[,2]))

        print(paste0('Starting new n: ', n))

        n_dirpath <- NULL
        if (!is.null(save_full_path)) {
            n_dirpath <- paste0(distpath, 'n', n, '/')
            if (!dir.exists(n_dirpath)) dir.create(n_dirpath, recursive = TRUE)
        }

        clust <- makeSOCKcluster(n_cores)
        registerDoSNOW(clust)
        clusterEvalQ(clust, { source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R') })
        pb <- txtProgressBar(min = 1, max = n_sims, style = 3)
        progress <- function(iter) setTxtProgressBar(pb, iter)
        opts <- list(progress = progress)

        parallelizedCode <- function(ind) {
            dat <- sampling_func(n)
            regions <- computeExtremeRegionBC_FractionalLogAsymmetric(
                dat=dat, alphas=alphas, p=p, B=B, gamma1=gamma1, gamma2=gamma2,
                ncoords=n_coords, lbs=lbs)
            is_covereds <- evaluateCoverage(regions, isoline)
            regions$is_covereds <- is_covereds
            if (!is.null(n_dirpath)) saveRDS(regions, file=paste0(n_dirpath, 'simulation', ind, '.RData'))

            alpha_names <- names(regions$is_covereds)
            size_df <- do.call(rbind, regions$size_estimates[alpha_names])
            data.frame(
                dist = dist, p = p, n = n,
                alpha = alpha_names,
                is_covered = unlist(regions$is_covereds, use.names=FALSE),
                c_plus_adj  = unlist(regions$c_plus_estimates[alpha_names]),
                c_minus_adj = unlist(regions$c_minus_estimates[alpha_names]),
                c_plus_raw  = unlist(regions$raw_c_plus_estimates[alpha_names]),
                c_minus_raw = unlist(regions$raw_c_minus_estimates[alpha_names]),
                beta_frac_pos = regions$beta_frac_pos,
                beta_frac_neg = regions$beta_frac_neg,
                xi_hat = regions$xi_hat,
                is_unbounded     = size_df$is_unbounded,
                width_abs_median = size_df$width_abs_median,
                width_abs_max    = size_df$width_abs_max,
                width_rel_median = size_df$width_rel_median,
                width_rel_max    = size_df$width_rel_max
            )
        }

        samp_df_list <- foreach(l = 1:n_sims, .options.snow = opts,
                .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt', 'matrixStats')) %dopar% parallelizedCode(l)
        close(pb); stopCluster(clust)

        if (!is.null(args$save_df_path)) {
            coverage_df <- do.call(rbind, samp_df_list)
            save_fname <- paste0(args$save_df_path, '0.0.1.2_bivt4copula_loggammamargins_gamma1_0.5_gamma2_0.67_bootxi.csv')
            file_exists <- file.exists(save_fname)
            write.table(coverage_df, file = save_fname, sep = ",", row.names = FALSE,
                        col.names = !file_exists, append = file_exists)
        }
    }
}