# Coverage for confidence tubes (Mammen-Polonik style), extreme regime.
# MARGINAL-TRANSFORMATION PROCEDURE, asymptotic independence (AI/HRV) version
# with Cooley et al. (2019) smoothing -- NATIVE on the standardized scale.
#
# DISTRIBUTION: bivariate GAUSSIAN, unit variances, correlation 0.7.
#   Margins are LIGHT-TAILED (xi = 0): the raw-scale HRV machinery could not
#   handle this at all (it required regularly varying raw margins). The marginal
#   transform removes that requirement -- this run is the demonstration.
#   Gaussian copula => asymptotically independent, eta = (1 + 0.7)/2 = 0.85.
#
# What the standardization buys in the AI case:
#   * axis/marginal index on the Pareto scale is EXACTLY 1 (known) -- the old
#     xi_marg Hill estimator is gone;
#   * only eta is estimated: Hill on the standardized minima min(Z1, Z2)
#     (the classical Ledford-Tawn / Draisma et al. estimator);
#   * Cooley's beta = 200 applies NATIVELY (plain weight max(z)/(z1+z2)) --
#     no raw-scale power adaptation;
#   * Cooley transition: xi_eff(theta) = m(theta)*eta_hat + (1 - m(theta))*1.
#
# Pipeline per simulation:
#   1. Fit each margin: empirical body + GPD tail above u_j = F_emp^{-1}(1 - q1),
#      plain unconstrained MLE; fitted CDF clamped to [1/(2n), 1 - 1/(2n)].
#   2. Standardize to unit Pareto: Z_ij = 1 / (1 - F~_j(X_ij)).
#   3. Standardized isoline: min-projection anchor at q1, projection
#      r_p(theta) = r_q1(theta) * (q1/p)^{xi_eff(theta)}; back-map through the
#      full-sample fits to FIXED original-scale monitoring points.
#      THETA GRID: hybrid -- uniform interior + log-spaced points hugging each
#      axis, so the beta = 200 transition zone (within ~0.2 deg of the axes on
#      the standardized scale) is actually monitored by the sup/inf band.
#   4. Variance (d+-): bootstrap with MARGINS REFIT PER DRAW and eta RE-ESTIMATED
#      PER DRAW on the refit-standardized minima; composed estimator re-evaluated
#      at the fixed monitoring points; log-scale centered sup/inf quantiles.
#   5. Bias pads: radial pad on the standardized scale (now the HIDDEN measure's
#      convergence -- the Gaussian copula's notoriously slow, near-rho = 0 HRV
#      approach; rho -> 0 worst case, M0 multiplier) + marginal pad
#      b_marg = (1/eta_hat) * M0 * max_j |B_j|  (kappa = 1/eta under AI).
#   6. Tube: p*exp(c-) <= Fhat(x) <= p*exp(c+), Fhat the composed clamped estimator.
#
# COST NOTE: margins refit + data re-standardized + eta re-estimated every draw;
# expect ~2x the common-marginal runtime at equal B.
#
# Jimmy Butler, July 2026

set.seed(45678)

library(dplyr); library(mvtnorm); library(mnormt)
library(foreach); library(doSNOW); library(parallel)
library(argparse); library(matrixStats)

source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')            # computeTubeSize + estimate_xi_hill live here
source('~/isolines_uq/scripts/R/confidence_regions/modules/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/coverageEvaluation.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')

parser <- ArgumentParser(description = "Marginal-transform AI/HRV tube sims: bivariate Gaussian (0.7), GPD margins, Pareto standardization, Cooley smoothing.")
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

# --- distribution parameters (data-generating) ---
rho_corr <- 0.7                 # Gaussian correlation; eta_true = (1 + rho_corr)/2 = 0.85
eta_true <- (1 + rho_corr) / 2  # for reference/diagnostics only (NOT used by the method)

# --- Cooley smoothing parameter, native on the standardized (Frechet/Pareto) scale ---
beta_cooley <- 200

dists <- c('bivgauss')
lbs <- c(0, 0)
ubs <- c(200, 200)

n_cores <- args$n_cores
save_full_path <- args$save_full_path

# sampler: bivariate Gaussian, unit variances, correlation 0.7
sample_bivgauss <- function(n) {
    Sigma <- matrix(c(1, rho_corr, rho_corr, 1), 2, 2)
    X <- mvtnorm::rmvnorm(n, mean = c(0, 0), sigma = Sigma)
    data.frame(X1 = X[,1], X2 = X[,2])
}

# hybrid theta grid: uniform interior + log-spaced points hugging each axis so the
# beta = 200 transition zone is monitored (uniform spacing ~1.8 deg would skip it)
make_theta_grid <- function(ncoords, n_axis = 8, ax_lo = 5e-4, ax_hi = 0.05) {
    ax <- exp(seq(log(ax_lo), log(ax_hi), length.out = n_axis))
    interior <- seq(0, pi/2, length.out = ncoords - 2 * n_axis)
    sort(unique(c(interior, ax, (pi/2) - ax)))
}

# ======================================================================
# Marginal machinery: unconstrained GPD MLE + blended (empirical/GPD) CDF
# with the [1/(2n), 1 - 1/(2n)] clamp, forward and inverse transforms.
# (Gaussian margins => expect xi_gpd near 0; the machinery does not care.)
# ======================================================================

gpd_nll <- function(par, exc) {
    sigma <- exp(par[1]); xi <- par[2]
    z <- 1 + xi * exc / sigma
    if (any(z <= 0) || !is.finite(sigma)) return(1e10)
    if (abs(xi) < 1e-8) {
        sum(log(sigma) + exc / sigma)
    } else {
        sum(log(sigma) + (1/xi + 1) * log(z))
    }
}

fit_gpd_mle <- function(exc) {
    m <- mean(exc); v <- var(exc)
    xi0 <- 0.5 * (1 - m^2 / v)
    if (!is.finite(xi0)) xi0 <- 0.1
    xi0 <- max(min(xi0, 0.9), -0.4)
    sigma0 <- m * (1 - xi0); if (!is.finite(sigma0) || sigma0 <= 0) sigma0 <- m
    out <- tryCatch(
        optim(c(log(sigma0), xi0), gpd_nll, exc = exc, method = "Nelder-Mead"),
        error = function(e) NULL)
    if (is.null(out) || out$convergence != 0 || !is.finite(out$par[1])) {
        return(list(sigma = sigma0, xi = xi0, failed = TRUE))
    }
    list(sigma = exp(out$par[1]), xi = out$par[2], failed = FALSE)
}

fit_marginal <- function(x, gamma1, fallback = NULL) {
    n <- length(x)
    q1 <- n^(-gamma1)
    xs <- sort(x)
    u <- quantile(x, probs = 1 - q1, type = 1, names = FALSE)
    exc <- x[x > u] - u
    g <- fit_gpd_mle(exc)
    used_fallback <- FALSE
    if (g$failed && !is.null(fallback)) {
        g <- list(sigma = fallback$sigma, xi = fallback$xi)
        used_fallback <- TRUE
    }
    list(u = u, sigma = g$sigma, xi = g$xi, q1 = q1, n = n,
         x_sorted = xs, used_fallback = used_fallback)
}

ptilde <- function(t, fit) {
    n <- fit$n
    Femp <- findInterval(t, fit$x_sorted) / n
    y <- pmax(t - fit$u, 0)
    if (abs(fit$xi) < 1e-8) {
        Fgpd <- 1 - fit$q1 * exp(-y / fit$sigma)
    } else {
        Fgpd <- 1 - fit$q1 * pmax(1 + fit$xi * y / fit$sigma, 0)^(-1/fit$xi)
    }
    Fv <- ifelse(t <= fit$u, Femp, Fgpd)
    pmin(pmax(Fv, 1/(2*n)), 1 - 1/(2*n))     # <- the clamp safeguard
}

qblend <- function(s, fit) {
    emp_idx <- pmax(1L, pmin(fit$n, as.integer(ceiling((1 - s) * fit$n))))
    x_emp <- fit$x_sorted[emp_idx]
    if (abs(fit$xi) < 1e-8) {
        x_gpd <- fit$u + fit$sigma * log(fit$q1 / s)
    } else {
        x_gpd <- fit$u + (fit$sigma / fit$xi) * ((s / fit$q1)^(-fit$xi) - 1)
    }
    ifelse(s >= fit$q1, x_emp, x_gpd)
}

standardize <- function(x1, x2, fit1, fit2) {
    cbind(1 / (1 - ptilde(x1, fit1)), 1 / (1 - ptilde(x2, fit2)))
}

backmap <- function(zpts, fit1, fit2) {
    s1 <- pmin(1 / pmax(zpts[, 1], 1), 1)
    s2 <- pmin(1 / pmax(zpts[, 2], 1), 1)
    cbind(qblend(s1, fit1), qblend(s2, fit2))
}

# ======================================================================
# Standardized-scale machinery, AI/HRV version:
#   axis index = 1 (KNOWN), interior index = eta (estimated),
#   Cooley transition with the NATIVE Frechet-scale weight, beta = 200.
# ======================================================================

# Cooley smoothed effective exponent: xi_marg = 1 on the standardized scale,
# so the weight is the plain max(z)/(z1+z2) -- beta = 200 as Cooley calibrated it
cooley_xi_eff <- function(thetas, xi_int, xi_marg, beta_cooley) {
    a   <- 1 / xi_marg                       # = 1 on the standardized scale
    cs  <- cos(thetas)^a
    sn  <- sin(thetas)^a
    wmx <- pmax(cs, sn) / (cs + sn)
    m   <- 1 - wmx^beta_cooley               # ~1 interior, -> 0 at the axes
    m * xi_int + (1 - m) * xi_marg
}

std_radial_anchor <- function(Z, thetas, qlev) {
    rq <- rep(NA_real_, length(thetas))
    ax0  <- which(thetas == 0)
    ax90 <- which(thetas == pi/2)
    if (length(ax0))  rq[ax0]  <- quantile(Z[,1], probs = 1 - qlev, type = 1, names = FALSE)
    if (length(ax90)) rq[ax90] <- quantile(Z[,2], probs = 1 - qlev, type = 1, names = FALSE)
    intr <- which(!(thetas == 0 | thetas == pi/2))
    if (length(intr)) {
        M <- pmin(Z[,1] %o% (1/cos(thetas[intr])), Z[,2] %o% (1/sin(thetas[intr])))
        rq[intr] <- colQuantiles(M, probs = 1 - qlev, type = 1, drop = TRUE)
    }
    rq
}

# blended standardized survival at query points: empirical above qn,
# Cooley-smoothed HRV projection (exponent xi_eff(theta)) below
std_blended_surv_hrv <- function(zpts, Z, qn, eta_hat, beta_cooley) {
    zpts <- as.matrix(zpts)
    emp <- colMeans(outer(Z[,1], zpts[,1], ">") & outer(Z[,2], zpts[,2], ">"))
    theta <- atan2(zpts[,2], zpts[,1])
    rp <- sqrt(rowSums(zpts^2))
    rq <- std_radial_anchor(Z, theta, qn)
    xi_eff <- cooley_xi_eff(theta, eta_hat, 1, beta_cooley)
    tail_probs <- ((rq / rp)^(1 / xi_eff)) * qn
    ifelse(emp >= qn, emp, tail_probs)
}

# ======================================================================
# Tube constructor: marginal transform + AI/HRV standardized tube,
# margins + eta in the bootstrap, radial + marginal bias pads
# ======================================================================
computeRegion_MargTransform_HRV <- function(dat, alphas, p, B, gamma1 = 1/2, gamma2 = 2/3,
                                            ncoords = 50, lbs = c(0, 0), beta_cooley = 200,
                                            progress_bar = FALSE, verbose = FALSE) {

    n_dat <- nrow(dat)
    q1 <- n_dat^(-gamma1); q2 <- n_dat^(-gamma2); qn <- q1
    prob_floor <- 0.5 / n_dat

    # --- Stage 1: full-sample marginal fits ---
    fit1 <- fit_marginal(dat[,1], gamma1)
    fit2 <- fit_marginal(dat[,2], gamma1)

    # --- Stage 2: standardize; estimate eta on the standardized minima ---
    Z <- standardize(dat[,1], dat[,2], fit1, fit2)
    eta_hat <- estimate_xi_hill(pmin(Z[,1], Z[,2]), gamma1)   # Ledford-Tawn eta

    # --- Stage 3: standardized isoline (Cooley xi_eff) + fixed monitoring points ---
    thetas <- make_theta_grid(ncoords)
    xi_eff_grid <- cooley_xi_eff(thetas, eta_hat, 1, beta_cooley)
    rq1_vec <- std_radial_anchor(Z, thetas, q1)
    rp_std  <- rq1_vec * (q1 / p)^xi_eff_grid
    ziso    <- cbind(rp_std * cos(thetas), rp_std * sin(thetas))
    x_iso   <- backmap(ziso, fit1, fit2)                # FIXED original-scale points

    # --- Stage 5a: radial pad (hidden-measure channel), standardized scale ---
    rq2_vec <- rq1_vec * (q1 / q2)^xi_eff_grid
    zq2 <- cbind(rq2_vec * cos(thetas), rq2_vec * sin(thetas))
    emp_q2 <- colMeans(outer(Z[,1], zq2[,1], ">") & outer(Z[,2], zq2[,2], ">"))
    delta_P <- emp_q2 - q2
    C_theta <- delta_P / (q2 * log(q1 / q2))
    beta_frac_theta <- C_theta * log(q1 / p)            # = M0 * fractional gap
    beta_frac_pos <- max(0, max(beta_frac_theta))
    beta_frac_neg <- max(0, -min(beta_frac_theta))

    # --- Stage 5b: marginal pad; kappa = 1/eta under AI (interior Euler sum) ---
    M0 <- log(q1 / p) / log(q1 / q2)
    xq2_1 <- qblend(q2, fit1); xq2_2 <- qblend(q2, fit2)
    B_marg1 <- mean(dat[,1] > xq2_1) / q2 - 1
    B_marg2 <- mean(dat[,2] > xq2_2) / q2 - 1
    b_marg  <- (1 / eta_hat) * M0 * max(abs(B_marg1), abs(B_marg2))

    # --- Stage 4: bootstrap with margins refit + eta re-estimated per draw ---
    boot_probs_mat <- matrix(NA_real_, nrow = B, ncol = length(thetas))
    n_fallback <- 0L
    n_clamp    <- 0L
    zmax <- 2 * n_dat

    if (progress_bar) pbb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (k in 1:B) {
        idx <- sample.int(n_dat, n_dat, replace = TRUE)
        d1 <- dat[idx, 1]; d2 <- dat[idx, 2]

        f1b <- fit_marginal(d1, gamma1, fallback = fit1)
        f2b <- fit_marginal(d2, gamma1, fallback = fit2)
        n_fallback <- n_fallback + f1b$used_fallback + f2b$used_fallback

        Zb <- standardize(d1, d2, f1b, f2b)
        if (any(Zb >= zmax - 1e-9)) n_clamp <- n_clamp + 1L

        # eta re-estimated on the refit-standardized minima (carries copula +
        # marginal-fit noise); transition weight function is fixed (geometry)
        eta_b <- estimate_xi_hill(pmin(Zb[,1], Zb[,2]), gamma1)

        zx <- standardize(x_iso[,1], x_iso[,2], f1b, f2b)
        probs_b <- std_blended_surv_hrv(zx, Zb, qn, eta_b, beta_cooley)
        boot_probs_mat[k, ] <- log(pmax(probs_b, prob_floor))
        if (progress_bar) utils::setTxtProgressBar(pbb, k)
    }
    if (progress_bar) close(pbb)

    # --- centered log deviations, sup/inf, quantiles ---
    mean_log_boot <- colMeans(boot_probs_mat)
    boot_log_devs <- sweep(boot_probs_mat, MARGIN = 2, STATS = mean_log_boot, FUN = "-")
    boot_draws_Wplus  <- apply(boot_log_devs, 1, max)
    boot_draws_Wminus <- apply(boot_log_devs, 1, min)

    c_plus_estimates <- list(); c_minus_estimates <- list()
    raw_c_plus_estimates <- list(); raw_c_minus_estimates <- list()
    size_estimates <- list()

    for (i in 1:length(alphas)) {
        alpha <- alphas[i]
        akey <- as.character(alpha)

        d_plus_log  <- as.numeric(quantile(boot_draws_Wplus,  probs = 1 - (alpha/2)))
        d_minus_log <- as.numeric(quantile(boot_draws_Wminus, probs = (alpha/2)))

        raw_c_plus_estimates[[akey]]  <- p * exp(d_plus_log)  - p
        raw_c_minus_estimates[[akey]] <- p * exp(d_minus_log) - p

        cplus  <- p * exp(d_plus_log  + beta_frac_neg + b_marg) - p
        cminus <- p * exp(d_minus_log - beta_frac_pos - b_marg) - p

        c_plus_estimates[[akey]]  <- cplus
        c_minus_estimates[[akey]] <- cminus

        # --- tube size (bias-corrected walls); wall inversion uses xi_eff(theta) ---
        s_lo <- p + cminus
        if (s_lo <= 0) {
            size_estimates[[akey]] <- data.frame(
                is_unbounded = TRUE,
                width_abs_median = NA_real_, width_abs_max = NA_real_,
                width_rel_median = NA_real_, width_rel_max = NA_real_)
        } else {
            s_hi <- min(p + cplus, 0.999 * qn)
            r_in_std  <- rp_std * (p / s_hi)^xi_eff_grid
            r_out_std <- rp_std * (p / s_lo)^xi_eff_grid

            inner_iso <- backmap(cbind(r_in_std * cos(thetas),  r_in_std * sin(thetas)),  fit1, fit2)
            outer_iso <- backmap(cbind(r_out_std * cos(thetas), r_out_std * sin(thetas)), fit1, fit2)

            rad <- function(iso) sqrt((pmax(iso[,1], lbs[1]) - lbs[1])^2 +
                                      (pmax(iso[,2], lbs[2]) - lbs[2])^2)
            r_inner <- rad(inner_iso)
            r_outer <- rad(outer_iso)
            r_ref   <- rad(x_iso)

            sz <- computeTubeSize(r_inner, r_outer, r_ref)
            size_estimates[[akey]] <- cbind(data.frame(is_unbounded = FALSE), sz)
        }
    }

    # composed clamped survival estimator (Cooley HRV on the standardized scale)
    survFunc <- function(pts) {
        pts <- as.matrix(pts)
        zq <- standardize(pts[,1], pts[,2], fit1, fit2)
        std_blended_surv_hrv(zq, Z, qn, eta_hat, beta_cooley)
    }

    res_lst <- list(
        dat = dat, p = p, gamma1 = gamma1, gamma2 = gamma2,
        fit1 = fit1[c("u","sigma","xi","q1")], fit2 = fit2[c("u","sigma","xi","q1")],
        eta_hat = eta_hat, beta_cooley = beta_cooley,
        thetas = thetas, xi_eff = xi_eff_grid, r_naive = rp_std, x_iso = x_iso,
        beta_frac_pos = beta_frac_pos, beta_frac_neg = beta_frac_neg,
        b_marg = b_marg, B_marg1 = B_marg1, B_marg2 = B_marg2,
        c_plus_estimates = c_plus_estimates, c_minus_estimates = c_minus_estimates,
        raw_c_plus_estimates = raw_c_plus_estimates, raw_c_minus_estimates = raw_c_minus_estimates,
        size_estimates = size_estimates,
        fallback_frac = n_fallback / (2 * B),
        clamp_frac    = n_clamp / B,
        survFunc = survFunc,
        transformed = FALSE
    )

    return(res_lst)
}

# ======================================================================
# Main simulation loop
# ======================================================================
for (i in 1:length(dists)) {
    dist <- dists[i]
    sampling_func <- sample_bivgauss

    if (!is.null(save_full_path)) {
        distpath <- paste0(save_full_path, dist, '/')
        if (!dir.exists(distpath)) dir.create(distpath, recursive = TRUE)
    }
    print(paste0('Starting new distribution: ', dist))

    for (k in 1:length(ns)) {
        n <- ns[k]
        p <- pn(n)

        # TRUE p-isoline: the distribution IS the bivariate Gaussian -- drawn directly
        isoline <- drawIsoline(dist = dist, numCoords = n_coords, gridUbs = ubs, gridLbs = lbs, prob = p)

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

            regions <- computeRegion_MargTransform_HRV(
                dat = dat, alphas = alphas, p = p, B = B,
                gamma1 = gamma1, gamma2 = gamma2,
                ncoords = n_coords, lbs = lbs, beta_cooley = beta_cooley)

            is_covereds <- evaluateCoverage(regions, isoline)
            regions$is_covereds <- is_covereds

            if (!is.null(n_dirpath)) {
                saveRDS(regions, file = paste0(n_dirpath, 'simulation', ind, '.RData'))
            }

            alpha_names <- names(regions$is_covereds)
            size_df <- do.call(rbind, regions$size_estimates[alpha_names])

            data.frame(
                dist = dist, p = p, n = n,
                alpha = alpha_names,
                is_covered = unlist(regions$is_covereds, use.names = FALSE),
                c_plus_adj  = unlist(regions$c_plus_estimates[alpha_names]),
                c_minus_adj = unlist(regions$c_minus_estimates[alpha_names]),
                c_plus_raw  = unlist(regions$raw_c_plus_estimates[alpha_names]),
                c_minus_raw = unlist(regions$raw_c_minus_estimates[alpha_names]),
                beta_frac_pos = regions$beta_frac_pos,
                beta_frac_neg = regions$beta_frac_neg,
                b_marg  = regions$b_marg,
                B_marg1 = regions$B_marg1,
                B_marg2 = regions$B_marg2,
                eta_hat = regions$eta_hat,
                xi_gpd1 = regions$fit1$xi,
                xi_gpd2 = regions$fit2$xi,
                fallback_frac = regions$fallback_frac,
                clamp_frac    = regions$clamp_frac,
                is_unbounded     = size_df$is_unbounded,
                width_abs_median = size_df$width_abs_median,
                width_abs_max    = size_df$width_abs_max,
                width_rel_median = size_df$width_rel_median,
                width_rel_max    = size_df$width_rel_max
            )
        }

        samp_df_list <- foreach(l = 1:n_sims, .options.snow = opts,
                .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt', 'matrixStats')) %dopar% parallelizedCode(l)

        close(pb)
        stopCluster(clust)

        if (!is.null(args$save_df_path)) {
            coverage_df <- do.call(rbind, samp_df_list)
            save_fname <- paste0(args$save_df_path, '1.1.1.2_bivgauss_gamma1_0.5_gamma2_0.67_bootmargs_cooleysmooth.csv')
            file_exists <- file.exists(save_fname)
            write.table(coverage_df, file = save_fname, sep = ",", row.names = FALSE,
                        col.names = !file_exists, append = file_exists)
        }
    }
}