# Coverage for confidence tubes (Mammen-Polonik style), extreme regime.
# ASYMMETRIC MARGINS: t_4 on margin 1 (xi=1/4), t_10 on margin 2 (xi=1/10),
# t_4-copula dependence. Tests the marginal-transform procedure when the two
# margins have DIFFERENT tail heaviness -- nonstandard on the raw scale, routine here.
# MARGINAL-TRANSFORMATION PROCEDURE, asymptotic dependence (AD) version.
#
# Pipeline per simulation:
#   1. Fit each margin: empirical body + GPD tail above u_j = F_emp^{-1}(1 - q1),
#      GPD by PLAIN UNCONSTRAINED MLE (Cooley-faithful). Fitted CDF clamped to
#      [1/(2n), 1 - 1/(2n)] (the sample's probability resolution).
#   2. Standardize to unit Pareto: Z_ij = 1 / (1 - F~_j(X_ij)). On this scale the
#      radial index is EXACTLY 1 under AD -- no Hill estimation anywhere.
#   3. Standardized-scale isoline via min-projection anchor at q1 and projection
#      r_p(theta) = r_q1(theta) * (q1/p). Back-map through the full-sample fits to
#      get FIXED original-scale monitoring points.
#   4. Variance (d+-): bootstrap with MARGINS REFIT PER DRAW (fallback to the
#      full-sample fit on MLE failure, counted); composed estimator re-evaluated
#      at the fixed monitoring points; log-scale centered sup/inf quantiles.
#   5. Bias pads: radial pad on the standardized scale (copula channel; rho -> 0
#      worst case, M0 multiplier) + marginal pad b_marg = M0 * max_j |B_j| from the
#      per-margin two-anchor gap at q2 (kappa = 1 under AD). Cross-mapped radial
#      signs; marginal pad symmetric on both walls.
#   6. Tube: p*exp(c-) <= Fhat(x) <= p*exp(c+), Fhat the composed clamped estimator.
#
# COST NOTE: margins are refit and the data re-standardized in every bootstrap
# draw, so expect ~2x the runtime of the common-marginal scripts at equal B.
#
# Jimmy Butler, July 2026

set.seed(45678)

library(dplyr); library(mvtnorm); library(mnormt)
library(foreach); library(doSNOW); library(parallel)
library(argparse); library(matrixStats)

source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')            # computeTubeSize lives here
source('~/isolines_uq/scripts/R/confidence_regions/modules/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/coverageEvaluation.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')

parser <- ArgumentParser(description = "Marginal-transform AD tube sims: bivariate t_4, GPD margins, Pareto standardization.")
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

# --- marginal degrees of freedom (ASYMMETRIC: heavy vs. light margin) + copula ---
df1 <- 4                       # margin 1: t_4  (heavier tail, xi_1 = 1/4)
df2 <- 10                      # margin 2: t_10 (lighter tail, xi_2 = 1/10)
copula_df <- 4                 # t-copula degrees of freedom (asymptotically DEPENDENT)
dists <- c('t4copula_t4t10margs')
lbs <- c(0, 0)
ubs <- c(200, 200)

n_cores <- args$n_cores
save_full_path <- args$save_full_path

# sampler: t_4-copula with t_4 and t_10 margins (different marginal tail heaviness).
# Build the copula from a bivariate t_4, map each coordinate to its copula uniform
# via the t_4 CDF, then push through the per-margin t quantile functions.
sample_t4copula_t4t10 <- function(n) {
    Sigma <- matrix(c(1, 0.7, 0.7, 1), 2, 2)
    Tm <- mvtnorm::rmvt(n, sigma = Sigma, df = copula_df)
    U  <- pt(Tm, df = copula_df)
    X1 <- qt(U[,1], df = df1)
    X2 <- qt(U[,2], df = df2)
    data.frame(X1 = X1, X2 = X2)
}

# ======================================================================
# Marginal machinery: unconstrained GPD MLE + blended (empirical/GPD) CDF
# with the [1/(2n), 1 - 1/(2n)] clamp, forward and inverse transforms.
# ======================================================================

# negative log-likelihood, parametrized (log sigma, xi); support violations
# penalized (that IS the likelihood being -Inf, not an added constraint)
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

# plain MLE; moment-estimator start; returns moment estimates with failed=TRUE
# if optim errors/doesn't converge (caller decides the fallback)
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

# fit one margin: empirical body + GPD tail above the (1 - q1) empirical quantile.
# On MLE failure, use `fallback` (the full-sample fit) if provided, else the
# moment estimates. used_fallback is counted by the caller.
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

# forward transform: blended fitted CDF, clamped to [1/(2n), 1 - 1/(2n)]
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

# inverse transform: survival prob s -> original-scale quantile
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

# Pareto standardization of a data matrix / point matrix
standardize <- function(x1, x2, fit1, fit2) {
    cbind(1 / (1 - ptilde(x1, fit1)), 1 / (1 - ptilde(x2, fit2)))
}

# back-map standardized points to the original scale (T^{-1} then qblend)
backmap <- function(zpts, fit1, fit2) {
    s1 <- pmin(1 / pmax(zpts[, 1], 1), 1)
    s2 <- pmin(1 / pmax(zpts[, 2], 1), 1)
    cbind(qblend(s1, fit1), qblend(s2, fit2))
}

# ======================================================================
# Standardized-scale machinery (xi = 1 known under AD)
# ======================================================================

# radial anchor at level qlev along arbitrary angles (axis-aware; Z > 0 always)
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
# xi = 1 radial projection below
std_blended_surv <- function(zpts, Z, qn) {
    zpts <- as.matrix(zpts)
    emp <- colMeans(outer(Z[,1], zpts[,1], ">") & outer(Z[,2], zpts[,2], ">"))
    theta <- atan2(zpts[,2], zpts[,1])
    rp <- sqrt(rowSums(zpts^2))
    rq <- std_radial_anchor(Z, theta, qn)
    tail_probs <- (rq / rp) * qn                 # xi = 1: (rq/rp)^{1/xi} = rq/rp
    ifelse(emp >= qn, emp, tail_probs)
}

# ======================================================================
# Tube constructor: marginal transform + AD standardized tube,
# margins-in-the-bootstrap, radial + marginal bias pads
# ======================================================================
computeRegion_MargTransform_AD <- function(dat, alphas, p, B, gamma1 = 1/2, gamma2 = 2/3,
                                           ncoords = 50, lbs = c(0, 0),
                                           progress_bar = FALSE, verbose = FALSE) {

    n_dat <- nrow(dat)
    q1 <- n_dat^(-gamma1); q2 <- n_dat^(-gamma2); qn <- q1
    prob_floor <- 0.5 / n_dat

    # --- Stage 1: full-sample marginal fits (plain MLE; clamp inside ptilde) ---
    fit1 <- fit_marginal(dat[,1], gamma1)
    fit2 <- fit_marginal(dat[,2], gamma1)

    # --- Stage 2: standardize the data ---
    Z <- standardize(dat[,1], dat[,2], fit1, fit2)

    # --- Stage 3: standardized isoline (xi = 1) + fixed monitoring points ---
    thetas <- seq(0, pi/2, length.out = ncoords)
    rq1_vec <- std_radial_anchor(Z, thetas, q1)
    rp_std  <- rq1_vec * (q1 / p)                       # projection, known xi = 1
    ziso    <- cbind(rp_std * cos(thetas), rp_std * sin(thetas))
    x_iso   <- backmap(ziso, fit1, fit2)                # FIXED original-scale points

    # --- Stage 5a: radial pad (copula channel), standardized scale ---
    rq2_vec <- rq1_vec * (q1 / q2)
    zq2 <- cbind(rq2_vec * cos(thetas), rq2_vec * sin(thetas))
    emp_q2 <- colMeans(outer(Z[,1], zq2[,1], ">") & outer(Z[,2], zq2[,2], ">"))
    delta_P <- emp_q2 - q2
    C_theta <- delta_P / (q2 * log(q1 / q2))
    beta_frac_theta <- C_theta * log(q1 / p)            # = M0 * fractional gap
    beta_frac_pos <- max(0, max(beta_frac_theta))
    beta_frac_neg <- max(0, -min(beta_frac_theta))

    # --- Stage 5b: marginal pad (transform channel), two-anchor per margin ---
    M0 <- log(q1 / p) / log(q1 / q2)
    xq2_1 <- qblend(q2, fit1); xq2_2 <- qblend(q2, fit2)
    B_marg1 <- mean(dat[,1] > xq2_1) / q2 - 1
    B_marg2 <- mean(dat[,2] > xq2_2) / q2 - 1
    b_marg  <- M0 * max(abs(B_marg1), abs(B_marg2))     # kappa = 1 under AD

    # --- Stage 4: bootstrap with margins refit per draw ---
    boot_probs_mat <- matrix(NA_real_, nrow = B, ncol = ncoords)
    n_fallback <- 0L
    n_clamp    <- 0L
    zmax <- 2 * n_dat                                    # clamped upper standardized value

    if (progress_bar) pbb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (k in 1:B) {
        idx <- sample.int(n_dat, n_dat, replace = TRUE)
        d1 <- dat[idx, 1]; d2 <- dat[idx, 2]

        f1b <- fit_marginal(d1, gamma1, fallback = fit1)
        f2b <- fit_marginal(d2, gamma1, fallback = fit2)
        n_fallback <- n_fallback + f1b$used_fallback + f2b$used_fallback

        Zb <- standardize(d1, d2, f1b, f2b)
        if (any(Zb >= zmax - 1e-9)) n_clamp <- n_clamp + 1L   # upper clamp fired this draw

        # standardized coordinates of the FIXED original-scale monitoring points
        zx <- standardize(x_iso[,1], x_iso[,2], f1b, f2b)

        probs_b <- std_blended_surv(zx, Zb, qn)
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

        # cross-mapped radial pads + symmetric marginal pad
        cplus  <- p * exp(d_plus_log  + beta_frac_neg + b_marg) - p
        cminus <- p * exp(d_minus_log - beta_frac_pos - b_marg) - p

        c_plus_estimates[[akey]]  <- cplus
        c_minus_estimates[[akey]] <- cminus

        # --- tube size (bias-corrected walls) ---
        s_lo <- p + cminus
        if (s_lo <= 0) {   # inert on the log scale; kept for interface parity
            size_estimates[[akey]] <- data.frame(
                is_unbounded = TRUE,
                width_abs_median = NA_real_, width_abs_max = NA_real_,
                width_rel_median = NA_real_, width_rel_max = NA_real_)
        } else {
            # closed-form standardized wall radii (xi = 1); cap the inner wall at qn
            # so the inversion stays in the projection branch
            s_hi <- min(p + cplus, 0.999 * qn)
            r_in_std  <- rp_std * (p / s_hi)
            r_out_std <- rp_std * (p / s_lo)

            inner_iso <- backmap(cbind(r_in_std * cos(thetas),  r_in_std * sin(thetas)),  fit1, fit2)
            outer_iso <- backmap(cbind(r_out_std * cos(thetas), r_out_std * sin(thetas)), fit1, fit2)

            # original-scale radii from lbs; clamp coordinates into the quadrant so
            # negative empirical minima at the axes don't pollute the radii
            rad <- function(iso) sqrt((pmax(iso[,1], lbs[1]) - lbs[1])^2 +
                                      (pmax(iso[,2], lbs[2]) - lbs[2])^2)
            r_inner <- rad(inner_iso)
            r_outer <- rad(outer_iso)
            r_ref   <- rad(x_iso)

            sz <- computeTubeSize(r_inner, r_outer, r_ref)
            size_estimates[[akey]] <- cbind(data.frame(is_unbounded = FALSE), sz)
        }
    }

    # composed clamped survival estimator (the estimator of record)
    survFunc <- function(pts) {
        pts <- as.matrix(pts)
        zq <- standardize(pts[,1], pts[,2], fit1, fit2)
        std_blended_surv(zq, Z, qn)
    }

    res_lst <- list(
        dat = dat, p = p, gamma1 = gamma1, gamma2 = gamma2,
        fit1 = fit1[c("u","sigma","xi","q1")], fit2 = fit2[c("u","sigma","xi","q1")],
        thetas = thetas, r_naive = rp_std, x_iso = x_iso,
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
    sampling_func <- sample_t4copula_t4t10

    if (!is.null(save_full_path)) {
        distpath <- paste0(save_full_path, dist, '/')
        if (!dir.exists(distpath)) dir.create(distpath, recursive = TRUE)
    }
    print(paste0('Starting new distribution: ', dist))

    for (k in 1:length(ns)) {
        n <- ns[k]
        p <- pn(n)

        # TRUE p-isoline for the bivariate t_4 (elliptical, drawn directly)
        # TRUE p-isoline: t_4-copula survival isoline -> pt(copula_df) -> per-margin t quantiles.
        # The copula is the same t_4 used to generate the data; margins are t_4 and t_10.
        tcop_isoline   <- drawIsoline(dist='bivt', numCoords=n_coords, gridUbs=ubs, gridLbs=lbs, prob=p)
        copula_isoline <- data.frame(X1 = pt(tcop_isoline[,1], df=copula_df),
                                     X2 = pt(tcop_isoline[,2], df=copula_df))
        isoline <- data.frame(X1 = qt(copula_isoline[,1], df=df1),
                              X2 = qt(copula_isoline[,2], df=df2))

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

            regions <- computeRegion_MargTransform_AD(
                dat = dat, alphas = alphas, p = p, B = B,
                gamma1 = gamma1, gamma2 = gamma2,
                ncoords = n_coords, lbs = lbs)

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
            save_fname <- paste0(args$save_df_path, '0.1.1.2_bivt4_t4t10margins_gamma1_0.5_gamma2_0.67_bootmargs.csv')
            file_exists <- file.exists(save_fname)
            write.table(coverage_df, file = save_fname, sep = ",", row.names = FALSE,
                        col.names = !file_exists, append = file_exists)
        }
    }
}