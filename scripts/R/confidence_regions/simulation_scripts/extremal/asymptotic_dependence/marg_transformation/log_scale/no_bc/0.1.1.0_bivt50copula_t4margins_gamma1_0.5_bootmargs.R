# Coverage evaluation for the NON-BIAS-CORRECTED (raw, variance-only) walls of
# the COPULA-LADDER slow-dependence cell: t_50 copula with t_4/t_4 margins.
#
# This is the structural-channel stress test read raw: at accessible levels the
# t_50 dependence is near-Gaussian (interior exponent ~ 0.85) while the AD
# estimator projects with frozen exponent 1 -- a dependence-channel bias no
# marginal refit absorbs, and exactly what the M0 pad exists to cover. The raw
# rows here quantify how much coverage that pad is buying.
#
# Loads each saved tube object from the 0.1.1.2 runs, reconstructs the survival
# estimator from the stored trimmed fits + data (the saved survFunc closure
# does not survive readRDS), swaps the RAW walls (raw_c_plus/minus_estimates)
# into the evaluation slots, and re-evaluates against the same true isoline
# (drawBivtIsoline with df = copula_df -> pt -> t_4 quantiles). Machinery
# inlined VERBATIM from the simulation script.
#
# Jimmy Butler, July 2026

library(dplyr); library(matrixStats); library(mvtnorm)   # pmvt used by drawBivtIsoline
source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')             # drawBivtIsoline
source('~/isolines_uq/scripts/R/confidence_regions/modules/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/coverageEvaluation.R')

# NOTE: point save_full_path at wherever the 0.1.1.2 t50-copula tubes were saved.
save_full_path <- '/pscratch/sd/j/jbbutler/isolines/simulations/extremal/0.1.1.2_bivt50copula_t4margins_gamma1_0.5_gamma2_0.67_bootmargs/'
save_df_path   <- '~/isolines_uq/outputs/simulations/extremal/'

ns <- c(1000, 5000, 10000, 50000, 100000)
alphas <- c(0.05, 0.1, 0.01)
C <- 5; pn <- function(n){ C/n }
n_coords <- 50

copula_df <- 50                # t-copula degrees of freedom (SLOW second order)
marg_df   <- 4                 # both margins: t_4 (validated heavy-tailed cell)
dists <- c('t50copula_t4t4margs')
lbs <- c(0, 0)
ubs <- c(200, 200)

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

# reconstruct a full marginal fit (add x_sorted, n) from a trimmed stored fit + data
rebuild_fit <- function(fit_trim, x) {
    fit_trim$n <- length(x)
    fit_trim$x_sorted <- sort(x)
    fit_trim
}

# ======================================================================
# Evaluation loop
# ======================================================================
results_list <- list()
print(paste0('Evaluating RAW tubes: ', dists[1]))

for (dist in dists) {
  for (n in ns) {
    p <- pn(n)
    print(paste0('Processing ', dist, ' | n = ', n))

        # TRUE p-isoline: t_50 elliptical survival isoline via the module's
        # drawBivtIsoline (same constructor as all existing bivt truth curves,
        # df passed explicitly) -> pt(copula_df) -> t_4 quantiles on both margins.
        tcop_isoline   <- drawBivtIsoline(numCoords = n_coords, gridUbs = ubs,
                                          gridLbs = lbs, prob = p, df = copula_df)
        copula_isoline <- data.frame(X1 = pt(tcop_isoline[,1], df = copula_df),
                                     X2 = pt(tcop_isoline[,2], df = copula_df))
        isoline <- data.frame(X1 = qt(copula_isoline[,1], df = marg_df),
                              X2 = qt(copula_isoline[,2], df = marg_df))

    sim_dir <- paste0(save_full_path, dist, '/n', n, '/')
    sim_files <- list.files(sim_dir, pattern = "\\.RData$", full.names = TRUE)
    if (length(sim_files) == 0) { warning(paste0("No .RData in ", sim_dir)); next }
    pb <- txtProgressBar(min = 0, max = length(sim_files), style = 3)

    for (i in seq_along(sim_files)) {
      regions <- readRDS(sim_files[i])
      dat <- regions$dat
      fit1 <- rebuild_fit(regions$fit1, dat[,1])
      fit2 <- rebuild_fit(regions$fit2, dat[,2])
      qn <- regions$fit1$q1
      Z  <- standardize(dat[,1], dat[,2], fit1, fit2)

      survFunc <- function(pts) {
          pts <- as.matrix(pts)
          zq <- standardize(pts[,1], pts[,2], fit1, fit2)
          std_blended_surv(zq, Z, qn)
      }

      obj <- list(survFunc = survFunc, p = p,
                  c_plus_estimates  = regions$raw_c_plus_estimates,
                  c_minus_estimates = regions$raw_c_minus_estimates,
                  transformed = FALSE)
      covs <- evaluateCoverage(obj, isoline)
      ak <- names(covs)

      results_list[[length(results_list) + 1]] <- data.frame(
          dist = dist, p = p, n = n, alpha = ak,
          is_covered = unlist(covs, use.names = FALSE),
          c_plus_raw  = unlist(regions$raw_c_plus_estimates[ak]),
          c_minus_raw = unlist(regions$raw_c_minus_estimates[ak]),
          b_marg = regions$b_marg, B_marg1 = regions$B_marg1, B_marg2 = regions$B_marg2,
          beta_frac_pos = regions$beta_frac_pos, beta_frac_neg = regions$beta_frac_neg,
          xi_gpd1 = regions$fit1$xi, xi_gpd2 = regions$fit2$xi,
          fallback_frac = regions$fallback_frac, clamp_frac = regions$clamp_frac)
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
}

if (length(results_list) > 0) {
    coverage_df <- do.call(rbind, results_list)
    save_fname <- paste0(save_df_path, '0.1.1.0_bivt50copula_t4margs_gamma1_0.5_bootmargs.csv')
    write.table(coverage_df, file = save_fname, sep = ",", row.names = FALSE, col.names = TRUE)
    print(paste0("Saved: ", save_fname))
} else print("No results; check .RData directories.")