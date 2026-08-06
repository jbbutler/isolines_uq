# 1.1.1.0 raw-walls AI tubes (bivgauss, n = 100000 ONLY): pure-HRV radial
# construction on the INTERIOR cone (constant exponent eta_hat) -- NO bias
# correction anywhere (no beta pad, no marginal pad): walls are centered
# bootstrap quantiles only.
#
# Method rows per simulation:
#   hrv_full      : raw HRV-interior tube vs the FULL true isoline
#   hrv_interior  : same tube vs the interior truth subset (true pseudo-angle
#                   in [w_int, 1 - w_int], w_int = 0.02) -- the cut the tube
#                   is designed for; hrv_full is the arm-inclusive stress row.
#
# Construction matches the dual-tube 1.1.1.2 sims' interior tube verbatim
# (theta grid on [atan(w_int/(1-w_int)), pi/2 - .], constant-eta radial
# projection, std_blended_surv_hrvconst), so rows are directly comparable.
#
# Jimmy Butler, July 2026

set.seed(45678)

library(dplyr); library(mvtnorm); library(mnormt)
library(foreach); library(doSNOW); library(parallel)
library(argparse); library(matrixStats)

source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')            # estimate_xi_hill
source('~/isolines_uq/scripts/R/confidence_regions/modules/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/coverageEvaluation.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')

parser <- ArgumentParser(description = "Percoord-only raw-wall AI tube, bivariate Gaussian, n = 100000.")
parser$add_argument("--save_full_path", type = "character", default = NULL, required = FALSE)
parser$add_argument("--n_cores", type = "integer", default = 64)
parser$add_argument("--save_df_path", type = "character", default = NULL)
args <- parser$parse_args()

if (is.null(args$save_full_path) && is.null(args$save_df_path)) {
    stop("Error: Both --save_full_path and --save_df_path are NULL.")
}

# --- simulation parameters ---
ns <- c(100000)
alphas <- c(0.05, 0.1, 0.01)
C <- 5
gamma1 <- 0.5
B <- 5000
n_sims <- 500
n_coords <- 50

pn <- function(n){ C/n }

# --- distribution / method parameters ---
rho_corr <- 0.7
w_int <- 0.02                      # interior pseudo-angle cut (truth subset only)

dists <- c('bivgauss')
lbs <- c(0, 0)
ubs <- c(200, 200)

n_cores <- args$n_cores
save_full_path <- args$save_full_path

sample_bivgauss <- function(n) {
    Sigma <- matrix(c(1, rho_corr, rho_corr, 1), 2, 2)
    X <- mvtnorm::rmvnorm(n, mean = c(0, 0), sigma = Sigma)
    data.frame(X1 = X[,1], X2 = X[,2])
}

### BEGIN CORE #########################################################

# hybrid theta grid (full-space tube): uniform interior + log-spaced points
# hugging each axis so the transition zone is monitored

# ---------- marginal machinery (unconstrained GPD MLE + clamped blend) -------
gpd_nll <- function(par, exc) {
    sigma <- exp(par[1]); xi <- par[2]
    z <- 1 + xi * exc / sigma
    if (any(z <= 0) || !is.finite(sigma)) return(1e10)
    if (abs(xi) < 1e-8) sum(log(sigma) + exc / sigma)
    else sum(log(sigma) + (1/xi + 1) * log(z))
}
fit_gpd_mle <- function(exc) {
    m <- mean(exc); v <- var(exc)
    xi0 <- 0.5 * (1 - m^2 / v)
    if (!is.finite(xi0)) xi0 <- 0.1
    xi0 <- max(min(xi0, 0.9), -0.4)
    sigma0 <- m * (1 - xi0); if (!is.finite(sigma0) || sigma0 <= 0) sigma0 <- m
    out <- tryCatch(optim(c(log(sigma0), xi0), gpd_nll, exc = exc, method = "Nelder-Mead"),
                    error = function(e) NULL)
    if (is.null(out) || out$convergence != 0 || !is.finite(out$par[1]))
        return(list(sigma = sigma0, xi = xi0, failed = TRUE))
    list(sigma = exp(out$par[1]), xi = out$par[2], failed = FALSE)
}
fit_marginal <- function(x, gamma1, fallback = NULL) {
    n <- length(x); q1 <- n^(-gamma1); xs <- sort(x)
    u <- quantile(x, probs = 1 - q1, type = 1, names = FALSE)
    exc <- x[x > u] - u
    g <- fit_gpd_mle(exc)
    used_fallback <- FALSE
    if (g$failed && !is.null(fallback)) { g <- list(sigma = fallback$sigma, xi = fallback$xi); used_fallback <- TRUE }
    list(u = u, sigma = g$sigma, xi = g$xi, q1 = q1, n = n, x_sorted = xs, used_fallback = used_fallback)
}
ptilde <- function(t, fit) {
    n <- fit$n
    Femp <- findInterval(t, fit$x_sorted) / n
    y <- pmax(t - fit$u, 0)
    if (abs(fit$xi) < 1e-8) Fgpd <- 1 - fit$q1 * exp(-y / fit$sigma)
    else Fgpd <- 1 - fit$q1 * pmax(1 + fit$xi * y / fit$sigma, 0)^(-1/fit$xi)
    Fv <- ifelse(t <= fit$u, Femp, Fgpd)
    pmin(pmax(Fv, 1/(2*n)), 1 - 1/(2*n))
}
qblend <- function(s, fit) {
    emp_idx <- pmax(1L, pmin(fit$n, as.integer(ceiling((1 - s) * fit$n))))
    x_emp <- fit$x_sorted[emp_idx]
    if (abs(fit$xi) < 1e-8) x_gpd <- fit$u + fit$sigma * log(fit$q1 / s)
    else x_gpd <- fit$u + (fit$sigma / fit$xi) * ((s / fit$q1)^(-fit$xi) - 1)
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
std_radial_anchor <- function(Z, thetas, qlev) {
    rq <- rep(NA_real_, length(thetas))
    ax0  <- which(thetas == 0); ax90 <- which(thetas == pi/2)
    if (length(ax0))  rq[ax0]  <- quantile(Z[,1], probs = 1 - qlev, type = 1, names = FALSE)
    if (length(ax90)) rq[ax90] <- quantile(Z[,2], probs = 1 - qlev, type = 1, names = FALSE)
    intr <- which(!(thetas == 0 | thetas == pi/2))
    if (length(intr)) {
        M <- pmin(Z[,1] %o% (1/cos(thetas[intr])), Z[,2] %o% (1/sin(thetas[intr])))
        rq[intr] <- colQuantiles(M, probs = 1 - qlev, type = 1, drop = TRUE)
    }
    rq
}


std_blended_surv_hrvconst <- function(zpts, Z, qn, eta_hat) {
    zpts <- as.matrix(zpts)
    emp <- colMeans(outer(Z[, 1], zpts[, 1], ">") & outer(Z[, 2], zpts[, 2], ">"))
    theta <- atan2(zpts[, 2], zpts[, 1])
    rp <- sqrt(rowSums(zpts^2))
    rq <- std_radial_anchor(Z, theta, qn)
    tail_probs <- ((rq / rp)^(1 / eta_hat)) * qn
    ifelse(emp >= qn, emp, tail_probs)
}

# ---------- raw-walls-only HRV-interior tube constructor ----------------------
computeRegion_HRV_Raw <- function(dat, alphas, p, B, gamma1 = 1/2,
                                  ncoords = 50, lbs = c(0, 0), w_int = 0.02) {
    n_dat <- nrow(dat)
    q1 <- n_dat^(-gamma1); qn <- q1
    prob_floor <- 0.5 / n_dat

    fit1 <- fit_marginal(dat[,1], gamma1)
    fit2 <- fit_marginal(dat[,2], gamma1)
    Z <- standardize(dat[,1], dat[,2], fit1, fit2)
    eta_hat <- estimate_xi_hill(pmin(Z[,1], Z[,2]), gamma1)

    ## ---- pure-HRV radial tube on the interior cone (constant eta) ----
    theta_cut <- atan(w_int / (1 - w_int))
    thetas_int <- seq(theta_cut, pi/2 - theta_cut, length.out = ncoords)
    rq1_int <- std_radial_anchor(Z, thetas_int, q1)
    rp_int  <- rq1_int * (q1 / p)^eta_hat
    ziso_int  <- cbind(rp_int * cos(thetas_int), rp_int * sin(thetas_int))
    x_iso_int <- backmap(ziso_int, fit1, fit2)

    ## ---- bootstrap: raw (variance-only) walls ----
    m_int <- nrow(ziso_int)
    boot_int <- matrix(NA_real_, B, m_int)
    n_fallback <- 0L; n_clamp <- 0L; zmax <- 2 * n_dat
    for (k in 1:B) {
        idx <- sample.int(n_dat, n_dat, replace = TRUE)
        d1 <- dat[idx, 1]; d2 <- dat[idx, 2]
        f1b <- fit_marginal(d1, gamma1, fallback = fit1)
        f2b <- fit_marginal(d2, gamma1, fallback = fit2)
        n_fallback <- n_fallback + f1b$used_fallback + f2b$used_fallback
        Zb <- standardize(d1, d2, f1b, f2b)
        if (any(Zb >= zmax - 1e-9)) n_clamp <- n_clamp + 1L
        eta_b <- estimate_xi_hill(pmin(Zb[,1], Zb[,2]), gamma1)
        zx_int <- standardize(x_iso_int[,1], x_iso_int[,2], f1b, f2b)
        boot_int[k, ] <- log(pmax(std_blended_surv_hrvconst(zx_int, Zb, qn, eta_b), prob_floor))
    }

    devs <- sweep(boot_int, 2, colMeans(boot_int), FUN = "-")
    Wp <- apply(devs, 1, max); Wm <- apply(devs, 1, min)
    rcps <- list(); rcms <- list()
    for (alpha in alphas) {
        ak <- as.character(alpha)
        dp <- as.numeric(quantile(Wp, probs = 1 - alpha/2))
        dm <- as.numeric(quantile(Wm, probs = alpha/2))
        rcps[[ak]] <- p * exp(dp) - p
        rcms[[ak]] <- p * exp(dm) - p
    }

    sf_hrv <- function(pts) { pts <- as.matrix(pts)
        zq <- standardize(pts[,1], pts[,2], fit1, fit2)
        std_blended_surv_hrvconst(zq, Z, qn, eta_hat) }

    list(dat = dat, p = p, gamma1 = gamma1,
         fit1 = fit1[c("u","sigma","xi","q1")], fit2 = fit2[c("u","sigma","xi","q1")],
         eta_hat = eta_hat, w_int = w_int,
         fallback_frac = n_fallback / (2 * B), clamp_frac = n_clamp / B,
         interior = list(raw_c_plus_estimates = rcps, raw_c_minus_estimates = rcms,
                         survFunc = sf_hrv, thetas = thetas_int, x_iso = x_iso_int,
                         projection = "hrv_const"),
         transformed = FALSE)
}
### END CORE ###########################################################

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

        # TRUE p-isoline (full) + interior subset via TRUE pseudo-angle
        isoline <- drawIsoline(dist = dist, numCoords = n_coords, gridUbs = ubs, gridLbs = lbs, prob = p)
        Zt <- cbind(1 / (1 - pnorm(isoline[,1])), 1 / (1 - pnorm(isoline[,2])))
        w_true <- Zt[,1] / (Zt[,1] + Zt[,2])
        keep_int <- w_true >= w_int & w_true <= 1 - w_int
        isoline_int <- isoline[keep_int, , drop = FALSE]
        print(paste0('n = ', n, ' : interior truth points = ', sum(keep_int), ' / ', n_coords))
        if (sum(keep_int) < 5) warning("Fewer than 5 interior truth points; widen n_coords or w_int.")

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
            regions <- computeRegion_HRV_Raw(dat = dat, alphas = alphas, p = p, B = B,
                                                  gamma1 = gamma1, ncoords = n_coords,
                                             lbs = lbs, w_int = w_int)

            # raw-wall coverage: same tube vs full truth and interior truth
            ev_raw <- function(iso) {
                obj <- list(survFunc = regions$interior$survFunc, p = p,
                            c_plus_estimates  = regions$interior$raw_c_plus_estimates,
                            c_minus_estimates = regions$interior$raw_c_minus_estimates,
                            transformed = FALSE)
                evaluateCoverage(obj, iso)
            }
            cov_hrv_full <- ev_raw(isoline)
            cov_hrv_int  <- ev_raw(isoline_int)

            if (!is.null(n_dirpath)) {
                regions_save <- regions
                regions_save$interior$survFunc <- NULL   # closures don't survive readRDS; dropping keeps files small
                saveRDS(regions_save, file = paste0(n_dirpath, 'simulation', ind, '.RData'))
            }

            mk_rows <- function(method, covs, tube, n_truth) {
                ak <- names(covs)
                data.frame(dist = dist, p = p, n = n, method = method, alpha = ak,
                    is_covered = unlist(covs, use.names = FALSE),
                    c_plus_raw  = unlist(tube$raw_c_plus_estimates[ak]),
                    c_minus_raw = unlist(tube$raw_c_minus_estimates[ak]),
                    b_marg = NA_real_, B_marg1 = NA_real_, B_marg2 = NA_real_,
                    eta_hat = regions$eta_hat, xi_gpd1 = regions$fit1$xi, xi_gpd2 = regions$fit2$xi,
                    fallback_frac = regions$fallback_frac, clamp_frac = regions$clamp_frac,
                    n_truth_pts = n_truth)
            }
            rbind(mk_rows('hrv_full',     cov_hrv_full, regions$interior, nrow(isoline)),
                  mk_rows('hrv_interior', cov_hrv_int,  regions$interior, nrow(isoline_int)))
        }

        samp_df_list <- foreach(l = 1:n_sims, .options.snow = opts,
                .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt', 'matrixStats')) %dopar% parallelizedCode(l)

        close(pb)
        stopCluster(clust)

        if (!is.null(args$save_df_path)) {
            coverage_df <- do.call(rbind, samp_df_list)
            save_fname <- paste0(args$save_df_path, '1.1.1.0_bivgauss_gamma1_0.5_bootmargs_hrv_n100000.csv')
            file_exists <- file.exists(save_fname)
            write.table(coverage_df, file = save_fname, sep = ",", row.names = FALSE,
                        col.names = !file_exists, append = file_exists)
        }
    }
}