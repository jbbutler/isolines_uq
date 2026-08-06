# Coverage for confidence tubes (Mammen-Polonik style), extreme regime.
# MARGINAL-TRANSFORMATION PROCEDURE, AI, KARACHI BETA-KDE distribution:
# PERCOORD-ONLY, RAW-WALL simulation at n = 100000.
#
# Purpose: the 1.1.1.2 dual-tube runs (percoord full + hrv interior, with bias
# correction) timed out at n = 100000. This script fills in that sample size by
# simulating ONLY the per-coordinate (Cooley-exact) tube and evaluating coverage
# directly on the RAW variance-only walls
#     c_plus_raw  = p * exp(dp) - p,   c_minus_raw = p * exp(dm) - p,
# i.e. NO bias correction: the q2 beta-fraction terms, the marginal pad b_marg,
# and gamma2 are never computed, and the hrv-interior tube (and its bootstrap
# column block) is dropped entirely.
#
# DISTRIBUTION: the beta-kernel KDE of the Karachi temperature/relative-humidity
# data (Chen 1999, first method; karachiTools.R), support [50,140] x [0,100].
#
# MARGINAL MODELS (marg_models): per-margin transform choice, as in the 1.1.1.2
# dual-tube script.
#   gpd:       empirical body + GPD tail above the (1 - q1) quantile
#   empirical: interpolated tie-aware empirical CDF (Hazen positions averaged
#              over ties), one coherent forward/inverse pair
# SET BELOW to match the configuration of the n <= 50000 runs being extended:
# the existing 1.1.1.0_bivkarachi CSV has non-NA xi_gpd2 throughout, i.e.
# c("gpd","gpd"). Change only if deliberately switching configurations.
#
# Coverage rows written per sim (method column):
#   percoord_full     : raw percoord tube vs the full true isoline
#   percoord_interior : same raw tube vs the interior truth subset
# (Both come from the one tube; the interior evaluation is free. Delete the
# percoord_interior line in mk_rows/rbind below if only the full row is wanted.)
#
# Output schema matches the existing
#   1.1.1.0_bivkarachi_gamma1_0.5_bootmargs_cooleysmoothpercoord.csv
# exactly (b_marg, B_marg1, B_marg2 are written as NA), so the n = 100000 rows
# can be rbind-ed / bind_rows-ed onto the n <= 50000 raw-wall results.
#
# INTERIOR DEFINITION: true standardized pseudo-angle from the KDE's TRUE
# marginal CDFs (computable exactly from the kernel mixture), kept in
# [w_int, 1 - w_int], w_int = 0.02 (Cooley weight there ~ 0.018).
#
# CLAMP DIAGNOSTIC (as redefined in the 1.1.1.2 Karachi script): clamp_frac is
# the average FRACTION of standardized points at the ceiling 2n per draw.
#
# All machinery is inlined verbatim from the 1.1.1.2 Karachi simulation script.
#
# Jimmy Butler, July 2026

set.seed(45678)

library(dplyr); library(mvtnorm); library(mnormt)
library(foreach); library(doSNOW); library(parallel)
library(argparse); library(matrixStats)

source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')            # estimate_xi_hill
source('~/isolines_uq/scripts/R/confidence_regions/modules/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/coverageEvaluation.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')      # rKarachiBetaKDE, pKarachiBetaKDE, karachiDat, *_CONST

parser <- ArgumentParser(description = "Percoord-only raw-wall AI tube on the Karachi beta-KDE, n = 100000.")
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

# --- method parameters ---
beta_cooley <- 200
w_int <- 0.02
marg_models <- c("gpd", "gpd")     # (temp, relHum); see header

dists <- c('karachiBetaKDE')
lbs <- LBS_CONST                          # c(50, 0)  -- polar origin = support corner
ubs <- UBS_CONST                          # c(140, 100)

n_cores <- args$n_cores
save_full_path <- args$save_full_path

# TRUE marginal survivals of the KDE (kernel mixture), for the interior cut
true_marg_surv <- function(x, j, dat = karachiDat, b = b_CONST,
                           lb = LBS_CONST[j], ub = UBS_CONST[j]) {
    y <- (dat[, j] - lb) / (ub - lb)
    s1 <- y / b + 1; s2 <- (1 - y) / b + 1
    x01 <- (x - lb) / (ub - lb)
    sapply(x01, function(u) mean(1 - pbeta(u, shape1 = s1, shape2 = s2)))
}

### BEGIN CORE #########################################################

make_theta_grid_ai <- function(ncoords, n_axis = 8, ax_lo = 5e-4, ax_hi = 0.05) {
    ax <- exp(seq(log(ax_lo), log(ax_hi), length.out = n_axis))
    interior <- seq(0, pi/2, length.out = ncoords - 2 * n_axis)
    sort(unique(c(interior, ax, (pi/2) - ax)))
}

# ---------- marginal machinery with PER-MARGIN MODEL DISPATCH ----------------
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

# fit one margin under the requested model.
#   gpd:       empirical body + GPD tail above the (1 - q1) quantile (as before)
#   empirical: interpolated tie-aware empirical CDF -- Hazen positions
#              (i - 0.5)/n averaged within tie groups; the SAME (value, position)
#              grid serves forward (ptilde) and inverse (qblend): one coherent
#              monotone pair, clamped to [1/(2n), 1 - 1/(2n)] by construction.
fit_marginal <- function(x, gamma1, model = "gpd", fallback = NULL) {
    n <- length(x); q1 <- n^(-gamma1); xs <- sort(x)
    if (model == "empirical") {
        pos_all <- (seq_len(n) - 0.5) / n
        ux <- unique(xs)
        pos <- as.numeric(tapply(pos_all, match(xs, ux), mean))
        return(list(model = "empirical", u = NA_real_, sigma = NA_real_, xi = NA_real_,
                    q1 = q1, n = n, x_sorted = xs, ux = ux, pos = pos,
                    used_fallback = FALSE))
    }
    u <- quantile(x, probs = 1 - q1, type = 1, names = FALSE)
    exc <- x[x > u] - u
    g <- fit_gpd_mle(exc)
    used_fallback <- FALSE
    if (g$failed && !is.null(fallback)) { g <- list(sigma = fallback$sigma, xi = fallback$xi); used_fallback <- TRUE }
    list(model = "gpd", u = u, sigma = g$sigma, xi = g$xi, q1 = q1, n = n,
         x_sorted = xs, used_fallback = used_fallback)
}

ptilde <- function(t, fit) {
    n <- fit$n
    if (identical(fit$model, "empirical")) {
        Fv <- approx(fit$ux, fit$pos, xout = t, rule = 2)$y
        return(pmin(pmax(Fv, 1/(2*n)), 1 - 1/(2*n)))
    }
    Femp <- findInterval(t, fit$x_sorted) / n
    y <- pmax(t - fit$u, 0)
    if (abs(fit$xi) < 1e-8) Fgpd <- 1 - fit$q1 * exp(-y / fit$sigma)
    else Fgpd <- 1 - fit$q1 * pmax(1 + fit$xi * y / fit$sigma, 0)^(-1/fit$xi)
    Fv <- ifelse(t <= fit$u, Femp, Fgpd)
    pmin(pmax(Fv, 1/(2*n)), 1 - 1/(2*n))
}
qblend <- function(s, fit) {
    if (identical(fit$model, "empirical")) {
        return(approx(fit$pos, fit$ux, xout = 1 - s, rule = 2)$y)
    }
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

# ---------- per-coordinate (Cooley-exact) machinery --------------------------
cooley_exponents <- function(w, eta_hat, beta_cooley) {
    cbind(eta_hat + (1 - eta_hat) * w^beta_cooley,
          eta_hat + (1 - eta_hat) * (1 - w)^beta_cooley)
}
build_base_curve <- function(z1, z2) {
    keep <- is.finite(z1) & is.finite(z2) & z1 > 0 & z2 > 0
    z1 <- z1[keep]; z2 <- z2[keep]
    w <- z1 / (z1 + z2)
    o <- order(w)
    w <- w[o]; lf1 <- log(z1[o]); lf2 <- log(z2[o])
    lf1 <- cummax(lf1); lf2 <- -cummax(-lf2)
    eps <- 1e-10 * seq_along(w)
    w <- w + 1e-12 * seq_along(w)
    list(w = w, lf1 = lf1 + eps, lf2 = lf2 - eps, w_lo = w[1], w_hi = w[length(w)])
}
percoord_project <- function(z1, z2, t, eta_hat, beta_cooley) {
    w <- ifelse(z1 + z2 > 0, z1 / (z1 + z2), 0.5)
    E <- cooley_exponents(w, eta_hat, beta_cooley)
    cbind(t^E[, 1] * z1, t^E[, 2] * z2)
}
percoord_level <- function(zq1, zq2, base, q1, eta_hat, beta_cooley, n_iter = 40) {
    m <- length(zq1)
    lz1 <- log(pmax(zq1, 1e-12)); lz2 <- log(pmax(zq2, 1e-12))
    lf1_at <- function(w) approx(base$w, base$lf1, xout = w, rule = 2)$y
    lf2_at <- function(w) approx(base$w, base$lf2, xout = w, rule = 2)$y
    wB <- approx(base$lf1, base$w, xout = lz1, rule = 2)$y
    wA <- approx(rev(base$lf2), rev(base$w), xout = lz2, rule = 2)$y
    wA <- pmax(wA, base$w_lo); wB <- pmin(wB, base$w_hi)
    ok <- wA < wB
    logt <- numeric(m)
    if (any(ok)) {
        lo <- wA[ok]; hi <- wB[ok]; l1 <- lz1[ok]; l2 <- lz2[ok]
        for (it in 1:n_iter) {
            mid <- 0.5 * (lo + hi)
            E <- cooley_exponents(mid, eta_hat, beta_cooley)
            g <- E[, 2] * (l1 - lf1_at(mid)) - E[, 1] * (l2 - lf2_at(mid))
            mv <- g > 0
            lo <- ifelse(mv, mid, lo); hi <- ifelse(mv, hi, mid)
        }
        ws <- 0.5 * (lo + hi)
        E <- cooley_exponents(ws, eta_hat, beta_cooley)
        logt[ok] <- 0.5 * ((l1 - lf1_at(ws)) / E[, 1] + (l2 - lf2_at(ws)) / E[, 2])
    }
    if (any(!ok)) {
        l1 <- lz1[!ok]; l2 <- lz2[!ok]
        z1dom <- zq1[!ok] / (zq1[!ok] + zq2[!ok]) >= 0.5
        E_hi <- cooley_exponents(base$w_hi, eta_hat, beta_cooley)
        E_lo <- cooley_exponents(base$w_lo, eta_hat, beta_cooley)
        logt[!ok] <- ifelse(z1dom, (l1 - lf1_at(base$w_hi)) / E_hi[, 1],
                                   (l2 - lf2_at(base$w_lo)) / E_lo[, 2])
    }
    logt <- pmax(logt, 0)
    q1 * exp(-logt)
}
std_blended_surv_pc <- function(zpts, Z, qn, thetas, eta_hat, beta_cooley) {
    zpts <- as.matrix(zpts)
    emp <- colMeans(outer(Z[, 1], zpts[, 1], ">") & outer(Z[, 2], zpts[, 2], ">"))
    rq <- std_radial_anchor(Z, thetas, qn)
    base <- build_base_curve(rq * cos(thetas), rq * sin(thetas))
    tail_probs <- percoord_level(zpts[, 1], zpts[, 2], base, qn, eta_hat, beta_cooley)
    ifelse(emp >= qn, emp, tail_probs)
}

# ---------- percoord-only, raw-wall tube constructor -------------------------
computeRegion_Percoord_Raw <- function(dat, alphas, p, B, gamma1 = 1/2,
                                       ncoords = 50, lbs = c(50, 0),
                                       beta_cooley = 200,
                                       marg_models = c("gpd", "gpd")) {
    n_dat <- nrow(dat)
    q1 <- n_dat^(-gamma1); qn <- q1
    prob_floor <- 0.5 / n_dat

    fit1 <- fit_marginal(dat[,1], gamma1, model = marg_models[1])
    fit2 <- fit_marginal(dat[,2], gamma1, model = marg_models[2])
    Z <- standardize(dat[,1], dat[,2], fit1, fit2)
    eta_hat <- estimate_xi_hill(pmin(Z[,1], Z[,2]), gamma1)

    ## ---- percoord tube, full space ----
    thetas_full <- make_theta_grid_ai(ncoords)
    rq1_full <- std_radial_anchor(Z, thetas_full, q1)
    base <- build_base_curve(rq1_full * cos(thetas_full), rq1_full * sin(thetas_full))
    fb1 <- exp(base$lf1); fb2 <- exp(base$lf2)
    rq_ax1 <- rq1_full[which.min(thetas_full)]
    rq_ax2 <- rq1_full[which.max(thetas_full)]
    pc_pts <- function(t) {
        P <- percoord_project(fb1, fb2, t, eta_hat, beta_cooley)
        rbind(c(0, t * rq_ax2), P, c(t * rq_ax1, 0))
    }
    ziso_pc  <- pc_pts(q1 / p)
    x_iso_pc <- backmap(ziso_pc, fit1, fit2)

    ## ---- bootstrap: percoord tube only, raw (variance-only) walls ----
    m_pc <- nrow(ziso_pc)
    boot_pc <- matrix(NA_real_, B, m_pc)
    n_fallback <- 0L; clamp_sum <- 0; zmax <- 2 * n_dat
    for (k in 1:B) {
        idx <- sample.int(n_dat, n_dat, replace = TRUE)
        d1 <- dat[idx, 1]; d2 <- dat[idx, 2]
        f1b <- fit_marginal(d1, gamma1, model = marg_models[1], fallback = fit1)
        f2b <- fit_marginal(d2, gamma1, model = marg_models[2], fallback = fit2)
        n_fallback <- n_fallback + f1b$used_fallback + f2b$used_fallback
        Zb <- standardize(d1, d2, f1b, f2b)
        clamp_sum <- clamp_sum + mean(Zb >= zmax - 1e-9)   # avg FRACTION at ceiling
        eta_b <- estimate_xi_hill(pmin(Zb[,1], Zb[,2]), gamma1)
        zx_pc <- standardize(x_iso_pc[,1], x_iso_pc[,2], f1b, f2b)
        boot_pc[k, ] <- log(pmax(std_blended_surv_pc(zx_pc, Zb, qn, thetas_full, eta_b, beta_cooley), prob_floor))
    }

    devs <- sweep(boot_pc, 2, colMeans(boot_pc), FUN = "-")
    Wp <- apply(devs, 1, max); Wm <- apply(devs, 1, min)
    rcps <- list(); rcms <- list()
    for (alpha in alphas) {
        ak <- as.character(alpha)
        dp <- as.numeric(quantile(Wp, probs = 1 - alpha/2))
        dm <- as.numeric(quantile(Wm, probs = alpha/2))
        rcps[[ak]] <- p * exp(dp) - p
        rcms[[ak]] <- p * exp(dm) - p
    }

    sf_pc <- function(pts) { pts <- as.matrix(pts)
        zq <- standardize(pts[,1], pts[,2], fit1, fit2)
        std_blended_surv_pc(zq, Z, qn, thetas_full, eta_hat, beta_cooley) }

    list(dat = dat, p = p, gamma1 = gamma1,
         marg_models = marg_models,
         fit1 = fit1[c("model","u","sigma","xi","q1")], fit2 = fit2[c("model","u","sigma","xi","q1")],
         eta_hat = eta_hat, beta_cooley = beta_cooley,
         fallback_frac = n_fallback / (2 * B), clamp_frac = clamp_sum / B,
         percoord = list(raw_c_plus_estimates = rcps, raw_c_minus_estimates = rcms,
                         survFunc = sf_pc, thetas = thetas_full, x_iso = x_iso_pc,
                         projection = "percoord"),
         transformed = FALSE)
}
### END CORE ###########################################################

# ======================================================================
# Main simulation loop
# ======================================================================
for (i in 1:length(dists)) {
    dist <- dists[i]
    sampling_func <- rKarachiBetaKDE
    if (!is.null(save_full_path)) {
        distpath <- paste0(save_full_path, dist, '/')
        if (!dir.exists(distpath)) dir.create(distpath, recursive = TRUE)
    }
    print(paste0('Starting new distribution: ', dist))

    for (k in 1:length(ns)) {
        n <- ns[k]
        p <- pn(n)

        # TRUE p-isoline (full, KDE truth) + interior subset via TRUE pseudo-angle
        isoline <- drawIsoline(dist = 'karachi', numCoords = n_coords,
                               gridUbs = ubs, gridLbs = lbs, prob = p)
        S1 <- true_marg_surv(isoline[,1], 1)
        S2 <- true_marg_surv(isoline[,2], 2)
        Zt <- cbind(1 / S1, 1 / S2)
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
            regions <- computeRegion_Percoord_Raw(dat = dat, alphas = alphas, p = p, B = B,
                                                  gamma1 = gamma1, ncoords = n_coords,
                                                  lbs = lbs, beta_cooley = beta_cooley,
                                                  marg_models = marg_models)

            # raw-wall coverage: same tube vs full truth and interior truth
            ev_raw <- function(iso) {
                obj <- list(survFunc = regions$percoord$survFunc, p = p,
                            c_plus_estimates  = regions$percoord$raw_c_plus_estimates,
                            c_minus_estimates = regions$percoord$raw_c_minus_estimates,
                            transformed = FALSE)
                evaluateCoverage(obj, iso)
            }
            cov_pc_full <- ev_raw(isoline)
            cov_pc_int  <- ev_raw(isoline_int)

            if (!is.null(n_dirpath)) {
                regions_save <- regions
                regions_save$percoord$survFunc <- NULL   # closures don't survive readRDS; dropping keeps files small
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
            rbind(mk_rows('percoord_full',     cov_pc_full, regions$percoord, nrow(isoline)),
                  mk_rows('percoord_interior', cov_pc_int,  regions$percoord, nrow(isoline_int)))
        }

        samp_df_list <- foreach(l = 1:n_sims, .options.snow = opts,
                .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt', 'matrixStats')) %dopar% parallelizedCode(l)

        close(pb)
        stopCluster(clust)

        if (!is.null(args$save_df_path)) {
            coverage_df <- do.call(rbind, samp_df_list)
            save_fname <- paste0(args$save_df_path, '1.1.1.0_bivkarachi_gamma1_0.5_bootmargs_cooleysmoothpercoord_n100000.csv')
            file_exists <- file.exists(save_fname)
            write.table(coverage_df, file = save_fname, sep = ",", row.names = FALSE,
                        col.names = !file_exists, append = file_exists)
        }
    }
}