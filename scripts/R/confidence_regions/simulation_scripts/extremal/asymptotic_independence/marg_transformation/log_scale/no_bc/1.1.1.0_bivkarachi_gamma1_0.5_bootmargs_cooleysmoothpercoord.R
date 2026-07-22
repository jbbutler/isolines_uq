# Coverage evaluation for the NON-BIAS-CORRECTED (raw, variance-only) walls of
# the DUAL-TUBE Karachi beta-KDE simulations (1.1.1.2 bivkarachi
# bootmargs_cooleysmoothpercoord): per-coordinate full-space tube + pure-HRV
# interior tube, with per-margin model dispatch (marg_models).
#
# Loads each saved 'regions' object, reconstructs BOTH survival estimators from
# the stored fits/eta (the survFunc closures do not survive readRDS), swaps each
# tube's RAW walls into the evaluation slots, and re-evaluates the three
# coverage rows (percoord_full / percoord_interior / hrv_interior) against the
# same true-pseudo-angle interior cut. Machinery inlined verbatim from the
# simulation script; the interior-truth cut uses the KDE's TRUE marginal
# survivals. marg_models is READ FROM the saved object so the reconstruction
# matches whatever the run used (gpd/gpd pathology or gpd/empirical fix).
#
# Jimmy Butler, July 2026

library(dplyr); library(mvtnorm); library(matrixStats)
source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/coverageEvaluation.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')      # karachiDat, *_CONST, drawIsoline('karachi')

# NOTE: point save_full_path at wherever the 1.1.1.2 bivkarachi dual tubes were saved.
save_full_path <- '/pscratch/sd/j/jbbutler/isolines/simulations/extremal/1.1.1.2_bivkarachi_gamma1_0.5_gamma2_0.67_bootmargs_cooleysmoothpercoord/'
save_df_path   <- '~/isolines_uq/outputs/simulations/extremal/'

dists <- c('karachiBetaKDE')
ns <- c(1000, 5000, 10000, 50000, 100000)
C <- 5; pn <- function(n){ C/n }
n_coords <- 50
lbs <- LBS_CONST; ubs <- UBS_CONST
w_int <- 0.02
beta_cooley <- 200

# TRUE marginal survivals of the KDE (kernel mixture), for the interior cut
true_marg_surv <- function(x, j, dat = karachiDat, b = b_CONST,
                           lb = LBS_CONST[j], ub = UBS_CONST[j]) {
    y <- (dat[, j] - lb) / (ub - lb)
    s1 <- y / b + 1; s2 <- (1 - y) / b + 1
    x01 <- (x - lb) / (ub - lb)
    sapply(x01, function(u) mean(1 - pbeta(u, shape1 = s1, shape2 = s2)))
}

# ======================================================================
# Machinery inlined verbatim from the simulation script's CORE block.
# ======================================================================


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
std_blended_surv_hrvconst <- function(zpts, Z, qn, eta_hat) {
    zpts <- as.matrix(zpts)
    emp <- colMeans(outer(Z[, 1], zpts[, 1], ">") & outer(Z[, 2], zpts[, 2], ">"))
    theta <- atan2(zpts[, 2], zpts[, 1])
    rp <- sqrt(rowSums(zpts^2))
    rq <- std_radial_anchor(Z, theta, qn)
    tail_probs <- ((rq / rp)^(1 / eta_hat)) * qn
    ifelse(emp >= qn, emp, tail_probs)
}

# ---------- dual tube constructor ---------------------------------------------

# reconstruct a full marginal fit from a trimmed stored fit + data, honoring model
rebuild_fit <- function(fit_trim, x, model) {
    n <- length(x); xs <- sort(x)
    fit_trim$n <- n; fit_trim$x_sorted <- xs; fit_trim$model <- model
    if (identical(model, "empirical")) {
        pos_all <- (seq_len(n) - 0.5) / n
        ux <- unique(xs)
        fit_trim$ux <- ux
        fit_trim$pos <- as.numeric(tapply(pos_all, match(xs, ux), mean))
    }
    fit_trim
}

# ======================================================================
# Evaluation loop
# ======================================================================
results_list <- list()
print("Evaluating RAW dual tubes: karachi percoord + hrv interior...")

for (dist in dists) {
  for (n in ns) {
    p <- pn(n)
    print(paste0("Processing ", dist, " | n = ", n))

    # TRUE p-isoline (full, KDE) + interior subset via TRUE pseudo-angle
    isoline <- drawIsoline(dist = 'karachi', numCoords = n_coords, gridUbs = ubs, gridLbs = lbs, prob = p)
    S1 <- true_marg_surv(isoline[,1], 1); S2 <- true_marg_surv(isoline[,2], 2)
    Zt <- cbind(1 / S1, 1 / S2)
    w_true <- Zt[,1] / (Zt[,1] + Zt[,2])
    keep_int <- w_true >= w_int & w_true <= 1 - w_int
    isoline_int <- isoline[keep_int, , drop = FALSE]

    sim_dir <- paste0(save_full_path, dist, '/n', n, '/')
    sim_files <- list.files(sim_dir, pattern = "\\.RData$", full.names = TRUE)
    if (length(sim_files) == 0) { warning(paste0("No .RData in ", sim_dir)); next }
    pb <- txtProgressBar(min = 0, max = length(sim_files), style = 3)

    for (i in seq_along(sim_files)) {
      regions <- readRDS(sim_files[i])
      dat <- regions$dat
      mm <- if (!is.null(regions$marg_models)) regions$marg_models else c("gpd","gpd")
      fit1 <- rebuild_fit(regions$fit1, dat[,1], mm[1])
      fit2 <- rebuild_fit(regions$fit2, dat[,2], mm[2])
      qn <- regions$fit1$q1
      Z  <- standardize(dat[,1], dat[,2], fit1, fit2)
      eta_hat <- regions$eta_hat

      sf_pc <- function(pts) { pts <- as.matrix(pts)
          zq <- standardize(pts[,1], pts[,2], fit1, fit2)
          std_blended_surv_pc(zq, Z, qn, regions$percoord$thetas, eta_hat, beta_cooley) }
      sf_int <- function(pts) { pts <- as.matrix(pts)
          zq <- standardize(pts[,1], pts[,2], fit1, fit2)
          std_blended_surv_hrvconst(zq, Z, qn, eta_hat) }

      ev_raw <- function(tube, survFunc, iso) {
          obj <- list(survFunc = survFunc, p = p,
                      c_plus_estimates  = tube$raw_c_plus_estimates,
                      c_minus_estimates = tube$raw_c_minus_estimates,
                      transformed = FALSE)
          evaluateCoverage(obj, iso)
      }
      cov_pc_full <- ev_raw(regions$percoord, sf_pc,  isoline)
      cov_pc_int  <- ev_raw(regions$percoord, sf_pc,  isoline_int)
      cov_hrv_int <- ev_raw(regions$interior, sf_int, isoline_int)

      mk_rows <- function(method, covs, tube, n_truth) {
          ak <- names(covs)
          data.frame(dist = dist, p = p, n = n, method = method, alpha = ak,
              is_covered = unlist(covs, use.names = FALSE),
              c_plus_raw  = unlist(tube$raw_c_plus_estimates[ak]),
              c_minus_raw = unlist(tube$raw_c_minus_estimates[ak]),
              b_marg = regions$b_marg, B_marg1 = regions$B_marg1, B_marg2 = regions$B_marg2,
              eta_hat = regions$eta_hat, xi_gpd1 = regions$fit1$xi, xi_gpd2 = regions$fit2$xi,
              fallback_frac = regions$fallback_frac, clamp_frac = regions$clamp_frac,
              n_truth_pts = n_truth)
      }
      res <- rbind(
          mk_rows('percoord_full',     cov_pc_full, regions$percoord, nrow(isoline)),
          mk_rows('percoord_interior', cov_pc_int,  regions$percoord, nrow(isoline_int)),
          mk_rows('hrv_interior',      cov_hrv_int, regions$interior, nrow(isoline_int)))
      results_list[[length(results_list) + 1]] <- res
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
}

if (length(results_list) > 0) {
    coverage_df <- do.call(rbind, results_list)
    save_fname <- paste0(save_df_path, '1.1.1.0_bivkarachi_gamma1_0.5_bootmargs_cooleysmoothpercoord.csv')
    write.table(coverage_df, file = save_fname, sep = ",", row.names = FALSE, col.names = TRUE)
    print(paste0("Saved: ", save_fname))
} else print("No results; check .RData directories.")