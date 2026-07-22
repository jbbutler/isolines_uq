# Coverage for confidence tubes (Mammen-Polonik style), extreme regime.
# MARGINAL-TRANSFORMATION PROCEDURE, AI/HRV, bivariate Gaussian (0.7):
# DUAL-TUBE simulation with the PER-COORDINATE (Cooley-exact) transition.
#
# Per simulated dataset, TWO tubes are built from the SAME data, marginal fits,
# eta_hat, and bootstrap resamples:
#
#   (1) PERCOORD (full space): the Cooley transition implemented exactly as in
#       Cooley et al. (2019) -- per-coordinate exponents
#           e1(w) = eta + (1-eta) w^beta,  e2(w) = eta + (1-eta)(1-w)^beta,
#           w = z1/(z1+z2),
#       applied to the q1 min-projection anchor (isotonized base curve). This
#       REPLACES the radial per-angle xi_eff of the 1.1.1.2 runs, whose level
#       sets can double back (fold) in the transition zone.
#   (2) HRV-INTERIOR: pure radial HRV scaling, exponent eta everywhere, NO
#       transition, NO angular dependence -- defined and monitored only on the
#       interior cone, and evaluated only on the INTERIOR portion of the truth.
#
# Coverage rows written per sim (method column):
#   percoord_full     : tube (1) vs the full true isoline
#   percoord_interior : tube (1) vs the interior truth subset (same tube!)
#   hrv_interior      : tube (2) vs the interior truth subset
# The three-way comparison decomposes any undercoverage: full-vs-interior rows
# for tube (1) localize misses to the arms/transition; hrv_interior vs
# percoord_interior isolates the cost of the transition machinery itself.
#
# INTERIOR DEFINITION: true standardized pseudo-angle w_true = Z1/(Z1+Z2),
# Z_i = 1/(1 - pnorm(x_i)) (exact for the Gaussian truth), kept in
# [w_int, 1 - w_int] with w_int = 0.02: the Cooley weight there is
# 0.98^200 ~ 0.018 (interior >98% pure), and the transition midpoint
# (~0.0035) is far outside the retained band.
#
# Jimmy Butler, July 2026

set.seed(45678)

library(dplyr); library(mvtnorm); library(mnormt)
library(foreach); library(doSNOW); library(parallel)
library(argparse); library(matrixStats)

source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')            # estimate_xi_hill, computeTubeSize
source('~/isolines_uq/scripts/R/confidence_regions/modules/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/coverageEvaluation.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')

parser <- ArgumentParser(description = "Dual AI tubes: per-coordinate Cooley transition (full space) + pure-HRV radial (interior), bivariate Gaussian.")
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

# --- distribution / method parameters ---
rho_corr <- 0.7
beta_cooley <- 200
w_int <- 0.02                      # interior pseudo-angle cut (see header)

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
make_theta_grid_ai <- function(ncoords, n_axis = 8, ax_lo = 5e-4, ax_hi = 0.05) {
    ax <- exp(seq(log(ax_lo), log(ax_hi), length.out = n_axis))
    interior <- seq(0, pi/2, length.out = ncoords - 2 * n_axis)
    sort(unique(c(interior, ax, (pi/2) - ax)))
}

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

# ---------- per-coordinate (Cooley-exact) machinery --------------------------
cooley_exponents <- function(w, eta_hat, beta_cooley) {
    e1 <- eta_hat + (1 - eta_hat) * w^beta_cooley
    e2 <- eta_hat + (1 - eta_hat) * (1 - w)^beta_cooley
    cbind(e1, e2)
}
build_base_curve <- function(z1, z2) {
    # prepare the base curve for use in isoline estimation
    # essentially, enforce nice regular properties of isolines (no duplicate pseudo-angles, etc.)
    # make the inversions work too
    keep <- is.finite(z1) & is.finite(z2) & z1 > 0 & z2 > 0
    z1 <- z1[keep]; z2 <- z2[keep]
    w <- z1 / (z1 + z2)
    o <- order(w)
    w <- w[o]; lf1 <- log(z1[o]); lf2 <- log(z2[o])
    lf1 <- cummax(lf1)                       # nondecreasing in w
    lf2 <- -cummax(-lf2)                     # nonincreasing in w
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

# ---------- pure-HRV (constant exponent eta) survival -------------------------
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
computeRegion_Dual_HRV <- function(dat, alphas, p, B, gamma1 = 1/2, gamma2 = 2/3,
                                   ncoords = 50, lbs = c(0, 0), beta_cooley = 200,
                                   w_int = 0.02) {
    n_dat <- nrow(dat)
    q1 <- n_dat^(-gamma1); q2 <- n_dat^(-gamma2); qn <- q1
    prob_floor <- 0.5 / n_dat

    fit1 <- fit_marginal(dat[,1], gamma1)
    fit2 <- fit_marginal(dat[,2], gamma1)
    Z <- standardize(dat[,1], dat[,2], fit1, fit2)
    eta_hat <- estimate_xi_hill(pmin(Z[,1], Z[,2]), gamma1)

    ## ---- tube 1: PERCOORD, full space ----
    thetas_full <- make_theta_grid_ai(ncoords)
    rq1_full <- std_radial_anchor(Z, thetas_full, q1)
    base <- build_base_curve(rq1_full * cos(thetas_full), rq1_full * sin(thetas_full))
    fb1 <- exp(base$lf1); fb2 <- exp(base$lf2)
    rq_ax1 <- rq1_full[which.min(thetas_full)]     # theta = 0    (Z1 axis)
    rq_ax2 <- rq1_full[which.max(thetas_full)]     # theta = pi/2 (Z2 axis)
    pc_pts <- function(t) {
        P <- percoord_project(fb1, fb2, t, eta_hat, beta_cooley)
        rbind(c(0, t * rq_ax2), P, c(t * rq_ax1, 0))
    }
    ziso_pc  <- pc_pts(q1 / p)
    x_iso_pc <- backmap(ziso_pc, fit1, fit2)

    zq2_pc <- pc_pts(q1 / q2)
    emp_q2 <- colMeans(outer(Z[,1], zq2_pc[,1], ">") & outer(Z[,2], zq2_pc[,2], ">"))
    bft <- (emp_q2 - q2) / (q2 * log(q1 / q2)) * log(q1 / p)
    bf_pos_pc <- max(0, max(bft)); bf_neg_pc <- max(0, -min(bft))

    ## ---- tube 2: pure HRV, interior cone only ----
    theta_cut <- atan(w_int / (1 - w_int))
    thetas_int <- seq(theta_cut, pi/2 - theta_cut, length.out = ncoords)
    rq1_int <- std_radial_anchor(Z, thetas_int, q1)
    rp_int  <- rq1_int * (q1 / p)^eta_hat
    ziso_int  <- cbind(rp_int * cos(thetas_int), rp_int * sin(thetas_int))
    x_iso_int <- backmap(ziso_int, fit1, fit2)

    rq2_int <- rq1_int * (q1 / q2)^eta_hat
    zq2_int <- cbind(rq2_int * cos(thetas_int), rq2_int * sin(thetas_int))
    emp_q2i <- colMeans(outer(Z[,1], zq2_int[,1], ">") & outer(Z[,2], zq2_int[,2], ">"))
    bfti <- (emp_q2i - q2) / (q2 * log(q1 / q2)) * log(q1 / p)
    bf_pos_int <- max(0, max(bfti)); bf_neg_int <- max(0, -min(bfti))

    ## ---- shared marginal pad ----
    M0 <- log(q1 / p) / log(q1 / q2)
    B_marg1 <- mean(dat[,1] > qblend(q2, fit1)) / q2 - 1
    B_marg2 <- mean(dat[,2] > qblend(q2, fit2)) / q2 - 1
    b_marg <- (1 / eta_hat) * M0 * max(abs(B_marg1), abs(B_marg2))

    ## ---- shared bootstrap: one resample loop feeds both tubes ----
    m_pc <- nrow(ziso_pc); m_int <- nrow(ziso_int)
    boot_pc  <- matrix(NA_real_, B, m_pc)
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
        zx_pc  <- standardize(x_iso_pc[,1],  x_iso_pc[,2],  f1b, f2b)
        zx_int <- standardize(x_iso_int[,1], x_iso_int[,2], f1b, f2b)
        boot_pc[k, ]  <- log(pmax(std_blended_surv_pc(zx_pc, Zb, qn, thetas_full, eta_b, beta_cooley), prob_floor))
        boot_int[k, ] <- log(pmax(std_blended_surv_hrvconst(zx_int, Zb, qn, eta_b), prob_floor))
    }

    finish_tube <- function(bootM, bf_pos, bf_neg, wall_pts, survFunc) {
        devs <- sweep(bootM, 2, colMeans(bootM), FUN = "-")
        Wp <- apply(devs, 1, max); Wm <- apply(devs, 1, min)
        cps <- list(); cms <- list(); rcps <- list(); rcms <- list(); sizes <- list()
        rad <- function(iso) sqrt((pmax(iso[,1], lbs[1]) - lbs[1])^2 + (pmax(iso[,2], lbs[2]) - lbs[2])^2)
        for (alpha in alphas) {
            ak <- as.character(alpha)
            dp <- as.numeric(quantile(Wp, probs = 1 - alpha/2))
            dm <- as.numeric(quantile(Wm, probs = alpha/2))
            rcps[[ak]] <- p * exp(dp) - p
            rcms[[ak]] <- p * exp(dm) - p
            cp <- p * exp(dp + bf_neg + b_marg) - p
            cm <- p * exp(dm - bf_pos - b_marg) - p
            cps[[ak]] <- cp; cms[[ak]] <- cm
            s_lo <- p + cm
            if (s_lo <= 0) {
                sizes[[ak]] <- data.frame(is_unbounded = TRUE, width_abs_median = NA_real_,
                    width_abs_max = NA_real_, width_rel_median = NA_real_, width_rel_max = NA_real_)
            } else {
                s_hi <- min(p + cp, 0.999 * qn)
                wi <- wall_pts(s_hi); wo <- wall_pts(s_lo)
                sizes[[ak]] <- cbind(data.frame(is_unbounded = FALSE),
                    computeTubeSize(rad(wi$orig), rad(wo$orig), rad(wall_pts(p)$orig)))
            }
        }
        list(c_plus_estimates = cps, c_minus_estimates = cms,
             raw_c_plus_estimates = rcps, raw_c_minus_estimates = rcms,
             size_estimates = sizes, survFunc = survFunc,
             beta_frac_pos = bf_pos, beta_frac_neg = bf_neg)
    }

    wall_pc  <- function(s) list(orig = backmap(pc_pts(q1 / s), fit1, fit2))
    wall_int <- function(s) list(orig = backmap(cbind(rq1_int * (q1/s)^eta_hat * cos(thetas_int),
                                                      rq1_int * (q1/s)^eta_hat * sin(thetas_int)), fit1, fit2))
    sf_pc  <- function(pts) { pts <- as.matrix(pts)
        zq <- standardize(pts[,1], pts[,2], fit1, fit2)
        std_blended_surv_pc(zq, Z, qn, thetas_full, eta_hat, beta_cooley) }
    sf_int <- function(pts) { pts <- as.matrix(pts)
        zq <- standardize(pts[,1], pts[,2], fit1, fit2)
        std_blended_surv_hrvconst(zq, Z, qn, eta_hat) }

    t_pc  <- finish_tube(boot_pc,  bf_pos_pc,  bf_neg_pc,  wall_pc,  sf_pc)
    t_int <- finish_tube(boot_int, bf_pos_int, bf_neg_int, wall_int, sf_int)

    t_pc$thetas <- thetas_full; t_pc$x_iso <- x_iso_pc; t_pc$projection <- "percoord"
    t_int$thetas <- thetas_int; t_int$x_iso <- x_iso_int; t_int$projection <- "hrv_const"

    list(dat = dat, p = p, gamma1 = gamma1, gamma2 = gamma2,
         fit1 = fit1[c("u","sigma","xi","q1")], fit2 = fit2[c("u","sigma","xi","q1")],
         eta_hat = eta_hat, beta_cooley = beta_cooley, w_int = w_int,
         b_marg = b_marg, B_marg1 = B_marg1, B_marg2 = B_marg2,
         fallback_frac = n_fallback / (2 * B), clamp_frac = n_clamp / B,
         percoord = t_pc, interior = t_int, transformed = FALSE)
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
            regions <- computeRegion_Dual_HRV(dat = dat, alphas = alphas, p = p, B = B,
                                              gamma1 = gamma1, gamma2 = gamma2,
                                              ncoords = n_coords, lbs = lbs,
                                              beta_cooley = beta_cooley, w_int = w_int)

            # three coverage evaluations
            ev <- function(tube, iso) {
                obj <- list(survFunc = tube$survFunc, p = p,
                            c_plus_estimates = tube$c_plus_estimates,
                            c_minus_estimates = tube$c_minus_estimates,
                            transformed = FALSE)
                evaluateCoverage(obj, iso)
            }
            cov_pc_full <- ev(regions$percoord, isoline)
            cov_pc_int  <- ev(regions$percoord, isoline_int)
            cov_hrv_int <- ev(regions$interior, isoline_int)

            if (!is.null(n_dirpath)) saveRDS(regions, file = paste0(n_dirpath, 'simulation', ind, '.RData'))

            mk_rows <- function(method, covs, tube, n_truth) {
                ak <- names(covs)
                szs <- do.call(rbind, tube$size_estimates[ak])
                data.frame(dist = dist, p = p, n = n, method = method, alpha = ak,
                    is_covered = unlist(covs, use.names = FALSE),
                    c_plus_adj  = unlist(tube$c_plus_estimates[ak]),
                    c_minus_adj = unlist(tube$c_minus_estimates[ak]),
                    c_plus_raw  = unlist(tube$raw_c_plus_estimates[ak]),
                    c_minus_raw = unlist(tube$raw_c_minus_estimates[ak]),
                    beta_frac_pos = tube$beta_frac_pos, beta_frac_neg = tube$beta_frac_neg,
                    b_marg = regions$b_marg, B_marg1 = regions$B_marg1, B_marg2 = regions$B_marg2,
                    eta_hat = regions$eta_hat,
                    xi_gpd1 = regions$fit1$xi, xi_gpd2 = regions$fit2$xi,
                    fallback_frac = regions$fallback_frac, clamp_frac = regions$clamp_frac,
                    n_truth_pts = n_truth,
                    is_unbounded = szs$is_unbounded,
                    width_abs_median = szs$width_abs_median, width_abs_max = szs$width_abs_max,
                    width_rel_median = szs$width_rel_median, width_rel_max = szs$width_rel_max)
            }
            rbind(mk_rows('percoord_full',     cov_pc_full, regions$percoord, nrow(isoline)),
                  mk_rows('percoord_interior', cov_pc_int,  regions$percoord, nrow(isoline_int)),
                  mk_rows('hrv_interior',      cov_hrv_int, regions$interior, nrow(isoline_int)))
        }

        samp_df_list <- foreach(l = 1:n_sims, .options.snow = opts,
                .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt', 'matrixStats')) %dopar% parallelizedCode(l)

        close(pb)
        stopCluster(clust)

        if (!is.null(args$save_df_path)) {
            coverage_df <- do.call(rbind, samp_df_list)
            save_fname <- paste0(args$save_df_path, '1.1.1.2_bivgauss_gamma1_0.5_gamma2_0.67_bootmargs_cooleysmoothpercoord.csv')
            file_exists <- file.exists(save_fname)
            write.table(coverage_df, file = save_fname, sep = ",", row.names = FALSE,
                        col.names = !file_exists, append = file_exists)
        }
    }
}