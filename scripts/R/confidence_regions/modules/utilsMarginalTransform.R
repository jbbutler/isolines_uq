# utils_margtransform.R
#
# Shared helper functions for the marginal-transformation confidence-tube
# procedure. These were previously duplicated across every marginal-transform
# simulation script, its raw-eval companion, and the plotting module; they are
# collected here so all of those can `source()` a single definition.
#
# Grouped as:
#   (A) Marginal machinery      -- GPD fit, blended CDF + clamp, inverse transform
#   (B) Pareto standardization  -- forward/back maps built on the fitted margins
#   (C) Standardized-scale HRV  -- min-projection anchor, Cooley exponent,
#                                  blended standardized survival (AD xi=1 and AI HRV)
#   (D) Reconstruction helpers  -- rebuild a full fit / survFunc from a saved tube
#
# Depends on: matrixStats (colQuantiles).
#
# Jimmy Butler

library(matrixStats)

# ======================================================================
# (A) Marginal machinery: unconstrained GPD MLE + blended (empirical/GPD)
#     CDF with the [1/(2n), 1 - 1/(2n)] clamp.
# ======================================================================

# GPD negative log-likelihood, parametrized (log sigma, xi).
# Support violations are penalized (that IS the likelihood being -Inf; not an
# added constraint -- the fit stays Cooley-faithful/unconstrained).
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

# Plain MLE of (sigma, xi) via Nelder-Mead with moment-estimator starts.
# Returns moment estimates with failed = TRUE if optim errors or fails to
# converge; the caller decides the fallback.
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

# Fit one margin: empirical body + GPD tail above the (1 - q1) empirical quantile.
# On MLE failure, use `fallback` (the full-sample fit) if provided, else the
# moment estimates. used_fallback is reported for the caller to count.
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

# Forward transform: blended fitted CDF, clamped to [1/(2n), 1 - 1/(2n)].
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
    pmin(pmax(Fv, 1/(2*n)), 1 - 1/(2*n))
}

# Inverse transform: survival prob s -> original-scale quantile.
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

# ======================================================================
# (B) Pareto standardization built on the fitted margins.
# ======================================================================

# Data / point matrix -> unit-Pareto standardized coordinates.
standardize <- function(x1, x2, fit1, fit2) {
    cbind(1 / (1 - ptilde(x1, fit1)), 1 / (1 - ptilde(x2, fit2)))
}

# Standardized points -> original scale (T^{-1} then qblend per margin).
backmap <- function(zpts, fit1, fit2) {
    s1 <- pmin(1 / pmax(zpts[, 1], 1), 1)
    s2 <- pmin(1 / pmax(zpts[, 2], 1), 1)
    cbind(qblend(s1, fit1), qblend(s2, fit2))
}

# ======================================================================
# (C) Standardized-scale machinery (min-projection, Cooley exponent, survival).
# ======================================================================

# Radial anchor at level qlev along arbitrary angles (axis-aware; Z > 0 always).
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

# Cooley smoothed effective exponent along rays. On the standardized scale the
# marginal (axis) index is xi_marg = 1, so the weight is the plain
# max(z)/(z1+z2); beta_cooley is used natively (Frechet-scale calibration).
#   AD case: pass xi_int = 1  -> xi_eff == 1 everywhere.
#   AI case: pass xi_int = eta_hat.
cooley_xi_eff <- function(thetas, xi_int, xi_marg, beta_cooley) {
    a   <- 1 / xi_marg
    cs  <- cos(thetas)^a
    sn  <- sin(thetas)^a
    wmx <- pmax(cs, sn) / (cs + sn)
    m   <- 1 - wmx^beta_cooley
    m * xi_int + (1 - m) * xi_marg
}

# Blended standardized survival at query points: empirical above qn, radial
# projection below with the (possibly angle-dependent) exponent xi_eff.
# For AD pass eta_hat = 1, beta_cooley = 1 (xi_eff collapses to 1); for AI pass
# the estimated eta_hat and the beta used at construction.
std_blended_surv <- function(zpts, Z, qn, eta_hat = 1, beta_cooley = 200) {
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
# (D) Reconstruction from a saved tube object.
#     Saved tube objects trim fit1/fit2 to (u, sigma, xi, q1); x_sorted and n
#     are rebuilt from the stored data. The stored survFunc closure does not
#     survive readRDS into a fresh session (it loses its helper environment),
#     so rebuild_survFunc reconstructs it from scratch.
# ======================================================================

# Reconstruct a full marginal fit from a trimmed stored fit + the margin's data.
rebuild_fit <- function(fit_trim, x) {
    fit_trim$n <- length(x)
    fit_trim$x_sorted <- sort(x)
    fit_trim
}

# Reconstruct the composed clamped survival estimator from a saved tube object.
# Auto-detects AI vs AD from the presence of eta_hat in the tube object.
rebuild_survFunc <- function(tube_obj) {
    fit1 <- rebuild_fit(tube_obj$fit1, tube_obj$dat[, 1])
    fit2 <- rebuild_fit(tube_obj$fit2, tube_obj$dat[, 2])
    Z    <- standardize(tube_obj$dat[, 1], tube_obj$dat[, 2], fit1, fit2)
    qn   <- tube_obj$fit1$q1
    eta_hat     <- if (!is.null(tube_obj$eta_hat))     tube_obj$eta_hat     else 1
    beta_cooley <- if (!is.null(tube_obj$beta_cooley)) tube_obj$beta_cooley else 200
    percoord <- !is.null(tube_obj$projection) && tube_obj$projection == "percoord"
    thetas <- tube_obj$thetas
    if (percoord) {
        function(pts) {
            pts <- as.matrix(pts)
            zq <- standardize(pts[, 1], pts[, 2], fit1, fit2)
            std_blended_surv_pc(zq, Z, qn, thetas, eta_hat, beta_cooley)
        }
    } else {
        function(pts) {
            pts <- as.matrix(pts)
            zq <- standardize(pts[, 1], pts[, 2], fit1, fit2)
            std_blended_surv(zq, Z, qn, eta_hat, beta_cooley)
        }
    }
}

# ======================================================================
# (E) PER-COORDINATE (Cooley-exact) projection machinery for the AI branch.
#
# The radial single-exponent adaptation of Cooley's transition can FOLD in the
# transition zone: the exponent surge t^{1 - eta} near an axis multiplies BOTH
# coordinates, including the shrinking one, and when
#       (1 - eta) * log(q1/p) > e
# (crude criterion; the anchor slope and deeper wall levels lower it), the
# subdominant coordinate's blow-up outraces its geometric collapse and the
# level set doubles back -- an illegal (two-valued) isoline. Cooley's original
# per-coordinate scaling cannot fold: each coordinate's exponent varies in the
# same direction as the coordinate itself, so monotonicity of the projected
# curve is unconditional. This section implements that scheme on our anchor.
#
#   Base curve  : the q1 min-projection anchor (interior points only),
#                 pseudo-angle w = z1/(z1+z2), isotonized in log-coords.
#   Projection  : (z1, z2) -> (t^{e1(w)} z1, t^{e2(w)} z2), t = q1/s,
#                 e1(w) = eta + (1-eta) w^beta, e2(w) = eta + (1-eta)(1-w)^beta.
#   Inversion   : level s of a query point = unique t s.t. the point lies on
#                 the t-projected base curve (isolines are nested), found by
#                 VECTORIZED BISECTION in w on the dominance bracket where the
#                 objective is provably monotone.
#
# AD limit: eta = 1 gives e1 = e2 = 1 and the scheme reduces exactly to the
# radial projection; the AD code path is unchanged and validated against this.
# ======================================================================

# Cooley per-coordinate exponents at pseudo-angle w = z1/(z1+z2).
cooley_exponents <- function(w, eta_hat, beta_cooley) {
    e1 <- eta_hat + (1 - eta_hat) * w^beta_cooley
    e2 <- eta_hat + (1 - eta_hat) * (1 - w)^beta_cooley
    cbind(e1, e2)
}

# Interior base curve from anchor points on the theta grid: drops exact-axis
# points (a zero coordinate has no log), sorts by pseudo-angle, and isotonizes
# the log-coordinates (log z1 nondecreasing, log z2 nonincreasing in w) --
# the true q1-isoline is monotone; this removes empirical quantile noise that
# would break the bracketing of the inversion.
build_base_curve <- function(z1, z2) {
    keep <- is.finite(z1) & is.finite(z2) & z1 > 0 & z2 > 0
    z1 <- z1[keep]; z2 <- z2[keep]
    w <- z1 / (z1 + z2)
    o <- order(w)
    w <- w[o]; lf1 <- log(z1[o]); lf2 <- log(z2[o])
    lf1 <- cummax(lf1)                     # nondecreasing in w
    lf2 <- -cummax(-lf2)                   # nonincreasing in w (running min)
    # strictify for interpolation (ties break approx)
    eps <- 1e-10 * seq_along(w)
    w <- w + 1e-12 * seq_along(w)
    list(w = w, lf1 = lf1 + eps, lf2 = lf2 - eps,
         w_lo = w[1], w_hi = w[length(w)])
}

# Project points (z1, z2) with pseudo-angles from the points themselves to a
# deeper level via per-coordinate multipliers; t = (shallow level)/(deep level).
# Zero coordinates (exact-axis anchor points) pass through as zero.
percoord_project <- function(z1, z2, t, eta_hat, beta_cooley) {
    w <- ifelse(z1 + z2 > 0, z1 / (z1 + z2), 0.5)
    E <- cooley_exponents(w, eta_hat, beta_cooley)
    cbind(t^E[, 1] * z1, t^E[, 2] * z2)
}

# Level (survival prob) of query points under the per-coordinate isoline
# family anchored at level q1 on the given base curve. Vectorized bisection:
# on the dominance bracket [wA, wB] (where the query dominates the base curve
# componentwise) the objective
#     g(w) = e2(w) (log z1 - lf1(w)) - e1(w) (log z2 - lf2(w))
# is strictly decreasing with g(wA) >= 0 >= g(wB), so the root is unique.
# Queries outside the base curve's angular span get the dominant-coordinate
# marginal assignment at the nearest boundary (exact on the axes, where the
# margins are unit Pareto).
percoord_level <- function(zq1, zq2, base, q1, eta_hat, beta_cooley, n_iter = 40) {
    m <- length(zq1)
    # floor only for log-safety; do NOT clamp at 1 (near-axis projected points
    # legitimately have a subdominant coordinate below 1)
    lz1 <- log(pmax(zq1, 1e-12)); lz2 <- log(pmax(zq2, 1e-12))

    lf1_at <- function(w) approx(base$w, base$lf1, xout = w, rule = 2)$y
    lf2_at <- function(w) approx(base$w, base$lf2, xout = w, rule = 2)$y

    # dominance bracket: f1(w) <= z1  =>  w <= wB ;  f2(w) <= z2  =>  w >= wA
    wB <- approx(base$lf1, base$w, xout = lz1, rule = 2)$y      # lf1 increasing
    wA <- approx(rev(base$lf2), rev(base$w), xout = lz2, rule = 2)$y  # lf2 decreasing
    wA <- pmax(wA, base$w_lo); wB <- pmin(wB, base$w_hi)

    ok <- wA < wB                       # interior queries with a genuine bracket
    logt <- numeric(m)

    if (any(ok)) {
        lo <- wA[ok]; hi <- wB[ok]
        l1 <- lz1[ok]; l2 <- lz2[ok]
        for (it in 1:n_iter) {
            mid <- 0.5 * (lo + hi)
            E <- cooley_exponents(mid, eta_hat, beta_cooley)
            g <- E[, 2] * (l1 - lf1_at(mid)) - E[, 1] * (l2 - lf2_at(mid))
            move_right <- g > 0        # g decreasing: positive => root to the right
            lo <- ifelse(move_right, mid, lo)
            hi <- ifelse(move_right, hi, mid)
        }
        wstar <- 0.5 * (lo + hi)
        E <- cooley_exponents(wstar, eta_hat, beta_cooley)
        # average the two (equal at the root) expressions for log t
        lt1 <- (l1 - lf1_at(wstar)) / E[, 1]
        lt2 <- (l2 - lf2_at(wstar)) / E[, 2]
        logt[ok] <- 0.5 * (lt1 + lt2)
    }
    if (any(!ok)) {                     # beyond angular span: dominant marginal
        z1n <- zq1[!ok]; z2n <- zq2[!ok]
        l1 <- lz1[!ok]; l2 <- lz2[!ok]
        z1dom <- z1n / (z1n + z2n) >= 0.5
        E_hi <- cooley_exponents(base$w_hi, eta_hat, beta_cooley)
        E_lo <- cooley_exponents(base$w_lo, eta_hat, beta_cooley)
        logt[!ok] <- ifelse(z1dom,
                            (l1 - lf1_at(base$w_hi)) / E_hi[, 1],
                            (l2 - lf2_at(base$w_lo)) / E_lo[, 2])
    }
    logt <- pmax(logt, 0)               # never assign a level above q1
    q1 * exp(-logt)
}

# Blended standardized survival, PER-COORDINATE tail branch: empirical above
# qn; below, the level of the query under the per-coordinate isoline family
# anchored at the q1 base curve built from Z on the supplied theta grid.
std_blended_surv_pc <- function(zpts, Z, qn, thetas, eta_hat, beta_cooley) {
    zpts <- as.matrix(zpts)
    emp <- colMeans(outer(Z[, 1], zpts[, 1], ">") & outer(Z[, 2], zpts[, 2], ">"))
    rq <- std_radial_anchor(Z, thetas, qn)
    base <- build_base_curve(rq * cos(thetas), rq * sin(thetas))
    tail_probs <- percoord_level(zpts[, 1], zpts[, 2], base, qn, eta_hat, beta_cooley)
    ifelse(emp >= qn, emp, tail_probs)
}