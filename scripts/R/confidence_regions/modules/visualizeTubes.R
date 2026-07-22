# visualizeTubes.R
#
# Plotting routines that turn a saved tube object (.RData) into bounding-curve
# data.frames on the ORIGINAL scale, ready to plot.
#
# There are now three entry points, one per construction:
#   drawEmpiricalTubes      -- empirical-survival tubes (additive walls)
#   drawExtremeTubes        -- raw-scale extreme tubes with a stored xi_hat and a
#                              directly usable survival estimator (drawExtremeIsoline)
#   drawMargTransformTubes  -- marginal-transformation tubes, where the walls are
#                              level sets of the COMPOSED estimator and must be
#                              inverted on the standardized scale then back-mapped
#                              (handles both the AD, xi = 1, and AI/HRV branches)
#
# The old `transformed` option inside drawExtremeTubes has been REMOVED. That
# path conflated two different estimators; the marginal-transformation case now
# has its own function with the correct standardized-scale inversion.
#
# Shared marginal / standardized helpers live in utils_margtransform.R; only the
# plotting-specific glue is defined here.
#
# Jimmy Butler
# January 2026 (restructured)

source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')                 # drawExtremeIsoline, drawEmpiricalIsoline, computeTubeSize, ...
source('~/isolines_uq/scripts/R/confidence_regions/modules/utilsMarginalTransform.R')   # GPD/standardization/HRV helpers + rebuild_*

# ======================================================================
# Empirical-survival tubes (additive walls; unchanged behavior)
# ======================================================================
drawEmpiricalTubes <- function(tube_obj, alpha, lbs = c(0,0), n_coords = 200) {
    # Bounding curves for empirical-survival confidence tubes.
    #
    # Inputs:
    #   tube_obj : list defining the tube (output of computeEmpiricalRegion, etc.)
    #   alpha    : which alpha's tube to draw
    #   lbs      : lower-left corner / polar origin
    #   n_coords : number of points per curve
    # Output: list(tube_upper, tube_lower) of data.frames (or NA if a wall is
    #   unattainable).

    alpha <- as.character(alpha)
    dat <- tube_obj$dat
    p <- tube_obj$p
    alphas <- names(tube_obj$c_estimates)
    c_estimates <- tube_obj$c_estimates

    if (!(alpha %in% alphas)) {
        stop('No confidence tube for desired alpha.')
    }

    tube_lower <- drawEmpiricalIsoline(dat = dat, n_coords = n_coords,
                                       grid_lbs = lbs, p = c_estimates[[alpha]] + p)

    if (p - c_estimates[[alpha]] > 0) {
        tube_upper <- drawEmpiricalIsoline(dat = dat, n_coords = n_coords,
                                           grid_lbs = lbs, p = -c_estimates[[alpha]] + p)
    } else {
        tube_upper <- NA
    }

    list(tube_upper = tube_upper, tube_lower = tube_lower)
}

# ======================================================================
# Extreme-survival tubes: stored xi_hat + directly usable estimator
# (common-marginal AD/AI runs, oracle-xi runs, HRV raw-scale runs).
# NO marginal transform is involved here.
# ======================================================================
drawExtremeTubes <- function(tube_obj, alpha, lbs = c(0,0), n_coords = 200) {
    # Bounding curves for extreme-survival confidence tubes whose survival
    # estimator is the raw-scale blended estimator, parametrized by a stored
    # xi_hat (scalar) and gamma. Use this for every construction EXCEPT the
    # marginal-transformation procedure (use drawMargTransformTubes for that).
    #
    # Inputs:
    #   tube_obj : list defining the tube (output of computeExtremeRegion, etc.);
    #              must contain dat, p, gamma, xi_hat, c_plus_estimates,
    #              c_minus_estimates.
    #   alpha    : which alpha's tube to draw
    #   lbs      : lower-left corner / polar origin
    #   n_coords : number of points per curve
    # Output: list(tube_upper, tube_lower) of data.frames (NA if a wall is
    #   unattainable: inner wall survival >= 1, or outer wall survival <= 0).

    alpha <- as.character(alpha)
    dat   <- tube_obj$dat
    p     <- tube_obj$p
    gamma <- tube_obj$gamma
    # accept either xi_hat (current scripts) or xi (older objects)
    xi <- if (!is.null(tube_obj$xi_hat)) tube_obj$xi_hat else tube_obj$xi

    alphas <- names(tube_obj$c_plus_estimates)
    if (!(alpha %in% alphas)) {
        stop('No confidence tube for desired alpha.')
    }

    c_plus_estimate  <- tube_obj$c_plus_estimates[[alpha]]
    c_minus_estimate <- tube_obj$c_minus_estimates[[alpha]]

    # inner wall: higher survival (p + c_plus), smaller radius
    if (p + c_plus_estimate >= 1) {
        tube_lower <- NA
    } else {
        tube_lower <- drawExtremeIsoline(dat = dat, p = p + c_plus_estimate,
                                         n_coords = n_coords, grid_lbs = lbs,
                                         gamma = gamma, xi = xi)
    }

    # outer wall: lower survival (p + c_minus), larger radius
    if (p + c_minus_estimate <= 0) {
        tube_upper <- NA
    } else {
        tube_upper <- drawExtremeIsoline(dat = dat, p = p + c_minus_estimate,
                                         n_coords = n_coords, grid_lbs = lbs,
                                         gamma = gamma, xi = xi)
    }

    list(tube_upper = tube_upper, tube_lower = tube_lower)
}

# ======================================================================
# Marginal-transformation tubes: walls are level sets of the COMPOSED
# estimator. Invert on the standardized scale (xi = 1 under AD; angle-dependent
# xi_eff under AI/HRV), then back-map through the single fitted transform.
# Handles both dependence branches; AI is auto-detected via eta_hat.
# ======================================================================
drawTransformedTubes <- function(tube_obj, alpha, lbs = c(0,0), n_coords = 200,
                                   use_raw = FALSE, clamp_quadrant = TRUE) {
    # Bounding curves for marginal-transformation confidence tubes.
    #
    # Inputs:
    #   tube_obj : saved marginal-transform tube object; must contain dat, p,
    #              fit1, fit2 (trimmed), thetas, r_naive (standardized p-isoline
    #              radii), and c_(plus/minus)_estimates. AI objects additionally
    #              carry eta_hat, beta_cooley, and xi_eff (angle-dependent).
    #   alpha    : which alpha's tube to draw
    #   lbs      : polar origin for the returned original-scale radii/curves
    #   n_coords : angular resolution of the drawn curves. NOTE: if the object
    #              stores an angle-dependent xi_eff (AI branch), the drawn grid
    #              is taken from tube_obj$thetas (so xi_eff aligns) and n_coords
    #              is ignored; under AD a fresh uniform grid of size n_coords is
    #              used and xi_eff = 1.
    #   use_raw  : FALSE -> bias-corrected walls; TRUE -> raw variance-only walls
    #   clamp_quadrant : clamp returned coordinates at lbs so axis-angle back-map
    #              (empirical minima) does not leave the quadrant
    # Output: list(tube_upper, tube_lower, estimate) of data.frames on the
    #   original scale (tube_upper = outer/c_minus wall; tube_lower = inner/
    #   c_plus wall; estimate = the composed estimator's p-isoline). NA for a
    #   wall whose outer survival level is <= 0 (does not occur on the log scale).

    alpha <- as.character(alpha)
    thetas <- seq(0, pi/2, length.out = n_coords)
    
    cps <- if (use_raw) tube_obj$raw_c_plus_estimates  else tube_obj$c_plus_estimates
    cms <- if (use_raw) tube_obj$raw_c_minus_estimates else tube_obj$c_minus_estimates
    if (!(alpha %in% names(cps))) {
        stop(sprintf("No confidence tube for alpha '%s'; stored: %s",
                     alpha, paste(names(cps), collapse = ", ")))
    }

    p  <- tube_obj$p
    q1 <- tube_obj$fit1$q1

    # rebuild the single full-sample transform and the standardized data
    fit1 <- rebuild_fit(tube_obj$fit1, tube_obj$dat[, 1])
    fit2 <- rebuild_fit(tube_obj$fit2, tube_obj$dat[, 2])
    Z    <- standardize(tube_obj$dat[, 1], tube_obj$dat[, 2], fit1, fit2)

    # dispatch: percoord (Cooley-exact AI) / radial-AI (legacy) / AD.
    is_pc <- !is.null(tube_obj$projection) && tube_obj$projection == "percoord"
    is_ai <- !is.null(tube_obj$eta_hat)
    if (is_pc) {
        rq1  <- std_radial_anchor(Z, thetas, q1)
        base <- build_base_curve(rq1 * cos(thetas), rq1 * sin(thetas))
        fb1 <- exp(base$lf1); fb2 <- exp(base$lf2)
        rq_ax1 <- rq1[which.min(thetas)]; rq_ax2 <- rq1[which.max(thetas)]
        eta <- tube_obj$eta_hat; bc <- tube_obj$beta_cooley
        s_hi <- min(p + as.numeric(cps[[alpha]]), 0.999 * q1)
        s_lo <- p + as.numeric(cms[[alpha]])
        pc_curve <- function(t) {
            P <- percoord_project(fb1, fb2, t, eta, bc)
            M <- rbind(c(0, t * rq_ax2), P[rev(seq_len(nrow(P))), , drop = FALSE],
                       c(t * rq_ax1, 0))
            xy <- backmap(M, fit1, fit2)
            if (clamp_quadrant) xy <- cbind(pmax(xy[,1], lbs[1]), pmax(xy[,2], lbs[2]))
            data.frame(X1 = xy[,1], X2 = xy[,2])
        }
        return(list(tube_upper = if (s_lo <= 0) NA else pc_curve(q1 / s_lo),
                    tube_lower = pc_curve(q1 / s_hi),
                    estimate   = pc_curve(q1 / p)))
    }
    if (is_ai) {
        xi_eff <- tube_obj$xi_eff
        rq1    <- std_radial_anchor(Z, thetas, q1)
    } else {
        xi_eff <- rep(1, length(thetas))
        rq1    <- std_radial_anchor(Z, thetas, q1)
    }
    rp_std <- rq1 * (q1 / p)^xi_eff

    # wall survival levels (walls stored as additive deviations from p)
    s_hi <- p + as.numeric(cps[[alpha]])            # inner wall (higher survival)
    s_lo <- p + as.numeric(cms[[alpha]])            # outer wall (lower  survival)
    s_hi <- min(s_hi, 0.999 * q1)                   # keep inner wall in projection branch

    backmap_curve <- function(r_std) {
        xy <- backmap(cbind(r_std * cos(thetas), r_std * sin(thetas)), fit1, fit2)
        if (clamp_quadrant) xy <- cbind(pmax(xy[,1], lbs[1]), pmax(xy[,2], lbs[2]))
        data.frame(X1 = xy[,1], X2 = xy[,2])
    }

    tube_lower <- backmap_curve(rp_std * (p / s_hi)^xi_eff)   # inner / c_plus wall
    estimate   <- backmap_curve(rp_std)                        # p-isoline
    if (s_lo <= 0) {
        tube_upper <- NA
    } else {
        tube_upper <- backmap_curve(rp_std * (p / s_lo)^xi_eff)   # outer / c_minus wall
    }

    list(tube_upper = tube_upper, tube_lower = tube_lower, estimate = estimate)
}
# ======================================================================
# Dual-object glue (Karachi percoord + hrv_interior)
# ----------------------------------------------------------------------
# The restructured dual RDS stores dat/fits/eta_hat/marg_models ONCE at top
# level, with $percoord and $interior carrying only walls/thetas/projection.
# These helpers reshape that object for the drawers above and add the one
# construction the dispatcher does not cover (constant-exponent hrv_const).
# karachiTools.R (karachiDat, LBS_CONST, UBS_CONST, b_CONST) must be sourced
# by the caller for drawHrvInteriorTube's defaults and true_marg_surv.
# ======================================================================

# Attach the dual object's top-level fields onto a chosen sub-tube, producing a
# single-tube object ready for drawTransformedTubes (percoord/radial) or
# drawHrvInteriorTube (interior).
assemble_tube <- function(full, which = c("percoord", "interior")) {
    which <- match.arg(which)
    tube <- full[[which]]
    tube$dat <- full$dat; tube$p <- full$p
    tube$gamma1 <- full$gamma1; tube$gamma2 <- full$gamma2
    tube$fit1 <- full$fit1; tube$fit2 <- full$fit2
    tube$eta_hat <- full$eta_hat; tube$beta_cooley <- full$beta_cooley
    tube$transformed <- FALSE
    tube
}

# Constant-exponent (hrv_const) interior tube: walls at
# r(theta) = rq1(theta) * (q1/s)^eta on the interior cone, back-mapped. This is
# neither the percoord nor the Cooley-blend radial path the dispatcher handles,
# so it gets its own drawer. Uses sourced rebuild_fit/standardize/
# std_radial_anchor/backmap from utilsMarginalTransform.R.
drawHrvInteriorTube <- function(tube_obj, alpha, use_raw = FALSE,
                                lbs = LBS_CONST, clamp_quadrant = TRUE) {
    akey <- as.character(alpha)
    cps <- if (use_raw) tube_obj$raw_c_plus_estimates  else tube_obj$c_plus_estimates
    cms <- if (use_raw) tube_obj$raw_c_minus_estimates else tube_obj$c_minus_estimates
    p <- tube_obj$p; q1 <- tube_obj$fit1$q1
    fit1 <- rebuild_fit(tube_obj$fit1, tube_obj$dat[, 1])
    fit2 <- rebuild_fit(tube_obj$fit2, tube_obj$dat[, 2])
    Z <- standardize(tube_obj$dat[, 1], tube_obj$dat[, 2], fit1, fit2)
    rq1 <- std_radial_anchor(Z, tube_obj$thetas, q1)
    eta <- tube_obj$eta_hat
    curve_at <- function(s) {
        r <- rq1 * (q1 / s)^eta
        xy <- backmap(cbind(r * cos(tube_obj$thetas), r * sin(tube_obj$thetas)), fit1, fit2)
        if (clamp_quadrant) xy <- cbind(pmax(xy[, 1], lbs[1]), pmax(xy[, 2], lbs[2]))
        data.frame(X1 = xy[, 1], X2 = xy[, 2])
    }
    s_hi <- min(p + as.numeric(cps[[akey]]), 0.999 * q1)
    s_lo <- p + as.numeric(cms[[akey]])
    list(estimate   = curve_at(p),
         tube_lower = curve_at(s_hi),
         tube_upper = if (s_lo <= 0) NA else curve_at(s_lo))
}

# TRUE marginal survivals of the Karachi beta-KDE, for the interior-truth cut.
true_marg_surv <- function(x, j, dat = karachiDat, b = b_CONST,
                           lb = LBS_CONST[j], ub = UBS_CONST[j]) {
    y <- (dat[, j] - lb) / (ub - lb)
    s1 <- y / b + 1; s2 <- (1 - y) / b + 1
    x01 <- (x - lb) / (ub - lb)
    sapply(x01, function(u) mean(1 - pbeta(u, shape1 = s1, shape2 = s2)))
}