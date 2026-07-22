# computeMargTransformTube.R
#
# End-to-end construction of an extreme marginal-transformation confidence tube
# for a NEW bivariate dataset.
#
#   computeMargTransformTube(dat, p, alphas, dependence = "AD" | "AI", ...)
#
#   AD : standardized radial index exactly 1; radial projection (correct: a
#        single global homogeneity order); no index estimation anywhere.
#   AI : eta estimated by Hill on the standardized minima; PER-COORDINATE
#        (Cooley-exact) projection through the AD/AI transition:
#            (z1, z2) -> (t^{e1(w)} z1, t^{e2(w)} z2),   w = z1/(z1+z2),
#            e1 = eta + (1-eta) w^beta,  e2 = eta + (1-eta)(1-w)^beta,
#        which is UNCONDITIONALLY monotone (each coordinate's exponent moves
#        with that coordinate). The legacy radial single-exponent transition
#        (projection = "radial") is retained for reproducing archived results;
#        it can FOLD when (1 - eta) log(q1/p) approaches e, and warns.
#
# Returns a tube object shape-compatible with the simulation scripts' saveRDS
# output; AI objects carry projection/eta_hat/beta_cooley, which downstream
# code (drawMargTransformTubes, rebuild_survFunc) uses to dispatch.
#
# Depends on: utils_margtransform.R; estimate_xi_hill + computeTubeSize from
# the project utils.R; matrixStats; foreach/doSNOW/parallel if parallel = TRUE.
#
# Jimmy Butler

source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/utilsMarginalTransform.R')

library(matrixStats)

# hybrid theta grid for the AI branch: uniform interior + log-spaced points
# hugging each axis so the transition zone is monitored.
.make_theta_grid_ai <- function(ncoords, n_axis = 8, ax_lo = 5e-4, ax_hi = 0.05) {
    ax <- exp(seq(log(ax_lo), log(ax_hi), length.out = n_axis))
    interior <- seq(0, pi/2, length.out = ncoords - 2 * n_axis)
    sort(unique(c(interior, ax, (pi/2) - ax)))
}

computeMargTransformTube <- function(dat, p, alphas,
                                     dependence = c("AD", "AI"),
                                     projection = c("percoord", "radial"),
                                     gamma1 = 0.5, gamma2 = 0.67,
                                     B = 5000, ncoords = 50, lbs = c(0, 0),
                                     beta_cooley = 200,
                                     parallel = FALSE, n_cores = 4,
                                     seed = NULL, progress_bar = FALSE) {

    dependence <- match.arg(dependence)
    projection <- match.arg(projection)
    if (!is.null(seed)) set.seed(seed)
    dat <- as.matrix(dat)
    n_dat <- nrow(dat)
    q1 <- n_dat^(-gamma1); q2 <- n_dat^(-gamma2); qn <- q1
    prob_floor <- 0.5 / n_dat
    is_ai <- (dependence == "AI")
    use_pc <- is_ai && (projection == "percoord")   # AD is radial (identical to percoord at eta = 1)

    # --- Stage 1: full-sample marginal fits ---
    fit1 <- fit_marginal(dat[, 1], gamma1)
    fit2 <- fit_marginal(dat[, 2], gamma1)

    # --- Stage 2: standardize; (AI) estimate eta on the standardized minima ---
    Z <- standardize(dat[, 1], dat[, 2], fit1, fit2)
    eta_hat <- if (is_ai) estimate_xi_hill(pmin(Z[, 1], Z[, 2]), gamma1) else 1

    # fold diagnostic for the legacy radial transition
    if (is_ai && projection == "radial") {
        fold_stat <- (1 - eta_hat) * log(q1 / p)
        if (fold_stat > 0.5 * exp(1)) {
            warning(sprintf(paste0(
                "Radial AI transition: (1 - eta_hat) * log(q1/p) = %.2f ",
                "(fold threshold ~ e = 2.72; anchor slope and outer-wall depth lower it). ",
                "Level sets may fold in the transition zone; use projection = 'percoord'."),
                fold_stat))
        }
    }

    thetas <- seq(0, pi/2, length.out = ncoords)
    rq1_vec <- std_radial_anchor(Z, thetas, q1)
    bz1 <- rq1_vec * cos(thetas)                 # base (q1 anchor) points
    bz2 <- rq1_vec * sin(thetas)

    if (use_pc) {
        # isotonized interior base curve: the SAME object the survFunc inversion
        # uses, so projection and level-assignment are exactly consistent
        base <- build_base_curve(bz1, bz2)
        fb1 <- exp(base$lf1); fb2 <- exp(base$lf2)
        rq_ax1 <- rq1_vec[which.min(thetas)]      # theta = 0    (Z1 axis)
        rq_ax2 <- rq1_vec[which.max(thetas)]      # theta = pi/2 (Z2 axis)
        # curve at level s = q1/t: Z2-axis arm tip, interior base projected
        # per-coordinate, Z1-axis arm tip (ordered as a continuous path)
        pc_pts <- function(t) {
            P <- percoord_project(fb1, fb2, t, eta_hat, beta_cooley)
            rbind(c(0, t * rq_ax2), P[rev(seq_len(nrow(P))), , drop = FALSE], c(t * rq_ax1, 0))
        }
        ziso <- pc_pts(q1 / p)
    } else {
        xi_eff <- if (is_ai) cooley_xi_eff(thetas, eta_hat, 1, beta_cooley) else rep(1, length(thetas))
        rp_std <- rq1_vec * (q1 / p)^xi_eff
        ziso <- cbind(rp_std * cos(thetas), rp_std * sin(thetas))
    }
    x_iso <- backmap(ziso, fit1, fit2)           # FIXED original-scale points

    # --- Stage 5a: radial pad (copula / hidden-measure channel) at q2 ---
    if (use_pc) {
        zq2 <- pc_pts(q1 / q2)
    } else {
        rq2_vec <- rq1_vec * (q1 / q2)^xi_eff
        zq2 <- cbind(rq2_vec * cos(thetas), rq2_vec * sin(thetas))
    }
    emp_q2 <- colMeans(outer(Z[, 1], zq2[, 1], ">") & outer(Z[, 2], zq2[, 2], ">"))
    C_theta <- (emp_q2 - q2) / (q2 * log(q1 / q2))
    beta_frac_theta <- C_theta * log(q1 / p)
    beta_frac_pos <- max(0, max(beta_frac_theta))
    beta_frac_neg <- max(0, -min(beta_frac_theta))

    # --- Stage 5b: marginal pad; kappa = 1 (AD) or 1/eta (AI) ---
    M0 <- log(q1 / p) / log(q1 / q2)
    xq2_1 <- qblend(q2, fit1); xq2_2 <- qblend(q2, fit2)
    B_marg1 <- mean(dat[, 1] > xq2_1) / q2 - 1
    B_marg2 <- mean(dat[, 2] > xq2_2) / q2 - 1
    kappa <- if (is_ai) 1 / eta_hat else 1
    b_marg <- kappa * M0 * max(abs(B_marg1), abs(B_marg2))

    # --- Stage 4: bootstrap (margins refit per draw; eta re-estimated in AI) ---
    zmax <- 2 * n_dat
    one_draw <- function(k) {
        idx <- sample.int(n_dat, n_dat, replace = TRUE)
        d1 <- dat[idx, 1]; d2 <- dat[idx, 2]
        f1b <- fit_marginal(d1, gamma1, fallback = fit1)
        f2b <- fit_marginal(d2, gamma1, fallback = fit2)
        Zb <- standardize(d1, d2, f1b, f2b)
        eta_b <- if (is_ai) estimate_xi_hill(pmin(Zb[, 1], Zb[, 2]), gamma1) else 1
        zx <- standardize(x_iso[, 1], x_iso[, 2], f1b, f2b)
        probs_b <- if (use_pc) {
            std_blended_surv_pc(zx, Zb, qn, thetas, eta_b, beta_cooley)
        } else {
            std_blended_surv(zx, Zb, qn, eta_b, beta_cooley)
        }
        list(logp = log(pmax(probs_b, prob_floor)),
             fb = f1b$used_fallback + f2b$used_fallback,
             clamp = as.integer(any(Zb >= zmax - 1e-9)))
    }

    if (parallel) {
        requireNamespace("foreach"); requireNamespace("doSNOW"); requireNamespace("parallel")
        clust <- parallel::makeSOCKcluster(n_cores)
        doSNOW::registerDoSNOW(clust)
        parallel::clusterCall(clust, function() {
            suppressMessages(source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R'))
            suppressMessages(source('~/isolines_uq/scripts/R/confidence_regions/modules/utils_margtransform.R'))
            NULL
        })
        draws <- foreach::foreach(k = 1:B, .packages = c("matrixStats")) %dopar% one_draw(k)
        parallel::stopCluster(clust)
    } else {
        if (progress_bar) pbb <- txtProgressBar(min = 0, max = B, style = 3)
        draws <- vector("list", B)
        for (k in 1:B) {
            draws[[k]] <- one_draw(k)
            if (progress_bar) setTxtProgressBar(pbb, k)
        }
        if (progress_bar) close(pbb)
    }

    boot_probs_mat <- do.call(rbind, lapply(draws, `[[`, "logp"))
    n_fallback <- sum(vapply(draws, `[[`, numeric(1), "fb"))
    n_clamp    <- sum(vapply(draws, `[[`, integer(1), "clamp"))

    # --- centered log deviations, sup/inf, per-alpha quantiles ---
    mean_log_boot <- colMeans(boot_probs_mat)
    boot_log_devs <- sweep(boot_probs_mat, 2, mean_log_boot, FUN = "-")
    W_plus  <- apply(boot_log_devs, 1, max)
    W_minus <- apply(boot_log_devs, 1, min)

    c_plus_estimates <- list(); c_minus_estimates <- list()
    raw_c_plus_estimates <- list(); raw_c_minus_estimates <- list()
    size_estimates <- list()

    rad <- function(iso) sqrt((pmax(iso[, 1], lbs[1]) - lbs[1])^2 +
                              (pmax(iso[, 2], lbs[2]) - lbs[2])^2)

    for (alpha in alphas) {
        akey <- as.character(alpha)
        d_plus_log  <- as.numeric(quantile(W_plus,  probs = 1 - alpha/2))
        d_minus_log <- as.numeric(quantile(W_minus, probs = alpha/2))

        raw_c_plus_estimates[[akey]]  <- p * exp(d_plus_log)  - p
        raw_c_minus_estimates[[akey]] <- p * exp(d_minus_log) - p

        cplus  <- p * exp(d_plus_log  + beta_frac_neg + b_marg) - p
        cminus <- p * exp(d_minus_log - beta_frac_pos - b_marg) - p
        c_plus_estimates[[akey]]  <- cplus
        c_minus_estimates[[akey]] <- cminus

        s_lo <- p + cminus
        if (s_lo <= 0) {
            size_estimates[[akey]] <- data.frame(
                is_unbounded = TRUE, width_abs_median = NA_real_, width_abs_max = NA_real_,
                width_rel_median = NA_real_, width_rel_max = NA_real_)
        } else {
            s_hi <- min(p + cplus, 0.999 * qn)
            if (use_pc) {
                inner_z <- pc_pts(q1 / s_hi)
                outer_z <- pc_pts(q1 / s_lo)
            } else {
                inner_z <- cbind(rp_std * (p / s_hi)^xi_eff * cos(thetas),
                                 rp_std * (p / s_hi)^xi_eff * sin(thetas))
                outer_z <- cbind(rp_std * (p / s_lo)^xi_eff * cos(thetas),
                                 rp_std * (p / s_lo)^xi_eff * sin(thetas))
            }
            inner_iso <- backmap(inner_z, fit1, fit2)
            outer_iso <- backmap(outer_z, fit1, fit2)
            size_estimates[[akey]] <- cbind(
                data.frame(is_unbounded = FALSE),
                computeTubeSize(rad(inner_iso), rad(outer_iso), rad(x_iso)))
        }
    }

    # --- composed clamped survival estimator (the estimator of record) ---
    survFunc <- if (use_pc) {
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

    res_lst <- list(
        dat = dat, p = p, gamma1 = gamma1, gamma2 = gamma2,
        dependence = dependence,
        fit1 = fit1[c("u", "sigma", "xi", "q1")], fit2 = fit2[c("u", "sigma", "xi", "q1")],
        thetas = thetas, anchor_rq1 = rq1_vec, x_iso = x_iso,
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
    if (!use_pc) res_lst$r_naive <- if (is_ai) rp_std else rq1_vec * (q1 / p)  # legacy field
    if (is_ai) {
        res_lst$eta_hat     <- eta_hat
        res_lst$beta_cooley <- beta_cooley
        res_lst$projection  <- projection
        if (!use_pc) res_lst$xi_eff <- xi_eff
    }
    res_lst
}