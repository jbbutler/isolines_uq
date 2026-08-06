library(mvtnorm)
library(dplyr)
library(purrr)
library(data.table)
library(ismev)
library(evd)
library(matrixStats)

source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')

#### Extreme Tube Helpers ####

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
std_blended_surv_ad <- function(zpts, Z, qn) {
    zpts <- as.matrix(zpts)
    emp <- colMeans(outer(Z[,1], zpts[,1], ">") & outer(Z[,2], zpts[,2], ">"))
    theta <- atan2(zpts[,2], zpts[,1])
    rp <- sqrt(rowSums(zpts^2))
    rq <- std_radial_anchor(Z, theta, qn)
    tail_probs <- (rq / rp) * qn                 # xi = 1: (rq/rp)^{1/xi} = rq/rp
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
        
##########

# Fraga Alves, Gomes & de Haan (2003) estimator of the second-order parameter rho.
# Works on the LOG-EXCESSES V_{ik} = log(X_{n-i+1:n}) - log(X_{n-k:n}), i=1..k.
# 
# tau is the tuning parameter: tau = 0 is the standard choice for |rho| <= 1
# (which covers your bivariate-t case). For |rho| > 1, use tau = 1.
estimate_rho_FAGH <- function(radii, gamma_rho = 0.001, tau = 0) {
    n <- length(radii)
    k <- min(n-1, floor(2*n/(log(log(n)))))
    
    sorted_radii <- sort(radii)
    # log-excesses over threshold X_{n-k:n}: V_i = log X_{n-i+1:n} - log X_{n-k:n}
    log_top <- log(sorted_radii[(n - k + 1):n])     # k largest log-values
    log_threshold <- log(sorted_radii[n - k])       # threshold
    V <- log_top - log_threshold                     # log-excesses
    
    # Moments of the log-excesses
    M1 <- mean(V)
    M2 <- mean(V^2)
    M3 <- mean(V^3)
    
    # T statistic (tau = 0 vs tau > 0 cases)
    if (tau == 0) {
        T_num <- log(M1) - 0.5 * log(M2 / 2)
        T_den <- 0.5 * log(M2 / 2) - (1/3) * log(M3 / 6)
    } else {
        T_num <- (M1)^tau - (M2 / 2)^(tau / 2)
        T_den <- (M2 / 2)^(tau / 2) - (M3 / 6)^(tau / 3)
    }
    
    T_stat <- T_num / T_den
    
    # Invert to get rho. The closed-form is rho = -|3(T-1)/(T-3)|;
    # the abs and negation enforce rho <= 0 which is the valid range.
    rho_hat <- -abs(3 * (T_stat - 1) / (T_stat - 3))
    
    return(rho_hat)
}


# Gomes-Martins (2002) estimator of the "scale" second-order parameter beta,
# given rho_hat. Estimated at the same high level k1 = floor(n^(1 - gamma_rho)).
estimate_beta_GM <- function(radii, rho_hat, gamma_rho = 0.001) {
    n <- length(radii)
    k <- min(n-1, floor(2*n/(log(log(n)))))
    
    sorted_radii <- sort(radii)
    log_top <- log(sorted_radii[(n - k + 1):n])
    log_threshold <- log(sorted_radii[n - k])
    V <- log_top - log_threshold     # log-excesses
    
    i_seq <- 1:k
    w_rho  <- (i_seq / k)^(-rho_hat)
    w_2rho <- (i_seq / k)^(-2 * rho_hat)
    
    d_k       <- mean(w_rho)
    D_k_0     <- mean(V)
    D_k_rho   <- mean(w_rho * V)
    D_k_2rho  <- mean(w_2rho)
    
    numerator   <- (k / n)^rho_hat * (d_k * D_k_0 - D_k_rho)
    denominator <- d_k * D_k_rho - D_k_0 * D_k_2rho
    
    beta_hat <- numerator / denominator
    return(beta_hat)
}


# Caeiro-Gomes-Pestana (2005) Corrected Hill (CH) estimator.
# xi_hat_CH = xi_hat_Hill * (1 - (beta_hat / (1 - rho_hat)) * (n/k)^rho_hat)
# 
# By default, beta and rho are estimated externally at a higher level
# (gamma_rho = 0.001 -> k1 ~ n) and treated as fixed inputs. You can also
# pass in pre-computed rho_hat / beta_hat (e.g. fixed at -1) to avoid the
# noise from FAGH or to do sensitivity analysis.
estimate_xi_CH <- function(radii, gamma, gamma_rho = 0.001,
                           rho_hat = NULL, beta_hat = NULL, tau = 0) {
    n <- length(radii)
    k <- floor(n^(1 - gamma))
    
    # Plain Hill at level k
    xi_hill <- estimate_xi_hill(radii, gamma)
    
    # External estimation of (rho, beta) at level k1 if not supplied
    if (is.null(rho_hat)) {
        rho_hat <- estimate_rho_FAGH(radii, gamma_rho = gamma_rho, tau = tau)
    }
    if (is.null(beta_hat)) {
        beta_hat <- estimate_beta_GM(radii, rho_hat = rho_hat, 
                                     gamma_rho = gamma_rho)
    }
    
    # CH correction
    correction <- 1 - (beta_hat / (1 - rho_hat)) * (n / k)^rho_hat
    xi_ch <- xi_hill * correction

    res_lst <- list()
    res_lst$xi_ch <- xi_ch
    res_lst$xi_hill <- xi_hill
    
    return(res_lst)
}

estimate_xi_hill <- function(radii, gamma) {
    n <- length(radii)
    
    k <- floor(n^(1 - gamma))
    sorted_radii <- sort(radii)
    threshold_val <- sorted_radii[n - k]
    top_k_values <- sorted_radii[(n - k + 1):n]

    xi_hat <- mean(log(top_k_values) - log(threshold_val))
    
    return(xi_hat)
}

addCIs <- function(alpha, res) {
    #res is a dataframe that has the number of successful covers as one of its columns
    #Uses Clopper-Pearson to construct confidence intervals

    confint_lb <- qbeta(alpha/2, res$num_covered, 500 - res$num_covered + 1)
    confint_ub <- qbeta(1-(alpha/2), res$num_covered + 1, 500 - res$num_covered)

    res$confint_lb <- round(confint_lb,3)
    res$confint_ub <- round(confint_ub,3)

    return(res)
}

estimate_bias_bounds <- function(dat, gamma1, gamma2, p_target, xi, lbs=c(0,0)) {
    # Function to estimate conservative bias bounds b1 and b2 for a confidence tube 
    # using a 3-gamma hold-out projection framework.
    #
    # Arguments:
    # dat: the original data, in the form of a 2-column R data.frame
    # gamma1: the shallower threshold used to fit the baseline MRV model
    # gamma2: the deeper test threshold used to evaluate the projection error
    # p_target: your final extreme target probability (e.g., 5/n)
    # xi: extreme value index
    # lbs: lower bounds of the grid
    #
    # Output:
    # A list containing the absolute bounds (b1, b2), the relative bounds (r1, r2), 
    # and the test probability (p2).

    if (gamma1 <= gamma2) stop('gamma1 must be strictly greater than gamma2')

    n_dat <- nrow(dat)
    p2 <- n_dat^(-gamma2)

    # get the MRV-estimated isoline for the test probability p2, anchored at gamma1
    iso_p2 <- drawExtremeIsoline(dat, p = p2, n_coords = 200, grid_lbs = lbs, gamma = gamma1, xi = xi)

    # evaluate the purely empirical survival function along this generated isoline
    match_x <- outer(dat[,1], iso_p2[,1], '>')
    match_y <- outer(dat[,2], iso_p2[,2], '>')
    hit_matrix <- (match_x & match_y) * 1
    
    # the empirical survival probability for each point on the isoline
    emp_surv_p2 <- colMeans(hit_matrix)

    # calculate the relative bias bounds at the test level p2
    rel_errors <- (emp_surv_p2 - p2) / p2
    r1 <- min(rel_errors)
    r2 <- max(rel_errors)

    # scale the relative bounds to the final extreme target probability
    b1 <- r1 * p_target
    b2 <- r2 * p_target

    # rxeturn the bounds alongside the intermediate metrics for diagnostics
    res_lst <- list(
        b1 = b1, 
        b2 = b2, 
        r1 = r1, 
        r2 = r2, 
        p2 = p2
    )
    
    return(res_lst)
}

est_cdf <- function(x, dat, gamma, method = "ecdf") {
    # Function to estimate a univariate CDF, with more standard central statistics
    # methods for non-extreme probabilities, and extremal methods for extreme
    # probabilities. Defaults to ecdf as the standard central statistical method.
    
  if (method == "ecdf") { 
    edf <- ecdf(dat)
    vals <- edf(x)
  } else if (method == "kde") {
      
    bw <- bw.nrd0(dat)

    vals <- vapply(x, function(val) {
      mean(pnorm(val, mean = dat, sd = bw))
    }, numeric(1))
  } else {
    stop("Method must be either 'ecdf' or 'kde'")
  }
  # 2. Determine the threshold for the tail
  threshold_prob <- length(dat)^(-gamma)
  threshold <- quantile(dat, 1 - threshold_prob)
  # 3. Fit GPD to the tail
  is_tail <- x > threshold
  if (any(is_tail)) {
    x_gpd <- x[is_tail]
    # Fit GPD to the original data above threshold
    gpdOut <- ismev::gpd.fit(dat, threshold = threshold, show = FALSE)
    # Calculate GPD tail probabilities
    tail_probs <- 1-(1-evd::pgpd(x_gpd, 
                                     loc = gpdOut$threshold, 
                                     scale = gpdOut$mle[1], 
                                     shape = gpdOut$mle[2])) * threshold_prob
    vals[is_tail] <- tail_probs
  }
  return(vals)
}

# the inverse of the above function
est_inv_cdf <- function(p, dat, gamma) {

    threshold_prob <- length(dat)^(-gamma)
    threshold <- as.numeric(quantile(dat, 1-threshold_prob))
    quantiles <- as.numeric(quantile(dat, p))
    gpd_probs <- p[p > 1-threshold_prob]
    gpdOut <- ismev::gpd.fit(dat,
					threshold = threshold, show = F)
    quantiles[p > 1-threshold_prob] <- evd::qgpd(1-((1-gpd_probs)/threshold_prob), loc = gpdOut$threshold, scale = gpdOut$mle[1], shape = gpdOut$mle[2])
    quantiles <- as.numeric(quantiles)
    return(quantiles)

}

# function to estimate the survival function using the blended empirical + regularly varying method
# function assumes data lives in the nonnegative orthant
blendedSurvivalFunc <- function(pt, dat, gamma, xi) {

    qn <- nrow(dat)^(-gamma)
    empsurv_prob <- mean((dat[,1] > pt[1]) & (dat[,2] > pt[2]))

    if (empsurv_prob >= qn) {
        return(empsurv_prob)
    }
    else {
        pos_orthant_mask <- (dat[,1] > 0) & (dat[,2] > 0)
        theta <- atan2(pt[2], pt[1])
        # handle theta = 0, pi/2 cases separately
        if (theta == 0) {
            M <- dat[,1]*pos_orthant_mask
            rq <- quantile(M, probs=1-qn, names=FALSE, type=1)
        } else if (theta == pi/2) {
            M <- dat[,2]*pos_orthant_mask
            rq <- quantile(M, probs=1-qn, names=FALSE, type=1)
        } else {
            proj_x <- dat[,1]/cos(theta)
            proj_y <- dat[,2]/sin(theta)
            M <- pmin(proj_x, proj_y)
            rq <- quantile(M, probs=1-qn, type=1, names=FALSE)
        }
        rp <- sqrt(sum(pt^2))
        return(((rq/rp)^(1/xi))*qn)
    }
}

loadSamplingFunction <- function(dist) {
    if (dist == 'bivt') {
        return(function(n) data.frame(rmvt(n, sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2), df = 4)))
    } else if (dist == 'bivgauss') {
        return(function(n) data.frame(rmvnorm(n, mean = rep(0, 2), sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2))))
    } else if (dist == 'karachi') {
        return(function(n) rKarachiBetaKDE(n))
    }
}
# function to to compute the empirical survival function (given observed data) on a region of points which are not
# a regular grid
# works by converting all of the points to a regular grid, doing a fast empirical survival function operation
# and then taking only those points we were interested in in the first place
vectorizedEmpSurv <- function(region, dat) {
    
    sup_reg_xs_unique <- sort(unique(region$X1))
    sup_reg_ys_unique <- sort(unique(region$X2))
    full_grid <- expand.grid(X1=sup_reg_xs_unique, X2=sup_reg_ys_unique)
    surv <- fastEmpSurv(full_grid, dat)

    full_surv_res <- data.table(X1=full_grid$X1, X2=full_grid$X2, surv=surv)
    sup_reg <- data.table(region)

    res <- full_surv_res[sup_reg, on=c('X1', 'X2')]
    
    return(res$surv)
}

vectorizedBlendedSurv <- function(pts, dat, gamma, xi) {
    # Function to compute the blended survival function in vectorized fashion.
    # Assuming that the origin is (0,0).

    if (any(pts < 0)) {
        stop('Error: some points outside the nonnegative orthant.')
    }

    qn <- nrow(dat)^(-gamma)
    # compute empsurv value for all points
    match_x <- outer(dat[,1], pts[,1], ">") 
    match_y <- outer(dat[,2], pts[,2], ">")
    empsurv_probs <- colMeans(match_x & match_y)

    pos_orthant_mask <- (dat[,1] > 0) & (dat[,2] > 0)
    # compute projections of data onto radii of desired points
    theta <- atan2(pts[,2], pts[,1])
    # preallocate the list of radii
    rq <- rep(NA, length(theta))
    # if any 0s or pi/2, handle those separately
    if (any(theta == 0)) {
        ind_0 <- which(theta == 0)
        M_0 <- dat[,1]*pos_orthant_mask
        rq_0 <- quantile(M_0, probs=1-qn, names=FALSE, type=1)
        rq[ind_0] <- rq_0
    }
    if (any(theta == pi/2)) {
        ind_90 <- which(theta == pi/2)
        M_90 <- dat[,2]*pos_orthant_mask
        rq_90 <- quantile(M_90, probs=1-qn, names=FALSE, type=1)
        rq[ind_90] <- rq_90
    }
    
    inds_no_ax <- which(!(theta == 0 | theta == pi/2))
    # otherwise, handle all angles which are neither 0 nor pi/2
    if (length(inds_no_ax) > 0) {
        angles_no_ax <- theta[!(theta == 0 | theta == pi/2)] 
        rp <- sqrt(rowSums(pts^2))
        inv_cos <- 1/cos(angles_no_ax)
        inv_sin <- 1/sin(angles_no_ax)
        proj_x <- dat[,1] %o% inv_cos
        proj_y <- dat[,2] %o% inv_sin

        M_int <- pmin(proj_x, proj_y)
        rq[inds_no_ax] <- colQuantiles(M_int, probs=1-qn, type=1, drop=TRUE)
    }

    tail_probs <- ((rq/rp)^(1/xi))*qn
    exceedance_probs <- ifelse(empsurv_probs >= qn, empsurv_probs, tail_probs)

    return(exceedance_probs)
}

# function that selects points from an empirical isoline (from the empirical survival function)
# function that draws a grid of points over which we wish to evaluate the sup difference between
# the true survival function (in bootstrap world) and the empirical estimate (bootstrap survival function in bootstrap world)
drawEmpiricalIsoline <- function(dat, n_coords, grid_lbs, p) {

    # recenter data coordinate system for polar transformation
    dat_shifted <- sweep(dat, 2, grid_lbs, "-")

    # choose angles non uniformly spaced in theta to account for different axis scales
    x_seq <- seq(max(dat_shifted[,1]), 0, length.out=n_coords)
    y_seq <- seq(0, max(dat_shifted[,2]), length.out=n_coords)
    angles <- atan2(y_seq, x_seq)
    # handle angles between 0 and pi/2 exclusive
    angles_no_ax <- angles[2:(n_coords-1)]
    inv_cos <- 1/cos(angles_no_ax)
    inv_sin <- 1/sin(angles_no_ax)
    
    proj_x <- dat_shifted[,1] %o% inv_cos
    proj_y <- dat_shifted[,2] %o% inv_sin
    
    M_int <- pmin(proj_x, proj_y)

    # handle 0 and pi/2
    pos_orthant_mask <- (dat_shifted[,1] > 0) & (dat_shifted[,2] > 0)
    M_0 <- dat_shifted[,1]*pos_orthant_mask
    M_90 <- dat_shifted[,2]*pos_orthant_mask

    # stitch together into the same matrix
    M <- cbind(M_0, M_int, M_90)
    colnames(M) <- NULL
    radii <- colQuantiles(M, probs=1-p, type=1, drop=TRUE)

    iso_pts_x <- radii*cos(angles) + grid_lbs[1]
    iso_pts_y <- radii*sin(angles) + grid_lbs[2]

    return(data.frame(X1=iso_pts_x, X2=iso_pts_y))
}

# function to draw the isoline estimate using the blended empirical + regularly varying method
# Arguments
# dat: the data
# n_coords: the number of coordates you want
# gridLbs: the lower lefthand corner of the grid
# gridUbs: the upper bounds of the grid that will contain the full isoline
# gamma: parameter controlling how far into the tail you will be using empirical df
# xi: 1/the index of regular variation, EV index
drawExtremeIsoline <- function(dat, p, n_coords, grid_lbs, gamma, xi) {
    q <- nrow(dat)^(-gamma)
    q_isoline <- drawEmpiricalIsoline(dat, n_coords, grid_lbs, q)
    p_isoline <- q_isoline*((q/p)^(xi))
    
    return(p_isoline)
}

# Size of a Mammen-Polonik tube. Assumes that r_inner, r_outer, and r_ref have been computed
# and defined on the same angles, angle by angle. r_inner is the list of radii from the inner tube
# r_outer is the list of radii from the outer tube, and r_ref is the list of radii from the best estimate
computeTubeSize <- function(r_inner, r_outer, r_ref) {
 
  stopifnot(length(r_outer) == length(r_inner),
            length(r_ref)   == length(r_inner))
 
  w_abs <- r_outer - r_inner                # data-space distance
  w_rel <- w_abs / r_ref                    # dimensionless
 
  data.frame(
    width_abs_median = median(w_abs),
    width_abs_max    = max(w_abs),
    width_rel_median = median(w_rel),
    width_rel_max    = max(w_rel)
  )
}