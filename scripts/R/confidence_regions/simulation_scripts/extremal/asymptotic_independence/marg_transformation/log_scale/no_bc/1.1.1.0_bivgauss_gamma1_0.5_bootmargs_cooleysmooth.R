# Script to loop through saved MARGINAL-TRANSFORM (AI/HRV, Cooley-smoothed)
# bivariate GAUSSIAN confidence tubes and evaluate un-bias-corrected coverage.
#
# Loads the saved tube objects from the marginal-transformation AI simulations
# (1.1.1.2 bivgauss bootmargs cooleysmooth), swaps the RAW (variance-only) walls
# into the evaluation slots, and re-evaluates coverage of the UNCORRECTED tubes
# against the true bivariate Gaussian isoline. Also records raw tube sizes.
#
# THE AI RAW COLUMN: the hidden channel here carries the Gaussian copula's slowly
# varying factor (effectively rho = 0 for the hidden-measure convergence), so the
# prediction is that RAW coverage degrades with n while the corrected column
# (already in the sim CSV) holds -- the radial pad load-bearing. Read eta_hat
# (truth 0.85, expected biased and slow to improve) and B_marg1/2 (expected small:
# Gaussian margins are GPD-friendly in amplitude) alongside to attribute the
# degradation to the hidden channel rather than the margins.
#
# IMPORTANT: the tube's survFunc is a closure over inline helpers defined in the
# SIMULATION script's global environment; after readRDS in a fresh session those
# are not found. This script REDEFINES them (including the HRV/Cooley pieces) and
# rebuilds survFunc explicitly from the stored dat + trimmed fits + eta_hat.
library(dplyr)
library(matrixStats)
source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')            # computeTubeSize lives here
source('~/isolines_uq/scripts/R/confidence_regions/modules/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/coverageEvaluation.R')

# NOTE: point save_full_path at wherever the 1.1.1.2 bivgauss bootmargs tubes were saved.
save_full_path <- '/pscratch/sd/j/jbbutler/isolines/simulations/extremal/1.1.1.2_bivgauss_gamma1_0.5_gamma2_0.67_bootmargs_cooleysmooth/'
save_df_path   <- '~/isolines_uq/outputs/simulations/extremal/'

dists <- c('bivgauss')
ns <- c(1000, 5000, 10000, 50000, 100000)
C <- 5
pn <- function(n){ C/n }
n_coords <- 50
lbs <- c(0, 0)
ubs <- c(200, 200)

# ======================================================================
# Marginal + standardized machinery (copied from the simulation script so the
# reconstructed survFunc and the size back-mapping resolve their dependencies).
# ======================================================================

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
    pmin(pmax(Fv, 1/(2*n)), 1 - 1/(2*n))
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

standardize <- function(x1, x2, fit1, fit2) {
    cbind(1 / (1 - ptilde(x1, fit1)), 1 / (1 - ptilde(x2, fit2)))
}

backmap <- function(zpts, fit1, fit2) {
    s1 <- pmin(1 / pmax(zpts[, 1], 1), 1)
    s2 <- pmin(1 / pmax(zpts[, 2], 1), 1)
    cbind(qblend(s1, fit1), qblend(s2, fit2))
}

# Cooley smoothed effective exponent: xi_marg = 1 on the standardized scale
cooley_xi_eff <- function(thetas, xi_int, xi_marg, beta_cooley) {
    a   <- 1 / xi_marg
    cs  <- cos(thetas)^a
    sn  <- sin(thetas)^a
    wmx <- pmax(cs, sn) / (cs + sn)
    m   <- 1 - wmx^beta_cooley
    m * xi_int + (1 - m) * xi_marg
}

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

# blended standardized survival: empirical above qn, Cooley-smoothed HRV
# projection (exponent xi_eff(theta)) below
std_blended_surv_hrv <- function(zpts, Z, qn, eta_hat, beta_cooley) {
    zpts <- as.matrix(zpts)
    emp <- colMeans(outer(Z[,1], zpts[,1], ">") & outer(Z[,2], zpts[,2], ">"))
    theta <- atan2(zpts[,2], zpts[,1])
    rp <- sqrt(rowSums(zpts^2))
    rq <- std_radial_anchor(Z, theta, qn)
    xi_eff <- cooley_xi_eff(theta, eta_hat, 1, beta_cooley)
    tail_probs <- ((rq / rp)^(1 / xi_eff)) * qn
    ifelse(emp >= qn, emp, tail_probs)
}

# reconstruct the full marginal fit (add x_sorted, n) from the trimmed stored fit + data
rebuild_fit <- function(fit_trim, x) {
    fit_trim$n <- length(x)
    fit_trim$x_sorted <- sort(x)
    fit_trim
}

# original-scale radius from lbs, coordinates clamped into the quadrant
rad_from_lbs <- function(iso) {
    sqrt((pmax(iso[,1], lbs[1]) - lbs[1])^2 + (pmax(iso[,2], lbs[2]) - lbs[2])^2)
}

# Initialize a list to hold the results of all simulations
results_list <- list()

print("Starting evaluation of raw (uncorrected) marginal-transform bivgauss HRV tubes...")

for (dist in dists) {
  for (n in ns) {
    p <- pn(n)
    print(paste0("\nProcessing dist: ", dist, " | n: ", n))

    # TRUE p-isoline: the distribution IS the bivariate Gaussian -- drawn directly
    # (truth on the ORIGINAL scale; survFunc composes the transform internally)
    true_isoline <- drawIsoline(dist = dist, numCoords = n_coords, gridUbs = ubs, gridLbs = lbs, prob = p)

    sim_dir <- paste0(save_full_path, dist, '/n', n, '/')
    sim_files <- list.files(sim_dir, pattern = "\\.RData$", full.names = TRUE)
    total_files <- length(sim_files)

    if (total_files == 0) {
        warning(paste0("No .RData files found in: ", sim_dir, ". Skipping..."))
        next
    }

    pb <- txtProgressBar(min = 0, max = total_files, style = 3)

    for (i in 1:total_files) {
      file <- sim_files[i]
      tube_obj <- readRDS(file)

      # --- reconstruct full marginal fits + fresh HRV survFunc (see header note) ---
      fit1 <- rebuild_fit(tube_obj$fit1, tube_obj$dat[, 1])
      fit2 <- rebuild_fit(tube_obj$fit2, tube_obj$dat[, 2])
      qn   <- tube_obj$fit1$q1                       # = n^(-gamma1)
      Z    <- standardize(tube_obj$dat[, 1], tube_obj$dat[, 2], fit1, fit2)
      eta_hat     <- tube_obj$eta_hat
      beta_cooley <- tube_obj$beta_cooley

      tube_obj$survFunc <- function(pts) {
          pts <- as.matrix(pts)
          zq <- standardize(pts[, 1], pts[, 2], fit1, fit2)
          std_blended_surv_hrv(zq, Z, qn, eta_hat, beta_cooley)
      }

      # swap the RAW bounds into the standard evaluation slots
      tube_obj$c_plus_estimates  <- tube_obj$raw_c_plus_estimates
      tube_obj$c_minus_estimates <- tube_obj$raw_c_minus_estimates

      # evaluate coverage using the raw bounds (transformed = FALSE: survFunc composes internally)
      is_covereds <- evaluateCoverage(tube_obj, true_isoline)
      alphas <- names(is_covereds)

      # --- tube size (RAW walls): angle-dependent xi_eff inversion -> backmap -> radii ---
      rp_std <- tube_obj$r_naive     # standardized p-isoline radius, per angle (hybrid grid)
      thetas <- tube_obj$thetas
      xi_eff <- tube_obj$xi_eff      # angle-dependent smoothed exponent on that grid
      r_ref  <- rad_from_lbs(tube_obj$x_iso)

      size_rows <- lapply(alphas, function(a) {
        cplus  <- as.numeric(tube_obj$raw_c_plus_estimates[[a]])
        cminus <- as.numeric(tube_obj$raw_c_minus_estimates[[a]])
        s_lo <- p + cminus
        if (s_lo <= 0) {               # inert on the log scale; kept for parity
          data.frame(is_unbounded = TRUE,
                     width_abs_median = NA_real_, width_abs_max = NA_real_,
                     width_rel_median = NA_real_, width_rel_max = NA_real_)
        } else {
          s_hi <- min(p + cplus, 0.999 * qn)   # keep inner wall in the projection branch
          r_in_std  <- rp_std * (p / s_hi)^xi_eff
          r_out_std <- rp_std * (p / s_lo)^xi_eff
          inner_iso <- backmap(cbind(r_in_std  * cos(thetas), r_in_std  * sin(thetas)), fit1, fit2)
          outer_iso <- backmap(cbind(r_out_std * cos(thetas), r_out_std * sin(thetas)), fit1, fit2)
          r_inner <- rad_from_lbs(inner_iso)
          r_outer <- rad_from_lbs(outer_iso)
          cbind(data.frame(is_unbounded = FALSE),
                computeTubeSize(r_inner, r_outer, r_ref))
        }
      })
      size_df <- do.call(rbind, size_rows)

      res_df <- data.frame(
        dist = dist,
        p = p,
        n = n,
        alpha = alphas,
        is_covered = unlist(is_covereds, use.names=FALSE),
        c_plus_raw  = unlist(tube_obj$raw_c_plus_estimates[alphas]),
        c_minus_raw = unlist(tube_obj$raw_c_minus_estimates[alphas]),
        b_marg  = tube_obj$b_marg,
        B_marg1 = tube_obj$B_marg1,
        B_marg2 = tube_obj$B_marg2,
        eta_hat = tube_obj$eta_hat,
        xi_gpd1 = tube_obj$fit1$xi,
        xi_gpd2 = tube_obj$fit2$xi,
        fallback_frac = tube_obj$fallback_frac,
        clamp_frac    = tube_obj$clamp_frac,
        is_unbounded     = size_df$is_unbounded,
        width_abs_median = size_df$width_abs_median,
        width_abs_max    = size_df$width_abs_max,
        width_rel_median = size_df$width_rel_median,
        width_rel_max    = size_df$width_rel_max
      )

      results_list[[length(results_list) + 1]] <- res_df
      setTxtProgressBar(pb, i)
    }

    close(pb)
  }
}

# combine all results and save the cleaned CSV
if (length(results_list) > 0) {
    coverage_df <- do.call(rbind, results_list)

    save_fname <- paste0(save_df_path, '1.1.1.0_bivgauss_gamma1_0.5_bootmargs_cooleysmooth.csv')

    write.table(coverage_df,
                file = save_fname,
                sep = ",",
                row.names = FALSE,
                col.names = TRUE)

    print(paste0("\nSuccessfully saved results to: ", save_fname))
} else {
    print("\nNo results to save. Check if the directories contain .RData files.")
}