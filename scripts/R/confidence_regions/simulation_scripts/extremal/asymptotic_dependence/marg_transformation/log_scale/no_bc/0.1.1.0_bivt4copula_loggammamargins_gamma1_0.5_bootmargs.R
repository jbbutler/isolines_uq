# Script to loop through saved MARGINAL-TRANSFORM (AD) t4-copula/LOG-GAMMA-margin
# confidence tubes and evaluate un-bias-corrected coverage.
#
# Loads the saved tube objects from the marginal-transformation AD simulations with
# log-gamma margins (0.1.1.2 bivt4copula_loggammamargins bootmargs), swaps the RAW
# (variance-only) walls into the evaluation slots, and re-evaluates coverage of the
# UNCORRECTED tubes against the true (original-scale) log-gamma isoline. Also
# records the size of the raw tubes.
#
# THE DECISIVE-TEST RAW COLUMN: the marginal channel here is rho = 0 (GPD misfit to
# the log-gamma tail decays only logarithmically), so the prediction is that RAW
# coverage degrades with n -- slowly, log-rate -- in contrast to the flat-at-nominal
# raw coverage of the bivt4 run, while the corrected column (already in the sim CSV)
# holds. Read B_marg1/2 and xi_gpd1/2 alongside to confirm the mechanism.
#
# IMPORTANT: the tube's survFunc is a closure over inline transform helpers defined
# in the SIMULATION script's global environment; after readRDS in a fresh session
# those are not found. This script REDEFINES them and rebuilds survFunc explicitly
# from the stored dat + trimmed marginal fits (x_sorted, n reconstructed from dat).
library(dplyr)
library(matrixStats)
source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')            # computeTubeSize lives here
source('~/isolines_uq/scripts/R/confidence_regions/modules/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/coverageEvaluation.R')

# NOTE: point save_full_path at wherever the 0.1.1.2 log-gamma bootmargs tubes were saved.
save_full_path <- '/pscratch/sd/j/jbbutler/isolines/simulations/extremal/0.1.1.2_bivt4copula_loggammamargins_gamma1_0.5_gamma2_0.67_bootmargs/'
save_df_path   <- '~/isolines_uq/outputs/simulations/extremal/'

dists <- c('t4copula_loggamma')
ns <- c(1000, 5000, 10000, 50000, 100000)
C <- 5
pn <- function(n){ C/n }
n_coords <- 50
lbs <- c(0, 0)
ubs <- c(200, 200)

# --- log-gamma / copula parameters (MUST match the 0.1.1.2 simulation) ---
lg_shape  <- 2      # a  (a != 1 => marginal rho = 0)
lg_rate   <- 4      # b = 1/xi
copula_df <- 4      # t-copula degrees of freedom

# marginal quantile for the log-gamma margins (for the TRUE isoline)
qloggamma <- function(u) exp(qgamma(u, shape = lg_shape, rate = lg_rate))

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

std_blended_surv <- function(zpts, Z, qn) {
    zpts <- as.matrix(zpts)
    emp <- colMeans(outer(Z[,1], zpts[,1], ">") & outer(Z[,2], zpts[,2], ">"))
    theta <- atan2(zpts[,2], zpts[,1])
    rp <- sqrt(rowSums(zpts^2))
    rq <- std_radial_anchor(Z, theta, qn)
    tail_probs <- (rq / rp) * qn
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

print("Starting evaluation of raw (uncorrected) marginal-transform log-gamma tubes...")

for (dist in dists) {
  for (n in ns) {
    p <- pn(n)
    print(paste0("\nProcessing dist: ", dist, " | n: ", n))

    # TRUE p-isoline: t-copula survival isoline -> pt -> log-gamma quantile
    # (truth on the ORIGINAL scale; survFunc composes the transform internally)
    tcop_isoline   <- drawIsoline(dist='bivt', numCoords=n_coords, gridUbs=ubs, gridLbs=lbs, prob=p)
    copula_isoline <- data.frame(X1 = pt(tcop_isoline[,1], df=copula_df),
                                 X2 = pt(tcop_isoline[,2], df=copula_df))
    true_isoline   <- data.frame(X1 = qloggamma(copula_isoline[,1]),
                                 X2 = qloggamma(copula_isoline[,2]))

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

      # --- reconstruct full marginal fits + fresh survFunc (see header note) ---
      fit1 <- rebuild_fit(tube_obj$fit1, tube_obj$dat[, 1])
      fit2 <- rebuild_fit(tube_obj$fit2, tube_obj$dat[, 2])
      qn   <- tube_obj$fit1$q1                       # = n^(-gamma1)
      Z    <- standardize(tube_obj$dat[, 1], tube_obj$dat[, 2], fit1, fit2)

      tube_obj$survFunc <- function(pts) {
          pts <- as.matrix(pts)
          zq <- standardize(pts[, 1], pts[, 2], fit1, fit2)
          std_blended_surv(zq, Z, qn)
      }

      # swap the RAW bounds into the standard evaluation slots
      tube_obj$c_plus_estimates  <- tube_obj$raw_c_plus_estimates
      tube_obj$c_minus_estimates <- tube_obj$raw_c_minus_estimates

      # evaluate coverage using the raw bounds (transformed = FALSE: survFunc composes internally)
      is_covereds <- evaluateCoverage(tube_obj, true_isoline)
      alphas <- names(is_covereds)

      # --- tube size (RAW walls): standardized-scale wall radii -> backmap -> original radii ---
      rp_std <- tube_obj$r_naive     # standardized p-isoline radius, per angle
      thetas <- tube_obj$thetas
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
          r_in_std  <- rp_std * (p / s_hi)     # xi = 1 on the standardized scale
          r_out_std <- rp_std * (p / s_lo)
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

    save_fname <- paste0(save_df_path, '0.1.1.0_bivt4copula_loggammamargins_gamma1_0.5_bootmargs.csv')

    write.table(coverage_df,
                file = save_fname,
                sep = ",",
                row.names = FALSE,
                col.names = TRUE)

    print(paste0("\nSuccessfully saved results to: ", save_fname))
} else {
    print("\nNo results to save. Check if the directories contain .RData files.")
}