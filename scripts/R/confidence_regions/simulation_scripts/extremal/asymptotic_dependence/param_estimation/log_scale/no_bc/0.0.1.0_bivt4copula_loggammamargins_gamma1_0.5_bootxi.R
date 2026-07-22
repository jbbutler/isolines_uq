# Script to loop through saved log-gamma confidence tubes and evaluate un-bias-corrected coverage.
#
# Loads the saved tube objects from the t4-copula / log-gamma-margin simulations
# (0.0.1.2), swaps the RAW (variance-only) walls into the evaluation slots, and
# re-evaluates coverage of the UNCORRECTED tubes against the true log-gamma isoline.
# Also records the size of the raw tubes.
library(dplyr)
library(matrixStats)
source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')            # computeTubeSize lives here
source('~/isolines_uq/scripts/R/confidence_regions/modules/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/coverageEvaluation.R')

# NOTE: point save_full_path at wherever the 0.0.1.2 log-gamma tubes were saved.
save_full_path <- '/pscratch/sd/j/jbbutler/isolines/simulations/extremal/0.0.1.2_bivt4copula_loggammamargins_gamma1_0.5_gamma2_0.67_bootxi/'
save_df_path   <- '~/isolines_uq/outputs/simulations/extremal/'

dists <- c('t4copula_loggamma')
ns <- c(1000, 5000, 10000, 50000, 100000)
C <- 5
pn <- function(n){ C/n }
n_coords <- 50

# --- log-gamma / copula parameters (MUST match the 0.0.1.2 simulation that made the tubes) ---
lg_shape  <- 2      # a  (a != 1 => rho = 0)
lg_rate   <- 4      # b = 1/xi
copula_df <- 4      # t-copula degrees of freedom (nu)

# marginal quantile for the log-gamma margins (for the TRUE isoline)
qloggamma <- function(u) exp(qgamma(u, shape = lg_shape, rate = lg_rate))

# Initialize a list to hold the results of all simulations
results_list <- list()

print("Starting evaluation of raw (uncorrected) log-gamma tubes...")

for (dist in dists) {
  for (n in ns) {
    p <- pn(n)
    print(paste0("\nProcessing dist: ", dist, " | n: ", n))

    # TRUE p-isoline for log-gamma margins: t-copula survival isoline -> pt -> log-gamma quantile
    # (same construction as the 0.0.1.2 simulation script; NOT drawIsoline('bivt') directly)
    tcop_isoline   <- drawIsoline(dist='bivt', numCoords=n_coords, gridUbs=c(200,200), gridLbs=c(0,0), prob=p)
    copula_isoline <- data.frame(X1 = pt(tcop_isoline[,1], df=copula_df),
                                 X2 = pt(tcop_isoline[,2], df=copula_df))
    true_isoline   <- data.frame(X1 = qloggamma(copula_isoline[,1]),
                                 X2 = qloggamma(copula_isoline[,2]))

    # locate the saved .RData files for this specific n and dist
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

      # load the saved tube object
      tube_obj <- readRDS(file)

      # swap the RAW bounds into the standard evaluation slots
      # evaluateCoverage assumes desired tube bounds are called c_plus_estimates / c_minus_estimates
      tube_obj$c_plus_estimates  <- tube_obj$raw_c_plus_estimates
      tube_obj$c_minus_estimates <- tube_obj$raw_c_minus_estimates

      # evaluate coverage using the newly swapped raw bounds
      is_covereds <- evaluateCoverage(tube_obj, true_isoline)
      alphas <- names(is_covereds)

      # --- tube size (RAW / uncorrected walls) ------------------------------
      # invert the walls to radii with the same radial projection the tube used:
      #   r(s) = r_naive * (p / s)^{xi_hat}. Log-scale => outer wall s_lo > 0 always.
      r_ref  <- tube_obj$r_naive     # estimated p-isoline radius, per angle
      xi_hat <- tube_obj$xi_hat

      size_rows <- lapply(alphas, function(a) {
        cplus  <- as.numeric(tube_obj$raw_c_plus_estimates[[a]])
        cminus <- as.numeric(tube_obj$raw_c_minus_estimates[[a]])
        s_hi <- p + cplus              # inner wall (higher survival, smaller radius)
        s_lo <- p + cminus             # outer wall (lower  survival, larger  radius)

        if (s_lo <= 0) {               # unbounded (won't trigger on log scale; kept for safety)
          data.frame(is_unbounded = TRUE,
                     width_abs_median = NA_real_, width_abs_max = NA_real_,
                     width_rel_median = NA_real_, width_rel_max = NA_real_)
        } else {
          r_inner <- r_ref * (p / s_hi)^xi_hat
          r_outer <- r_ref * (p / s_lo)^xi_hat
          cbind(data.frame(is_unbounded = FALSE),
                computeTubeSize(r_inner, r_outer, r_ref))
        }
      })
      size_df <- do.call(rbind, size_rows)

      # build the dataframe row (raw walls + raw tube sizes; no beta bias-fraction columns)
      res_df <- data.frame(
        dist = dist,
        p = p,
        n = n,
        alpha = alphas,
        is_covered = unlist(is_covereds, use.names=FALSE),
        c_plus_raw  = unlist(tube_obj$raw_c_plus_estimates[alphas]),
        c_minus_raw = unlist(tube_obj$raw_c_minus_estimates[alphas]),
        xi_hat = tube_obj$xi_hat,
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

    save_fname <- paste0(save_df_path, '0.0.1.0_bivt4copula_loggammamargins_gamma1_0.5_bootxi.csv')

    write.table(coverage_df,
                file = save_fname,
                sep = ",",
                row.names = FALSE,
                col.names = TRUE)

    print(paste0("\nSuccessfully saved results to: ", save_fname))
} else {
    print("\nNo results to save. Check if the directories contain .RData files.")
}