# Script to loop through saved confidence tubes and evaluate un-bias corrected coverage
library(dplyr)
source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/coverageEvaluation.R')

save_full_path <- '/pscratch/sd/j/jbbutler/isolines/simulations/extremal/1.0.1.2_gausscopula_t4margs_gamma1_0.5_gamma2_0.67_bootxi_cooleysmooth/'
save_df_path <- '~/isolines_uq/outputs/simulations/extremal/'

dists <- c('gauss_copula_t4')
ns <- c(1000, 5000, 10000, 50000, 100000)
C <- 5
pn <- function(n){ C/n } 
n_coords <- 50

# Initialize a list to hold the results of all simulations
results_list <- list()

print("Starting evaluation of raw tubes...")

for (dist in dists) {
  for (n in ns) {
    p <- pn(n)
    print(paste0("\nProcessing dist: ", dist, " | n: ", n))
    
    gauss_isoline <- drawIsoline(dist='bivgauss', numCoords=n_coords, gridUbs=c(200,200), gridLbs=c(0,0), prob=p)
    copula_isoline <- data.frame(X1=pnorm(gauss_isoline[,1]), X2=pnorm(gauss_isoline[,2]))
    true_isoline <- data.frame(X1=qt(copula_isoline[,1], df=4), X2=qt(copula_isoline[,2], df=4))
    
    # locate the saved .RData files for this specific n and dist
    sim_dir <- paste0(save_full_path, dist, '/n', n, '/')
    
    # grab all .RData files in this folder
    sim_files <- list.files(sim_dir, pattern = "\\.RData$", full.names = TRUE)
    total_files <- length(sim_files)
    
    if (total_files == 0) {
        warning(paste0("No .RData files found in: ", sim_dir, ". Skipping..."))
        next
    }
    
    pb <- txtProgressBar(min = 0, max = total_files, style = 3)
    
    # loop through each saved simulation sequentially using an index
    for (i in 1:total_files) {
      file <- sim_files[i]
      tube_obj <- readRDS(file)

      tube_obj$c_plus_estimates  <- tube_obj$raw_c_plus_estimates
      tube_obj$c_minus_estimates <- tube_obj$raw_c_minus_estimates

      is_covereds <- evaluateCoverage(tube_obj, true_isoline)
      alphas <- names(is_covereds)

      r_ref  <- tube_obj$r_naive          # estimated p-isoline radius, per angle
      xi_eff <- tube_obj$xi_eff           # angle-dependent smoothed exponent

      size_rows <- lapply(alphas, function(a) {
        cplus  <- as.numeric(tube_obj$raw_c_plus_estimates[[a]])
        cminus <- as.numeric(tube_obj$raw_c_minus_estimates[[a]])
        s_hi <- p + cplus                 # inner wall (higher survival, smaller radius)
        s_lo <- p + cminus                # outer wall (lower  survival, larger  radius)

        if (s_lo <= 0) {                  # unbounded (won't trigger on log scale, kept for safety)
          data.frame(is_unbounded = TRUE,
                     width_abs_median = NA_real_, width_abs_max = NA_real_,
                     width_rel_median = NA_real_, width_rel_max = NA_real_)
        } else {
          # re-project the walls from the same anchor the tube used:
          #   r_wall = r_ref * (p / s_wall)^{xi_eff}   (radial HRV/Cooley projection)
          r_inner <- r_ref * (p / s_hi)^xi_eff
          r_outer <- r_ref * (p / s_lo)^xi_eff
          cbind(data.frame(is_unbounded = FALSE),
                computeTubeSize(r_inner, r_outer, r_ref))
        }
      })
      size_df <- do.call(rbind, size_rows)
      # ----------------------------------------------------------------------

      res_df <- data.frame(
        dist = dist, p = p, n = n,
        alpha = alphas,
        is_covered = unlist(is_covereds, use.names=FALSE),
        c_plus_raw  = unlist(tube_obj$raw_c_plus_estimates[alphas]),
        c_minus_raw = unlist(tube_obj$raw_c_minus_estimates[alphas]),
        xi_hat      = tube_obj$xi_hat,
        xi_marg_hat = tube_obj$xi_marg_hat,
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
      
    save_fname <- paste0(save_df_path, '1.0.1.0_gausscopula_t4margs_gamma0.5_bootxi_cooleysmooth.csv')
      
    write.table(coverage_df, 
                file = save_fname, 
                sep = ",", 
                row.names = FALSE, 
                col.names = TRUE)
                
    print(paste0("\nSuccessfully saved results to: ", save_fname))
} else {
    print("\nNo results to save. Check if the directories contain .RData files.")
}