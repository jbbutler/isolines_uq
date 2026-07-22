# Script to loop through saved confidence tubes and evaluate un-bias corrected coverage
library(dplyr)
source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/coverageEvaluation.R')

save_full_path <- '/pscratch/sd/j/jbbutler/isolines/simulations/extremal_logprobability_bias_correction/exp2.2.2/gamma1_0.5_gamma2_0.67/'
save_df_path <- '~/isolines_uq/outputs/simulations/bias_correction/'

dists <- c('bivt')
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
    
    true_isoline <- drawIsoline(dist=dist, numCoords=n_coords, gridUbs=c(200,200), gridLbs=c(0,0), prob=p)
    
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
      
      # load the saved tube object
      tube_obj <- readRDS(file)
      
      # swap the raw bounds into the standard evaluation slots
      # evaluteCoverage assumes desired tube bounds are called c_plus_estimates and c_minus_estimates
      tube_obj$c_plus_estimates <- tube_obj$raw_c_plus_estimates
      tube_obj$c_minus_estimates <- tube_obj$raw_c_minus_estimates
      
      # evaluate coverage using the newly swapped raw bounds
      is_covereds <- evaluateCoverage(tube_obj, true_isoline)
      
      # extract the alphas from the list names
      alphas <- names(is_covereds)
      
      # build the dataframe row (omitting the beta bias fraction columns)
      res_df <- data.frame(
        dist = dist,
        p = p,
        n = n,
        alpha = alphas,
        is_covered = unlist(is_covereds, use.names=FALSE),
        c_plus_raw = unlist(tube_obj$raw_c_plus_estimates[alphas]),
        c_minus_raw = unlist(tube_obj$raw_c_minus_estimates[alphas]),
        xi_hat = tube_obj$xi_hat
      )
      
      # append to our master list
      results_list[[length(results_list) + 1]] <- res_df
      
      setTxtProgressBar(pb, i)
    }
    
    close(pb)
  }
}

# combine all results and save the cleaned CSV
if (length(results_list) > 0) {
    coverage_df <- do.call(rbind, results_list)
      
    save_fname <- paste0(save_df_path, 'extremal_coverage_exp2.2.2_ctrl.csv')
      
    write.table(coverage_df, 
                file = save_fname, 
                sep = ",", 
                row.names = FALSE, 
                col.names = TRUE)
                
    print(paste0("\nSuccessfully saved results to: ", save_fname))
} else {
    print("\nNo results to save. Check if the directories contain .RData files.")
}