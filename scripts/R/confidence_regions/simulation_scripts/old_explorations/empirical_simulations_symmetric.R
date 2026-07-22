# Script to evaluate the coverage properties of the confidence tube method using the empirical survival function,
# by simulating data from known distributions with different distributional characteristics.
#
# Jimmy Butler
# February 2025

set.seed(45678)

library(dplyr)
library(mvtnorm)
library(mnormt)
library(foreach)
library(doSNOW)
library(parallel)
library(argparse)

source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/confidenceRegions.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/coverageEvaluation.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')

# parsing command line specified parameters
parser <- ArgumentParser(description = "Args for empirical survival function simulations.")
parser$add_argument("--save_full_path", type = "character", default = NULL, required = FALSE,
                    help = "Optional file path to save full details about simulation (data sample, etc.).")
parser$add_argument("--n_cores", type = "integer", default = 64, 
                    help = "Number of cores for parallel processing.")
parser$add_argument("--save_df_path", type = "character", default = NULL, 
                    help = "Optional file path to save just the csv of coverage results for each simulation.")

args <- parser$parse_args()

if (is.null(args$save_full_path) && is.null(args$save_df_path)) {
  stop("Error: Both --save_full_path and --save_df_path are NULL. 
        \n You must provide at least one output destination so results are saved.")
}

# simulation parameters
ns <- c(1000, 5000, 10000, 50000, 100000) # sample sizes to test
alphas <- c(0.05, 0.1, 0.01) # alphas for 1-alpha tubes
C <- 5
p_funcs <- list()
# add in central functions
p_funcs[['central']] <- list()
p_funcs[['central']][['0.05']] <- function(n) 0.05
# add in intermediate functions
p_funcs[['intermediate']] <- list()
p_funcs[['intermediate']][['1_3']] <- function(n) n^(-1/3)
p_funcs[['intermediate']][['1_2']] <- function(n) n^(-1/2)
p_funcs[['intermediate']][['2_3']] <- function(n) n^(-2/3)
p_funcs[['extremal']] <- list()
# add in extremal functions
p_funcs[['extremal']][[as.character(C)]] <- function(n) C/n

dists <- c('bivt', 'bivgauss', 'karachi') # distributions to simulate from
# list of lower left hand corners of space, depending on distribution
lbs_lst <- list('bivt'=c(-2,-2), 'bivgauss'=c(-2,-2), 'karachi'=c(50, 0))
# list of upper right hand corners of space, depending on distribution
# note: this is for determining which points are on the isoline, not actually used to draw tubes
# need to be big enough for root-finding algorithm
ubs_lst <- list('bivt'=c(200, 200), 'bivgauss'=c(200,200), 'karachi'=c(140,100))

n_sims <- 500 # number of simulations for coverage statistics
B <- 1000 # number of bootstrap iterations for each confidence tube

# command line parameters
n_cores <- args$n_cores
save_full_path <- args$save_full_path

# part of loop that will be parallelized: the computation of c for each simulation
parallelizedCode <- function(ind) {

    # sample dataset
    dat <- sampling_func(n)
    # compute confidence regions for different alpha, for a particular p
    regions <- computeEmpiricalRegion(dat=dat,
                                 alphas=alphas,
                                 p=p,
                                 B=B,
                                 lbs=lbs,
                                 verbose=FALSE)

    # evaluate coverage for each region
    is_covereds <- evaluateCoverage(regions, isoline)
    # add into the list
    regions$is_covereds <- is_covereds
    
    return(regions)
}

# loop over the distributions
for (i in 1:length(dists)) {
    dist <- dists[i]
    lbs <- lbs_lst[[dist]]
    ubs <- ubs_lst[[dist]]
    
    sampling_func <- loadSamplingFunction(dist)
    print(paste0('Starting new distribution: ', dist))

    # loop over regimes
    for (regime in names(p_funcs)) {
        
        for (label in names(p_funcs[[regime]])) {
            
            current_p_func <- p_funcs[[regime]][[label]]
            # define path specific to Dist -> Regime -> Label
            if (!is.null(args$save_full_path)) {
                save_path_full <- paste0(args$save_full_path, dist, '/', regime, '/', label, '/')
                if (!dir.exists(save_path_full)) dir.create(save_path_full, recursive = TRUE)
            }

            print(paste0('  Running Regime: ', regime, ' | Label: ', label))
            
            for (k in 1:length(ns)) {
                n <- ns[k]
                p <- current_p_func(n)
                isoline <- drawIsoline(dist=dist, numCoords=200, gridUbs=ubs, gridLbs=lbs, prob=p)

                print(paste0('Starting new n: ', n))
    
                clust <- makeSOCKcluster(n_cores)
                registerDoSNOW(clust)
                clusterEvalQ(clust, {
                    source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')
                })
                
                pb <- txtProgressBar(min = 1, max = n_sims, style = 3)
                progress <- function(n) setTxtProgressBar(pb ,n)
                opts <- list(progress = progress)
    
                samp_regions <- foreach(l = 1:n_sims, 
                        .options.snow = opts, 
                        .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt', 'matrixStats')) %dopar% parallelizedCode(l)
    
                close(pb)
                stopCluster(clust)
    
                # only save .RData if save_full_path is provided
                if (!is.null(args$save_full_path)) {
                    saveRDS(samp_regions, file=paste0(save_path_full, 'n', n, '_tubes.RData'))
                }
                # save CSV if save_df_path argument is provided
                if (!is.null(args$save_df_path)) {
                    
                    coverage_df <- do.call(rbind, lapply(samp_regions, function(res) {
            
                    data.frame(
                      dist = dist,
                      regime = regime,
                      label = label,
                      n = n,
                      p = p,
                      alpha = names(res$is_covereds),
                      is_covered = unlist(res$is_covereds, use.names=FALSE),
                      c_estimate = unlist(res$c_estimates[names(res$is_covereds)])
                    )
                  }))
          
                  save_fname <- paste0(args$save_df_path, 'empirical_coverage.csv')
                  file_exists <- file.exists(save_fname)
          
                  write.table(coverage_df, 
                          file = save_fname, 
                          sep = ",", 
                          row.names = FALSE, 
                          col.names = !file_exists, 
                          append = file_exists)
                }
            }
        }
    }
}