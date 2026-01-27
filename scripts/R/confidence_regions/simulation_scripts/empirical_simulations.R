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
ps <- c(0.1, 0.01, 0.001, 0.0001) # exceedance probabilities
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
    
    # only create dist directory if save_full_path is provided
    if (!is.null(save_full_path)) {
        distpath <- paste0(save_full_path, dist, '/')
        if (!dir.exists(distpath)) dir.create(distpath, recursive = TRUE)
    }
    print(paste0('Starting new distribution: ', dist))
    
    for (j in 1:length(ps)) {
        p <- ps[j]
        # only create p directory if save_full_path is provided
        if (!is.null(save_full_path)) {
            ppath <- paste0(distpath, 'p', p, '/')
            if (!dir.exists(ppath)) dir.create(ppath, recursive = TRUE)
        }
        
        # create the true p-isoline for this distribution
        isoline <- drawIsoline(dist=dist, numCoords=200, gridUbs=ubs, gridLbs=lbs, prob=p)

        print(paste0('Starting new p: ', p))
        
        for (k in 1:length(ns)) {
            n <- ns[k]

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
                    .packages = c('mvtnorm', 'data.table', 'dplyr', 'mnormt')) %dopar% parallelizedCode(l)

            close(pb)
            stopCluster(clust)

            # only save .RData if save_full_path is provided
            if (!is.null(save_full_path)) {
                saveRDS(samp_regions, file=paste0(ppath, 'n', n, '_tubes.RData'))
            }

            # save CSV if save_df_path argument is provided
            if (!is.null(args$save_df_path)) {
                
                # combine results from the parallel list into a single data frame
                coverage_df <- do.call(rbind, lapply(samp_regions, function(res) {
                    data.frame(
                        dist = dist,
                        p = p,
                        n = n,
                        alpha = names(res$is_covereds), 
                        is_covered = unlist(res$is_covereds, use.names=FALSE)
                    )
                }))

                # determine if we need to write headers (only if file doesn't exist)
                save_fname <- paste0(args$save_df_path, 'empirical_coverage.csv')
                file_exists <- file.exists(save_fname)
                
                # write to CSV using write.table for append support
                write.table(coverage_df, 
                            file = save_fname, 
                            sep = ",", 
                            row.names = FALSE, 
                            col.names = !file_exists, # Write headers only if new file
                            append = file_exists)     # Append if file exists
            }
        }
    }
}