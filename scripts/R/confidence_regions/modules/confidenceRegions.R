# Module with functions to draw confidence tubes for isolines, using both the extreme
# and non-extreme methodologies.
#
# Jimmy Butler

library(ks)
library(dplyr)
library(ismev)

source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')
computeExtremeRegion <- function(dat, alphas, p, B, gamma, xi, lbs=c(0,0), verbose=FALSE) {
    # Function that constructs confidence tube(s) given a particular dataset and a desired isoline
    # exceedance probability p. Multiple confidence tubes will be returned if multiple alphas are
    # supplied, one for each alpha. NOTE: this function computes c_hat, an estimate of the worst-case difference
    # between the estimated survival function and the true survival function over the set of points corresponding
    # to the desired true isoline. This quantity + the dataset defines the confidence tube: to actually draw it
    # visually, use the drawExtremeRegion function.
    # 
    # Arguments:
    # dat: the data, in the form of a 2-column R data.frame
    # alphas: a vector of alphas indicating desired probability of miscoverage
    # (0.01 for 99% CI, 0.05 for 95% CI, etc.). NOTE: even if only one alpha desired, must put it in a vector.
    # p: the desired p-isoline you wish to capture
    # B: the number of bootstrap replicates in determining c_hat
    # lbs: the lower lefthand corner of your space
    # gamma: 1/n^(gamma), controls the smallest probability to use empirical df
    # xi: extreme value index
    # verbose: progress bar for bootstrap loop?
    #
    # Output:
    # A list with (1) the original data, (2) a vector of c_hat estimates, one for each alpha supplied.
    # (3) the desired exceedance probability

    ext_isoline <- drawExtremeIsoline(dat, p, n_coords=200, gridLbs=lbs, gamma, xi)
    boot_draws <- rep(NA, B)

    if (verbose) {
        pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
        for (k in 1:B) {
            boot_dat <- dat %>% sample_frac(1, replace = TRUE)
            boot_survfunc <- apply(ext_isoline, 1, 
                                   blendedSurvivalFunc, dat=boot_dat, 
                                   gamma=gamma, xi=xi)
            boot_draws[k] <- max(abs(boot_survfunc - p))
            utils::setTxtProgressBar(pb, k)
        }
        close(pb)
    } else {
        for (k in 1:B) {
            boot_dat <- dat %>% sample_frac(1, replace = TRUE)
            boot_survfunc <- apply(ext_isoline, 1, 
                                   blendedSurvivalFunc, dat=boot_dat, 
                                   gamma=gamma, xi=xi)
            boot_draws[k] <- max(abs(boot_survfunc - p))
        }
    }

    c_estimates <- list()
    for (i in 1:length(alphas)) {
        alpha <- alphas[i]       
        c_estimate <- as.numeric(quantile(boot_draws, probs = 1-alpha))
        c_estimates[as.character(alpha)] <- c_estimate
    }

    survFunc <- function(x) {
        survProb <- blendedSurvivalFunc(x, dat, gamma, xi)
        return(survProb)
    }

    res_lst <- list()
    res_lst$dat <- dat
    res_lst$c_estimates <- c_estimates
    res_lst$p <- p
    res_lst$gamma <- gamma
    res_lst$xi <- xi
    res_lst$survFunc <- survFunc

    return(res_lst)

}

computeEmpiricalRegion <- function(dat, alphas, p, B, lbs, verbose=FALSE) {
    # Function that constructs confidence tube(s) for a p-isoline of interest, given some data.
    # Will return multiple tube constructions if multiple alphas are passed as an argument (NOTE:
    # even if only one alpha is desired, you still must put it in an R vector). This function
    # uses the empirical survival function, and thus is more appropriate for use with non-extreme
    # isolines. As above, this function only computes the quantities needed to define
    # a confidence tube (i.e., c_hat). If you would like to obtain points in the tube so as to visualize
    # it, use the drawEmpiricalRegion function.
    #
    # Arguments:
    # dat: the data, in the form of a 2-column data.frame
    # alphas: a vector of alphas (even if one alpha desired, must put it in a vector as well)
    # p: the desired p-isoline to capture
    # B: the number of bootstrap replicates
    # lbs: the lower lefthand corner of the space on which the tube is drawn
    # verbose: boostrap loop progress bar?
    #
    # Output:
    # A list with (1) the original data, (2) a list of the c_hat estimates, one for each supplied alpha,
    # and (3) the desired exceedance probability.

    emp_isoline <- drawEmpiricalIsoline(dat=dat, n_coords=200, gridLbs=lbs, p)
    boot_draws <- rep(NA, B)

    if (verbose) {
        pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
        for (k in 1:B) {
            boot_dat <- dat %>% sample_frac(1, replace = TRUE)
            boot_survfunc <- computeEmpSurvIrregular(emp_isoline, boot_dat)
            boot_draws[k] <- max(abs(boot_survfunc-p))
            utils::setTxtProgressBar(pb, k)
        }
        close(pb)
    } else {
       for (k in 1:B) {
            boot_dat <- dat %>% sample_frac(1, replace = TRUE)
            boot_survfunc <- computeEmpSurvIrregular(emp_isoline, boot_dat)
            boot_draws[k] <- max(abs(boot_survfunc-p))
        }
    }
    c_estimates <- list()
    for (i in 1:length(alphas)) {
        alpha <- alphas[i]       
        c_estimate <- as.numeric(quantile(boot_draws, probs = 1-alpha))
        c_estimates[as.character(alpha)] <- c_estimate
    }

    survFunc <- function(x) {
        survProb <- mean((dat[,1] > x[1]) & (dat[,2] > x[2]))
        return(survProb)
    }

    res_lst <- list()
    res_lst$dat <- dat
    res_lst$c_estimates <- c_estimates
    res_lst$p <- p
    res_lst$survFunc <- survFunc

    return(res_lst)
}