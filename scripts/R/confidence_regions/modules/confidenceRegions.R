# Module with functions to draw confidence tubes for isolines, using both the extreme
# and non-extreme methodologies.
#
# Jimmy Butler

library(dplyr)
library(ismev)

source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')
computeExtremeRegion <- function(dat, alphas, p, B, gamma, gamma1=2/3, gamma2=1/2, ncoords=100, n_lines=10, xi=1, lbs=c(0,0), sup_region=FALSE, transformed=FALSE, progress_bar=FALSE) {
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
    # b1, b2: lower and upper thresholds around probability p to define region to take sup over
    # B: the number of bootstrap replicates in determining c_hat
    # lbs: the lower lefthand corner of your space
    # gamma: 1/n^(gamma), controls the smallest probability to use empirical df
    # xi: extreme value index of the data that will be used to construct tube (transformed if not already MRV)
    # progress_bar: progress bar for bootstrap loop?
    #
    # Output:
    # A list with (1) the original data, (2) a vector of c_hat estimates, one for each alpha supplied.
    # (3) the desired exceedance probability

    if (transformed) {

        # define a function for the transform
        transform <- function(pts, marginal_dat) {
            transformed_pts <- 1/(1-est_cdf(pts, marginal_dat, gamma)) - 1
        return(transformed_pts)
        }

        # define a function for the inverse transform
        inv_transform <- function(pts, marginal_dat) {
            orig_pts <- est_inv_cdf(pts/(1+pts), dat=marginal_dat, gamma=gamma)
        return(orig_pts)
        }

        # transform the data. We will now only work with the transformed data.
        transformed_dat_X1 <- transform(dat[,1], dat[,1])
        transformed_dat_X2 <- transform(dat[,2], dat[,2])
        # construction_dat is name for data used to construct tube
        construction_dat <- data.frame(X1=transformed_dat_X1, X2=transformed_dat_X2)

    } else {
        construction_dat <- dat
    }

    # if we're going to be constructing a sup region over which to take maxima instead of just the
    # estimated isoline
    if (sup_region) {
        bounds <- estimate_bias_bounds(construction_dat, gamma1, gamma2, p, xi)
        b1 <- min(bounds$b1, 0)
        b2 <- max(bounds$b2, 0)
        p_seq <- seq(p + b1, p + b2, length.out = n_lines)
        # generate points for all isolines in the sequence and stack them
        region_pts_list <- lapply(p_seq, function(prob) {
            drawExtremeIsoline(construction_dat, prob, n_coords=ncoords, grid_lbs=lbs, gamma, xi)
        })
        ext_region_pts <- do.call(rbind, region_pts_list)
        
    } else {
        ext_region_pts <- drawExtremeIsoline(construction_dat, p, n_coords=ncoords, grid_lbs=lbs, gamma, xi)
    }

    theta <- atan2(ext_region_pts[,2], ext_region_pts[,1])
    rp <- sqrt(rowSums(ext_region_pts^2))
    inv_cos <- 1/cos(theta)
    inv_sin <- 1/sin(theta)

    proj_x <- construction_dat[,1] %o% inv_cos
    proj_y <- construction_dat[,2] %o% inv_sin
    Z_static <- pmin(proj_x, proj_y)

    match_x <- outer(construction_dat[,1], ext_region_pts[,1], '>')
    match_y <- outer(construction_dat[,2], ext_region_pts[,2], '>')
    hit_static <- (match_x & match_y)*1

    n_dat <- nrow(construction_dat)
    qn <- n_dat^(-gamma)
    
    boot_draws_Wplus <- rep(NA, B)
    boot_draws_Wminus <- rep(NA, B)

    if (progress_bar) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (k in 1:B) {
        idx <- sample.int(n_dat, n_dat, replace=TRUE)
        Z_boot <- Z_static[idx, , drop=FALSE]
        hit_boot <- hit_static[idx, , drop=FALSE]
        emp_probs <- colMeans(hit_boot)
        
        rq <- colQuantiles(Z_boot, probs=1-qn, type=1, drop=TRUE)
        tail_probs <- ((rq/rp)^(1/xi))*qn
        final_probs <- ifelse(emp_probs >= qn, emp_probs, tail_probs)
        
        # Taking the max/min over the entire defined region minus the base p
        boot_draws_Wplus[k] <- max(final_probs - p)
        boot_draws_Wminus[k] <- min(final_probs - p)
        
        if (progress_bar) utils::setTxtProgressBar(pb, k)
    }
    if (progress_bar) close(pb)

    c_plus_estimates <- list()
    c_minus_estimates <- list()
    for (i in 1:length(alphas)) {
        alpha <- alphas[i]       
        c_plus_estimate <- as.numeric(quantile(boot_draws_Wplus, probs = 1-(alpha/2)))
        c_minus_estimate <- as.numeric(quantile(boot_draws_Wminus, probs = (alpha/2)))
        c_plus_estimates[as.character(alpha)] <- c_plus_estimate
        c_minus_estimates[as.character(alpha)] <- c_minus_estimate
    }

    survFunc <- function(pts) {
        survProb <- vectorizedBlendedSurv(pts, construction_dat, gamma, xi)
        return(survProb)
    }

    res_lst <- list()
    res_lst$dat <- dat
    res_lst$transformed <- transformed
    res_lst$c_plus_estimates <- c_plus_estimates
    res_lst$c_minus_estimates <- c_minus_estimates
    res_lst$p <- p
    res_lst$gamma <- gamma
    res_lst$xi <- xi
    res_lst$survFunc <- survFunc

    if (sup_region) {
        res_lst$gamma1 <- gamma1
        res_lst$gamma2 <- gamma2
        res_lst$b1 <- b1
        res_lst$b2 <- b2
    }

    if (transformed) {
        res_lst$transform_func1 <- function(pts) transform(pts, dat[,1])
        res_lst$transform_func2 <- function(pts) transform(pts, dat[,2])
        res_lst$inv_transform_func1 <- function(pts) inv_transform(pts, dat[,1])
        res_lst$inv_transform_func2 <- function(pts) inv_transform(pts, dat[,2])  
    }

    return(res_lst)
}

computeExtremeRegionBC <- function(dat, alphas, p, B, B_bias, gamma, q1, q2, ncoords=100, n_lines=10, xi=1, lbs=c(0,0), 
                                   sup_region=FALSE, transformed=FALSE, progress_bar=FALSE, verbose=FALSE) {
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
    # b1, b2: lower and upper thresholds around probability p to define region to take sup over
    # B: the number of bootstrap replicates in determining c_hat
    # lbs: the lower lefthand corner of your space
    # gamma: 1/n^(gamma), controls the smallest probability to use empirical df
    # xi: extreme value index of the data that will be used to construct tube (transformed if not already MRV)
    # progress_bar: progress bar for bootstrap loop?
    #
    # Output:
    # A list with (1) the original data, (2) a vector of c_hat estimates, one for each alpha supplied.
    # (3) the desired exceedance probability

    if (transformed) {

        # define a function for the transform
        transform <- function(pts, marginal_dat) {
            transformed_pts <- 1/(1-est_cdf(pts, marginal_dat, gamma)) - 1
        return(transformed_pts)
        }

        # define a function for the inverse transform
        inv_transform <- function(pts, marginal_dat) {
            orig_pts <- est_inv_cdf(pts/(1+pts), dat=marginal_dat, gamma=gamma)
        return(orig_pts)
        }

        # transform the data. We will now only work with the transformed data.
        transformed_dat_X1 <- transform(dat[,1], dat[,1])
        transformed_dat_X2 <- transform(dat[,2], dat[,2])
        # construction_dat is name for data used to construct tube
        construction_dat <- data.frame(X1=transformed_dat_X1, X2=transformed_dat_X2)

    } else {
        construction_dat <- dat
    }

    # if we're going to be constructing a sup region over which to take maxima instead of just the
    # estimated isoline
    if (sup_region) {
        if (verbose) print('Constructing bias corrected region via bootstrap.')
        # beginning bootstrap bias estimation #
        # construct the empirical isoline and evaluate extremal survival function over it via bootstrap
        emp_isoline <- drawEmpiricalIsoline(dat=construction_dat, n_coords=ncoords, grid_lbs=c(0,0), p=q2)

        # set up hit matrix for bootstrap sampling
        theta <- atan2(emp_isoline[,2], emp_isoline[,1])
        rp <- sqrt(rowSums(emp_isoline^2))
        inv_cos <- 1/cos(theta)
        inv_sin <- 1/sin(theta)

        proj_x <- construction_dat[,1] %o% inv_cos
        proj_y <- construction_dat[,2] %o% inv_sin
        Z_static <- pmin(proj_x, proj_y)

        match_x <- outer(construction_dat[,1], emp_isoline[,1], '>')
        match_y <- outer(construction_dat[,2], emp_isoline[,2], '>')
        hit_static <- (match_x & match_y)*1

        n_dat <- nrow(construction_dat)
    
        boot_ext_survfunc <- rep(NA, B_bias)
        boot_prob_matrix <- matrix(NA, nrow=B_bias, ncol=length(theta))

        if (progress_bar) pb <- utils::txtProgressBar(min = 0, max = B_bias, style = 3)
        for (k in 1:B_bias) {
            idx <- sample.int(n_dat, n_dat, replace=TRUE)
            Z_boot <- Z_static[idx, , drop=FALSE]
            hit_boot <- hit_static[idx, , drop=FALSE]
            emp_probs <- colMeans(hit_boot)
        
            rq <- colQuantiles(Z_boot, probs=1-q1, type=1, drop=TRUE)
            tail_probs <- ((rq/rp)^(1/xi))*q1
            boot_prob_matrix[k,] <- ifelse(emp_probs >= q1, emp_probs, tail_probs)

            if (progress_bar) utils::setTxtProgressBar(pb, k)
        }
        if (progress_bar) close(pb)
        
        biases <- (q2/colMeans(boot_prob_matrix)) - 1
        b1 <- p/(1+max(max(biases), 0))
        b2 <- p/(1+min(min(biases), 0))
        
        p_seq <- seq(b1, b2, length.out = n_lines)
        # generate points for all isolines in the sequence and stack them
        region_pts_list <- lapply(p_seq, function(prob) {
            drawExtremeIsoline(construction_dat, prob, n_coords=ncoords, grid_lbs=lbs, gamma, xi)
        })
        ext_region_pts <- do.call(rbind, region_pts_list)
        
    } else {
        ext_region_pts <- drawExtremeIsoline(construction_dat, p, n_coords=ncoords, grid_lbs=lbs, gamma, xi)
    }

    if (verbose) print('Estimating tube size parameters via bootstrap.')

    theta <- atan2(ext_region_pts[,2], ext_region_pts[,1])
    rp <- sqrt(rowSums(ext_region_pts^2))
    inv_cos <- 1/cos(theta)
    inv_sin <- 1/sin(theta)

    proj_x <- construction_dat[,1] %o% inv_cos
    proj_y <- construction_dat[,2] %o% inv_sin
    Z_static <- pmin(proj_x, proj_y)

    match_x <- outer(construction_dat[,1], ext_region_pts[,1], '>')
    match_y <- outer(construction_dat[,2], ext_region_pts[,2], '>')
    hit_static <- (match_x & match_y)*1

    qn <- n_dat^(-gamma)
    
    boot_draws_Wplus <- rep(NA, B)
    boot_draws_Wminus <- rep(NA, B)

    if (progress_bar) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (k in 1:B) {
        idx <- sample.int(n_dat, n_dat, replace=TRUE)
        Z_boot <- Z_static[idx, , drop=FALSE]
        hit_boot <- hit_static[idx, , drop=FALSE]
        emp_probs <- colMeans(hit_boot)
        
        rq <- colQuantiles(Z_boot, probs=1-qn, type=1, drop=TRUE)
        tail_probs <- ((rq/rp)^(1/xi))*qn
        final_probs <- ifelse(emp_probs >= qn, emp_probs, tail_probs)
        
        # Taking the max/min over the entire defined region minus the base p
        boot_draws_Wplus[k] <- max(final_probs - p)
        boot_draws_Wminus[k] <- min(final_probs - p)
        
        if (progress_bar) utils::setTxtProgressBar(pb, k)
    }
    if (progress_bar) close(pb)

    c_plus_estimates <- list()
    c_minus_estimates <- list()
    for (i in 1:length(alphas)) {
        alpha <- alphas[i]       
        c_plus_estimate <- as.numeric(quantile(boot_draws_Wplus, probs = 1-(alpha/2)))
        c_minus_estimate <- as.numeric(quantile(boot_draws_Wminus, probs = (alpha/2)))
        c_plus_estimates[as.character(alpha)] <- c_plus_estimate
        c_minus_estimates[as.character(alpha)] <- c_minus_estimate
    }

    survFunc <- function(pts) {
        survProb <- vectorizedBlendedSurv(pts, construction_dat, gamma, xi)
        return(survProb)
    }

    res_lst <- list()
    res_lst$dat <- dat
    res_lst$transformed <- transformed
    res_lst$c_plus_estimates <- c_plus_estimates
    res_lst$c_minus_estimates <- c_minus_estimates
    res_lst$p <- p
    res_lst$gamma <- gamma
    res_lst$xi <- xi
    res_lst$survFunc <- survFunc

    if (sup_region) {
        res_lst$b1 <- b1
        res_lst$b2 <- b2
    }

    if (transformed) {
        res_lst$transform_func1 <- function(pts) transform(pts, dat[,1])
        res_lst$transform_func2 <- function(pts) transform(pts, dat[,2])
        res_lst$inv_transform_func1 <- function(pts) inv_transform(pts, dat[,1])
        res_lst$inv_transform_func2 <- function(pts) inv_transform(pts, dat[,2])  
    }

    return(res_lst)
}

computeExtremeRegionSplit <- function(dat, alphas, p, B, gamma, xi, lbs=c(0,0), verbose=FALSE) {
    # Function that constructs confidence tube(s) given a particular dataset and a desired isoline
    # exceedance probability p. Multiple confidence tubes will be returned if multiple alphas are
    # supplied, one for each alpha. NOTE: this function computes c_hat, an estimate of the worst-case difference
    # between the estimated survival function and the true survival function over the set of points corresponding
    # to the desired true isoline. This quantity + the dataset defines the confidence tube: to actually draw it
    # visually, use the drawExtremeRegion function.
    # 
    # HOW THIS DIFFERS FROM OTHER FUNCTIONS: we split the sample in half: half used to estimate isoline, half for bootstrap
    # resampling
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

    n_dat <- nrow(dat)
    split_size <- floor(n_dat/2)
    split1 <- sample(seq_len(n_dat), size=split_size, replace=FALSE)
    iso_dat <- dat[split1,]
    boot_dat <- dat[-split1,]
    
    ext_isoline <- drawExtremeIsoline(iso_dat, p, n_coords=200, grid_lbs=lbs, gamma, xi)
    boot_draws <- rep(NA, B)

    theta <- atan2(ext_isoline[,2], ext_isoline[,1])
    rp <- sqrt(rowSums(ext_isoline^2))
    inv_cos <- 1/cos(theta)
    inv_sin <- 1/sin(theta)

    proj_x <- boot_dat[,1] %o% inv_cos
    proj_y <- boot_dat[,2] %o% inv_sin
    Z_static <- pmin(proj_x, proj_y)

    match_x <- outer(boot_dat[,1], ext_isoline[,1], '>')
    match_y <- outer(boot_dat[,2], ext_isoline[,2], '>')
    hit_static <- (match_x & match_y)*1

    n_boot <- nrow(boot_dat)
    qn <- (n_boot)^(-gamma)
    if (verbose) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (k in 1:B) {
        idx <- sample.int(n_boot, n_boot, replace=TRUE)
        Z_boot <- Z_static[idx, , drop=FALSE]
        hit_boot <- hit_static[idx, , drop=FALSE]
        emp_probs <- colMeans(hit_boot)
        rq <- colQuantiles(Z_boot, probs=1-qn, type=1, drop=TRUE)
        tail_probs <- (rq/rp)*qn
        final_probs <- ifelse(emp_probs >= qn, emp_probs, tail_probs)
        boot_draws[k] <- max(abs(final_probs - p))
        if (verbose) utils::setTxtProgressBar(pb, k)
    }
    if (verbose) close(pb)

    c_estimates <- list()
    for (i in 1:length(alphas)) {
        alpha <- alphas[i]       
        c_estimate <- as.numeric(quantile(boot_draws, probs = 1-alpha))
        c_estimates[as.character(alpha)] <- c_estimate
    }

    # store the vectorized survival function computed on your dataset
    # to be used to easily compute coverage of a desired isoline
    survFunc <- function(pts) {
        survProb <- vectorizedBlendedSurv(pts, iso_dat, gamma, xi)
        return(survProb)
    }

    res_lst <- list()
    res_lst$trans_dat <- dat
    res_lst$c_estimates <- c_estimates
    res_lst$p <- p
    res_lst$gamma <- gamma
    res_lst$xi <- xi
    res_lst$survFunc <- survFunc

    return(res_lst)

}

computeExtremeRegionThreewayIsoBootSurvTwoCs <- function(dat, alphas, p, B, gamma, xi, lbs=c(0,0), verbose=FALSE) {
    # Function that constructs confidence tube(s) given a particular dataset and a desired isoline
    # exceedance probability p. Multiple confidence tubes will be returned if multiple alphas are
    # supplied, one for each alpha.
    # 
    # HOW THIS DIFFERS FROM OTHER FUNCTIONS: we split the sample in thirds: one for estimating the isoline, one for
    # bootstrap sampling, and one for constructing the survival function. First two used for chat, last used for survfunc
    # resampling. Also, instead of computing a single thickness paramater, we compute two, a c_minus and c_plus
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

    # split sample into three groups: one for transform, one for isoline estimation, one for bootstrapping
    group_assignments <- dat %>%
        mutate(group_assignment = sample(rep_len(1:3, length.out = n())))

    survfunc_dat <- group_assignments %>% filter(group_assignment==1) %>% select(-group_assignment)
    iso_dat <- group_assignments %>% filter(group_assignment==2) %>% select(-group_assignment)
    boot_dat <- group_assignments %>% filter(group_assignment==3) %>% select(-group_assignment)

    # define the transforms here
    # the full dataset is transformed using the dataset itself..
    transform <- function(pts, marginal_dat) {
        transformed_pts <- 1/(1-est_cdf(pts, marginal_dat, gamma)) - 1
        return(transformed_pts)
    }
    transformed_dat_X1 <- transform(dat[,1], dat[,1])
    transformed_dat_X2 <- transform(dat[,2], dat[,2])
    transformed_dat <- data.frame(X1=transformed_dat_X1, X2=transformed_dat_X2)

    inv_transform <- function(pts, marginal_dat) {
        orig_pts <- est_inv_cdf(pts/(1+pts), dat=marginal_dat, gamma=gamma)
        return(orig_pts)
    }

    # transform marginals of isoline, bootstrap, and survfunc datasets
    transformed_iso_dat <- data.frame(X1=transform(iso_dat[,1], dat[,1]), 
                                      X2=transform(iso_dat[,2], dat[,2]))
    transformed_boot_dat <- data.frame(X1=transform(boot_dat[,1], dat[,1]), 
                                      X2=transform(boot_dat[,2], dat[,2]))
    transformed_survfunc_dat <- data.frame(X1=transform(survfunc_dat[,1], dat[,1]), 
                                      X2=transform(survfunc_dat[,2], dat[,2]))
    
    ext_isoline <- drawExtremeIsoline(dat=transformed_iso_dat, 
                                      p=p,
                                      n_coords=200, 
                                      grid_lbs=lbs, 
                                      gamma=gamma, 
                                      xi=1)

    boot_draws_Wplus <- rep(NA, B)
    boot_draws_Wminus <- rep(NA, B)

    theta <- atan2(ext_isoline[,2], ext_isoline[,1])
    rp <- sqrt(rowSums(ext_isoline^2))

    n_boot <- nrow(transformed_boot_dat)
    # preconstruct M matrix for the bootstrap data
    M_static <- matrix(NA, nrow = n_boot, ncol = length(theta))
    pos_orthant_mask <- (transformed_boot_dat[,1] > 0) & (transformed_boot_dat[,2] > 0)

    # handle cases of theta = 0, pi/2 separately
    ind_0 <- which(theta == 0)
    if (length(ind_0) > 0) {
        M_static[, ind_0] <- transformed_boot_dat[,1] * pos_orthant_mask
    }
    ind_90 <- which(theta == pi/2)
    if (length(ind_90) > 0) {
        M_static[, ind_90] <- transformed_boot_dat[,2] * pos_orthant_mask
    }
    # handle all other angles now
    inds_no_ax <- which(!(theta == 0 | theta == pi/2))
    if (length(inds_no_ax) > 0) {
        angles_no_ax <- theta[inds_no_ax]
        inv_cos <- 1/cos(angles_no_ax)
        inv_sin <- 1/sin(angles_no_ax)
        
        proj_x <- transformed_boot_dat[,1] %o% inv_cos
        proj_y <- transformed_boot_dat[,2] %o% inv_sin
        M_static[, inds_no_ax] <- pmin(proj_x, proj_y)
    }

    # for each point (row), is it strictly northeast of each point on isoline (column)
    match_x <- outer(transformed_boot_dat[,1], ext_isoline[,1], '>')
    match_y <- outer(transformed_boot_dat[,2], ext_isoline[,2], '>')
    hit_static <- (match_x & match_y)*1

    qn <- (n_boot)^(-gamma)
    
    if (verbose) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (k in 1:B) {
        idx <- sample.int(n_boot, n_boot, replace=TRUE)
        # sample rows of M_boot and hit_boot in lieu of the actual data
        M_boot <- M_static[idx, , drop=FALSE]
        hit_boot <- hit_static[idx, , drop=FALSE]
        # fraction of strict exceedances over each point on isoline
        emp_probs <- colMeans(hit_boot)
        # radii for which 1-qn strict exceedances for each ray on isoline
        rq <- colQuantiles(M_boot, probs=1-qn, type=1, drop=TRUE)
        tail_probs <- (rq/rp)*qn
        final_probs <- ifelse(emp_probs >= qn, emp_probs, tail_probs)
        boot_draws_Wplus[k] <- max(final_probs-p)
        boot_draws_Wminus[k] <- min(final_probs-p)
        if (verbose) utils::setTxtProgressBar(pb, k)
    }
    if (verbose) close(pb)

    c_plus_estimates <- list()
    c_minus_estimates <- list()
    for (i in 1:length(alphas)) {
        alpha <- alphas[i]
        c_plus_estimate <- as.numeric(quantile(boot_draws_Wplus, probs = 1-(alpha/2), type=1))
        c_minus_estimate <- as.numeric(quantile(boot_draws_Wminus, probs=(alpha/2), type=1))
        c_plus_estimates[as.character(alpha)] <- c_plus_estimate
        c_minus_estimates[as.character(alpha)] <- c_minus_estimate
    }

    # store the vectorized survival function computed on your dataset
    # to be used to easily compute coverage of a desired isoline
    survFunc <- function(pts) {
        survProb <- vectorizedBlendedSurv(pts, transformed_survfunc_dat, gamma, xi)
        return(survProb)
    }

    res_lst <- list()
    res_lst$dat <- survfunc_dat
    res_lst$full_dat <- group_assignments
    res_lst$c_plus_estimates <- c_plus_estimates
    res_lst$c_minus_estimates <- c_minus_estimates
    res_lst$p <- p
    res_lst$gamma <- gamma
    res_lst$xi <- xi
    res_lst$survFunc <- survFunc
    res_lst$transform_func1 <- function(pts) transform(pts, dat[,1])
    res_lst$transform_func2 <- function(pts) transform(pts, dat[,2])
    res_lst$inv_transform_func1 <- function(pts) inv_transform(pts, dat[,1])
    res_lst$inv_transform_func2 <- function(pts) inv_transform(pts, dat[,2])

    return(res_lst)

}

computeExtremeRegionThreewayIsoBootSurv <- function(dat, alphas, p, B, gamma, xi, lbs=c(0,0), verbose=FALSE) {
    # Function that constructs confidence tube(s) given a particular dataset and a desired isoline
    # exceedance probability p. Multiple confidence tubes will be returned if multiple alphas are
    # supplied, one for each alpha. NOTE: this function computes c_hat, an estimate of the worst-case difference
    # between the estimated survival function and the true survival function over the set of points corresponding
    # to the desired true isoline. This quantity + the dataset defines the confidence tube: to actually draw it
    # visually, use the drawExtremeRegion function.
    # 
    # HOW THIS DIFFERS FROM OTHER FUNCTIONS: we split the sample in thirds: one for estimating the isoline, one for
    # bootstrap sampling, and one for constructing the survival function. First two used for chat, last used for survfunc
    # resampling
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

    # split sample into three groups: one for transform, one for isoline estimation, one for bootstrapping
    group_assignments <- dat %>%
        mutate(group_assignment = sample(rep_len(1:3, length.out = n())))

    survfunc_dat <- group_assignments %>% filter(group_assignment==1) %>% select(-group_assignment)
    iso_dat <- group_assignments %>% filter(group_assignment==2) %>% select(-group_assignment)
    boot_dat <- group_assignments %>% filter(group_assignment==3) %>% select(-group_assignment)

    # define the transforms here
    # the full dataset is transformed using the dataset itself..
    transform <- function(pts, marginal_dat) {
        transformed_pts <- 1/(1-est_cdf(pts, marginal_dat, gamma)) - 1
        return(transformed_pts)
    }
    transformed_dat_X1 <- transform(dat[,1], dat[,1])
    transformed_dat_X2 <- transform(dat[,2], dat[,2])
    transformed_dat <- data.frame(X1=transformed_dat_X1, X2=transformed_dat_X2)

    inv_transform <- function(pts, marginal_dat) {
        orig_pts <- est_inv_cdf(pts/(1+pts), dat=marginal_dat, gamma=gamma)
        return(orig_pts)
    }

    # transform marginals of isoline, bootstrap, and survfunc datasets
    transformed_iso_dat <- data.frame(X1=transform(iso_dat[,1], dat[,1]), 
                                      X2=transform(iso_dat[,2], dat[,2]))
    transformed_boot_dat <- data.frame(X1=transform(boot_dat[,1], dat[,1]), 
                                      X2=transform(boot_dat[,2], dat[,2]))
    transformed_survfunc_dat <- data.frame(X1=transform(survfunc_dat[,1], dat[,1]), 
                                      X2=transform(survfunc_dat[,2], dat[,2]))
    
    ext_isoline <- drawExtremeIsoline(dat=transformed_iso_dat, 
                                      p=p,
                                      n_coords=200, 
                                      grid_lbs=lbs, 
                                      gamma=gamma, 
                                      xi=1)

    boot_draws <- rep(NA, B)

    theta <- atan2(ext_isoline[,2], ext_isoline[,1])
    rp <- sqrt(rowSums(ext_isoline^2))

    n_boot <- nrow(transformed_boot_dat)
    # preconstruct M matrix for the bootstrap data
    M_static <- matrix(NA, nrow = n_boot, ncol = length(theta))
    pos_orthant_mask <- (transformed_boot_dat[,1] > 0) & (transformed_boot_dat[,2] > 0)

    # handle cases of theta = 0, pi/2 separately
    ind_0 <- which(theta == 0)
    if (length(ind_0) > 0) {
        M_static[, ind_0] <- transformed_boot_dat[,1] * pos_orthant_mask
    }
    ind_90 <- which(theta == pi/2)
    if (length(ind_90) > 0) {
        M_static[, ind_90] <- transformed_boot_dat[,2] * pos_orthant_mask
    }
    # handle all other angles now
    inds_no_ax <- which(!(theta == 0 | theta == pi/2))
    if (length(inds_no_ax) > 0) {
        angles_no_ax <- theta[inds_no_ax]
        inv_cos <- 1/cos(angles_no_ax)
        inv_sin <- 1/sin(angles_no_ax)
        
        proj_x <- transformed_boot_dat[,1] %o% inv_cos
        proj_y <- transformed_boot_dat[,2] %o% inv_sin
        M_static[, inds_no_ax] <- pmin(proj_x, proj_y)
    }

    # for each point (row), is it strictly northeast of each point on isoline (column)
    match_x <- outer(transformed_boot_dat[,1], ext_isoline[,1], '>')
    match_y <- outer(transformed_boot_dat[,2], ext_isoline[,2], '>')
    hit_static <- (match_x & match_y)*1

    qn <- (n_boot)^(-gamma)
    
    if (verbose) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (k in 1:B) {
        idx <- sample.int(n_boot, n_boot, replace=TRUE)
        # sample rows of M_boot and hit_boot in lieu of the actual data
        M_boot <- M_static[idx, , drop=FALSE]
        hit_boot <- hit_static[idx, , drop=FALSE]
        # fraction of strict exceedances over each point on isoline
        emp_probs <- colMeans(hit_boot)
        # radii for which 1-qn strict exceedances for each ray on isoline
        rq <- colQuantiles(M_boot, probs=1-qn, type=1, drop=TRUE)
        tail_probs <- (rq/rp)*qn
        final_probs <- ifelse(emp_probs >= qn, emp_probs, tail_probs)
        boot_draws[k] <- max(abs(final_probs - p))
        if (verbose) utils::setTxtProgressBar(pb, k)
    }
    if (verbose) close(pb)

    c_estimates <- list()
    for (i in 1:length(alphas)) {
        alpha <- alphas[i]
        c_estimate <- as.numeric(quantile(boot_draws, probs = 1-alpha))
        c_estimates[as.character(alpha)] <- c_estimate
    }

    # store the vectorized survival function computed on your dataset
    # to be used to easily compute coverage of a desired isoline
    survFunc <- function(pts) {
        survProb <- vectorizedBlendedSurv(pts, transformed_survfunc_dat, gamma, xi)
        return(survProb)
    }

    res_lst <- list()
    res_lst$dat <- dat
    res_lst$full_dat <- group_assignments
    res_lst$c_estimates <- c_estimates
    res_lst$p <- p
    res_lst$gamma <- gamma
    res_lst$xi <- xi
    res_lst$survFunc <- survFunc
    res_lst$transform_func1 <- function(pts) transform(pts, dat[,1])
    res_lst$transform_func2 <- function(pts) transform(pts, dat[,2])
    res_lst$inv_transform_func1 <- function(pts) inv_transform(pts, dat[,1])
    res_lst$inv_transform_func2 <- function(pts) inv_transform(pts, dat[,2])

    return(res_lst)

}

computeExtremeRegionThreewaySplit <- function(dat, alphas, p, B, gamma, xi, lbs=c(0,0), verbose=FALSE) {
    # Function that constructs confidence tube(s) given a particular dataset and a desired isoline
    # exceedance probability p. Multiple confidence tubes will be returned if multiple alphas are
    # supplied, one for each alpha. NOTE: this function computes c_hat, an estimate of the worst-case difference
    # between the estimated survival function and the true survival function over the set of points corresponding
    # to the desired true isoline. This quantity + the dataset defines the confidence tube: to actually draw it
    # visually, use the drawExtremeRegion function.
    # 
    # HOW THIS DIFFERS FROM OTHER FUNCTIONS: we split the sample in thirds: one for estimating marginal distributions for the marginal transformation step, one for estimating the isoline, and one for generating boostrap draws
    # resampling
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

    # split sample into three groups: one for transform, one for isoline estimation, one for bootstrapping
    group_assignments <- dat %>%
        mutate(group_assignment = sample(rep_len(1:3, length.out = n())))

    trans_dat <- group_assignments %>% filter(group_assignment==1) %>% select(-group_assignment)
    iso_dat <- group_assignments %>% filter(group_assignment==2) %>% select(-group_assignment)
    boot_dat <- group_assignments %>% filter(group_assignment==3) %>% select(-group_assignment)

    # the transform: inverse probability transform...
    # this is no longer using the GPD
    transform <- function(pts, marginal_dat) {
        edf <- ecdf(marginal_dat)
        n_marg <- length(marginal_dat)
        uniform_dat <- (n_marg/(n_marg+1))*edf(pts) # avoid boundary issues with plugging 1 into a quantile function
        transformed_pts <- 1/(1-uniform_dat) - 1
        return(transformed_pts)
    }

    # inverse transform, using just (n/(n+1))*ecdf as fitted marginal
    inv_transform <- function(pts, marginal_dat) {
        n_marg <- length(marginal_dat)
        unif_scale <- (1-(1/(pts+1)))*((n_marg+1)/n_marg)
        orig_pts <- quantile(marginal_dat, probs=unif_scale, type=2, names=FALSE)
        return(orig_pts)
    }

    # transform isoline drawing and bootstrapping marginals, using trans_dat to fit marginals
    transformed_iso_dat <- data.frame(X1=transform(iso_dat[,1], trans_dat[,1]), 
                                      X2=transform(iso_dat[,2], trans_dat[,2]))
    transformed_boot_dat <- data.frame(X1=transform(boot_dat[,1], trans_dat[,1]), 
                                      X2=transform(boot_dat[,2], trans_dat[,2]))
    
    ext_isoline <- drawExtremeIsoline(dat=transformed_iso_dat, 
                                      p=p,
                                      n_coords=200, 
                                      grid_lbs=lbs, 
                                      gamma=gamma, 
                                      xi=1)

    boot_draws <- rep(NA, B)

    theta <- atan2(ext_isoline[,2], ext_isoline[,1])
    rp <- sqrt(rowSums(ext_isoline^2))

    n_boot <- nrow(transformed_boot_dat)
    # preconstruct M matrix for the bootstrap data
    M_static <- matrix(NA, nrow = n_boot, ncol = length(theta))
    pos_orthant_mask <- (transformed_boot_dat[,1] > 0) & (transformed_boot_dat[,2] > 0)

    # handle cases of theta = 0, pi/2 separately
    ind_0 <- which(theta == 0)
    if (length(ind_0) > 0) {
        M_static[, ind_0] <- transformed_boot_dat[,1] * pos_orthant_mask
    }
    ind_90 <- which(theta == pi/2)
    if (length(ind_90) > 0) {
        M_static[, ind_90] <- transformed_boot_dat[,2] * pos_orthant_mask
    }
    # handle all other angles now
    inds_no_ax <- which(!(theta == 0 | theta == pi/2))
    if (length(inds_no_ax) > 0) {
        angles_no_ax <- theta[inds_no_ax]
        inv_cos <- 1/cos(angles_no_ax)
        inv_sin <- 1/sin(angles_no_ax)
        
        proj_x <- transformed_boot_dat[,1] %o% inv_cos
        proj_y <- transformed_boot_dat[,2] %o% inv_sin
        M_static[, inds_no_ax] <- pmin(proj_x, proj_y)
    }

    # for each point (row), is it strictly northeast of each point on isoline (column)
    match_x <- outer(transformed_boot_dat[,1], ext_isoline[,1], '>')
    match_y <- outer(transformed_boot_dat[,2], ext_isoline[,2], '>')
    hit_static <- (match_x & match_y)*1

    qn <- (n_boot)^(-gamma)
    
    if (verbose) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (k in 1:B) {
        idx <- sample.int(n_boot, n_boot, replace=TRUE)
        # sample rows of M_boot and hit_boot in lieu of the actual data
        M_boot <- M_static[idx, , drop=FALSE]
        hit_boot <- hit_static[idx, , drop=FALSE]
        # fraction of strict exceedances over each point on isoline
        emp_probs <- colMeans(hit_boot)
        # radii for which 1-qn strict exceedances for each ray on isoline
        rq <- colQuantiles(M_boot, probs=1-qn, type=1, drop=TRUE)
        tail_probs <- (rq/rp)*qn
        final_probs <- ifelse(emp_probs >= qn, emp_probs, tail_probs)
        boot_draws[k] <- max(abs(final_probs - p))
        if (verbose) utils::setTxtProgressBar(pb, k)
    }
    if (verbose) close(pb)

    c_estimates <- list()
    for (i in 1:length(alphas)) {
        alpha <- alphas[i]
        c_estimate <- as.numeric(quantile(boot_draws, probs = 1-alpha))
        c_estimates[as.character(alpha)] <- c_estimate
    }

    # store the vectorized survival function computed on your dataset
    # to be used to easily compute coverage of a desired isoline
    survFunc <- function(pts) {
        survProb <- vectorizedBlendedSurv(pts, transformed_iso_dat, gamma, xi)
        return(survProb)
    }

    res_lst <- list()
    res_lst$full_dat <- group_assignments
    res_lst$c_estimates <- c_estimates
    res_lst$p <- p
    res_lst$gamma <- gamma
    res_lst$xi <- xi
    res_lst$survFunc <- survFunc
    res_lst$transform_func1 <- function(pts) transform(pts, trans_dat[,1])
    res_lst$transform_func2 <- function(pts) transform(pts, trans_dat[,2])
    res_lst$inv_transform_func1 <- function(pts) inv_transform(pts, trans_dat[,1])
    res_lst$inv_transform_func2 <- function(pts) inv_transform(pts, trans_dat[,2])

    return(res_lst)

}

computeExtremeRegionPretransformGPD <- function(dat, alphas, p, B, gamma, verbose=FALSE) {
    # Function to compute an extremal confidence region, but as opposed to the method above,
    # we input the data on it's original scale, perform the marginal transformation in here,
    # and then re-fit the marginals within the bootstrap. The idea is that we are not accounting
    # for the variability in fitting the marginals in the above function when determining c, so
    # we will refit the marginals in each iteration of the bootstrap loop to account for this variability.
    # 
    # Arguments:
    # dat: the data, in the form of a 2-column R data.frame
    # alphas: a vector of alphas indicating desired probability of miscoverage
    # (0.01 for 99% CI, 0.05 for 95% CI, etc.). NOTE: even if only one alpha desired, must put it in a vector.
    # p: the desired p-isoline you wish to capture
    # B: the number of bootstrap replicates in determining c_hat
    # gamma: 1/n^(gamma), controls the smallest probability to use empirical df
    # verbose: print bootstrap loop progress bar and checkpoint messages

    # the transform: inverse probability transform, with CDF estimated via ecdf + GPD beyond n^-gamma quantile
    transform <- function(pts, marginal_dat) {
        transformed_pts <- 1/(1-est_cdf(pts, marginal_dat, gamma)) - 1
        return(transformed_pts)
    }

    # the inverse transform: not explicitly used in the following computations, but useful for
    # drawing on the original scale
    inv_transform <- function(pts, marginal_dat) {
        orig_pts <- est_inv_cdf(pts/(1+pts), dat=marginal_dat, gamma=gamma)
        return(orig_pts)
    }
    
    transformed_dat <- data.frame(X1=transform(dat[,1], dat[,1]), X2=transform(dat[,2], dat[,2]))
    n_dat <- nrow(transformed_dat)
    
    ext_isoline <- drawExtremeIsoline(transformed_dat, p, n_coords=200, grid_lbs=c(0,0), gamma=gamma, xi=1)
    boot_draws <- rep(NA, B)

    # bootstrap loop
    if (verbose) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (k in 1:B) {
        idx <- sample.int(n_dat, n_dat, replace=TRUE)
        
        boot_dat <- dat[idx,]
        # transform bootstrapped dataset
        boot_dat_transformed <- data.frame(X1=transform(boot_dat[,1], boot_dat[,1]), 
                                           X2=transform(boot_dat[,2], boot_dat[,2]))
        # compute the survival function using the transformed dataset, on your isoline of interest
        final_probs <- vectorizedBlendedSurv(pts=ext_isoline, dat=boot_dat_transformed, gamma=gamma, xi=1)
        boot_draws[k] <- max(abs(final_probs - p))
        if (verbose) utils::setTxtProgressBar(pb, k)
    }
    if (verbose) close(pb)

    c_estimates <- list()
    for (i in 1:length(alphas)) {
        alpha <- alphas[i]       
        c_estimate <- as.numeric(quantile(boot_draws, probs = 1-alpha))
        c_estimates[as.character(alpha)] <- c_estimate
    }

    # store the vectorized survival function computed on your dataset
    # to be used to easily compute coverage of a desired isoline
    survFunc <- function(pts) {
        survProb <- vectorizedBlendedSurv(pts, transformed_dat, gamma=gamma, xi=1)
        return(survProb)
    }

    res_lst <- list()
    res_lst$dat <- dat
    res_lst$trans_dat <- transformed_dat
    res_lst$c_estimates <- c_estimates
    res_lst$p <- p
    res_lst$gamma <- gamma
    res_lst$xi <- 1
    res_lst$survFunc <- survFunc
    res_lst$transform_func1 <- function(pts) transform(pts, dat[,1])
    res_lst$transform_func2 <- function(pts) transform(pts, dat[,2])
    res_lst$inv_transform_func1 <- function(pts) inv_transform(pts, dat[,1])
    res_lst$inv_transform_func2 <- function(pts) inv_transform(pts, dat[,2])

    return(res_lst)

}

computeExtremeRegionPretransform <- function(dat, alphas, p, B, gamma, verbose=FALSE) {
    
    n_dat <- nrow(dat)
    
    transform_gpd <- function(pts, marginal_dat) {
        probs <- est_cdf(pts, marginal_dat, gamma)
        transformed_pts <- (1 / (1 - probs)) - 1
        return(transformed_pts)
    }

    inv_transform_gpd <- function(pts, marginal_dat) {
        u_pts <- pts / (pts + 1)
        orig_pts <- est_inv_cdf(u_pts, dat=marginal_dat, gamma=gamma)
        return(orig_pts)
    }

    # draw the best estimate of the isoline using the hybrid ecdf/gpd transform
    transformed_dat <- data.frame(X1=transform_gpd(dat[,1], dat[,1]), 
                                  X2=transform_gpd(dat[,2], dat[,2]))
    ext_isoline <- drawExtremeIsoline(transformed_dat, p, n_coords=200, 
                                      grid_lbs=c(0,0), gamma=gamma, xi=1)
    
    boot_draws <- rep(NA, B)

    if (verbose) pb <- utils::txtProgressBar(min=0, max=B, style=3)    
    for (k in 1:B) {
        idx <- sample.int(n_dat, n_dat, replace=TRUE)
        boot_dat <- dat[idx,]   
        
        u1 <- rank(boot_dat[,1], ties.method = "average")/(n_dat+1)
        u2 <- rank(boot_dat[,2], ties.method = "average")/(n_dat+1)
    
        boot_dat_transformed <- cbind(
            (1/(1-u1))-1, 
            (1/(1-u2))-1
        )

        final_probs <- vectorizedBlendedSurv(pts=ext_isoline, 
                                             dat=boot_dat_transformed, 
                                             gamma=gamma, 
                                             xi=1)
        
        boot_draws[k] <- max(abs(final_probs - p))
        
        if (verbose) utils::setTxtProgressBar(pb, k)
    }
    if (verbose) close(pb)

    c_estimates <- list()
    for (i in 1:length(alphas)) {
        alpha <- alphas[i]        
        c_estimate <- as.numeric(quantile(boot_draws, probs = 1-alpha))
        c_estimates[as.character(alpha)] <- c_estimate
    }

    # store the survival function for later usage (evaluating coverage of a known isoline, etc.)
    survFunc <- function(pts) {
        survProb <- vectorizedBlendedSurv(pts, transformed_dat, gamma=gamma, xi=1)
        return(survProb)
    }

    res_lst <- list()
    res_lst$dat <- dat
    res_lst$trans_dat <- transformed_dat
    res_lst$c_estimates <- c_estimates
    res_lst$p <- p
    res_lst$gamma <- gamma
    res_lst$xi <- 1
    res_lst$survFunc <- survFunc
    
    # save the transforms, to be used to evaluate coverage of isolines + plot bounds of tubes
    res_lst$transform_func1 <- function(pts) transform_gpd(pts, dat[,1])
    res_lst$transform_func2 <- function(pts) transform_gpd(pts, dat[,2])
    res_lst$inv_transform_func1 <- function(pts) inv_transform_gpd(pts, dat[,1])
    res_lst$inv_transform_func2 <- function(pts) inv_transform_gpd(pts, dat[,2])

    return(res_lst)
}

computeEmpiricalRegionSymmetric <- function(dat, alphas, p, B, lbs, verbose=FALSE) {
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

    emp_isoline <- drawEmpiricalIsoline(dat=dat, n_coords=200, grid_lbs=lbs, p)
    boot_draws <- rep(NA, B)

    if (verbose) pb <- utils::txtProgressBar(min=0, max=B, style=3) 
    for (k in 1:B) {
        boot_dat <- dat %>% sample_frac(1, replace = TRUE)
        boot_survfunc <- vectorizedEmpSurv(emp_isoline, boot_dat)
        boot_draws[k] <- max(abs(boot_survfunc-p))
        if (verbose) utils::setTxtProgressBar(pb, k)
        }
    if (verbose) close(pb)

    c_estimates <- list()
    for (i in 1:length(alphas)) {
        alpha <- alphas[i]       
        c_estimate <- as.numeric(quantile(boot_draws, probs = 1-alpha))
        c_estimates[as.character(alpha)] <- c_estimate
    }

    survFunc <- function(pts) {
        survProb <- vectorizedEmpSurv(pts, dat)
        return(survProb)
    }

    res_lst <- list()
    res_lst$dat <- dat
    res_lst$c_estimates <- c_estimates
    res_lst$p <- p
    res_lst$survFunc <- survFunc

    return(res_lst)
}

computeEmpiricalRegion <- function(dat, alphas, p, B, lbs, verbose=FALSE) {
    # Function that constructs confidence tube(s) for a p-isoline of interest, given some data.
    # Will return multiple tube constructions if multiple alphas are passed as an argument (NOTE:
    # even if only one alpha is desired, you still must put it in an R vector). This function
    # uses the empirical survival function, and thus is more appropriate for use with non-extreme
    # isolines. As above, this function only computes the quantities needed to define
    # a confidence tube (i.e., cminus and cplus). 
    # it, use the drawEmpiricalRegion function.
    # NOTE: this function differs from the previous version in that we compute both
    # a separate lower and upper vertical bound, instead of using the same one.
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

    emp_isoline <- drawEmpiricalIsoline(dat=dat, n_coords=200, grid_lbs=lbs, p)
    boot_draws_Wplus <- rep(NA, B)
    boot_draws_Wminus <- rep(NA, B)

    if (verbose) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (k in 1:B) {
        boot_dat <- dat %>% sample_frac(1, replace = TRUE)
        boot_survfunc <- vectorizedEmpSurv(emp_isoline, boot_dat)
        boot_draws_Wplus[k] <- max(boot_survfunc-p)
        boot_draws_Wminus[k] <- min(boot_survfunc-p)
        if (verbose) utils::setTxtProgressBar(pb, k)
    }
    if (verbose) close(pb)
    
    c_plus_estimates <- list()
    c_minus_estimates <- list()
    for (i in 1:length(alphas)) {
        alpha <- alphas[i]       
        c_plus_estimate <- as.numeric(quantile(boot_draws_Wplus, probs=1-(alpha/2), type=1))
        c_minus_estimate <- as.numeric(quantile(boot_draws_Wminus, probs=alpha/2, type=1))
        c_plus_estimates[as.character(alpha)] <- c_plus_estimate
        c_minus_estimates[as.character(alpha)] <- c_minus_estimate
    }

    survFunc <- function(pts) {
        survProb <- vectorizedEmpSurv(pts, dat)
        return(survProb)
    }

    res_lst <- list()
    res_lst$dat <- dat
    res_lst$c_plus_estimates <- c_plus_estimates
    res_lst$c_minus_estimates <- c_minus_estimates
    res_lst$p <- p
    res_lst$survFunc <- survFunc

    return(res_lst)
}