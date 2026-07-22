source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')

computeExtremeRegionBC_GivenBounds <- function(dat, alphas, p, B, gamma, b1, b2, ncoords=100, n_lines=10, xi=1, lbs=c(0,0), 
                                               sup_region=TRUE, transformed=TRUE, progress_bar=FALSE, verbose=FALSE) {
    # function that constructs confidence tube(s) given a particular dataset and bounds b1 and b2.
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
        construction_dat <- data.frame(X1=transformed_dat_X1, X2=transformed_dat_X2)

    } else {
        construction_dat <- dat
    }

    # Construct the supremum region B using the provided b1 and b2 bounds
    if (sup_region) {
        if (verbose) print('Constructing region using provided true expected bounds.')
        
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

    # Set up hit matrix for the bootstrap loop over the constructed region
    theta <- atan2(ext_region_pts[,2], ext_region_pts[,1])
    rp <- sqrt(rowSums(ext_region_pts^2))
    
    n_dat <- nrow(construction_dat)
    
    # SAFELY construct Z_static avoiding NaN generation on the axes
    Z_static <- matrix(NA, nrow = n_dat, ncol = length(theta))
    pos_orthant_mask <- (construction_dat[,1] > 0) & (construction_dat[,2] > 0)
    
    ind_0 <- which(theta == 0)
    if (length(ind_0) > 0) Z_static[, ind_0] <- construction_dat[,1] * pos_orthant_mask
    
    ind_90 <- which(theta == pi/2)
    if (length(ind_90) > 0) Z_static[, ind_90] <- construction_dat[,2] * pos_orthant_mask
    
    inds_no_ax <- which(!(theta == 0 | theta == pi/2))
    if (length(inds_no_ax) > 0) {
        Z_static[, inds_no_ax] <- pmin(
            construction_dat[,1] %o% (1/cos(theta[inds_no_ax])), 
            construction_dat[,2] %o% (1/sin(theta[inds_no_ax]))
        )
    }

    hit_static <- (outer(construction_dat[,1], ext_region_pts[,1], '>') & 
                   outer(construction_dat[,2], ext_region_pts[,2], '>')) * 1.0

    qn <- n_dat^(-gamma)
    
    boot_draws_Wplus <- rep(NA, B)
    boot_draws_Wminus <- rep(NA, B)

    if (progress_bar) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (k in 1:B) {
        idx <- sample.int(n_dat, n_dat, replace=TRUE)
        
        # Tabulate trick for massive memory savings
        w <- tabulate(idx, nbins = n_dat)
        emp_probs <- as.numeric(w %*% hit_static) / n_dat
        
        Z_boot <- Z_static[idx, , drop=FALSE]
        rq <- colQuantiles(Z_boot, probs=1-qn, type=1, drop=TRUE)

        rm(Z_boot) # explicit garbage collection
        
        tail_probs <- ((rq/rp)^(1/xi))*qn
        final_probs <- ifelse(emp_probs >= qn, emp_probs, tail_probs)
        
        # Taking the max/min over the entire defined region minus the base p
        boot_draws_Wplus[k] <- max(final_probs - p)
        boot_draws_Wminus[k] <- min(final_probs - p)
        
        if (progress_bar) utils::setTxtProgressBar(pb, k)
    }
    if (progress_bar) close(pb)
    
    rm(Z_static, hit_static) # Final cleanup

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

computeExtremeRegionBCThreeway_GivenBounds <- function(dat, alphas, p, B, gamma, b1, b2, ncoords=100, n_lines=10, xi=1, lbs=c(0,0), sup_region=TRUE, transformed=TRUE, verbose=FALSE) {
    
    group_assignments <- dat %>%
        mutate(group_assignment = sample(rep_len(1:3, length.out = n())))

    survfunc_dat <- group_assignments %>% filter(group_assignment==1) %>% select(-group_assignment)
    iso_dat <- group_assignments %>% filter(group_assignment==2) %>% select(-group_assignment)
    boot_dat <- group_assignments %>% filter(group_assignment==3) %>% select(-group_assignment)

    apply_transform <- function(pts, marginal_dat) {
        1/(1-est_cdf(pts, marginal_dat, gamma)) - 1
    }
    inv_transform <- function(pts, marginal_dat) {
        est_inv_cdf(pts/(1+pts), dat=marginal_dat, gamma=gamma)
    }

    transformed_iso_dat <- data.frame(X1=apply_transform(iso_dat[,1], dat[,1]), 
                                      X2=apply_transform(iso_dat[,2], dat[,2]))
    transformed_boot_dat <- data.frame(X1=apply_transform(boot_dat[,1], dat[,1]), 
                                       X2=apply_transform(boot_dat[,2], dat[,2]))
    transformed_survfunc_dat <- data.frame(X1=apply_transform(survfunc_dat[,1], dat[,1]), 
                                           X2=apply_transform(survfunc_dat[,2], dat[,2]))
    
    if (sup_region) {
        p_seq <- seq(b1, b2, length.out = n_lines)
        region_pts_list <- lapply(p_seq, function(prob) {
            drawExtremeIsoline(transformed_iso_dat, prob, n_coords=ncoords, grid_lbs=lbs, gamma, xi)
        })
        ext_region_pts <- do.call(rbind, region_pts_list)
    } else {
        ext_region_pts <- drawExtremeIsoline(transformed_iso_dat, p, n_coords=ncoords, grid_lbs=lbs, gamma, xi)
    }

    theta <- atan2(ext_region_pts[,2], ext_region_pts[,1])
    rp <- sqrt(rowSums(ext_region_pts^2))
    n_boot <- nrow(transformed_boot_dat)

    # projection matrix setup
    M_static <- matrix(NA, nrow = n_boot, ncol = length(theta))
    pos_orthant_mask <- (transformed_boot_dat[,1] > 0) & (transformed_boot_dat[,2] > 0)
    
    ind_0 <- which(theta == 0)
    if (length(ind_0) > 0) M_static[, ind_0] <- transformed_boot_dat[,1] * pos_orthant_mask
    
    ind_90 <- which(theta == pi/2)
    if (length(ind_90) > 0) M_static[, ind_90] <- transformed_boot_dat[,2] * pos_orthant_mask
    
    inds_no_ax <- which(!(theta == 0 | theta == pi/2))
    if (length(inds_no_ax) > 0) {
        M_static[, inds_no_ax] <- pmin(
            transformed_boot_dat[,1] %o% (1/cos(theta[inds_no_ax])), 
            transformed_boot_dat[,2] %o% (1/sin(theta[inds_no_ax]))
        )
    }

    hit_static <- (outer(transformed_boot_dat[,1], ext_region_pts[,1], '>') & 
                   outer(transformed_boot_dat[,2], ext_region_pts[,2], '>')) * 1.0

    qn <- n_boot^(-gamma)
    boot_draws_Wplus <- rep(NA, B)
    boot_draws_Wminus <- rep(NA, B)

    if (verbose) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (k in 1:B) {
        idx <- sample.int(n_boot, n_boot, replace=TRUE)
        
        w <- tabulate(idx, nbins = n_boot)
        emp_probs <- as.numeric(w %*% hit_static) / n_boot
        
        M_boot <- M_static[idx, , drop=FALSE]
        rq <- colQuantiles(M_boot, probs=1-qn, type=1, drop=TRUE)
        rm(M_boot)
        
        tail_probs <- (rq/rp)*qn
        final_probs <- ifelse(emp_probs >= qn, emp_probs, tail_probs)
        
        boot_draws_Wplus[k] <- max(final_probs-p)
        boot_draws_Wminus[k] <- min(final_probs-p)
        if (verbose) utils::setTxtProgressBar(pb, k)
    }
    if (verbose) close(pb)
    
    rm(M_static, hit_static)

    c_plus_estimates <- list()
    c_minus_estimates <- list()
    for (i in 1:length(alphas)) {
        alpha <- alphas[i]
        c_plus_estimates[as.character(alpha)] <- as.numeric(quantile(boot_draws_Wplus, probs = 1-(alpha/2), type=1))
        c_minus_estimates[as.character(alpha)] <- as.numeric(quantile(boot_draws_Wminus, probs=(alpha/2), type=1))
    }

    survFunc <- function(pts) {
        survProb <- vectorizedBlendedSurv(pts, transformed_survfunc_dat, gamma, xi)
        return(survProb)
    }

    res_lst <- list(
        dat = survfunc_dat,
        full_dat = group_assignments,
        c_plus_estimates = c_plus_estimates,
        c_minus_estimates = c_minus_estimates,
        p = p,
        gamma = gamma,
        xi = xi,
        survFunc = survFunc,
        transformed = transformed,
        transform_func1 = function(pts) apply_transform(pts, dat[,1]),
        transform_func2 = function(pts) apply_transform(pts, dat[,2]),
        inv_transform_func1 = function(pts) inv_transform(pts, dat[,1]),
        inv_transform_func2 = function(pts) inv_transform(pts, dat[,2]),
        b1 = b1, 
        b2 = b2
    )

    return(res_lst)
}