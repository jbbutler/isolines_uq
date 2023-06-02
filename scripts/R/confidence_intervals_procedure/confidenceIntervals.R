# R script containing functions to make bootstrap confidence intervals for given datasets

library(bootstrap)

source('/global/u1/j/jbbutler/isolines_uq/scripts/R/orig_isolines.R')

interpEstCoord <- function(iso_dat, fixed_coord, ax) {

    less_mask <- iso_dat[, ax] < fixed_coord
    greater_mask <- iso_dat[, ax] > fixed_coord

    if (any(less_mask) & any(greater_mask)) {

        less_vals <- iso_dat[less_mask,]
        greater_vals <- iso_dat[greater_mask,]
        
	low_mask <- which.min(abs(less_vals[, ax] - fixed_coord))
	greater_mask <- which.min(abs(greater_vals[, ax] - fixed_coord))
        xl <- less_vals[,ax][low_mask]
	yl <- less_vals[,-ax][low_mask]

	xh <- greater_vals[, ax][greater_mask]
        yh <- greater_vals[, -ax][greater_mask]

        coord <- yl + ((yh - yl)/(xh - xl))*(fixed_coord - xl)

	return(coord)

    }
    return(NA)
}

percentileMethod <- function(bootstrap_samps, alpha, isoline_num, fixed_coord, ax) {
    
    n_boots <- length(bootstrap_samps)
    coords <- rep(NA, n_boots)

    for (k in 1:length(bootstrap_samps)) {

	iso_dat <- data.frame(bootstrap_samps[[k]][[1]][[isoline_num]])
        coords[k] <- interpEstCoord(iso_dat, fixed_coord, ax)
    }

    cI <- as.numeric(quantile(coords, probs = c(alpha/2, 1-(alpha/2)), na.rm = TRUE))
    percentile_out <- list()
    percentile_out$cI <- cI
    percentile_out$boot_samps <- coords
    return(percentile_out)
}

# Helper Function to do jackknife estimation "by-hand"
jackknifeHelper <- function(entry, dat, isoline_num, asympIndep, fixed_coord, ax) {

    jack_dat <- dat[-entry,]
    boot_out_full <- suppressWarnings(xContours(dat = jack_dat, faster = TRUE, asympIndep = asympIndep))
    isolines <- c(list(boot_out_full$contourOrig), boot_out_full$projContours)
    iso_dat <- as.matrix(isolines[[isoline_num]])

    coord <- interpEstCoord(iso_dat, fixed_coord, ax)

    return(coord)

}

# Function to do the jackknife estimation "by-hand"
jackknifeAcceleration_alt <- function(dat, isoline_num, asympIndep, fixed_coord, ax, n_cores) {

    clust <- makeCluster(n_cores)
    registerDoParallel(cl = clust)

    coords <- foreach (i = 1:nrow(dat), .packages = c('mvtnorm', 'dplyr', 'MASS', 'ismev', 'evd', 'ks')) %dopar% {

	jackknifeHelper(i, dat, isoline_num, asympIndep, fixed_coord, ax)

    }

    stopCluster(clust)

    mcoords <- mean(unlist(coords))

    return(sum((mcoords - coords)**3)/(6*(sum((mcoords - coords)**2)**(3/2))))

}

jackknifeAcceleration_setup <- function(inds, dat, isoline_num, asympIndep, fixed_coord, ax) {

    jack_dat <- cbind(dat[inds,1], dat[inds, 2])
    boot_out_full <- suppressWarnings(xContours(dat = jack_dat, faster = TRUE, asympIndep = asympIndep, projContourLevels = c(0.001)))
    isolines <- c(list(boot_out_full$contourOrig), boot_out_full$projContours)
    iso_dat <- as.matrix(isolines[[isoline_num]])

    coord <- interpEstCoord(iso_dat, fixed_coord, ax)

    return(coord)

}



biasCorrection <- function(bootstrap_ests, orig_est) {
    return(qnorm(sum(bootstrap_ests < orig_est)/length(bootstrap_ests)))
}

bcaMethod <- function(bootstrap_samps, orig_samp, alpha, isoline_num, fixed_coord, ax, asympIndep) {


    len_start <- proc.time()
    n_boots <- length(bootstrap_samps)
    coords <- rep(NA, n_boots)
    len_elapsed <- proc.time() - len_start
    #print(len_elapsed)
    coord_start <- proc.time()
    for (k in 1:n_boots) {
k
        iso_dat <- as.matrix(bootstrap_samps[[k]][[1]][[isoline_num]])
        coords[k] <- interpEstCoord(iso_dat, fixed_coord, ax)
    }
    coord_elapsed <- proc.time() - coord_start
    #print(coord_elapsed)
    na_mask <- is.na(coords)
    num_na_boot <- sum(na_mask)
    coords <- coords[!na_mask]

    orig_dat_out <- suppressWarnings(xContours(dat = orig_samp, faster = TRUE, asympIndep = asympIndep, projContourLevels = c(0.001)))
    isolines <- c(list(orig_dat_out$contourOrig), orig_dat_out$projContours)
    iso_dat <- as.matrix(isolines[[isoline_num]])
    orig_coord <- interpEstCoord(iso_dat, fixed_coord, ax)

    bias_start <- proc.time()
    bc <- biasCorrection(coords, orig_coord)
    bias_elapsed <- proc.time() - bias_start
    #print(bias_elapsed)
    jack_start <- proc.time()
    jack_res <- jackknife(1:nrow(orig_samp), jackknifeAcceleration_setup, 
			  dat = orig_samp, isoline_num = isoline_num, 
			  asympIndep = asympIndep, fixed_coord = fixed_coord, ax = ax)
    jack_elapsed <- proc.time() - jack_start
    #print(jack_elapsed)
    jack_res <- jack_res[[3]]
    na_mask <- is.na(jack_res)
    num_na_jack <- sum(na_mask)
    jack_res <- jack_res[!na_mask]
    mjack_res <- mean(jack_res)

    rest_start <- proc.time()

    ac <- sum((mjack_res - jack_res)**3)/(6*(sum((mjack_res - jack_res)**2)**(3/2)))


    alpha1 <- pnorm(bc + ((bc + qnorm(alpha/2))/(1 - ac*(bc + qnorm(alpha/2)))))
    alpha2 <- pnorm(bc + ((bc + qnorm(1 - alpha/2))/(1 - ac*(bc + qnorm(1 - alpha/2)))))

    bca_cI <- as.numeric(quantile(coords, probs = c(alpha1, alpha2)))
    per_cI <- as.numeric(quantile(coords, probs = c(alpha/2, 1-(alpha/2))))

    rest_elapsed <- proc.time() - rest_start
    #print(rest_elapsed)

    bca_out <- list()
    bca_out$bca_cI <- bca_cI
    bca_out$per_cI <- per_cI
    bca_out$num_na_boot <- num_na_boot
    bca_out$num_na_jack <- num_na_jack
    bca_out$boot_samps <- coords
    bca_out$ac <- ac
    bca_out$bc <- bc

    return(bca_out)

}

parseConfints <- function(out, true_coord, n, B, na_tol) {
    
    error_inds <- c()

    for (i in 1:length(out)) {

	if (length(out[[i]]) <= 2) {
            error_inds <- c(error_inds, i)
	}
    }

    if (length(error_inds) > 0) {
        out <- out[-error_inds]
        print(paste('Indices of excluded pathological simulations:', toString(error_inds)))
     }
    confints_bca <- lapply(seq_along(out), function(i) out[[i]][[1]])
    confints_bca <- as.data.frame(do.call(rbind, confints_bca))
    na_coord <- unlist(lapply(seq_along(out), function(i) out[[i]][[3]]))
    prop_na_coord <- na_coord/B
    na_jack <- unlist(lapply(seq_along(out), function(i) out[[i]][[4]]))
    prop_na_jack <- na_jack/n
    confints_per <- lapply(seq_along(out), function(i) out[[i]][[2]])
    confints_per <- as.data.frame(do.call(rbind, confints_per))

    bca_na <- rowSums(is.na(confints_bca)) > 0
    per_na <- rowSums(is.na(confints_per)) > 0
    confints_bca <- confints_bca[!bca_na, ]
    confints_per <- confints_per[!per_na, ]

    cvg_bca <- sum(confints_bca$V1 < true_coord & confints_bca$V2 > true_coord)/nrow(confints_bca)
    cvg_per <- sum(confints_per$V1 < true_coord & confints_per$V2 > true_coord)/nrow(confints_per)
    avg_bca <- mean(confints_bca$V2 - confints_bca$V1)
    avg_per <- mean(confints_per$V2 - confints_per$V1)

    res <- list()
    res$confints_bca <- confints_bca
    res$cvg_bca <- cvg_bca
    res$avg_bca <- avg_bca
    res$na_bca <- sum(bca_na) 
    res$confints_per <- confints_per
    res$cvg_per <- cvg_per
    res$avg_per <- avg_per
    res$na_per <- sum(per_na)

    res$num_na_coord <- sum(prop_na_coord >= na_tol)
    res$num_na_jack <- sum(prop_na_jack >= na_tol)

    return(res)

}

