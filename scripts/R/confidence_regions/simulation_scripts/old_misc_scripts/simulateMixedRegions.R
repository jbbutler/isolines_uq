# load packages
library(ggplot2)
library(dplyr)
library(mvtnorm)
library(data.table)
library(gridExtra)
library(FNN)
library(foreach)
library(doSNOW)
library(parallel)

# load perlmutter scratch filepath signature
perlpath <- Sys.getenv('PSCRATCH')

# load self-made functions
source('~/isolines_uq/scripts/R/auxiliary_scripts/distributionIsolines.R')
source('~/isolines_uq/scripts/R/auxiliary_scripts/utils.R')


##################### USEFUL HELPER FUNCTIONS (workshop/put in new script later) #######################

# function to find the angle corresponding to a specific point in the cartesian grid
drawEmpiricalIsoline <- function(dat, angles, p) {

    survivalDiff <- function(r, theta, p, dat) {

        xPt <- r*cos(theta)
        yPt <- r*sin(theta)

        actual_survival <- mean((dat[,1] > xPt) & dat[,2] > yPt)

        return(actual_survival - p)
        }

    maxX <- max(dat[,1])
    maxY <- max(dat[,2])
    maxRad <- sqrt(maxX**2 + maxY**2)

    num_radii <- length(angles)

    radii <- rep(NA, num_radii)

    for (i in 1:num_radii) {
        radii[i] <- uniroot(survivalDiff, interval=c(0, maxRad), theta=angles[i], dat=dat, p=p)$root
    }

    xs <- radii*cos(angles)
    ys <- radii*sin(angles)

    return(data.frame(X1=xs, X2=ys, r=radii, theta=angles))
}

# function to evaluate survival function across a grid, incorporating both empirical survfunc and regular variation
# seems like a strategy would be use nearest-neighbor interpolation, a la https://stackoverflow.com/questions/30262434/interpolate-values-from-a-grid-efficiently-in-r
# whether nearest neighbor interpolation makes sense is another story...
# returns the q-isoline as well if you are not interpolating on another grid
# rv_ub is the upper extent of grid for drawing rv survfunc

evaluateSurvivalGrid <- function(dat, q, emp_ticks, num_angles, num_radii, rv_ub, interp_grid=NULL) {

    # first, find the q_isoline
    angles <- seq(0, pi/2, length.out=num_angles)
    q_isoline <- drawEmpiricalIsoline(dat, angles, q)
    # get upper bounds for grid on which you will evaluate empirical survival function
    # add 3 tick marks worth of padding onto it
    ub <- max(q_isoline[1,1], tail(q_isoline, 1)[1,2])
    ub <- ub + 3*(ub/emp_ticks)

    # make the grid and evaluate empirical survival function on it
    # subsetting the points for which exceedance probability is at least q
    grid <- expand.grid(X1=seq(0, ub, length.out=emp_ticks),
                   X2=seq(0, ub, length.out=emp_ticks))

    surv_func <- fastEmpSurv(grid, dat)
    empirical_mask <- surv_func >= q
    empirical_gridpts <- grid[empirical_mask,]
    survfunc <- surv_func[empirical_mask]
    survgrid <- cbind(empirical_gridpts, survfunc)

    # now specify how far out you wish to draw the regularly-varying portion of the survival function
    max_r <- sqrt(rv_ub**2 + rv_ub**2)

    # next, for each angle/intersecting ray of the q_isoline, estimate survfunc using RV for several radii along the ray
    surv_rays <- vector(mode='list', length=length(angles))

    for (i in 1:nrow(q_isoline)) {
        theta <- q_isoline$theta[i]
        starting_r <- q_isoline$r[i]
        ray_radii <- seq(starting_r, max_r, length.out=num_radii+1)[-1]
        ray_probs <- (starting_r/ray_radii)*q
        surv_rays[[i]] <- data.frame(X1=ray_radii*cos(theta), X2=ray_radii*sin(theta), survfunc=ray_probs)

    }

    # concatenate the empirical and RV parts into single dataframe
    full_survgrid <- rbind(survgrid, bind_rows(surv_rays))
    full_res <- list()

    # use nearest neighbor interpolation to interpolate onto grid, if desired, and return that interpolation
    if (!is.null(interp_grid)) {
        interp_res <- get.knnx(data=full_survgrid[,c(1,2)], query=interp_grid, k=1)
        survfunc <- full_survgrid[as.vector(interp_res$nn.index),3]
        full_interp_res <- cbind(interp_grid, survfunc)
        full_res$full_survgrid_interpolated <- full_interp_res
    }

    full_res$q_isoline <- q_isoline
    full_res$full_survgrid <- full_survgrid

    return(full_res)
}

# function to evaluate the survival probability of a point in polar coords using this blended method
# must input a radius for a corresponding point on the desired q_isoline, which shares the same theta
### questionable if i even need this, since desired level set radii can be computed analytically
evaluateSurvivalPointDiff <- function(r, dat, q, theta, q_iso_r, lvlset) {
    
    x <- r*cos(theta)
    y <- r*sin(theta)
    empirical_prob <- mean(dat[,1] > x & dat[,2] > y)
    
    if (empirical_prob >= q) { 
        return(empirical_prob - lvlset) 
    }
    
    rv_prob <- (q_iso_r/r)*q
    return(rv_prob - lvlset)
    
}

# now, rewrite how to draw the tube bound given a bhat and a means of computing the survival function
drawExtremeTubeBound <- function(dat, q, q_isoline, lvlset) {

    numCoords <- length(q_isoline$r)
    radii <- rep(NA, numCoords)
    angles <- q_isoline$theta
    # 50 is arbitrary, just some padding
    maxRad <- max((q_isoline$r)*(q/lvlset)) + 50

    for (i in 1:numCoords) {
        # note: may need to pay attention to the fact that this is not a continuous function..
        angle <- angles[i]
        radii[i] <- uniroot(evaluateSurvivalPointDiff, interval=c(0, maxRad), dat=dat, q=q, theta=angle, q_iso_r=q_isoline$r[i], lvlset=lvlset)$root

    }

    xs <- radii*cos(angles)
    ys <- radii*sin(angles)

    tube <- data.frame(X1=xs, X2=ys, r=radii, theta=angles)

    return(tube)

}


drawExtremeTubeBounds <- function(dat, bnhat, p, q, q_isoline) {

    lvlset_top <- p - bnhat
    lvlset_bottom <- p + bnhat

    # by default, tube is unbounded up and to the right
    tube_top <- 'unbounded'

    # if not, actually find the upper bound of the tube
    if (lvlset_top > 0) {
        tube_top <- drawExtremeTubeBound(dat, q, q_isoline, lvlset_top)
    }

    tube_bottom <- drawExtremeTubeBound(dat, q, q_isoline, lvlset_bottom)

    tube_bounds <- list()
    tube_bounds$top <- tube_top
    tube_bounds$bottom <- tube_bottom

    return(tube_bounds)

}

# Big function to draw all the necessary regions!

drawExtremeRegions <- function(dat, emp_grid_ticks, alphas, ps, beta_funcs_dict, B, num_angles, num_radii, rv_ub) {
    # Function to draw one or more uncertainty regions for some base isoline of a specified
    # exceedence probability (p). Uses the methods in Mammen and Polonik (2013)
    # and also blends a bit of the code from Cooley.
    #
    # Arguments:
    # rv_ubs is a dictionary of desired exceedance probabilities as keys and the upper bound of rv grid as values

    # get all combinations of p and beta functions
    beta_func_labs <- names(beta_funcs_dict)
    combos <- expand.grid(beta_func_labs=beta_func_labs, ps=ps)
    p_beta_lst <- vector(mode='list', length=nrow(combos))
    # start list of names
    names <- c()

    # choose the starting q, to be a function of n
    # for now, chosen such that we have on average 50 joint exceedances per point on isoline
    q <- 100/nrow(dat)

    # evaluate survival function on whole dataset
    main_surv_res <- evaluateSurvivalGrid(dat, q, emp_grid_ticks, num_angles, num_radii, rv_ub)
    main_surv_grid <- main_surv_res$full_survgrid
    q_isoline <- main_surv_res$q_isoline

    # for each p-beta function combo, get the hhat values, the boolean mask for
    # the points in Delta, and preallocate the bootstrap list of Zs
    # (all in the language of Mammen and Polonik (2013))
    for (i in 1:nrow(combos)) {
        p <- combos$ps[i]
        beta_func_lab <- combos$beta_func_labs[i]
        hhat_vals <- -main_surv_grid$survfunc + p
        # boolean mask to extract points in blown up boundary around est isoline
        deltamask <- abs(hhat_vals) <= beta_funcs_dict[[beta_func_lab]](nrow(dat))
        Zs <- rep(0, B)
        lst <- list(hhat_vals=hhat_vals, deltamask=deltamask, Zs=Zs)
        names <- c(names, paste0('beta', beta_func_lab, '_p', p))
        p_beta_lst[[i]] <- lst
    }
    names(p_beta_lst) <- names

    # for every bootstrap resample
    for (i in 1:B) {

        boot_dat <- dat %>% sample_frac(1, replace = TRUE)
        boot_surv_grid <- evaluateSurvivalGrid(boot_dat, q, emp_grid_ticks, num_angles, num_radii, rv_ub, main_surv_grid[,c(1,2)])

        # for every p-beta combo, get the bootstrap distribution of Z(beta)
        for (j in 1:nrow(combos)) {
            p <- combos$ps[j]
            deltamask <- p_beta_lst[[j]]$deltamask
            hhat_vals <- p_beta_lst[[j]]$hhat_vals
            boot_hhat_vals <- -boot_surv_grid$full_survgrid_interpolated$survfunc + p

            p_beta_lst[[j]]$Zs[i] <- max(abs((boot_hhat_vals - hhat_vals)[deltamask]))
        }
    }
    # create lists of output confidence regions
    output_lst <- list()
    # labels for the output list
    newnames <- c()
    # index counter variable
    k <- 1

    # for every p-beta combo
    for (i in 1:length(p_beta_lst)) {

        name <- names(p_beta_lst)[i]
        Zs <- p_beta_lst[[i]]$Zs
        hhat_vals <- p_beta_lst[[i]]$hhat_vals

        # for every p-beta-alpha triple, draw confidence regions and store results
        for (j in 1:length(alphas)) {
            alpha <- alphas[j]
            bhat <- as.numeric(quantile(Zs, probs = 1 - alpha))
            #tube_bounds <- drawExtremeTubeBounds(dat, bhat, combos$ps[i], q, q_isoline)
            out <- list()
            #out$tube_top <- tube_bounds$top
            #out$tube_bottom <- tube_bounds$bottom
            out$alpha <- alpha
            out$p <- combos$ps[i]
            out$bhat <- bhat
            out$B <- B
            out$beta_func <- combos$beta_func_labs[i]
            out$beta_function <- beta_funcs_dict[[combos$beta_func_labs[i]]]
	    out$tube_bounds <- drawExtremeTubeBounds(dat, bhat, combos$ps[i], q, q_isoline)

            newnames <- c(newnames, paste0(name, '_alpha', alpha))
            output_lst[[k]] <- out
            # update counter
            k <- k+1
        }
    }
    # label the list elements for easier bookkeeping
    names(output_lst) <- newnames

    output_lst$data <- dat

    return(output_lst)
}


##################################################################################################
######################################## Running the Code ########################################
ns <- c(10000, 50000, 100000)
alphas <- c(0.01, 0.05, 0.1)
ps <- c(0.0001)
rv_ub <- 10000
gticks <- 200
num_angles <- 100
num_radii <- 100
beta_funcs_dict <- list()
beta_funcs_dict[[as.character(1/2)]] <- function(n) {return((1/n)^(1/2))}
B <- 500

########################

n_cores <- 64
n_iter <- 500

path <- '/pscratch/sd/j/jbbutler/sims/regions/mixed_confregs/'
distribution <- 'bivtcopula_frecmargins'

# store simulation results for particular grid settings
subpath <- paste0(path, distribution, '/')
# create the directory
dir.create(subpath)

parallelizedCode <- function(ind) {

    bivt_dat <- data.frame(rmvt(n, sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2), df = 4))
    unif_dat <- data.frame(X1=pt(bivt_dat[,1], df=4), X2=pt(bivt_dat[,2], df=4))
    frec_dat <- -1/log(unif_dat)
    
    regions <- drawExtremeRegions(frec_dat, gticks, alphas, ps, beta_funcs_dict, B, num_angles, num_radii, rv_ub)
    iso_prob_tag <- paste0(ps[1], '_prob')
    dir.create(path = paste0(subpath, '/', iso_prob_tag), showWarnings=FALSE)
    fname <- paste0(n, '_samples')
    full_path <- paste0(subpath, '/', iso_prob_tag, '/', fname)
    dir.create(path = full_path, showWarnings = FALSE)

    saveRDS(regions, file = paste0(full_path, '/', 'simulation_', ind, '.RData'))

}

for (i in 1:length(ns)) {

    clust <- makeSOCKcluster(n_cores)
    registerDoSNOW(clust)
    pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
    progress <- function(n) setTxtProgressBar(pb ,n)
    opts <- list(progress = progress)
    
    n <- ns[i]
    foreach(i = 1:n_iter,
                   .packages = c('mvtnorm', 'purrr', 'data.table', 'dplyr', 'ismev', 'ks', 'FNN'),
                   .options.snow = opts) %dopar% {

        parallelizedCode(i)

    }
    close(pb)
    stopCluster(clust)
}



