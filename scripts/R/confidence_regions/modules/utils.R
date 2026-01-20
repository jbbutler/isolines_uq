library(mvtnorm)
library(dplyr)
library(purrr)
library(data.table)
library(ismev)
library(evd)

source('~/isolines_uq/scripts/R/confidence_regions/modules/karachiTools.R')

# function to evaluate an estimate of a univariate cdf on a set of points
# uses ecdf for points below the 1-n^(-gamma) quantile, and uses GPD for points above
est_cdf <- function(x, dat, gamma) {
    edf <- ecdf(dat)
    edf_vals <- edf(x)
    threshold_prob <- length(dat)^(-gamma)
    threshold <- quantile(dat, 1-threshold_prob)
    x_gpd <- x[x > threshold]
    gpdOut <- ismev::gpd.fit(dat,
					threshold = threshold, show = F)
    gpd_probs <- 1-(1-evd::pgpd(x_gpd, loc = gpdOut$threshold, scale = gpdOut$mle[1], shape = gpdOut$mle[2]))*(threshold_prob)
    edf_vals[x > threshold] <- gpd_probs
    return(edf_vals)   
}

# the inverse of the above function
est_inv_cdf <- function(p, dat, gamma) {

    threshold_prob <- length(dat)^(-gamma)
    threshold <- as.numeric(quantile(dat, 1-threshold_prob))
    quantiles <- as.numeric(quantile(dat, p))
    gpd_probs <- p[p > 1-threshold_prob]
    gpdOut <- ismev::gpd.fit(dat,
					threshold = threshold, show = F)
    quantiles[p > 1-threshold_prob] <- evd::qgpd(1-((1-gpd_probs)/threshold_prob), loc = gpdOut$threshold, scale = gpdOut$mle[1], shape = gpdOut$mle[2])
    quantiles <- as.numeric(quantiles)
    return(quantiles)

}

# function to estimate the survival function using the blended empirical + regularly varying method
blendedSurvivalFunc <- function(pt, dat, gamma, xi) {

    qn <- nrow(dat)^(-gamma)
    empsurv_prob <- mean((dat[,1] > pt[1]) & (dat[,2] > pt[2]))

    if (empsurv_prob >= qn) {
        return(empsurv_prob)
    }
    else {
        theta <- atan(pt[2]/pt[1])
        empsurvDiff <- function(r) { return(mean((dat[,1] > r*cos(theta)) & (dat[,2] > r*sin(theta)))-qn) }
        rp <- sqrt(pt[2]**2 + pt[1]**2)
        rq <- uniroot(empsurvDiff, interval=c(0, rp))$root

        return(((rq/rp)^(1/xi))*qn)
    }
}

# function to load up the different distributions (sampling functions)
loadSamplingFunction <- function(dist) {
    
    if (dist=='bivt') {
        samplingFunction <- function(n) return(data.frame(rmvt(n, sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2), df = 4)))
    }
    if (dist=='bivgauss') {
        samplingFunction <- function(n) return(data.frame(rmvnorm(n, mean = rep(0, 2), sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2))))
    }
    if (dist=='karachi') {
        samplingFunction <- function(n) return(rKarachiBetaKDE(n))
    }

    return(samplingFunction)
}

# function to to compute the empirical survival function (given observed data) on a region of points which are not
# a regular grid
# works by converting all of the points to a regular grid, doing a fast empirical survival function operation
# and then taking only those points we were interested in in the first place
computeEmpSurvIrregular <- function(region, dat) {
    
    sup_reg_xs_unique <- sort(unique(region$X1))
    sup_reg_ys_unique <- sort(unique(region$X2))
    full_grid <- expand.grid(X1=sup_reg_xs_unique, X2=sup_reg_ys_unique)
    surv <- fastEmpSurv(full_grid, dat)

    full_surv_res <- data.table(X1=full_grid$X1, X2=full_grid$X2, surv=surv)
    sup_reg <- data.table(region)

    res <- full_surv_res[sup_reg, on=c('X1', 'X2')]
    
    return(res$surv)
}

# function that selects points from an empirical isoline (from the empirical survival function)
# function that draws a grid of points over which we wish to evaluate the sup difference between
# the true survival function (in bootstrap world) and the empirical estimate (bootstrap survival function in bootstrap world)
drawEmpiricalIsoline <- function(dat, n_coords, gridLbs, p) {
    
    # function that we will use our root-finding algorithm for
    empSurvFixedCoord <- function(radius, angle, prob, xCenter, yCenter) {
        xCoord <- xCenter + radius*cos(angle)
        yCoord <- yCenter + radius*sin(angle)
        exceedanceProb <- mean((dat[,1] > xCoord) & (dat[,2] > yCoord))
        return(exceedanceProb - prob)
    }
    
    radii <- rep(NA, n_coords)
    angles <- seq(0, pi/2, length.out=n_coords)
    max_dat_r <- max(10 + sqrt((dat[,1]**2) + (dat[,2]**2)))
    gridUbs <- c(max_dat_r, max_dat_r)
    maxRad <- sqrt(sum((max_dat_r-gridLbs)**2))

    # for each angle in first quadrant, find radius that gives a point with desired exceedance probability
    # by finding roots of pmvtFixedCoord given the angle
    # find both the lower and upper radius, and then store a mesh of them
    for (i in 1:n_coords) {
        radii[i] <- uniroot(empSurvFixedCoord, interval=c(0, maxRad), angle=angles[i],
                             prob=p, xCenter=gridLbs[1], yCenter=gridLbs[2])$root
    }

    # convert back to cartesian coordinates
    X1 <- radii*cos(angles) + gridLbs[1]
    X2 <- radii*sin(angles) + gridLbs[2]
    return(data.frame(cbind(X1, X2)))
    
}

# function to draw the isoline estimate using the blended empirical + regularly varying method
# Arguments
# dat: the data
# n_coords: the number of coordates you want
# gridLbs: the lower lefthand corner of the grid
# gridUbs: the upper bounds of the grid that will contain the full isoline
# gamma: parameter controlling how far into the tail you will be using empirical df
# xi: 1/the index of regular variation, EV index
drawExtremeIsoline <- function(dat, p, n_coords, gridLbs, gamma, xi) {
    
    q <- nrow(dat)^(-gamma)
    q_isoline <- drawEmpiricalIsoline(dat, n_coords, gridLbs, q)
    p_isoline <- q_isoline*((q/p)^(xi))
    
    return(p_isoline)
    
}

# function that draws a grid of points over which we wish to evaluate the sup difference between
# the true survival function (in bootstrap world) and the empirical estimate (bootstrap survival function in bootstrap world)
drawSupRegion <- function(dat, n_boundary_coords, n_interior_coords, gridLbs, gridUbs, p, beta) {
    
    # function that we will use our root-finding algorithm for
    empSurvFixedCoord <- function(radius, angle, prob, xCenter, yCenter) {
        xCoord <- xCenter + radius*cos(angle)
        yCoord <- yCenter + radius*sin(angle)
        exceedanceProb <- mean((dat[,1] > xCoord) & (dat[,2] > yCoord))
        return(exceedanceProb - prob)
    }
    
    radii <- rep(NA, n_boundary_coords*n_interior_coords)
    full_angles <- rep(NA, n_boundary_coords*n_interior_coords)
    angles <- seq(0, pi/2, length.out=n_boundary_coords)
    maxRad <- sqrt(sum((gridUbs-gridLbs)**2))

    # for each angle in first quadrant, find radius that gives a point with desired exceedance probability
    # by finding roots of pmvtFixedCoord given the angle
    # find both the lower and upper radius, and then store a mesh of them
    for (i in 1:n_boundary_coords) {
        angle <- angles[i]
        lower_radius <- uniroot(empSurvFixedCoord, interval=c(0, maxRad), angle=angle,
                            prob=p+beta, xCenter=gridLbs[1], yCenter=gridLbs[2])$root
        upper_radius <- uniroot(empSurvFixedCoord, interval=c(0, maxRad), angle=angle,
                             prob=p-beta, xCenter=gridLbs[1], yCenter=gridLbs[2])$root
        
        radii_per_angle <- seq(lower_radius, upper_radius, length.out=n_interior_coords)
        angles_per_angle <- rep(angle, length.out=n_interior_coords)
        
        radii[((i-1)*n_interior_coords + 1):(i*n_interior_coords)] <- radii_per_angle
        full_angles[((i-1)*n_interior_coords + 1):(i*n_interior_coords)] <- angles_per_angle
        
        
    }

    # convert back to cartesian coordinates
    X1 <- radii*cos(full_angles) + gridLbs[1]
    X2 <- radii*sin(full_angles) + gridLbs[2]
    return(data.frame(cbind(X1, X2)))
    
}

# function to evaluate the exact Gaussian KDE at a single grid point
# credit to Cooley et al. source code
kernSurv <- function(loc, dat, bw)
{
        p1 <- 1 - pnorm(loc[1], mean = dat[,1], sd = bw[1]/4)
        p2 <- 1 - pnorm(loc[2], mean = dat[,2], sd = bw[2]/4)
        return(mean(p1*p2))
}

empSurv <- function(point, dat) {
    points <- as.numeric(point)
	val <- mean(dat[,1] > point[1] & dat[,2] > point[2])
        return(val)
}

# function to draw a single bounding curve for the uncertainty tube, given a particular sample size and
# value for bnhat
# Note that this function exactly draws the tube for, assuming your survival function is the empirical survival function
drawTubeBoundExact <- function(dat, exceedances, lbs, ubs) {
    
    subdat <- dat %>% filter(X1 >= lbs[1], X2 >= lbs[2])
    y_ordered <- sort(subdat$X2, decreasing=TRUE)

    # maximum number of points in each isoline is:
    # 1 for the starting point, and 2 for each remaining point in subdat
    ys <- rep(0, nrow(subdat)*2 + 1)
    xs <- rep(0, nrow(subdat)*2 + 1)

    ys[1] <- y_ordered[exceedances]
    xs[1] <- lbs[1]

    # use for faster lookups/boolean indexing
    arranged_xs <- as.data.table(subdat %>% arrange(X1))
    arranged_ys <- as.data.table(subdat %>% arrange(desc(X2)))

    ys[2] <- ys[1]
    xs[2] <- as.numeric(arranged_xs[X2 >= ys[1]][1,1])

    end_alg <- nrow(arranged_xs[X1 >= xs[2]]) == exceedances

    # will start adding points to the isoline lists in the third index
    ind <- 3

    while (!end_alg) {
        # grab the last point you are stepping down and to the right from
        last_pt_x <- xs[ind-1]
        last_pt_y <- ys[ind-1]
    
        # y you step down to is highest y among all points down and to the right
        poss_ys <- arranged_ys[X1 > last_pt_x & X2 < last_pt_y]
    
        new_y <- as.numeric(poss_ys[1,2])
    
        poss_xs <- arranged_xs[X1 > last_pt_x & X2 >= new_y]
    
        new_x <- as.numeric(poss_xs[1,1])

        if (is.na(new_x) | is.na(new_y)) {
            break
        }

        xs[ind] <- last_pt_x
        xs[ind+1] <- new_x
        ys[ind] <- new_y
        ys[ind+1] <- new_y
    
        ind <- ind + 2
    
        end_alg <- nrow(arranged_xs[X1 >= new_x]) == exceedances
    }

    xs[ind] <- xs[ind-1]
    ys[ind] <- lbs[2]

    points_x <- xs[1:ind]
    points_y <- ys[1:ind]
    
    # subset to have isoline contained in upper bound of grid
    curve <- data.frame(X1=points_x, X2=points_y)
    curve <- curve[curve$X1 <= ubs[1] & curve$X2 <= ubs[2],]
    
    return(curve)
    
}

# function to draw the bounds of the p-isoline region for certain sample with a certain bnhat
# Note this draws it exactly, assuming your survival function is the empirical survival function
drawTubeBoundsExact <- function(dat, bnhat, p, lbs, ubs) {

    n <- nrow(dat)
    exceedances_top <- ceiling(n*(p-bnhat))
    exceedances_bottom <- floor(n*(p+bnhat)) + 1

    if (exceedances_top <= 0) {
	tube_top <- data.frame(X1=c(lbs[1], ubs[1], ubs[1]), 
			       X2=c(ubs[2], ubs[2], lbs[2]))
    } else {
        tube_top <- drawTubeBound(dat, exceedances_top, lbs, ubs)
    }

    tube_bottom <- drawTubeBound(dat, exceedances_bottom, lbs, ubs)
    
    tube_bounds <- list()
    tube_bounds$top <- tube_top
    tube_bounds$bottom <- tube_bottom
    
    return(tube_bounds)
    
}

# function to draw one of the bounding curves of your confidence tube, assuming use of the empirical survival
# function as your survival function (needed to be handled separately since this function has flat parts)
# Note: obtains some finite number of points on the curve, so it is approximate but faster than the previous
# exact methods of drawing the tube bounds!
drawTubeBoundES <- function(dat, numCoords, gridLbs, exceedances) {
    
    subdat <- dat %>% filter(X1 >= gridLbs[1], X2 >= gridLbs[2])    
    y_ordered <- sort(subdat$X2, decreasing=TRUE)
    x_ordered <- sort(subdat$X1, decreasing=TRUE)
    # find largest y and x coordinates of any points on isoline
    # use these to find upper bound of search window for root-finding algorithm
    # for good measure, add 1 to make it definitively larger!
    ymax <- y_ordered[exceedances] + 1
    xmax <- x_ordered[exceedances] + 1
    ubs <- c(xmax, ymax)
    
    numExceedances <- function(radius, angle, dat, desired_exceedances, xCenter, yCenter) {
        xCoord <- xCenter + radius*cos(angle)
        yCoord <- yCenter + radius*sin(angle)
        actual_exceedances <- sum((dat[,1] > xCoord) & (dat[,2] > yCoord))
        diff <- actual_exceedances - (desired_exceedances - 0.5)   
        return(diff)
    }
    
    radii <- rep(NA, numCoords)
    angles <- seq(0, pi/2, length.out=numCoords)
    maxRad <- sqrt(sum((ubs - gridLbs)**2))
    
    for (i in 1:numCoords) {
        # note: may need to pay attention to the fact that this is not a continuous function..
        angle <- angles[i]
        radii[i] <- uniroot(numExceedances, interval=c(0, maxRad), angle=angle, dat=dat, desired_exceedances=exceedances, xCenter=gridLbs[1], yCenter=gridLbs[2])$root
        
    }
    
    xs <- radii*cos(angles) + gridLbs[1]
    ys <- radii*sin(angles) + gridLbs[2]
    
    tube <- data.frame(xs, ys)
    colnames(tube) <- c('X1', 'X2')
    return(tube)
    
}

# Function that draws both of the tube bounds for the confidence tubes approximately, using
# the empirical survival function as the function estimate
drawTubeBoundsES <- function(dat, numCoords, bnhat, p, gridLbs, gridUbs) {

    n <- nrow(dat)
    exceedances_top <- ceiling(n*(p-bnhat))
    exceedances_bottom <- floor(n*(p+bnhat)) + 1

    if (exceedances_top <= 0) {
        tube_top <- data.frame(X1=c(gridLbs[1], gridUbs[1], gridUbs[1]),
                               X2=c(gridUbs[2], gridUbs[2], gridLbs[2]))
    } else {
        tube_top <- drawTubeBoundES(dat, numCoords, gridLbs, exceedances_top)
    }

    tube_bottom <- drawTubeBoundES(dat, numCoords, gridLbs, exceedances_bottom)

    tube_bounds <- list()
    tube_bounds$top <- tube_top
    tube_bounds$bottom <- tube_bottom

    return(tube_bounds)

}


    
fastEmpSurv <- function(grid, dat) {
  # Note: code is adapted from the mltools package source code

  dat <- data.table(dat)
  # flip the grid points around
  grid <- data.table(grid[nrow(grid):1,])
  
  dat_copy <- copy(dat[, names(grid), with=FALSE])
  uboundDT <- unique(data.table(grid[['X1']], grid[['X1']]))
  setnames(uboundDT, c('X1', paste0("Bound.X1")))
  dat_copy <- uboundDT[dat_copy, on='X1', roll=Inf, nomatch=0]
  uboundDT <- unique(data.table(grid[['X2']], grid[['X2']]))
  setnames(uboundDT, c('X2', paste0("Bound.X2")))
  dat_copy <- uboundDT[dat_copy, on='X2', roll=Inf, nomatch=0]
  
  binned.uniques <- dat_copy[, .N, keyby=eval(paste0("Bound.", names(grid)))]
  setnames(binned.uniques, paste0("Bound.", names(grid)), names(grid))
  grid <- binned.uniques[grid, on=names(grid)]
  grid[is.na(N), N := 0]
  
  grid[, N.cum := cumsum(N), by='X2']
  grid[, N.cum := cumsum(N.cum), by='X1']
  grid[, `:=`(N = NULL, Surv = N.cum/nrow(dat))]
  
  # reflip survival function vals to be in same order that grid points came in
  vals <- grid$Surv[nrow(grid):1]
  return(vals)
}

gatherNaRes <- function(output_lst, samp_sizes) {

    na_pers <- c()
    na_bcas <- c()
    num_na_coords <- c()
    num_na_jacks <- c()

    for (i in 1:length(output_lst)) {

	na_pers <- c(na_pers, output_lst[[i]]$na_per)
        na_bcas <- c(na_bcas, output_lst[[i]]$na_bca)
        num_na_coords <- c(num_na_coords, output_lst[[i]]$num_na_coord)
	num_na_jacks <- c(num_na_jacks, output_lst[[i]]$num_na_jack)

    }

    na_df <- data.frame(samp_sizes, na_pers, na_bcas, num_na_coords, num_na_jacks)

    colnames(na_df) <- c('n', 'NA Percentile CIs', 'NA BCA CIs', 'NA Coords', 'NA Jackknifes')

    return(na_df)
}

gatherNaResComparison <- function(output_lst, coords) {

    na_pers <- c()
    na_bcas <- c()
    num_na_coords <- c()
    num_na_jacks <- c()

    for (i in 1:length(output_lst)) {

        na_pers <- c(na_pers, output_lst[[i]]$na_per)
        na_bcas <- c(na_bcas, output_lst[[i]]$na_bca)
        num_na_coords <- c(num_na_coords, output_lst[[i]]$num_na_coord)
        num_na_jacks <- c(num_na_jacks, output_lst[[i]]$num_na_jack)

    }

    na_df <- data.frame(coords, na_pers, na_bcas, num_na_coords, num_na_jacks)

    colnames(na_df) <- c('Coord.', 'NA Percentile CIs', 'NA BCA CIs', 'NA Coords', 'NA Jackknifes')

    return(na_df)


}

gatherConfintRes <- function(output_lst, samp_sizes) {

    cvg_pers <- c()
    avg_pers <- c()
    cvg_bcas <- c()
    avg_bcas <- c()

    for (i in 1:length(output_lst)) {

        cvg_pers <- c(cvg_pers, output_lst[[i]]$cvg_per)
        avg_pers <- c(avg_pers, output_lst[[i]]$avg_per)
	cvg_bcas <- c(cvg_bcas, output_lst[[i]]$cvg_bca)
	avg_bcas <- c(avg_bcas, output_lst[[i]]$avg_bca)

    }

    res_df <- data.frame(samp_sizes, cvg_pers, avg_pers, cvg_bcas, avg_bcas)

    colnames(res_df) <- c('n', 'Percentile Coverage', 'Percentile Avg. Width', 'BCA Coverage', 'BCA Avg. Width')

    return(res_df)

}

extractBoundary <- function(region) {
    # region: dataframe of points corresponding to confidence set, restricted to 1st quadrant
    
    colnames(region) <- c('Var1', 'Var2')
    xvals <- sort(unique(region$Var1), decreasing = FALSE)
    num_xs <- length(xvals)

    left_line <- filter(region, Var1 == min(region$Var1))
    right_line <- filter(region, Var1 == max(region$Var1))

    upper_pts <- vector(mode = 'list')
    lower_pts <- vector(mode = 'list')

    for (i in 2:(num_xs-1)) {
    
        x <- xvals[i]
        next_x <- xvals[i+1]
        last_x <- xvals[i-1]
        upper_pts_x <- filter(region, Var1 == x, Var2 >= max(filter(region, Var1 == next_x)$Var2))
        lower_pts_x <- filter(region, Var1 == x, Var2 <= min(filter(region, Var1 == last_x)$Var2))
    
        upper_pts[[i]] <- upper_pts_x
        lower_pts[[i]] <- lower_pts_x
    
    }

    upper_pts[[1]] <- left_line
    upper_pts[[num_xs]] <- right_line

    boundary <- bind_rows(upper_pts, lower_pts)
    
    return(boundary)
    
}

extendIsoline <- function(isoline, lb) {
    # function to extend the sides of a pre-drawn isoline to the edges of a larger grid

    left <- isoline[which.min(isoline$X1),]
    bottom <- isoline[which.min(isoline$X2),]

    left <- sfg_linestring(rbind(left, c(lb, left$X2)))
    bottom <- sfg_linestring(rbind(bottom, c(bottom$X1, lb)))

    isoline <- sfg_multilinestring(isoline)

    isoline_aug <- st_union(st_union(left, isoline), bottom)

    return(isoline_aug)
}

extendRegion <- function(region_pts, polyregion_pts, bds_x1, bds_x2) {

    # extends the projection region to edges of the full grid if the extent of points does not match full grid size
    # to do: investigate wy the projection region works this way?


    mins <- apply(region_pts, 2, min)
    lefts <- c(min((region_pts %>% filter(X1 == mins[[1]]))$X2), max((region_pts %>% filter(X1 == mins[[1]]))$X2))
    rights <- c(min((region_pts %>% filter(X2 == mins[[2]]))$X1), max((region_pts %>% filter(X2 == mins[[2]]))$X1))
    mins <- ceiling(mins*10^(4))/(10^4)
    left_topoff <- data.frame(X1 = c(rep(bds_x1[1], 2), rep(mins[[1]], 2)), X2 = c(lefts, lefts[2], lefts[1]))
    right_topoff <- data.frame(X1 = c(rights, rights[2], rights[1]), X2 = c(rep(bds_x2[1], 2), rep(mins[[2]], 2)))

    left_topoff <- sfg_polygon(left_topoff)
    right_topoff <- sfg_polygon(right_topoff)
    polyregion_aug <- st_union(st_union(polyregion_pts, left_topoff), right_topoff)

    return(polyregion_aug)
}
computeUpperPadding <- function(beta, beta_zero, surv_func, boot_diffs, p, alpha) {
    
    # function to compute the padding on the upper bound, for use in `optimize` function
    # when minimizing with respect to beta

    getZ <- function(beta, surv_func, p, boot_diffs) {
    deltamask <- abs(-surv_func$estimate + p) <= beta
    return(max(boot_diffs[deltamask]))
    }
    
    Z0s <- map(boot_diffs, getZ, beta = beta_zero, surv_func = surv_func, p = p)
    
    Zbetas <- map_dbl(boot_diffs, getZ, beta = beta, surv_func = surv_func, p = p)
    bhat <- quantile(Zbetas, prob = 1 - alpha)[[1]]
    
    ts <- seq(min(Zbetas), max(Zbetas), length.out = 1000)
    e <- ecdf(Zbetas)

    lower_bool <- e(ts) <= 1-alpha
    upper_bool <- e(ts) >= 1-alpha

    bminus <- ts[upper_bool][which.min(e(ts)[upper_bool])]
    
    pn0 <- mean(Z0s <= bhat)
    FZn_bminus <- mean(Zbetas <= bminus)
    
    ub_padding <- pn0 - FZn_bminus
    
    return(abs(ub_padding))
    
}
