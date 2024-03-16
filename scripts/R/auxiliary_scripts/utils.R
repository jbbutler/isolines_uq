library(mvtnorm)
library(dplyr)
library(purrr)
library(data.table)

# function to evaluate the exact Gaussian KDE at a single grid point
# credit to Cooley et al. source code
kernSurv <- function(loc, dat, bw)
{
        p1 <- 1 - pnorm(loc[1], mean = dat[,1], sd = bw[1]/4)
        p2 <- 1 - pnorm(loc[2], mean = dat[,2], sd = bw[2]/4)
        return(mean(p1*p2))
}

empSurv <- function(point, dat) {

	val <- mean(dat[,1] > point[1] & dat[,2] > point[2])
        return(val)
}

# function to draw a single bounding curve for the uncertainty tube, given a particular sample size an
drawTubeBound <- function(dat, exceedances, lbs, ubs) {
    
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
drawTubeBounds <- function(dat, bnhat, p, lbs, ubs) {

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
