# Script containing code to simulate the upper bound of Lemma 2.1 in Mammen and Polonik (2013)
#
# Jimmy Butler

library(mvtnorm)
library(ggplot2)
library(mnormt)
library(tidyr)
library(foreach)
library(doSNOW)
library(parallel)

source('/global/u1/j/jbbutler/isolines_uq/scripts/R/auxiliary_scripts/utils.R')
source('~/isolines_uq/scripts/R/auxiliary_scripts/distributionIsolines.R')

##### Helper Functions #####
evalbns <- function(alpha, Zn) {

    n <- length(Zn)
    ecdf_vals <- 1:length(Zn)/length(Zn)
    sorted <- sort(Zn)

    bminus <- sorted[which(ecdf_vals >= alpha)[1]]
    bplus <- sorted[which(ecdf_vals > alpha)[1]]

    return(c(bminus, bplus))
}

# function to find the angle corresponding to a specific point in the cartesian grid
drawEmpiricalIsoline <- function(dat, angles, p, lbs) {

    survivalDiff <- function(r, theta, p, dat) {

        xPt <- lbs[1] + r*cos(theta)
        yPt <- lbs[2] + r*sin(theta)

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

    xs <- lbs[1] + radii*cos(angles)
    ys <- lbs[2] + radii*sin(angles)

    return(data.frame(X1=xs, X2=ys, r=radii, theta=angles))
}

empSurv <- function(pt, dat) {
    point <- as.numeric(pt)
    return(mean((point[1] < dat[,1]) & (point[2] < dat[,2])))
}

rawSupRegion <- function(n_boundary_coords, n_interior_coords, gridLbs, gridUbs, p, beta) {

    # function that we will use our root-finding algorithm for
    pmvtFixedCoord <- function(radius, angle, sigma, df, prob, xCenter, yCenter) {
        xCoord <- xCenter + radius*cos(angle)
        yCoord <- yCenter + radius*sin(angle)
        exceedanceProb <- pmvt(lower=c(xCoord, yCoord), upper=c(Inf, Inf), df=df, sigma=sigma)
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
        lower_radius <- uniroot(pmvtFixedCoord, interval=c(0, maxRad), angle=angle,
                             sigma=matrix(c(1, 0.7, 0.7, 1), nrow=2), df=4, prob=p+beta, xCenter=gridLbs[1], yCenter=gridLbs[2])$root
        upper_radius <- uniroot(pmvtFixedCoord, interval=c(0, maxRad), angle=angle,
                             sigma=matrix(c(1, 0.7, 0.7, 1), nrow=2), df=4, prob=p-beta, xCenter=gridLbs[1], yCenter=gridLbs[2])$root

        radii_per_angle <- seq(lower_radius, upper_radius, length.out=n_interior_coords)
        angles_per_angle <- rep(angle, length.out=n_interior_coords)

        radii[((i-1)*n_interior_coords + 1):(i*n_interior_coords)] <- radii_per_angle
        full_angles[((i-1)*n_interior_coords + 1):(i*n_interior_coords)] <- angles_per_angle
    }

    # convert back to cartesian coordinates
    xs <- radii*cos(full_angles) + gridLbs[1]
    ys <- radii*sin(full_angles) + gridLbs[2]

    return(cbind(xs, ys))
}

#######################################################################

### simulation specs ###
n_cores <- 64
# specifications for the survival function grid
lbs <- c(-2,-2)
ubs <- c(10, 10)
n_boundary_coords <- 100
n_interior_coords <- 5

# specifications for the confidence regions
p <- 0.1
alpha <- 0.05
n <- 10000

# list of betas to choose from (ranging from small to large)
# highest beta chosen as 1+p because this beta covers all of R2, any higher beta will still cover this
n_betas <- 128
betas <- seq(0.05, 0.000001, length.out=n_betas)

n_monte_carlo <- 10000

#original_dat <- data.frame(rmvt(1000000, df=4))
#true_survfunc <- fastEmpSurv(grid, original_dat)
#isoline <- drawEmpiricalIsoline(original_dat, seq(0, pi/2, length.out=n_iso_pts), p, c(-2,-2))
true_survfunc <- pmt(x=-as.matrix(grid), mean=rep(0,2), S=matrix(c(1,0.7,0.7,1), 2,2), df=4)
isoline <- drawBivtIsoline(numCoords=n_iso_pts, gridUbs=c(10,10), gridLbs=c(-2,-2), prob=p)

results <- list()

parallelizedCode <- function(k) {

    # sample_dat <- original_dat[sample(1:n, replace=TRUE),]
    # sample_dat <- data.frame(rmvnorm(n, mean=rep(0,2), sigma=matrix(c(1,0.7,0.7,1), 2,2)))
    sample_dat <- data.frame(rmvt(n, sigma = matrix(c(1, 0.7, 0.7, 1), nrow = 2), df = 4))
    sample_survfunc <- fastEmpSurv(grid, sample_dat)

    Zn <- max(abs((sample_survfunc - true_survfunc)[deltamask]))
    isoline_survfunc <- apply(isoline, 1, empSurv, dat=sample_dat)
    Z0 <- max(abs(isoline_survfunc - p))

    return(c(Zn, Z0))
}

for (i in 1:length(beta_funcs)) {

    print(paste0('Starting betafunc ', i))

    bound <- rep(NA, length(ns))

    for (j in 1:length(ns)) {

        deltamask <- abs(true_survfunc-p) <= beta_funcs[[i]](ns[j])
	n <- ns[j]

	clust <- makeSOCKcluster(n_cores)
    	registerDoSNOW(clust)
    	pb <- txtProgressBar(min = 1, max = n_reps, style = 3)
    	progress <- function(n) setTxtProgressBar(pb ,n)
    	opts <- list(progress = progress)

        replicates <- foreach(k = 1:n_reps, 
			      .options.snow = opts, 
			      .packages = c('mvtnorm', 'data.table', 'dplyr'), 
			      .combine='rbind') %dopar% parallelizedCode(k)

	close(pb)
    	stopCluster(clust)

	#Zn <- as.numeric(replicates[,1])
        #Z0 <- as.numeric(replicates[,2])

        #bns <- evalbns(1-alpha, Zn)
        #bnminus <- bns[1]
        #bnplus <- bns[2]
        #bhat <- as.numeric(quantile(x=Zn, probs=1-alpha))

        #pminus <- mean(Zn <= bnminus)
        #p0 <- mean(Z0 <= bhat)

        #bound[j] <- 1-alpha+p0-pminus
	results[[as.character(n)]] <- replicates
        } 
}

saveRDS(results, file='/pscratch/sd/j/jbbutler/lemma_bounds/no_bootstrap/bivt.RData')
