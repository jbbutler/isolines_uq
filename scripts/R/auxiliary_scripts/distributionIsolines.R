# Script with functions to draw isolines of distributions used in simulation analyses
#
# Jimmy Butler

library(mvtnorm)
source('~/isolines_uq/scripts/R/auxiliary_scripts/karachiTools.R')

drawBivtIsoline <- function(numCoords, gridUbs, gridLbs, prob) {
    # Draw an isoline of a bivariate t distribution with a particular covariance matrix (unit variances, 0.7 covariance)
    # and 4 degrees of freedom
    #
    # Inputs:
    #    + numCoords (integer): the number of points you want in your isoline
    #    + gridUbs (2-vector of floats): the upper bounds in x and y respectively of the window containing the isoline
    #    + gridLbs (2-vector of floats): same as gridUbs, but for lower bounds
    #    + prob (float): the desired isoline's exceedance probability
    #
    # Returns:
    #    + A data.frame of points on the desired isoline (coordinate axis labels X1, X2)

    # instantiate function to compute difference between exceedance probability past a point and desired exceedance prob.
    # note: origin in polar coordinates is given by gridLbs, point is specified in polar coordinates in this system	
    pmvtFixedCoord <- function(radius, angle, sigma, df, prob, xCenter, yCenter) {
        xCoord <- xCenter + radius*cos(angle)
        yCoord <- yCenter + radius*sin(angle)
        exceedanceProb <- pmvt(lower=c(xCoord, yCoord), upper=c(Inf, Inf), df=df, sigma=sigma)
        return(exceedanceProb - prob)
    }

    radii <- rep(NA, numCoords)
    angles <- seq(0, pi/2, length.out=numCoords)
    maxRad <- sqrt(sum((gridUbs-gridLbs)**2))

    # for each angle in first quadrant, find radius that gives a point with desired exceedance probability
    # by finding roots of pmvtFixedCoord given the angle
    for (i in 1:numCoords) {
        angle <- angles[i]
        radii[i] <- uniroot(pmvtFixedCoord, interval=c(0, maxRad), angle=angle,
                             sigma=matrix(c(1, 0.7, 0.7, 1), nrow=2), df=4, prob=prob, xCenter=gridLbs[1], yCenter=gridLbs[2])$root
    }

    # convert back to cartesian coordinates
    xs <- radii*cos(angles) + gridLbs[1]
    ys <- radii*sin(angles) + gridLbs[2]

    return(data.frame(X1=xs, X2=ys))

}

drawBivGaussIsoline <- function(numCoords, gridUbs, gridLbs, prob) {
    # Draw an isoline of a zero-mean bivariate normal distribution with covariance matrix with unit variances and
    # covariance of 0.7
    #
    # Inputs:
    #    + numCoords (integer): the number of points you want in your isoline
    #    + gridUbs (2-vector of floats): the upper bounds in x and y respectively of the window containing the isoline
    #    + gridLbs (2-vector of floats): same as gridUbs, but for lower bounds
    #    + prob (float): the desired isoline's exceedance probability
    #
    # Returns:
    #    + A data.frame of points on the desired isoline (coordinate axis labels X1, X2)

    # instantiate function to compute difference between exceedance probability past a point and desired exceedance prob.
    # note: origin in polar coordinates is given by gridLbs, point is specified in polar coordinates in this system
    pmvnormFixedCoord <- function(radius, angle, sigma, mean, prob, xCenter, yCenter) {
        xCoord <- xCenter + radius*cos(angle)
        yCoord <- yCenter + radius*sin(angle)
        exceedanceProb <- pmvnorm(lower=c(xCoord, yCoord), upper=c(Inf, Inf), mean=mean, sigma=sigma)
        return(exceedanceProb - prob)
    }

    radii <- rep(NA, numCoords)
    angles <- seq(0, pi/2, length.out=numCoords)
    maxRad <- sqrt(sum((gridUbs-gridLbs)**2))

    # for each angle in first quadrant, find radius that gives a point with desired exceedance probability
    # by finding roots of pmvtFixedCoord given the angle
    for (i in 1:numCoords) {
        angle <- angles[i]
        radii[i] <- uniroot(pmvnormFixedCoord, interval=c(0, maxRad), angle=angle, mean=rep(0,2),
                             sigma=matrix(c(1, 0.7, 0.7, 1), nrow=2), prob=prob, xCenter=gridLbs[1], yCenter=gridLbs[2])$root
    }

    # convert back to cartesian coordinates
    xs <- radii*cos(angles) + gridLbs[1]
    ys <- radii*sin(angles) + gridLbs[2]

    return(data.frame(X1=xs, X2=ys))

}

drawKarachiIsoline <- function(numCoords, gridUbs, gridLbs, prob) {
    # Draw an isoline of the fitted beta-kernel Karachi data distribution using data from Cooley et al. (2019) 
    # with bandwidth parameter chosen via 5-fold CV in runKarachiCV.R 
    #
    # Inputs:
    #    + numCoords (integer): the number of points you want in your isoline
    #    + gridUbs (2-vector of floats): the upper bounds in x and y respectively of the window containing the isoline
    #    + gridLbs (2-vector of floats): same as gridUbs, but for lower bounds
    #    + prob (float): the desired isoline's exceedance probability
    #
    # Returns:
    #    + A data.frame of points on the desired isoline (coordinate axis labels X1, X2)

    # instantiate function to compute difference between exceedance probability past a point and desired exceedance prob.
    # note: origin in polar coordinates is given by gridLbs, point is specified in polar coordinates in this system 
    pKarachiFixedCoord <- function(radius, angle, prob, xCenter, yCenter) {
        xCoord <- xCenter + radius*cos(angle)
        yCoord <- yCenter + radius*sin(angle)
	# optimal bandwidth chosen via 5-fold CV
        exceedanceProb <- pKarachiBetaKDE(point=c(xCoord, yCoord),lbs=c(50, 0), ubs=c(140, 100), b=0.00073)
        return(exceedanceProb - prob)
    }

    radii <- rep(NA, numCoords)
    angles <- seq(0, pi/2, length.out=numCoords)
    maxRad <- sqrt(sum((gridUbs-gridLbs)**2))

    # for each angle in first quadrant, find radius that gives a point with desired exceedance probability
    # by finding roots of pmvtFixedCoord given the angle
    for (i in 1:numCoords) {
        angle <- angles[i]
        radii[i] <- uniroot(pKarachiFixedCoord, interval=c(0, maxRad), angle=angle,
                            prob=prob, xCenter=gridLbs[1], yCenter=gridLbs[2])$root
    }

    # convert back to cartesian coordinates
    xs <- radii*cos(angles) + gridLbs[1]
    ys <- radii*sin(angles) + gridLbs[2]

    return(data.frame(X1=xs, X2=ys))

}

