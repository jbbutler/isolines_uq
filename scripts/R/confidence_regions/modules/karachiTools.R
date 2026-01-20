# Functions to work with the Karachi temperature-humidity distribution, built from data from Cooley et al.
# Distribution is a KDE with beta-distributed kernels, according to the first method presented
# in "Beta kernel estimators for density functions" by Song Xi Chen (1999)

# loading the Karachi temperature-humidity data, used to build the distribution
load('~/isolines_uq/data/cooley_data/karachiDatDaily.Rdata')
karachiDat <- karachiDatDaily[c('temp', 'relHum')]

# constants to set the lower and upper bounds of the distribution support
LBS_CONST <- c(50, 0)
UBS_CONST <- c(140, 100)
b_CONST <- 0.00073

dKarachiBetaKDE <- function(point, dat=karachiDat, b=b_CONST, lbs=LBS_CONST, ubs=UBS_CONST) {
    # Find the value of the bivariate density of Karachi data at desired point

    dat <- data.frame(X=(dat[,1]-lbs[1])/(ubs[1]-lbs[1]),
                      Y=(dat[,2]-lbs[2])/(ubs[2]-lbs[2]))

    shape1sX <- dat[,1]/b + 1
    shape2sX <- (1-dat[,1])/b + 1

    shape1sY <- dat[,2]/b + 1
    shape2sY <- (1-dat[,2])/b + 1

    margXs <- dbeta((point[1]-lbs[1])/(ubs[1]-lbs[1]),
                    shape1=shape1sX, shape2=shape2sX)*(1/(ubs[1]-lbs[1]))
    margYs <- dbeta((point[2]-lbs[2])/(ubs[2]-lbs[2]),
                    shape1=shape1sY, shape2=shape2sY)*(1/(ubs[2]-lbs[2]))

    return(mean(margXs*margYs))
}

rKarachiBetaKDE <- function(n, b=b_CONST, dat=karachiDat, ubs=UBS_CONST, lbs=LBS_CONST) {
    # Draw a random sample from the distriubtion of Karachi data
    # First column is temperature, second is relative humidity

    kernel_inds <- sample(1:nrow(dat), n, replace=TRUE)
    dat <- dat[kernel_inds,]
    dat <- data.frame(X=(dat[,1]-lbs[1])/(ubs[1]-lbs[1]), 
                      Y=(dat[,2]-lbs[2])/(ubs[2]-lbs[2]))
    
    shape1sX <- dat[,1]/b + 1
    shape2sX <- (1-dat[,1])/b + 1
    
    shape1sY <- dat[,2]/b + 1
    shape2sY <- (1-dat[,2])/b + 1
    
    X <- rbeta(n, shape1=shape1sX, shape2=shape2sX)*(ubs[1]-lbs[1]) + lbs[1]
    Y <- rbeta(n, shape1=shape1sY, shape2=shape2sY)*(ubs[2]-lbs[2]) + lbs[2]
    
    draws <- data.frame(X1=X, X2=Y)
    
    return(draws)
}

pKarachiBetaKDE <- function(point, dat=karachiDat, b=b_CONST, lbs=LBS_CONST, ubs=UBS_CONST) {
    # Find the probability of jointly exceeding a particular point
    # with a particular set of parameters (b, lbs, ubs)

    dat <- data.frame(X=(dat[,1]-lbs[1])/(ubs[1]-lbs[1]),
                      Y=(dat[,2]-lbs[2])/(ubs[2]-lbs[2]))

    shape1sX <- dat[,1]/b + 1
    shape2sX <- (1-dat[,1])/b + 1

    shape1sY <- dat[,2]/b + 1
    shape2sY <- (1-dat[,2])/b + 1

    margXs <- 1-(pbeta((point[1]-lbs[1])/(ubs[1]-lbs[1]),
                    shape1=shape1sX, shape2=shape2sX))
    margYs <- 1-(pbeta((point[2]-lbs[2])/(ubs[2]-lbs[2]),
                    shape1=shape1sY, shape2=shape2sY))

    return(mean(margXs*margYs))

}
