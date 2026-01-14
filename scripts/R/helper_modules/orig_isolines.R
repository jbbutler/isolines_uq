# Functins taken from Dan Cooley's R repo for this project

library(MASS)
library(ismev)
library(evd)
library(chron)
library(ks)


#############
# secondary functions
#############

#functions to transform marginals (nonparametric below threshold, gpd above, blended)
#note:  blending is not "consistent" for both trans and invTrans--should it be?
#no guarantee that blending will result in non-decreasing cdf
transMar <- function(x, dat, marginal = 1, marProp, marWidth, gpdOut)
{
	n <- length(dat[,marginal])
	u <- (rank(dat[,marginal])-1)/(n)
	empirical <- approx(x=dat[,marginal], y=u, xout = x)$y
	tooLow <- which(x <= min(dat[,marginal]))
	tooHigh <- which(x > max(dat[,marginal]))
	empirical[tooLow] <- 0
	empirical[tooHigh] <- 1
	pgpdOut <- evd::pgpd(x, loc = gpdOut$threshold, scale = gpdOut$mle[1], shape = gpdOut$mle[2])
	gpdModel <- (1 - marProp) + marProp * pgpdOut
	width <- diff(quantile(dat[,marginal], probs = c(1 - marProp, 1 - marProp + marWidth)))
	weight <- numeric(length(x))
	below <- which(x < gpdOut$threshold)
	above <- which(x > gpdOut$threshold + width)
	between <- which(x >= gpdOut$threshold & x <= gpdOut$threshold+width)
	weight[below] <- 0
	weight[above] <- 1
	weight[between] <- (sin(pi/width*(x[between] - gpdOut$threshold) - pi/2) + 1)/2
	return(weight*gpdModel + (1 - weight)*empirical)
}

invTransMar <- function(u, dat, marginal = 1, marProp, marWidth, gpdOut)
{
	if(any(u < 0) | any(u >= 1)){stop("u must be =< 0 and > 1")}
	n <- length(dat[,marginal])
	uDat <- (rank(dat[,marginal])-1)/(n)
	empirical <- approx(x = uDat, y=dat[,marginal], xout = u)$y
	tooLow <- which(u == 0)
	tooHigh <- which(u > max(uDat))
	empirical[tooLow] <- min(dat[,marginal])
	empirical[tooHigh] <- max(dat[,marginal])
	below <- which(u <= 1 - marProp)
	above <- which(u > 1 - marProp + marWidth)
	between <- which(u > 1 - marProp & u <= 1 - marProp + marWidth)
	qgpdOut <- numeric(length(u))

	# bug I fixed..
	inds <- seq(1, length(u))

	qgpdOut[setdiff(inds, below)] <- evd::qgpd((u[setdiff(inds, below)]-(1 - marProp))/marProp, 
			loc = gpdOut$threshold, scale = gpdOut$mle[1], shape = gpdOut$mle[2])
	qgpdOut[below] <- 0
	weight <- numeric(length(u))
	weight[below] <- 0
	weight[above] <- 1
	weight[between] <- (u[between] - (1 - marProp))/marWidth
	weighted <- weight*qgpdOut + (1 - weight)*empirical
	return(weighted)
}

#kernel survival est function
kernSurv <- function(loc, dat, bw)
{
	p1 <- 1 - pnorm(loc[1], mean = dat[,1], sd = bw[1]/4)
	p2 <- 1 - pnorm(loc[2], mean = dat[,2], sd = bw[2]/4)
	return(mean(p1*p2))
}

#function to actually estimate the cdf
estCDF <- function(dat, bw, gridSize, faster)
{
	surv_lst <- list()
	if (faster) {
		H <- (diag(bw)/4)**2
		res <- kcde(dat, H = H, gridsize = gridSize, tail.flag = 'upper')
		surv <- as.vector(res$estimate)
		xGrid <- res$eval.points[[1]]
		yGrid <- res$eval.points[[2]]
	}

	else {
		xGrid <- seq(min(dat[,1]), max(dat[,1]), l = gridSize)
        	yGrid <- seq(min(dat[,2]), max(dat[,2]), l = gridSize)
        	grid <- expand.grid(xGrid, yGrid)
        	surv <- apply(grid, 1, kernSurv, dat = dat, bw = bw)
	}

	surv_lst$xGrid <- xGrid
	surv_lst$yGrid <- yGrid
	surv_lst$surv <- surv

	return(surv_lst)

}




#################
# primary function
#################

xContours <- function(dat, mar1Prop = .03, mar2Prop = .03, mar1Width = .01, mar2Width = .01,
		baseContourLevel = .01, projContourLevels = c(.005, .001, .0005), bw = NULL,
		gridSize = 100, asympIndep = F, etaProp = .02, beta = NULL, faster = FALSE)
{

	#mar1Prop is proportion of data above threshold for GPD fitting of marginal 1
	#mar1Width is proportion of data over which to smooth GPD and empirical fit (must be < mar1Prop)
	#baseContourLevel level for nonparametric estimation of base contour
	#projContourLevels levels for projection of higher contours
	#bw bandwidth for kernel density survival function (2 dimensional)
	#gridSize number of grid cells used in evaluating smoothed survival function
	#asympIndep tag for asymptotic independence
	#etaProp proportion of data above threshold used for estimating eta
#beta eta-smoothing tuning parameter
	#faster is whether or not to use an FFT-based computation of the kernel survival function

	if(is.null(bw))
	{
		bw <- c(bandwidth.nrd(dat[,1]), bandwidth.nrd(dat[,2]))  #sets bandwidth by kde2d
	}

	n <- dim(dat)[1]

	#estimates tail
        start_gpdfit <- proc.time()[[1]]
	gpdOut1 <- ismev::gpd.fit(dat[,1],
					threshold = quantile(dat[,1], probs = 1 - mar1Prop), show = F)
	gpdOut2 <- ismev::gpd.fit(dat[,2],
					threshold = quantile(dat[,2], probs = 1 - mar2Prop), show = F)
	elapsed_gpdfit <- proc.time()[[1]] - start_gpdfit

	#kernel cdf estimation
	start_cdfest <- proc.time()[[1]]
	surv_lst <- estCDF(dat, bw, gridSize, faster)
	elapsed_cdfest <- proc.time()[[1]] - start_cdfest
	surv <- surv_lst$surv
	xGrid <- surv_lst$xGrid
	yGrid <- surv_lst$yGrid
	survMtx <- matrix(surv, nrow = gridSize, ncol = gridSize)
	start_contour <- proc.time()[[1]]
	contourList <- contourLines(x = xGrid, y = yGrid, z = survMtx, levels = baseContourLevel)
	elapsed_contour <- proc.time()[[1]] - start_contour
	start_cbinding <- proc.time()[[1]]
	contourOrig <- cbind(contourList[[1]]$x, contourList[[1]]$y)
	elapsed_cbinding <- proc.time()[[1]] - start_cbinding
	start_martrans <- proc.time()[[1]]
	contourFrec <- cbind(	-1/log(transMar(contourOrig[,1], dat, 1, mar1Prop, mar1Width, gpdOut1)),
						     -1/log(transMar(contourOrig[,2], dat, 2, mar2Prop, mar2Width, gpdOut2)))
	contourAng <- contourFrec[,1]/(contourFrec[,1] + contourFrec[,2])
	elapsed_martrans <- proc.time()[[1]] - start_martrans

        u <- f <- matrix(ncol = ncol(dat), nrow = nrow(dat))
        u[,1] <- transMar(dat[,1], dat, 1, mar1Prop, mar1Width, gpdOut1)
        u[,2] <- transMar(dat[,2], dat, 2, mar2Prop, mar2Width, gpdOut2)
        f[,1] <- -1/log(u[,1])
        f[,2] <- -1/log(u[,2])

	if(asympIndep)
	{

		if(is.null(beta))
		{
			beta <- 200
		}

	        u <- f <- matrix(ncol = ncol(dat), nrow = nrow(dat))
	        start_martrans <- proc.time()[[1]]
		u[,1] <- transMar(dat[,1], dat, 1, mar1Prop, mar1Width, gpdOut1)
		u[,2] <- transMar(dat[,2], dat, 2, mar2Prop, mar2Width, gpdOut2)
		elapsed_martrans <- proc.time()[[1]] - start_martrans
		f[,1] <- -1/log(u[,1])
		f[,2] <- -1/log(u[,2])

		#estimate eta
		minVec <- apply(f, 1, min)
		sl <- sort.list(minVec, decreasing = T)
		k <- round(n * etaProp)
		keep <- sl[1:k]
		etaHat <- sum( log(minVec[keep]) - log(minVec[keep[k]]) )/k

		#smooth eta
		weight1 <- 1 - contourAng^beta
		smoothEta1 <- weight1*etaHat + (1 - weight1)*1
		weight2 <- 1 - (1 - contourAng)^beta
		smoothEta2 <- weight2*etaHat + (1 - weight2)*1

		#project to desired levels
		projContoursFrec <- list()
		for(i in seq(1, length(projContourLevels)))
		{
			t <- baseContourLevel/projContourLevels[i]
			multiplier <- cbind(t^smoothEta1, t^smoothEta2)
			projContoursFrec[[i]] <- contourFrec * multiplier
		}
	}
	else
	{
		projContoursFrec <- list()
		for(i in seq(1, length(projContourLevels)))
		{
			t <- baseContourLevel/projContourLevels[i]
			projContoursFrec[[i]] <- contourFrec * t
		}
	}

	#transform to uniform[0,1] scale

	projContoursU <- list()
	for(i in seq(1, length(projContourLevels)))
	{
		projContoursU[[i]] <- exp(-1/projContoursFrec[[i]])
	}

	#transform to original scale
	projContours <- list()
	for(i in seq(1, length(projContourLevels)))
	{
		projContours[[i]] <- cbind(
						invTransMar(projContoursU[[i]][,1], dat, 1, mar1Prop, mar1Width, gpdOut1),
						invTransMar(projContoursU[[i]][,2], dat, 2, mar2Prop, mar2Width, gpdOut2) )
	}

	out <- list()

	out$elapsed_contour <- elapsed_contour
	out$elapsed_cdfest <- elapsed_cdfest
	out$elapsed_martrans <- elapsed_martrans
	out$elapsed_gpdfit <- elapsed_gpdfit
	out$elapsed_cbinding <- elapsed_cbinding

	out$data <- dat
	out$f <- f  #frechet transformed data
	out$contourOrig <- contourOrig
	out$projContours <- projContours
	out$contourFrec <- contourFrec
	out$projContoursFrec <- projContoursFrec
	out$gpdOut1 <- gpdOut1
	out$gpdOut2 <- gpdOut2
	if(asympIndep)
	{
		out$minVec <- minVec
		out$etaHat <- etaHat
		out$k <- k
	}
	out$setup <- list()
	out$setup$mar1Prop <- mar1Prop
	out$setup$mar2Prop <- mar2Prop
	out$setup$mar1Width <- mar1Width
	out$setup$mar2Width <- mar2Width
	out$setup$baseContourLevel <- baseContourLevel
	out$setup$projContourLevels <- projContourLevels
	out$setup$bw <- bw
	out$setup$gridSize = gridSize
	out$setup$asympIndep <- asympIndep
	out$setup$etaProp <- etaProp
	out$setup$beta <- beta

	return(out)

}  #end of primary function


###############
# plotting functions
###############


plotOrig <- function(out, xlab = NULL, ylab = NULL, cex = NULL, col = NULL, lwd = NULL,
				xlim = NULL, ylim = NULL)
{
	if(is.null(col)){col <- 1}
	if(is.null(xlab)){xlab <- colnames(out$data[1])}
	if(is.null(ylab)){ylab <- colnames(out$data[2])}
	plot(out$data, xlab = xlab, ylab = ylab, cex = cex, col = col)
	lines(out$contourOrig, col = 2, lwd = lwd)
	for(i in seq(1, length(out$projContours)))
	{
		lines(out$projContours[[i]], col = i + 2, lwd = lwd)
	}
}

plotFrec <- function(out, xlim = c(0,50), ylim = c(0,50), xlab = NULL, ylab = NULL, cex = NULL, col = NULL, lwd = NULL)
{
	if(is.null(col)){col <- 1}
	if(is.null(xlab)){xlab <- colnames(out$data[1])}
	if(is.null(ylab)){ylab <- colnames(out$data[2])}
	plot(out$f, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, cex = cex, col = col)
	lines(out$contourFrec, col = 2, lwd = lwd)
	for(i in seq(1, length(out$projContoursFrec)))
	{
		lines(out$projContoursFrec[[i]], col = i + 2, lwd = lwd)
	}
}

etaDiag <- function(out, end = NULL)
{
	#requires library
	n <- dim(out$data)[1]
	hillOut <- evir::hill(out$minVec, end = n/10, option = "xi")
	abline(v = out$k, col = 3, lty = 2)
	abline(h = out$etaHat, col = 3, lty = 2)
}

marginalDiag <- function(out)
{
	#not sure about this one.  It would seem like we want to make some sort of plot which
	#assesses how well we have modeled each marginal (for transforming), but this is not so
	#useful
	empirical1 <- sort(out$data[,1])
	n <- dim(out$data)[1]
	model1 <- invTransMar(seq(1:n)/(n+1), out$data, 1, out$setup$mar1Prop, out$setup$mar1Width,
					out$gpdOut1)
	plot(model1, empirical1)
}

betaDiag <- function(out, level = 1, ylab = "Empirical Probabilities")
{
	p <- out$setup$projContourLevels[[level]]
	n <- dim(out$data)[1]
	upLevel <- ceiling(n*(p + 3*sqrt(p*(1-p)/n)))
	probVec <- dbinom(seq(0, upLevel), size = n, prob = p)
	sl <- sort.list(probVec, decreasing = T)
	num2Keep <- sum(cumsum(probVec[sl]) < .95) + 1
	keepers <- sl[1:num2Keep]
	minLine <- (min(keepers)-1)/n
	maxLine <- (max(keepers)-1)/n

	numContourPts <- dim(out$contourOrig)[1]
	empProbs <- numeric(numContourPts)
	for(i in seq(1, numContourPts))
	{
		empProbs[i] <- mean( out$projContours[[level]][i,1] < out$data[,1] &
							 out$projContours[[level]][i,2] < out$data[,2] )
	}
	yMax <- max(c(maxLine, empProbs))
	par(mar = c(2,4,2,2))
	plot(empProbs, ylim = c(0, yMax), ylab = ylab, xlab = "", axes = F)
	axis(1, labels = F)
	axis(2)
	box()
	abline(h = out$setup$projContourLevels[level], col = 2)
	abline(h = maxLine, lty = 2, col = 2)
	abline(h = minLine, lty = 2, col = 2)

}

bootstrapBvpot <- function(outObj, blockSize, iter, bootContourLevel)
{
	bootContour <- list()

	for(k in 1:iter)
	{
		print(k)
	  	#obtain bootstrap sample
	  	n <- dim(outObj$dat)[1]
	  	sampLocs <- sample(n, size = ceiling(n/blockSize), replace = T)
	  	sampVals <- NULL
	  	for(i in 0:(blockSize-1))
	  	{
	  		sampVals <- rbind(sampVals, outObj$dat[sampLocs+0,])
	  	}
		sampVals <- sampVals[1:n,]

		bootContour[[k]] <- xContours(dat = sampVals,
				mar1Prop = outObj$setup$mar1Prop, mar2Prop = outObj$setup$mar2Prop,
				mar1Width = outObj$setup$mar1Width, mar2Width = outObj$setup$mar2Width,
				baseContourLevel = outObj$setup$baseContourLevel,
				projContourLevel = outObj$setup$projContourLevels,
				bw = outObj$setup$bw,
				gridSize = outObj$setup$gridSize,
				asympIndep = outObj$setup$asympIndep,
				etaProp = outObj$setup$etaProp,
				beta = outObj$setup$beta)
	}

	return(bootContour)
}
