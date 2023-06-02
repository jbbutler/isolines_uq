# Helper functions to project bounding curves of uncertainty regions into the tails
# Mostly taken from source code of Cooley et al. (2019)
# 
# Jimmy Butler

source('/global/u1/j/jbbutler/isolines_uq/scripts/R/auxiliary_scripts/orig_isolines.R')

projectContour <- function(contour, dat, base_p, proj_p, mar1Prop, mar2Prop, 
			   mar1Width, mar2Width, asympIndep, etaProp, beta) {

    gpdOut1 <- gpd.fit(dat[,1], threshold = quantile(dat[,1], probs = 1 - mar1Prop), show = F)
    gpdOut2 <- gpd.fit(dat[,2], threshold = quantile(dat[,2], probs = 1 - mar2Prop), show = F)

    contourFrec <- cbind(-1/log(transMar(contour[,1], dat, 1, mar1Prop, mar1Width, gpdOut1)),
                                                     -1/log(transMar(contour[,2], dat, 2, mar2Prop, mar2Width, gpdOut2)))
    
    if (asympIndep) {
        if(is.null(beta)) {
            beta <- 200
        }
        u <- f <- matrix(ncol = ncol(dat), nrow = nrow(dat))
        u[,1] <- transMar(dat[,1], dat, 1, mar1Prop, mar1Width, gpdOut1)
        u[,2] <- transMar(dat[,2], dat, 2, mar2Prop, mar2Width, gpdOut2)

        f[,1] <- -1/log(u[,1])
        f[,2] <- -1/log(u[,2])

        #estimate eta
        minVec <- apply(f, 1, min)
        sl <- sort.list(minVec, decreasing = T)
        k <- round(nrow(dat) * etaProp)
        keep <- sl[1:k]
        etaHat <- sum( log(minVec[keep]) - log(minVec[keep[k]]) )/k

        contourAng <- contourFrec[,1]/(contourFrec[,1] + contourFrec[,2])

        #smooth eta
        weight1 <- 1 - contourAng^beta
        smoothEta1 <- weight1*etaHat + (1 - weight1)*1
        weight2 <- 1 - (1 - contourAng)^beta
        smoothEta2 <- weight2*etaHat + (1 - weight2)*1

        #project to desired level
        multiplier <- cbind((base_p/proj_p)^smoothEta1, (base_p/proj_p)^smoothEta2)
        projcontourFrec <- contourFrec*multiplier
    }

    else {
    projcontourFrec <- contourFrec*(base_p/proj_p)
    }

    projcontourU <- exp(-1/projcontourFrec)
    projcontour <- data.frame(cbind(invTransMar(projcontourU[,1], dat, 1, mar1Prop, mar1Width, gpdOut1),
                                                invTransMar(projcontourU[,2], dat, 2, mar2Prop, mar2Width, gpdOut2) ))

    return(projcontour)
}
