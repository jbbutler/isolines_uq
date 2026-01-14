# File to run a CV to choose the bandwidth parameter in fitting a beta KDE
# to the Karachi data

# CHANGED THIS FILE: RERUN TO ENSURE NO BUGS

library(ggplot2)
library(dplyr)
library(caret)
library(foreach)
library(doSNOW)
library(parallel)

source('~/isolines_uq/scripts/R/auxiliary_scripts/karachiTools.R')

load('~/isolines_uq/data/dans_data/karachiDatDaily.Rdata')
karachiDat <- karachiDatDaily %>% select(relHum, temp)

lower_bw <- 0.01
upper_bw <- 0.0001
n_bws <- 1000
bws <- seq(lower_bw, upper_bw, length.out=n_bws)

lbs <- c(0, 50)
ubs <- c(100, 140)
nfolds <- 5
n_cores <- 15

parallelizedCode <- function(ind){

    b <- bws[ind]
    inds <- seq(1, nrow(karachiDat))
    flds <- createFolds(inds, k = 5, list=TRUE, returnTrain=FALSE)
    loglikelihoods <- rep(NA, length(flds))

    for (j in 1:length(flds)) {

        train <- karachiDat[Reduce(c, flds[-j]),]
        test <- karachiDat[flds[[j]],]
        
        train <- data.frame(X=(train[,1]-lbs[1])/(ubs[1]-lbs[1]), 
                      Y=(train[,2]-lbs[2])/(ubs[2]-lbs[2]))
    
        loglikelihoods[j] <- sum(log(Reduce(c, apply(test, 1, dKarachiBetaKDE, 
						     dat=train, 
						     b=b, 
						     lbs=lbs, ubs=ubs))))
    
    }
    
    return(mean(loglikelihoods))
}

clust <- makeSOCKcluster(n_cores)
registerDoSNOW(clust)
pb <- txtProgressBar(min = 1, max = n_bws, style = 3)
progress <- function(n) setTxtProgressBar(pb ,n)
opts <- list(progress = progress)

results <- foreach(i = 1:n_bws,
               .packages = c('caret'),
               .options.snow = opts,
               .combine = 'c') %dopar% {

    parallelizedCode(i)
}

close(pb)
stopCluster(clust)

fname <- paste0('karachiKDECV_lbw', lower_bw, '_ubw', 
		upper_bw, '_nbw', n_bws, '.RData')

saveRDS(results, paste0('~/isolines_uq/results/karachiCVResults/', fname))


