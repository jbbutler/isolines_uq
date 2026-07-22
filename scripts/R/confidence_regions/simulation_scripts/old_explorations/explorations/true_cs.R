source('~/isolines_uq/scripts/R/confidence_regions/modules/distributionIsolines.R')
source('~/isolines_uq/scripts/R/confidence_regions/modules/utils.R')

sampling_func <- loadSamplingFunction(dist='bivt')

ns <- 3*c(1000, 5000, 10000, 50000)
C <- 15
pn <- function(n) C/n

for (n in ns) {

    

}



n_monte_carlo <- 10000
Wplus <- rep(NA, n_monte_carlo)
Wminus <- rep(NA, n_monte_carlo)