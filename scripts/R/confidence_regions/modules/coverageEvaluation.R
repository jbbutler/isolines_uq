# Module script with functions to evaluate coverage of isolines
#
# Jimmy Butler
# January 2026

evaluateCoverage <- function(tube_obj, true_iso) {
    # Function to evaluate whether a supplied true isoline is
    # coverged by a confidence tube passed in.
    #
    # Inputs:
    #     tube_obj: a list containing components to define the tube (outputs of computeExtremeRegion
    #         and computeEmpiricalRegion)
    #     true_iso: a data.frame of points on the true_isoline you wish to cover
    #
    # Outputs:
    #     whether the isoline is covered (a list for each alpha supplied)

    c_estimates <- tube_obj$c_estimates
    alphas <- names(tube_obj$c_estimates)
    p <- tube_obj$p
    survFunc <- tube_obj$survFunc

    # estimate survival function at points of true isoline
    est_survfunc <- apply(true_iso, 1, survFunc)
    
    isCovereds <- list()

    # evaluate whether estimates are within threshold for each alpha
    for (alpha in alphas) {
        c_estimate <- c_estimates[[alpha]]
        isCovered <- all((est_survfunc <= p + c_estimate) & (est_survfunc >= p - c_estimate))
        isCovereds[alpha] <- isCovered
    }

    return(isCovereds)
}