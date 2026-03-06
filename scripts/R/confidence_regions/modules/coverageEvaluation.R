# Module script with functions to evaluate coverage of isolines
#
# Jimmy Butler
# January 2026

evaluateCoverage <- function(tube_obj, true_iso, extreme=FALSE) {
    # Function to evaluate whether a supplied true isoline is
    # coverged by a confidence tube passed in.
    #
    # Inputs:
    #     tube_obj: a list containing components to define the tube (outputs of computeExtremeRegion
    #         and computeEmpiricalRegion)
    #     true_iso: a data.frame of points on the true_isoline you wish to cover
    #     extreme: whether we are evaluating coverage for a tube drawn using extreme method. If so, we will need
    #         to transform the isoline and evaluate coverage in the transformed space.
    #
    # Outputs:
    #     whether the isoline is covered (a list for each alpha supplied)

    c_estimates <- tube_obj$c_estimates
    alphas <- names(tube_obj$c_estimates)
    p <- tube_obj$p
    survFunc <- tube_obj$survFunc

    # if evaluating coverage for an extreme tube, need to transform isoline
    if (extreme) {
        true_iso_X1 <- tube_obj$transform_func1(true_iso[,1])
        true_iso_X2 <- tube_obj$transform_func2(true_iso[,2])
        true_iso <- data.frame(X1=true_iso_X1, X2=true_iso_X2)
    }

    # estimate survival function at points of true isoline
    est_survfunc <- survFunc(true_iso)
    
    isCovereds <- list()

    # evaluate whether estimates are within threshold for each alpha
    for (alpha in alphas) {
        c_estimate <- c_estimates[[alpha]]
        isCovered <- all((est_survfunc <= p + c_estimate) & (est_survfunc >= p - c_estimate))
        isCovereds[alpha] <- isCovered
    }

    return(isCovereds)
}