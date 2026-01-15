# A module script with functions that aid in the visualization of the confidence tubes.
#
# Jimmy Butler
# January 2026

drawEmpiricalTubes <- function(tube_obj, alpha, lbs=c(0,0), n_coords=200) {
    # Function to draw bounding curves for confidence tubes using empirical
    # survival function.
    #
    # Inputs:
    #     tube_obj: the list containing pieces to define tube (output of computeExtremeRegion,
    #         computeEmpiricalRegion, etc.)
    #     alpha: which alpha confidence tube you'd like to display
    #     lbs: the lower lefthand corner of the space to draw tube (origin of polar coord system)
    #     n_coords: the number of points you'd like to display on either curve
    #
    # Outputs:
    #     A list containing data.frames of points on either bounding curve

    alpha <- as.character(alpha)
    dat <- tube_obj$dat
    p <- tube_obj$p
    alphas <- names(tube_obj$c_estimates)
    c_estimates <- tube_obj$c_estimates
    
    if (!(alpha %in% alphas)) {
        stop('No confidence tube for desired alpha.')
    } else {
        c_est <- c_estimates[[alpha]]
        tube_upper <- drawEmpiricalIsoline(dat=dat, 
                              n_coords=n_coords, 
                              gridLbs=lbs, 
                              p=c_estimates[[alpha]] + p)
        tube_lower <- drawEmpiricalIsoline(dat=dat, 
                              n_coords=n_coords, 
                              gridLbs=lbs, 
                              p= -c_estimates[[alpha]] + p)
    }
    lst <- list()
    lst$tube_upper <- tube_upper
    lst$tube_lower <- tube_lower

    return(lst)
}

drawExtremeTubes <- function(tube_obj, alpha, lbs=c(0,0), n_coords=200) {
    # Function to draw bounding curves for confidence tubes using the extreme
    # survival function.
    #
    # Inputs:
    #     tube_obj: the list containing pieces to define tube (output of computeExtremeRegion,
    #         computeEmpiricalRegion, etc.)
    #     alpha: which alpha confidence tube you'd like to display
    #     lbs: the lower lefthand corner of the space to draw tube (origin of polar coord system)
    #     n_coords: the number of points you'd like to display on either curve
    #
    # Outputs:
    #     A list containing data.frames of points on either bounding curve

    alpha <- as.character(alpha)
    dat <- tube_obj$dat
    p <- tube_obj$p
    gamma <- tube_obj$gamma
    xi <- tube_obj$xi
    
    alphas <- names(tube_obj$c_estimates)
    c_estimates <- tube_obj$c_estimates
    
    if (!(alpha %in% alphas)) {
        stop('No confidence tube for desired alpha.')
    } else {
        c_est <- c_estimates[[alpha]]
        tube_upper <- drawExtremeIsoline(dat=dat, 
                              p=p + c_est,
                              n_coords=n_coords,
                              gridLbs=lbs, 
                              gamma=gamma,
                              xi=xi)

        tube_lower <- drawExtremeIsoline(dat=dat, 
                              p= p - c_est,
                              n_coords=n_coords,
                              gridLbs=lbs, 
                              gamma=gamma,
                              xi=xi)
    }
    lst <- list()
    lst$tube_upper <- tube_upper
    lst$tube_lower <- tube_lower

    return(lst)

}