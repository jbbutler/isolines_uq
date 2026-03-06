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
        tube_lower <- drawEmpiricalIsoline(dat=dat, 
                              n_coords=n_coords, 
                              grid_lbs=lbs, 
                              p=c_estimates[[alpha]] + p)

        if (p-c_estimates[[alpha]] > 0) { 
            tube_upper <- drawEmpiricalIsoline(dat=dat, 
                              n_coords=n_coords, 
                              grid_lbs=lbs, 
                              p= -c_estimates[[alpha]] + p)
        } else {

            tube_upper <- NA

        }
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
    transformed_dat_X1 <- tube_obj$transform_func1(dat[,1])
    transformed_dat_X2 <- tube_obj$transform_func2(dat[,2])
    trans_dat <- data.frame(X1=transformed_dat_X1, X2=transformed_dat_X2)
    p <- tube_obj$p
    gamma <- tube_obj$gamma
    xi <- tube_obj$xi
    
    alphas <- names(tube_obj$c_estimates)
    c_estimates <- tube_obj$c_estimates
    
    if (!(alpha %in% alphas)) {
        stop('No confidence tube for desired alpha.')
    } else {
        c_est <- c_estimates[[alpha]]
        tube_lower <- drawExtremeIsoline(dat=trans_dat, 
                              p=p + c_est,
                              n_coords=n_coords,
                              grid_lbs=lbs, 
                              gamma=gamma,
                              xi=xi)

        if (p-c_est > 0) {
            tube_upper <- drawExtremeIsoline(dat=trans_dat, 
                              p=p-c_est,
                              n_coords=n_coords,
                              grid_lbs=lbs, 
                              gamma=gamma,
                              xi=xi)
        } else {
            tube_upper <- NA

        }
    }

    # inverse transform marginals of the lower bound of the tube
    tube_lower <- data.frame(X1=tube_obj$inv_transform_func1(tube_lower[,1]),
                             X2=tube_obj$inv_transform_func2(tube_lower[,2]))

    if (!any(is.na(tube_upper))) {
        tube_upper <- data.frame(X1=tube_obj$inv_transform_func1(tube_upper[,1]),
                             X2=tube_obj$inv_transform_func2(tube_upper[,2]))
    }
    
    lst <- list()
    lst$tube_upper <- tube_upper
    lst$tube_lower <- tube_lower

    return(lst)

}