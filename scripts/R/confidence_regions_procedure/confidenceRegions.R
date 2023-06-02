library(ks)
library(dplyr)
library(ismev)
#library(sf, warn.conflicts=FALSE)
#library(sfheaders)
#library(concaveman)

source('/global/u1/j/jbbutler/isolines_uq/scripts/R/auxiliary_scripts/orig_isolines.R')
source('/global/u1/j/jbbutler/isolines_uq/scripts/R/auxiliary_scripts/utils.R')
source('/global/u1/j/jbbutler/isolines_uq/scripts/R/auxiliary_scripts/projection.R')

drawBaseRegions <- function(dat, grid_obj, alphas, ps, B, beta_funcs_dict, emp_surv=NULL, boot_survs=NULL) {
    # Function to draw one or more uncertainty regions for some base isoline of a specified
    # exceedence probability (p). Uses the methods in Mammen and Polonik (2013)
    # and also blends a bit of the code from Cooley.
    #
    # Arguments:
    #
    #    

    # unpacking the grid specifications
    grid <- grid_obj$grid
    ubs <- grid_obj$ubs
    lbs <- grid_obj$lbs

    # whether or not you are estimating the survival function or plugging in
    if (is.null(emp_surv)) {
        surv_func <- fastEmpSurv(grid, dat)
    } else {
    	surv_func <- emp_surv
    }

    # get all combinations of p and beta functions
    beta_func_labs <- names(beta_funcs_dict)
    combos <- expand.grid(beta_func_labs=beta_func_labs, ps=ps)
    p_beta_lst <- vector(mode='list', length=nrow(combos))
    # start list of names
    names <- c()

    # for each p-beta function combo, get the hhat values, the boolean mask for
    # the points in Delta, and preallocate the bootstrap list of Zs
    # (all in the language of Mammen and Polonik (2013))
    for (i in 1:nrow(combos)) {
        p <- combos$ps[i]
        beta_func_lab <- combos$beta_func_labs[i]

	hhat_vals <- -surv_func + p
	# boolean mask to extract points in blown up boundary around est isoline
	deltamask <- abs(hhat_vals) <= beta_funcs_dict[[beta_func_lab]](nrow(dat))
	Zs <- rep(0, B)
	lst <- list(hhat_vals=hhat_vals, deltamask=deltamask, Zs=Zs)
	names <- c(names, paste0('beta', beta_func_lab, '_p', p))
        p_beta_lst[[i]] <- lst
    }
    names(p_beta_lst) <- names

    # for every bootstrap resample
    for (i in 1:B) {
	# whether or not you using pre-made bootstrap survfuncs
        if (is.null(boot_survs)) {
            boot_samp <- dat %>% sample_frac(1, replace = TRUE)
            boot_surv_func <- fastEmpSurv(grid, boot_samp)
	} else {
            boot_surv_func <- boot_survs[[i]]
	}
        
        # for every p-beta combo, get the bootstrap distribution of Z(beta)
        for (j in 1:nrow(combos)) {
	    p <- combos$ps[j]
	    deltamask <- p_beta_lst[[j]]$deltamask
	    hhat_vals <- p_beta_lst[[j]]$hhat_vals 
            boot_hhat_vals <- -boot_surv_func + p

            p_beta_lst[[j]]$Zs[i] <- max(abs((boot_hhat_vals - hhat_vals)[deltamask]))
        }
    }
    # create lists of output confidence regions
    output_lst <- list()
    # labels for the output list
    newnames <- c()
    # index counter variable
    k <- 1

    # for every p-beta combo
    for (i in 1:length(p_beta_lst)) {

        name <- names(p_beta_lst)[i]
        Zs <- p_beta_lst[[i]]$Zs
	hhat_vals <- p_beta_lst[[i]]$hhat_vals

	# for every p-beta-alpha triple, draw confidence regions and store results
	for (j in 1:length(alphas)) {
            alpha <- alphas[j]
	    bhat <- as.numeric(quantile(Zs, probs = 1 - alpha))
            tube_bounds <- drawTubeBounds(dat, bhat, combos$ps[i], lbs, ubs)
            out <- list()
            out$tube_top <- tube_bounds$top
	    out$tube_bottom <- tube_bounds$bottom
            out$alpha <- alpha
            out$data <- dat
            out$p <- combos$ps[i]
            out$B <- B
            out$Zs <- Zs
            out$grid <- grid
            out$bhat <- bhat
            out$beta_func <- combos$beta_func_labs[i]
            out$beta_function <- beta_funcs_dict[[combos$beta_func_labs[i]]]
            
	    newnames <- c(newnames, paste0(name, '_alpha', alpha))
	    output_lst[[k]] <- out
	    # update counter
	    k <- k+1
	}
    }
    # label the list elements for easier bookkeeping
    names(output_lst) <- newnames

    return(output_lst)
}

projectBaseRegion <- function(base_out, proj_p, mar1Prop=0.03, mar2Prop=0.03, mar1Width=0.01, mar2Width=0.01, 
			      asympIndep, etaProp = 0.02, beta=NULL) {

    dat <- base_out$data
    n <- nrow(dat)
    baseTubeTop <- base_out$tube_top
    baseTubeBottom <- base_out$tube_bottom
    base_p <- base_out$p
    
    proj_out <- list()
    proj_out$base_out <- base_out
    proj_out$proj_p <- proj_p
    proj_out$projTubeTop <- NULL

    # if the confidence interval is bounded, project the top; otherwise, make no attempt
    # to project
    if (nrow(baseTubeTop) > 3) {

        projTubeTop <- projectContour(baseTubeTop, dat, base_p, 
				  proj_p, mar1Prop, mar2Prop, mar1Width, 
				  mar2Width, asympIndep, etaProp, beta)
        proj_out$projTubeTop <- projTubeTop
    }

    projTubeBottom <- projectContour(baseTubeBottom, dat, base_p,
                                  proj_p, mar1Prop, mar2Prop, mar1Width,
                                  mar2Width, asympIndep, etaProp, beta)

    proj_out$projTubeBottom <- projTubeBottom

    return(proj_out)
}
