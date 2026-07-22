library(foreach)
library(doSNOW)
library(parallel)
library(mvtnorm)
library(data.table)
library(SimilarityMeasures)

load_path <- '/pscratch/sd/j/jbbutler/sims/regions/mixed_confregs/bivtcopula_frecmargins/'
save_path <- '/pscratch/sd/j/jbbutler/sims/regions/mixed_coverage/bivtcopula_frecmargins/'
source('~/isolines_uq/scripts/R/auxiliary_scripts/distributionIsolines.R')

n_cores <- 15

p_list <- strsplit(list.files(load_path), '_')
ps <- c()
for (i in 1:length(p_list)) {
    ps <- c(ps, p_list[[i]][1])
}

lst <- vector(mode='list')

total_results <- vector(mode='list', length=length(ps))

for (i in 1:length(ps)){

    isoline <- drawBivtIsoline(300, c(100,100), c(-2, -2), as.numeric(ps[i]))
    iso_xs <- isoline$X1
    iso_ys <- isoline$X2
    frec_isoline <- data.frame(X1=-1/log(pt(iso_xs, df=4)), X2=-1/log(pt(iso_ys, df=4)))
    lst[[ps[i]]] = frec_isoline

    n_list <- list.files(paste0(load_path, ps[i], '_prob'))

    p_results <- vector(mode='list', length=length(n_list))

    for (j in 1:length(n_list)){

        n <- as.numeric(strsplit(n_list, '_')[[j]][1])
        sims_path <- paste0(load_path, ps[i], '_prob/', n_list[j], '/')
        n_iter <- length(list.files(sims_path))


	parallelizedCode <- function(ind) {
    
    		sim_confregs <- readRDS(paste0(sims_path, '/simulation_', ind, '.RData')) 
    		confreg_lst <- names(sim_confregs)
    		res <- vector(mode='list', length=length(confreg_lst))
		dat <- sim_confregs$data
    
    		for (k in 1:(length(confreg_lst)-1)) {
        
        		sim_lst <- list()
        		confreg <- sim_confregs[[k]]
        		sim_lst$n <- n
        		sim_lst$alpha <- confreg$alpha
        		sim_lst$p <- confreg$p
        		sim_lst$beta_func <- confreg$beta_func
        		sim_lst$sim_num <- ind
        
        		#tube_top <- confreg$tube_bounds$top
        		#tube_bottom <- confreg$tube_bounds$bottom
        
        		#comp2xs <- CJ(iso_x=frec_isoline[,1], bott_x=tube_bottom[,1], sorted=FALSE)
        		#comp2ys <- CJ(iso_y=frec_isoline[,2], bott_y=tube_bottom[,2], sorted=FALSE)
        		#lower_cov <- !any((comp2xs[[1]] <= comp2xs[[2]]) & (comp2ys[[1]] <= comp2ys[[2]]))

			#unbounded <- confreg$p <= confreg$bhat

        		#if (unbounded) {
            			#sim_lst$frec_distance <- Inf
            		#}
        
        		#else {
            			#comp1xs <- CJ(top_x=tube_top[,1], iso_x=frec_isoline[,1], sorted=FALSE)
            			#comp1ys <- CJ(top_y=tube_top[,2], iso_y=frec_isoline[,2], sorted=FALSE)
            			#upper_cov <- !any((comp1xs[[1]] <= comp1xs[[2]]) & (comp1ys[[1]] <= comp1ys[[2]]))
            
            			#sim_lst$frec_distance <- as.double(Frechet(as.matrix(tube_top), as.matrix(tube_bottom)))

			evaluateSurvivalFunction <- function(row) {

    				empirical_prob <- mean(dat[,1] > row[[1]] & dat[,2] > row[[2]])
                                q <- 100/nrow(dat)

    				if (empirical_prob >= q) {
        				return(empirical_prob)
    				}

				theta = atan(row[[2]]/row[[1]])
				r = sqrt(row[[1]]**2 + row[[2]]**2)

				# find the radius that correpsonds to this theta on the isoline

				survivalDiff <- function(r, theta, p, dat) {
        				xPt <- r*cos(theta)
        				yPt <- r*sin(theta)
        				actual_survival <- mean((dat[,1] > xPt) & dat[,2] > yPt)

        				return(actual_survival - p)
        			}

    				maxX <- max(dat[,1])
    				maxY <- max(dat[,2])
    				maxRad <- sqrt(maxX**2 + maxY**2)

        			iso_r <- uniroot(survivalDiff, interval=c(0, maxRad), theta=theta, dat=dat, p=q)$root

    				rv_prob <- (iso_r/r)*q
    				return(rv_prob)
			}

			# evaluate if each point on the isoline satisfies conditions of being in confidence set
                        survs <- apply(frec_isoline, 1, evaluateSurvivalFunction)
        		covered <- all((survs <= confreg$p + confreg$bhat) & (survs >= confreg$p - confreg$bhat))

        		sim_lst$covered <- covered
        		res[[k]] <- sim_lst
        	}
    		results_sim <- do.call(rbind.data.frame, res)
    		return(results_sim)
	}


        clust <- makeSOCKcluster(n_cores)
        registerDoSNOW(clust)
        pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
        progress <- function(n) setTxtProgressBar(pb ,n)
        opts <- list(progress = progress)

        nresult <- foreach(m = 1:n_iter,
                   .packages = c('mvtnorm', 'purrr', 'data.table', 'dplyr', 'ismev', 'ks', 'SimilarityMeasures'),
                   .options.snow = opts,
                   .combine = 'rbind') %dopar% {

            parallelizedCode(m)

            }
        p_results[[j]] <- nresult

        close(pb)
        stopCluster(clust)

        }

    p_results <- do.call(rbind.data.frame, p_results)
    total_results[[i]] <- p_results

}

total_results <- do.call(rbind.data.frame, total_results)
saveRDS(total_results, paste0(save_path, 'full_coverage_res.RData'))
