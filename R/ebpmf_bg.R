#' @title Empirical Bayes Poisson Matrix Factorization (Background Model)
#' @import ebpm
#' @import Matrix

#' @param X count matrix (dim(X) = c(n, p)).
#' @param k number of topics
#' @param pm_func functions for solving the \code{ebpm} subproblem for \code{L} and \code{F}; 
#' It is a list \code{list(l, f)};
#' For our purpose we use `mle_pm` or `ebpm_point_gammma`for \code{L}, and \code{ebpm_gamma_mixture} for \code{F}
#' @param pm_control control parameters for pm_func function
#' @param init Either \code{NULL} or \code{list(qg, l0, f0)}
#' @param fix_g list(l, f) where l, f are either TRUE or FALSE
#' @param maxiter maximum number of iterations
#' @param tol stopping tolerance for ELBO
#' @param verbose print progress if set TRUE
#'
#' @return A list containing elements:
#'     \describe{
#'       \item{\code{l0}}{sample-wise mean}
#'       \item{\code{f0}}{feature-wise mean}
#'       \item{\code{qg}}{list(ql, gl,qf, gf)}
#'       \item{\code{ELBO}}{ELBO objective for this VEB algorithm}
#'      }
#' @examples
#' To add
#' @export  ebpmf_bg

ebpmf_bg <- function(X, K, 
										 pm_func = list(f = ebpm::ebpm_gamma_mixture, 
																		l = ebpm::ebpm_gamma_mixture),
										 init = NULL, pm_control = NULL,
										 #fix_g = list(l = FALSE, f = FALSE), 
										 fix_option = list(l0 = FALSE, f0 = FALSE,
																				gl = FALSE, ql = FALSE,
																				gf = FALSE, qf = FALSE), 
										 maxiter = 100, tol = 1e-8, verbose = FALSE){
	## TODO: input check, require X_rs, X_cs to be nonzero

  ## transform to sparse matrix 
  X <- as(X, "sparseMatrix") 
	X_rs = Matrix::rowSums(X)
	X_cs = Matrix::colSums(X)
  d = summary(X)
	const = sum(apply.nonzeros(X = X, f = function(x) lgamma(x + 1)))
  ## initialization
	init_tmp <- init_ebpmf_bg(X = X, K = K, init = init, d = d)
  qg <- init_tmp$qg
  B <- init_tmp$B
	l0 <- init_tmp$l0
	f0 <- init_tmp$f0
  rm(init_tmp)
  ## update iteratively
  ELBOs <- c()
  for(i in 1:maxiter){
    for(k in 1:K){
			## store B_k
 			B_k = exp(qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k]) 
    	## compute q(Z)
      Ez <- compute_EZ(d = d,B = B,B_k = B_k)
      ## update (qL, gL, qF, gF) 
      #### TODO: change initialization`
			init_r1 = list(sf = sum(l0 * qg$qls_mean[,k]),
                     gl = qg$gls[[k]], gf = qg$gfs[[k]])
      rank1_qg <- rank1_bg(d = d, X_rs = Ez$rs, X_cs = Ez$cs,
														l0 = l0, f0 = f0, 
														pm_func = pm_func, pm_control = pm_control,
														ql = list(mean = qg$qls_mean[,k], 
																			mean_log = qg$qls_mean_log[,k]),
														qf = list(mean = qg$qfs_mean[,k], 
																			mean_log = qg$qfs_mean_log[,k]),
													 	gl = qg$gls[[k]],	
													 	gf = qg$gfs[[k]],	
														kl_l = qg$kl_l[k],
														kl_f = qg$kl_f[k],
														fix_option = fix_option)
      rm(Ez)
      qg = update_qg(tmp = rank1_qg, qg = qg, k = k)
      rm(rank1_qg)
			## update B, Lam
			B = B - B_k + exp(qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k])
		 #	B[B < 0] <- 1e-20 ## numerical issue gets - epsilon	
    }
		## update l0, f0
		if(!fix_option$l0){
			denom <- colSums( t(qg$qls_mean) * colSums(f0 * qg$qfs_mean)) 
			l0 <- X_rs/denom
		}
		if(!fix_option$f0){
			denom <- colSums( t(qg$qfs_mean) * colSums(l0 * qg$qls_mean)) 
    	f0 <- X_cs/denom
		}
    ## compute ELBO
		KL = sum(qg$kl_l) + sum(qg$kl_f)
    ELBO = - sum( colSums(l0 * qg$qls_mean) * colSums(f0 * qg$qfs_mean) ) +
						sum(d$x * (log(l0[d$i]) + log(f0[d$j]) + log(B)) ) - KL - const 	
		ELBOs <- c(ELBOs, ELBO)
		## verbose
    if(verbose){
      print("iter         ELBO")
      print(sprintf("%d:    %f", i, ELBO))
    }
    ## check convergence
    # diff = ifelse(i > 2, ELBOs[i] - ELBOs[i-1], Inf)
    # if(diff < tol){
    #   if(verbose){print(sprintf("reaches tol %f in %d iterations", tol, i))}
    #   break
    # }
  }
  return(list(l0 = l0, f0 = f0, qg = qg, ELBO = ELBOs))
}



#' @title Empirical Bayes Poisson Matrix Factorization, Background Model (rank 1)
#' @import ebpm

#' @export  rank1_bg

## TODO: uupdate KL, Lam
rank1_bg <- function(d, X_rs, X_cs, l0, f0, 
										 pm_func,pm_control, 
										 ql, gl, kl_l, 
										 qf, gf, kl_f, 
										 fix_option){
  p = length(X_cs)
  n = length(X_rs)
  ## fit for f, and compute kl_f
	if(!fix_option$qf){
		s <- sum(l0 * ql$mean) * f0
		fit = do.call(pm_func$f, 
									c(list(x = X_cs, s = s, g_init = gf, fix_g = fix_option$gf), pm_control))
  	qf = fit$posterior
		gf = fit$fitted_g
		kl_f = compute_kl_ebpm(y = X_cs, s = s, posterior = qf, ll = fit$log_likelihood)
		rm(fit)
	}
  ## fit for l, and compute kl_l
	if(!fix_option$ql){
  	s = sum(f0 * qf$mean) * l0
  	fit = do.call(pm_func$l, 
									c(list(x = X_rs, s = s, g_init = gl, fix_g = fix_option$gl), pm_control))
		ql = fit$posterior
		gl = fit$fitted_g
  	kl_l = compute_kl_ebpm(y = X_rs, s = s, posterior = ql, ll = fit$log_likelihood)
		rm(fit)
	}
	## list to return
  qg = list(ql = ql, gl = gl, kl_l = kl_l, qf = qf, gf = gf, kl_f = kl_f)
  return(qg)
}






