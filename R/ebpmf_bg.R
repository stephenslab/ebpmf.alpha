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
										 fix_g = list(l = FALSE, f = FALSE), 
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
    KL <- 0
    for(k in 1:K){
			## store B_k
 			B_k = exp(qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k]) 
    	## compute q(Z)
      Ez <- compute_EZ(d = d,B = B,B_k = B_k)
      ## update (qL, gL, qF, gF) 
      init_r1 = list(sf = sum(l0 * qg$qls_mean[,k]),
                     gl = qg$gls[[k]], gf = qg$gfs[[k]])
      rank1_tmp <- rank1_bg(d = d, X_rs = Ez$rs, X_cs = Ez$cs,
														l0 = l0, f0 = f0, 
														pm_func = pm_func, pm_control = pm_control,
														init = init_r1, fix_g = fix_g)
      rm(Ez)
      KL = KL + rank1_tmp$kl_l + rank1_tmp$kl_f
      qg = update_qg(rank1_tmp$qg, qg, k)
      rm(rank1_tmp)
			## update B, Lam
			B = B - B_k + exp(qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k])
		 #	B[B < 0] <- 1e-20 ## numerical issue gets - epsilon	
    }
		## update l0, f0
		denom <- colSums( t(qg$qls_mean) * colSums(f0 * qg$qfs_mean)) 
		l0 <- X_rs/denom
		denom <- colSums( t(qg$qfs_mean) * colSums(l0 * qg$qls_mean)) 
    f0 <- X_cs/denom
    ## compute ELBO
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

#' @param d: summary(X), so it has `i`, `j` and `x` as attributes
#' @param X_rs, rowSums of X,  count matrix (dim(X) = c(n, p)).
#' @param X_cs, colSums of X,  count matrix (dim(X) = c(n, p)).
#' @param pm_func function for solving the \code{ebpm} subproblem
#' @param pm_control control parameters for pm_func function
#' @param init list(sf, gl, gf), where sf is the scale for ebpm problem for F,
#' gl, gf are initialization for gl, gf; all can be set to NULL
#' @param fix_g list(l, f) where l, f are either TRUE or FALSE
#'
#' @return A list containing elements:
#'     \describe{
#'       \item{\code{ql}}{approximate posterior for l}
#'       \item{\code{gl}}{fitted g for l}
#'       \item{\code{kl_l}}{kl divergence between q and g for l}
#'       \item{\code{qf}}{approximate posterior for f}
#'       \item{\code{gf}}{fitted g for f}
#'       \item{\code{kl_f}}{kl divergence between q and g for f}
#'      }
#' @examples
#' To add
#' @export  rank1

## TODO: uupdate KL, Lam
rank1_bg <- function(d, X_rs, X_cs, l0, f0, pm_func,pm_control, init, fix_g){
  p = length(X_cs)
  n = length(X_rs)
  ## initialization (in fact, any non-negative number well do)
  sf = init$sf
  if(is.null(sf)){
		sf = 1
  }
  ## fit for f, and compute kl_f
	fit_f = do.call(pm_func$f, c(list(x = X_cs, s = sf*f0, g_init = init$gf, fix_g = fix_g$f), pm_control))
  kl_f = compute_kl_ebpm(y = X_cs, s = sf*f0, posterior = fit_f$posterior, ll = fit_f$log_likelihood)
  ## fit for l, and compute kl_l
  sl = sum(f0 * fit_f$posterior$mean)
  fit_l = do.call(pm_func$l, c(list(x = X_rs, s = sl*l0, g_init = init$gl, fix_g = fix_g$l), pm_control))
  kl_l = compute_kl_ebpm(y = X_rs, s = sl*l0, posterior = fit_l$posterior, ll = fit_l$log_likelihood)
  ## list to return
  qg = list(ql = fit_l$posterior, gl = fit_l$fitted_g, qf = fit_f$posterior, gf = fit_f$fitted_g)
  out = list(qg = qg, kl_l = kl_l, kl_f = kl_f)
  return(out)
}






