#' @title Empirical Bayes Poisson Matrix Factorization (Background Model with weights)
#' @import ebpm
#' @import Matrix

#' @param X count matrix (dim(X) = c(n, p)).
#' @param k number of topics
#' @param pm_func functions for solving the \code{ebpm} subproblem for \code{L} and \code{F}; 
#' It is a list \code{list(l, f)};
#' For our purpose we use `mle_pm` or `ebpm_point_gammma`for \code{L}, and \code{ebpm_gamma_mixture} for \code{F}
#' @param pm_control control parameters for pm_func function
#' @param init Either \code{NULL} or \code{list(qg, l0, f0, w)}
#' @param fix_g list(l, f) where l, f are either TRUE or FALSE
#' @param maxiter maximum number of iterations
#' @param tol stopping tolerance for ELBO
#' @param seed used when init is NULL
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
#' @export  ebpmf_wbg

ebpmf_wbg <- function(X, K, 
										 pm_func = list(f = ebpm::ebpm_gamma_mixture, 
																		l = ebpm::ebpm_gamma_mixture),
										 init = NULL, pm_control = NULL,
										 fix_option = list(l0 = FALSE, f0 = FALSE,
																				gl = FALSE, ql = FALSE,
																				gf = FALSE, qf = FALSE), 
										 maxiter = 100, tol = 1e-8, 
										 verbose = FALSE, seed = 123){
	## TODO: input check, require X_rs, X_cs to be nonzero

  ## transform to sparse matrix 
  X <- as(X, "sparseMatrix") 
	X_rs = Matrix::rowSums(X)
	X_cs = Matrix::colSums(X)
  d = summary(X)
	const = sum(apply.nonzeros(X = X, f = function(x) lgamma(x + 1)))
  ## initialization
	init_tmp <- init_ebpmf_wbg(X = X, K = K, init = init, d = d, seed = seed)
  qg <- init_tmp$qg
  b <- init_tmp$b
  a <- init_tmp$a
	l0 <- init_tmp$l0
	f0 <- init_tmp$f0
	w <- init_tmp$w
	w_log = log(w)
  rm(init_tmp)
  ## update iteratively
  ELBOs <- c()
	KLs <- c()
  for(i in 1:maxiter){
		b_k_max = replicate(length(d$x),0) ## max b_k
		for(k in 1:K){
			## store B_k
 			b_k = w_log[k] + qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k] - a
    	## compute q(Z)
      Ez <- compute_EZ(d = d,b = b,b_k = b_k)
      ## update (qL, gL, qF, gF) 
      rank1_qg <- rank1_wbg(d = d, X_rs = Ez$rs, X_cs = Ez$cs,
														l0 = l0, f0 = f0, w_log_k = w_log[k], 
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
			## update b
      b_k0 = b_k
      b_k = w_log[k] + qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k] - a
			b = log( exp(b) - exp(b_k0) + exp(b_k)  )
    }
		## update w[1:K], and b accordingly
		b_new = b
		for(k in 1:K){
			b_k = w_log[k] + qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k] - a
			w_log[k] = log( sum(d$x * exp(b_k - b)) ) - log( sum(l0 * qg$qls_mean[,k]) * sum(f0 * qg$qfs_mean[,k]) )
			b_k_new = w_log[k] + qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k] - a
			b_new =  log( exp(b_new) - exp(b_k) + exp(b_k_new)  )
			b_k_max = pmax(b_k, b_k_max)
		}
		b = b_new	
	## update l0, f0
    if(!fix_option$l0){
      denom <- colSums(w * t(qg$qls_mean) * colSums(f0 * qg$qfs_mean))
      l0 <- X_rs/denom
    }
    if(!fix_option$f0){
      denom <- colSums(w * t(qg$qfs_mean) * colSums(l0 * qg$qls_mean))
      f0 <- X_cs/denom
    }
		## compute ELBO
		w = exp(w_log)
		KL = sum(qg$kl_l) + sum(qg$kl_f)
		ELBO = compute_elbo_wbg(w = w, l0 = l0, f0 = f0, qg = qg, 
														b = b, a = a, d = d, const = const)
		ELBOs <- c(ELBOs, ELBO)
		KLs <- c(KLs, KL)
		## update a & b
    #a0 = a
    #a = b_k_max + a0
    #b = (b + a0) - a
    b = b - b_k_max
    a = b_k_max + a
		## verbose
    if(verbose){
      print("iter         ELBO")
      print(sprintf("%d:    %f", i, ELBO))
    }
  }
  return(list(w = w, l0 = l0, f0 = f0, qg = qg, ELBO = ELBOs, KL = KLs))
}









