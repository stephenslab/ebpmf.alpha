#' @title Nonparametric Empirical Bayes Poisson Matrix Factorization (Background Model)
#' @import ebpm
#' @import Matrix

#' @param X count matrix (dim(X) = c(n, p)).
#' @param k number of topics
#' @param pm_func functions for solving the \code{ebpm} subproblem for \code{L} and \code{F}; 
#' It is a list \code{list(l, f)};
#' For our purpose we use `mle_pm` or `ebpm_point_gammma`for \code{L}, and \code{ebpm_gamma_mixture} for \code{F}
#' @param pm_control control parameters for pm_func function
#' @param init Either \code{NULL} or \code{list(qg, l0_bar, f0)}
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
#' @export  np_ebpmf_bg

np_ebpmf_bg <- function(X, K, alpha = 1, beta = 1,
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
  fix_option$ql = TRUE ## hasty fix
  X <- as(X, "sparseMatrix") 
	X_rs = Matrix::rowSums(X)
	X_cs = Matrix::colSums(X)
  d = summary(X)
	c_alpha_log = digamma(1) - log(exp(digamma(1 + alpha)) - exp(digamma(alpha)))
	const = sum(apply.nonzeros(X = X, f = function(x) lgamma(x + 1)))
  ## initialization
	init_tmp <- init_np_ebpmf_bg(X = X, K = K, beta = beta, alpha = alpha, c_alpha_log = c_alpha_log,
															init = init, d = d, seed = seed)
  qg <- init_tmp$qg
  b <- init_tmp$b
  a <- init_tmp$a
	l0_bar <- init_tmp$l0_bar
	f0 <- init_tmp$f0
	tau <- init_tmp$tau
	eps_bar <- init_tmp$eps_bar
	eps_hat <- init_tmp$eps_hat	
	b_res <- init_tmp$b_res
  rm(init_tmp)
  ## update iteratively
  ELBOs <- c()
	KLs <- c()
  for(i in 1:maxiter){
		expb_sum =  replicate(length(d$x),0) 
		b_k_max = replicate(length(d$x),0) 
		for(k in 1:K){
			#print(k)
			#if(i == 9 && k == 7){browser()}
			## store b_k
 			b_k =  qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k] - a
			## compute q(Z)
      Ez <- compute_EZ(d = d,b = b, b_k = b_k)
      ## update (qF, gF) 
			rank1_qg <- rank1_bg(d = d, X_rs = Ez$rs, X_cs = Ez$cs,
													 l0 = l0_bar, f0 = f0, 
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
      ## update tau and q(L)
			zeta = exp(b_k - b)
			zeta_sum = expb_sum/exp(b) + zeta
			zeta_sum[zeta_sum >= 1] <- 1 ## just in case
			tau[,k] = optim_tau_k_vec(alpha = alpha, tau = tau, k = k, 
													 zeta_sum = zeta_sum, zeta = zeta, 
													 d = d, l0 = l0_bar, f0 = f0, qg = qg, eps_bar = eps_bar)
			qg$qls_mean_log[,k] = tau2L_k(tau, k, log = TRUE)
			qg$qls_mean[,k] = exp(qg$qls_mean_log[,k])
			## update b
      b_k0 = b_k
      b_k = qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k] - a
			b_k_max = pmax(b_k, b_k_max)
			b = log( exp(b) - exp(b_k0) + exp(b_k)  )
			expb_sum = expb_sum + exp(b_k)
    }
		## update b
		b_res0 = b_res
		b_res = c_alpha_log + rowSums( log(1-tau) ) + log(eps_hat) ## a vector with length I
		b = log( exp(b) - exp(b_res0[d$i] - a) + exp(b_res[d$i] - a) )
		b_k_max = pmax(b_res[d$i], b_k_max)
    #w_bar_res = exp(sum( log(1-tau) ))
    Lam_res = exp( rowSums( log(1-tau) ) + log(eps_hat) ) ## a vector with length I
		if(!fix_option$l0){
      alpha_l = alpha + X_rs
      #beta_l = c + qg$qls_mean * sum(f0 * qg$qfs_mean) + sum(f0)*Lam_res
			beta_l = beta + colSums(t(l0_bar * qg$qls_mean) * colSums(f0 * qg$qfs_mean)) + sum(f0)*Lam_res
      l0_bar <- (alpha_l)/(beta_l)
    }
    if(!fix_option$f0){
      denom <- colSums(t(qg$qfs_mean) * colSums(l0_bar * qg$qls_mean)) + sum(l0_bar * Lam_res)
      f0 <- X_cs/denom
    }
		## compute ELBO
		KL = sum(qg$kl_f)
		ELBO = compute_elbo_np_bg(alpha = alpha, beta = beta, 
															alpha_l = alpha_l, beta_l = beta_l,
															tau = tau, l0_bar = l0_bar, f0 = f0, qg = qg, 
															b = b, a = a, d = d, Lam_res = Lam_res, const = const)
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
  return(list(tau = tau, l0_bar = l0_bar, f0 = f0, qg = qg, ELBO = ELBOs, KL = KLs))
}




