#' @title Nonparametric Empirical Bayes Poisson Matrix Factorization (Background Model with weights)
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
#' @export  np_ebpmf_wbg

np_ebpmf_wbg <- function(X, K, alpha = 1, 
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
	c_alpha_log = digamma(1) - log(exp(digamma(1 + alpha)) - exp(digamma(alpha)))
	const = sum(apply.nonzeros(X = X, f = function(x) lgamma(x + 1)))
  ## initialization
	init_tmp <- init_np_ebpmf_wbg(X = X, K = K, alpha = alpha, c_alpha_log = c_alpha_log,
																init = init, d = d, seed = seed)
  qg <- init_tmp$qg
  b <- init_tmp$b
  a <- init_tmp$a
	l0 <- init_tmp$l0
	f0 <- init_tmp$f0
	tau <- init_tmp$tau
	eps_bar <- init_tmp$eps_bar
	eps_hat <- init_tmp$eps_hat	
	b_res <- init_tmp$b_res
	w_ <- tau2w(tau)
	w_bar_log <- w_$w_bar_log 
	w_hat_log <- w_$w_hat_log 
  rm(init_tmp)
  ## update iteratively
  ELBOs <- c()
	KLs <- c()
	#browser()
  for(i in 1:maxiter){
		#browser()
		expb_sum =  replicate(length(d$x),0) 
		b_k_max = replicate(length(d$x),0) ## max b_k
		for(k in 1:K){
			print(k)
			#if(i == 9 && k == 7){browser()}
			## store b_k
 			b_k = w_hat_log[k] + qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k] - a
    	## update sum_{t = 1: k} exp(b_ijt)
			## compute q(Z)
      Ez <- compute_EZ(d = d,b = b,b_k = b_k)
      ## update (qL, gL, qF, gF) 
      rank1_qg <- rank1_wbg(d = d, X_rs = Ez$rs, X_cs = Ez$cs,
														l0 = l0, f0 = f0, w_log_k = w_bar_log[k], 
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

			## update tau
			zeta = exp(b_k - b)
			zeta_sum = expb_sum/exp(b) + zeta
			zeta_sum[zeta_sum >= 1] <- 1 - 1e-6 ## TODO: how to prevent it
			#tau[k] = optim_tau_k(alpha = alpha, tau = tau, k = k, 
			#										 zeta_sum = zeta_sum, zeta = zeta, 
			#										 d= d, l0 = l0, f0 = f0, qg = qg, eps_bar = eps_bar)
			w_bar_log[k] = log(tau[k]) + ifelse(k==1,0,sum( log(1 - tau[1:(k-1)]) ))
			w_hat_log[k] = w_bar_log[k]

			## update b
      b_k0 = b_k
      b_k = w_hat_log[k] + qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k] - a
			b_k_max = pmax(b_k, b_k_max)
			b = log( exp(b) - exp(b_k0) + exp(b_k)  )
			expb_sum = expb_sum + exp(b_k)
    }
		## update l0, f0
		w_bar = exp(w_bar_log)
		w_hat = exp(w_hat_log)
		## update b
		b_res0 = b_res
		b_res = c_alpha_log + sum( log(1-tau) ) + log(eps_hat) - a
		b = log( exp(b) - exp(b_res0) + exp(b_res) )
		b_k_max = pmax(b_res, b_k_max)
    w_bar_res = exp(sum( log(1-tau) ))
    Lam_res = exp(sum( log(1-tau) ) + log(eps_bar))
		#browser()
		if(!fix_option$l0){
      denom <- colSums(w_bar * t(qg$qls_mean) * colSums(f0 * qg$qfs_mean)) + sum(f0)*w_bar_res*eps_bar
      l0 <- X_rs/denom
    }
    if(!fix_option$f0){
      denom <- colSums(w_bar * t(qg$qfs_mean) * colSums(l0 * qg$qls_mean)) + sum(l0)*w_bar_res*eps_bar
      f0 <- X_cs/denom
    }
		## compute ELBO
		KL = sum(qg$kl_l) + sum(qg$kl_f)
		ELBO = compute_elbo_np_wbg(alpha = alpha, tau = tau,  w_bar = w_bar, 
															 l0 = l0, f0 = f0, qg = qg, 
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
  return(list(tau = tau, w_bar = w_bar, w_hat = w_hat, l0 = l0, f0 = f0, qg = qg, ELBO = ELBOs, KL = KLs))
}


tau2w <- function(tau){
	K = length(tau)
	w_bar_log = log(tau) + cumsum(c(0, log(1-tau[1:(K-1)])))	
	w_hat_log = w_bar_log
	return(list(w_bar_log = w_bar_log, w_hat_log = w_hat_log))
}

w2tau <- function(w){
	K = length(w)
	tau = replicate(K, NA)
	tau[1] = w[1]
	for(k in 2:K){
		tau[k] = w[k]/( exp(sum ( log(1-tau[1:(k-1)]))))
	}
	return(tau)
}

optim_tau_k <- function(alpha, tau, k, zeta_sum, zeta, d, l0, f0, qg, eps_bar){
	K = ncol(qg$qls_mean)
	mu_sum = colSums(l0 * qg$qls_mean) * colSums(f0 * qg$qfs_mean)
	eps_sum = sum(l0) * sum(f0) * eps_bar
	A = ifelse(k == 1, 1, exp(cumsum( log(1 - tau[1:(k-1)]))))
	if (k == K){ tmp = 0}
	if (k == (K - 1)){ tmp = tau[K] * mu_sum[K] + eps_sum * (1 - tau[K])}
	if (k < (K-1)){
		tmp = sum( (tau[(k+1):K] * mu_sum[(k+1):K]) * 
							exp(cumsum( c(0, log(1-tau[(k+1):(K-1)])) ))) + 
					eps_sum * ( exp( sum( log(1- tau[(k+1):K])) ) )
	}
	A = A * (tmp - mu_sum[k])
	B = sum( d$x * (1 - zeta_sum) ) + alpha -1
	C = sum( d$x * zeta )
	tau_k = solve_quadratic(A = A, B = B, C = C)
	return(tau_k)
}


solve_quadratic <- function(A, B, C){
	a = A
	b = B - A +C
	c = -C
	s1 = (-b + sqrt(b^2 - 4*a*c))/(2*a)
	s2 = (-b - sqrt(b^2 - 4*a*c))/(2*a)
	if(s1*(1 - s1) >= 0){
		return(max(s1, 1e-6))
	}
	if(s2*(1 - s2) >= 0){
    return(max(s2, 1e-6))
  }
}



