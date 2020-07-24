#' @export init_np_ebpmf_bg
init_np_ebpmf_bg <- function(X, K, alpha, beta, c_alpha_log, init, d, seed = 123){
	set.seed(seed)
	n = nrow(X)
  p = ncol(X)
  if(is.null(init)){
		nnmf_fit = NNLM::nnmf(A = as.matrix(X), k = K,
                        loss = "mkl", method = "lee",
                        max.iter = 50, verbose = FALSE,
                        show.warning = FALSE)
		L = nnmf_fit$W
		F = t(nnmf_fit$H)
		init = ebpmf.alpha::initialize_np_qgl0f0_from_LF(L = L, F = F)
  }
  l0_bar = init$l0
	f0 = init$f0
	qg = init$qg
	rm(init)
	## get tau from qg$qls_mean
	tau = L2tau(qg$qls_mean)
	eps_bar = 1
	eps_hat = 1
	## compute `a`
	a = (c_alpha_log + rowSums( log(1-tau) ) + log(eps_hat))[d$i] ## start with residual
	for(k in 1:K){
    b_k_tmp <- qg$qls_mean_log[d$i, k] + qg$qfs_mean_log[d$j, k]
		a <- pmax(a, b_k_tmp)
  }
	## compute b
	b_res = c_alpha_log + rowSums( log(1-tau) ) + log(eps_hat) 
  b = b_res[d$i] - a
	for(k in 1:K){
    b_k = qg$qls_mean_log[d$i, k] + qg$qfs_mean_log[d$j, k] - a
    b <- log( exp(b) + exp(b_k)  )
  }
  return(list(qg = qg, l0_bar = l0_bar, f0 = f0, tau = tau,
  	eps_bar = eps_bar, eps_hat = eps_hat, a = a, b = b, b_res = b_res))
}

#' @export initialize_np_qgl0f0_from_LF
initialize_np_qgl0f0_from_LF <- function(L, F){
	L[L <  1e-8] <- 1e-8
  F[F <  1e-8] <- 1e-8
	l0 = apply(L, 1, sum) ## need to make sum_k l_ik = 1
  #l0[l0 == 0] <- 1e-8
  f0 = apply(F, 1, mean)
  #f0[f0 == 0] <- 1e-8
  L = L/l0
  F = F/f0
  qg = ebpmf.alpha::initialize_qg_from_LF(L0 = L, F0 = F)
  ## replace g with mixture of gamma
  K = ncol(L)
  qg$gfs = replicate(K, list(bg_prior()))
  return(list(qg = qg, l0 = l0, f0 = f0))
}

optim_tau_k_vec <- function(alpha, tau, k, zeta_sum, zeta, d, l0, f0, qg, eps_bar){
	K = ncol(qg$qls_mean)
	I = nrow(qg$qls_mean)
	mu_sum = l0 %o% colSums(f0 * qg$qfs_mean) ## I by K
	eps_sum = l0 * sum(f0) * eps_bar ## I	

	B = rowSums( sparseMatrix(i = d$i, j = d$j, x = d$x * (1 - zeta_sum)) ) + alpha -1
	C = rowSums( sparseMatrix(i = d$i, j = d$j, x = d$x * zeta) )

	if(k == 1){A = 1}
	if(k == 2){A = exp( log(1 - tau[,1:(k-1)])) }
	if(k > 2){A = exp( rowSums( log(1 - tau[,1:(k-1)])) )} ## I

	if (k == K){ tmp = eps_sum }
	if (k == (K - 1)){ tmp = tau[,K] * mu_sum[,K] + eps_sum * (1 - tau[,K])}
	if (k < (K-1)){
		tmp = rowSums( 
					(tau[,(k+1):K] * mu_sum[,(k+1):K]) * 
					#exp(cumsum( c(0, log(1-tau[(k+1):(K-1)])) ))) + 
					exp(cumsum_row(cbind(replicate(I, 0), log(1-tau[,(k+1):(K-1)])))) ) + 
					eps_sum * ( exp( rowSums( log(1- tau[,(k+1):K])) ) )
	}
	A = A * (tmp - mu_sum[, k])
	tau_k = solve_quadratic_vec(A = A, B = B, C = C)
	return(tau_k)
}


solve_quadratic_vec <- function(A, B, C){
	I = length(A)
	tau_k = replicate(I, NA)
	a = A
	b = B - A +C
	c = -C
	s1 = (-b + sqrt(b^2 - 4*a*c))/(2*a)
	s2 = (-b - sqrt(b^2 - 4*a*c))/(2*a)
	mask <- (s1*(1 - s1) >= 0)
	tau_k[mask] = s1[mask]
	mask <- (s2*(1 - s2) >= 0)
	tau_k[mask] = s2[mask]
	tau_k[tau_k == 0] = 1e-20
	return(tau_k)
}

compute_elbo_np_bg <- function(alpha, beta, alpha_l, beta_l, 
															 tau, l0_bar, f0, qg, b, a,
															 d, Lam_res, const){
	K = ncol(qg$qls_mean)
	n = nrow(qg$qls_mean)
  KL = sum(qg$kl_f) + 
  			sum(KL.gamma(c = alpha_l, d = beta_l, a = replicate(n, alpha), b = replicate(n, beta)))
  elbo = - sum(colSums(l0_bar * qg$qls_mean) * colSums(f0 * qg$qfs_mean) ) - sum(l0_bar * Lam_res)*sum(f0) +
  				sum(d$x * (log(l0_bar[d$i]) + log(f0[d$j]) + b + a) ) +
  				(alpha - 1)*(sum(log(1 - tau))) - n * K * lbeta(1, alpha)  - KL - const
  return(elbo)
}




L2tau <- function(L){
	K = ncol(L)
	n = nrow(L)
	tau = matrix(, nrow = n, ncol = K)
	tau[,1] = L[,1]
	if(K > 1){
		tau[,2] = L[,2]/( exp( log(1-tau[,1])))
	}
	for(k in 3:K){
		tau[,k] = L[,k]/( exp(rowSums ( log(1-tau[,1:(k-1)]))))
	}
	tau[tau >= (1-1e-5)] = 1 - 1e-5
	return(tau)
}

tau2L_k <- function(tau, k, log = TRUE){
	if(k == 1){
		l_log = log(tau[,k])
	}
	if(k == 2){
		l_log = log(tau[,k]) + log(1 - tau[,1:(k-1)]) 
	}
	if(k > 2){
		l_log = log(tau[,k]) + rowSums( log(1 - tau[,1:(k-1)]) )
	}
	if(log){
		return(l_log)
	}else{
		return(exp(l_log))
	}
}







