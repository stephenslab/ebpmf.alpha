#' @export init_np_ebpmf_wbg
init_np_ebpmf_wbg <- function(X, K, alpha, c_alpha_log, init, d, seed = 123){
	set.seed(seed)
	tmp = init_ebpmf_wbg(X = X, K = K, init = init, d = d, seed = seed)
	l0 = tmp$l0 * sum(tmp$w)
	f0 = tmp$f0
	qg = tmp$qg

#	tau = rbeta(n = K, shape1 = 1, shape2 = alpha) ## TOTHINK: may need to re-adjust q_k
	w = replicate(K, 0.9/K)
	tau = replicate(K, 0)
	tau[1] = w[1]
	for(k in 2:K){
		tau[k] = w[k]/(exp(sum( log(1 - tau[1:(k-1)]) )))
	}
	eps_bar = 1
	eps_hat = 1
	## compute `a`
  a = replicate(length(d$x), 0)
  for(k in 1:K){
    b_k_tmp <- log(w[k]) + qg$qls_mean_log[d$i, k] + qg$qfs_mean_log[d$j, k]
    a <- pmax(a, b_k_tmp)
  }
  ## compute b
	b_res = c_alpha_log + sum( log(1-tau) ) + log(eps_hat) - a
  b = b_res
	for(k in 1:K){
    b_k = log(w[k]) + qg$qls_mean_log[d$i, k] + qg$qfs_mean_log[d$j, k] - a
    b <- log( exp(b) + exp(b_k)  )
  }
	return(list(l0 = l0, f0 = f0, qg = qg, b = b, a = a,
							tau = tau, eps_bar = eps_bar, eps_hat = eps_hat, b_res = b_res))
}

#' @export compute_elbo_np_wbg
compute_elbo_np_wbg <- function(alpha, tau, w_bar, l0, f0, qg, b, a, d, Lam_res,const){
	K = ncol(qg$qls_mean)
  KL = sum(qg$kl_l) + sum(qg$kl_f)
  elbo = - sum(w_bar * colSums(l0 * qg$qls_mean) * colSums(f0 * qg$qfs_mean) ) - sum(l0)*sum(f0)*Lam_res +
					sum(d$x * (log(l0[d$i]) + log(f0[d$j]) + b + a) ) +
				 	(alpha - 1)*(sum(log(1 - tau))) - K * lbeta(1, alpha)  - KL - const
  return(elbo)
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
	A = ifelse(k == 1, 1, exp(sum( log(1 - tau[1:(k-1)]))))
	if (k == K){ tmp = eps_sum }
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
		return(max(s1, 1e-20))
	}
	if(s2*(1 - s2) >= 0){
    return(max(s2, 1e-20))
  }
}

