#' @export initialize_wbg2_subject_from_LFm
initialize_wbg2_subject_from_LFm <- function(L, F, m){
  E = ncol(m)
  K = ncol(L)
  L[L <  1e-8] <- 1e-8
  F[F <  1e-8] <- 1e-8
  m[m <  1e-8] <- 1e-8
  f0 = apply(F, 1, mean)
  F = F/f0
  m = f0 * m
  qg = ebpmf.alpha::initialize_qg_from_LF(L0 = L, F0 = F)
  w = replicate(K, 1)
  m = list(mean = m, mean_log = log(m), 
    g = replicate(E, list(NULL)), kl = replicate(E, list(NULL)))
  ## replace g with mixture of gamma
  qg$gfs = replicate(K, list(bg_prior()))
  return(list(qg = qg,m = m, w = w))
}


#' @export init_ebpmf_wbg2_subject
init_ebpmf_wbg2_subject <- function(X, K, u, init, d, seed = 123){
	n = nrow(X)
  p = ncol(X)
  if(is.null(init)){
 		init = ebpmf.alpha::pmf_subject(X = X, u = u, K = K, maxiter = 50, seed = seed)
		init = ebpmf.alpha::initialize_wbg2_subject_from_LFm(L = init$L, F = init$F, m = init$m)
	}
	m = init$m
	qg = init$qg
	w = init$w

	## compute `a`
  a = replicate(length(d$x), 0)
  for(k in 1:K){
    b_k_tmp <- log(w[k]) + qg$qls_mean_log[d$i, k] + qg$qfs_mean_log[d$j, k]
    a <- pmax(a, b_k_tmp)
  }
  ## compute b
  b = log(w[1]) + qg$qls_mean_log[d$i, 1] + qg$qfs_mean_log[d$j, 1] - a
  for(k in 2:K){
    b_k = log(w[k]) + qg$qls_mean_log[d$i, k] + qg$qfs_mean_log[d$j, k] - a
    b <- log( exp(b) + exp(b_k)  )
  }
	return(list(m =m, w = w, qg = qg, b = b, a = a))
}

#' @export compute_elbo_wbg2_subject
compute_elbo_wbg2_subject <- function(w, m, qg, b, a, d, X_cs_subject, S, const){
  KL = sum(qg$kl_l) + sum(qg$kl_f) + sum(unlist(m$kl))
  L_cs_subject = compute_cs_by_subject(qg$qls_mean, S)
  elbo = - sum(colSums((m$mean %*% L_cs_subject) * qg$qfs_mean) * w)
  elbo = elbo + sum(d$x * (b + a)) + sum(t(m$mean_log) * X_cs_subject) - KL - const
  return(elbo)
}

#' @title Empirical Bayes Poisson Matrix Factorization, Background Model (rank 1)
#' @import ebpm
#' @export  rank1_wbg2
## TODO: uupdate KL, Lam
rank1_wbg2_subject <- function(d, X_rs, X_cs, S, u, m, w_log_k,
                     pm_func,pm_control,
                     ql, gl, kl_l,
                     qf, gf, kl_f,
                     fix_option){
  p = length(X_cs)
  n = length(X_rs)
  ## fit for f, and compute kl_f
  w_k = exp(w_log_k)
  #browser()
  if(!fix_option$qf){
    l_sum_subject = compute_sum_by_subject(ql$mean, S)
    s <- as.vector(w_k * (m$mean %*% l_sum_subject)) 
    fit = do.call(pm_func$f,
                  c(list(x = X_cs, s = s, g_init = gf, fix_g = fix_option$gf), pm_control))
    qf = fit$posterior
    gf = fit$fitted_g
    kl_f = compute_kl_ebpm(y = X_cs, s = s, posterior = qf, ll = fit$log_likelihood)
    rm(fit)
  }
  ## fit for l, and compute kl_l
  if(!fix_option$ql){
    s = as.vector(w_k * (t(m$mean) %*% qf$mean)[u])
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




