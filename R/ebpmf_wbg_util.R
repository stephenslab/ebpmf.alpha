#' @export init_ebpmf_wbg
init_ebpmf_wbg <- function(X, K, init, d, seed = 123){
	n = nrow(X)
  p = ncol(X)
  if(is.null(init)){
 		nnmf_fit = NNLM::nnmf(A = as.matrix(X), k = K,
                        loss = "mkl", method = "scd",
                        max.iter = 50, verbose = FALSE,
                        show.warning = FALSE)
    L = nnmf_fit$W
    F = t(nnmf_fit$H)
		init = ebpmf.alpha::initialize_qgl0f0w_from_LF(L = L, F = F)
	}
	l0 = init$l0
	f0 = init$f0
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
	return(list(l0 = l0, f0 = f0, w = w, qg = qg, b = b, a = a))
}

#' @export compute_elbo_wbg
compute_elbo_wbg <- function(w, l0, f0, qg, b, a, d, const){
  KL = sum(qg$kl_l) + sum(qg$kl_f)
  elbo = - sum(w * colSums(l0 * qg$qls_mean) * colSums(f0 * qg$qfs_mean) ) +
            sum(d$x * (log(l0[d$i]) + log(f0[d$j]) + b + a) ) - KL - const
  return(elbo)
}

#' @title Empirical Bayes Poisson Matrix Factorization, Background Model (rank 1)
#' @import ebpm
#' @export  rank1_wbg
## TODO: uupdate KL, Lam
rank1_wbg <- function(d, X_rs, X_cs, l0, f0, w_log_k, 
                     pm_func,pm_control, 
                     ql, gl, kl_l, 
                     qf, gf, kl_f, 
                     fix_option){
  p = length(X_cs)
  n = length(X_rs)
  ## fit for f, and compute kl_f
  w_k = exp(w_log_k)
  if(!fix_option$qf){
    s <- sum(l0 * ql$mean) * f0 * w_k
    fit = do.call(pm_func$f, 
                  c(list(x = X_cs, s = s, g_init = gf, fix_g = fix_option$gf), pm_control))
    qf = fit$posterior
    gf = fit$fitted_g
    kl_f = compute_kl_ebpm(y = X_cs, s = s, posterior = qf, ll = fit$log_likelihood)
    rm(fit)
  }
  ## fit for l, and compute kl_l
  if(!fix_option$ql){
    s = sum(f0 * qf$mean) * l0 * w_k
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

check_progress_wbg <- function(elbo_prev, which_part,k, 
                               w_log, l0, f0, qg, b, a, d, const){
  elbo = compute_elbo_wbg(w = exp(w_log), l0 = l0, f0 = f0, qg = qg,
                          b = b, a = a, d = d, const = const)
  if(elbo < elbo_prev){
          print(sprintf("k = %d, %s updated, elbo_diff = %f", 
                        k, which_part, elbo - elbo_prev))
  }
  return(elbo)
}




