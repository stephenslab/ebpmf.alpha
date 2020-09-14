#' @export init_ebpmf_bg
init_ebpmf_bg <- function(X,K, init, d, seed = 123){
	n = nrow(X)
  p = ncol(X)
  if(is.null(init)){
		nnmf_fit = NNLM::nnmf(A = as.matrix(X), k = K,
                        loss = "mkl", method = "lee",
                        max.iter = 50, verbose = FALSE,
                        show.warning = FALSE)
		L = nnmf_fit$W
		F = t(nnmf_fit$H)
		init = ebpmf.alpha::initialize_qgl0f0_from_LF(L = L, F = F)
  }
	l0 = init$l0
	f0 = init$f0
	qg = init$qg
	## compute `a`
	a = replicate(length(d$x), 0)
	for(k in 1:K){
    b_k_tmp <- qg$qls_mean_log[d$i, k] + qg$qfs_mean_log[d$j, k]
		a <- pmax(a, b_k_tmp)
  }
	## compute b
  b = qg$qls_mean_log[d$i, 1] + qg$qfs_mean_log[d$j, 1] - a
  for(k in 2:K){
		b_k = qg$qls_mean_log[d$i, k] + qg$qfs_mean_log[d$j, k] - a
    b <- log( exp(b) + exp(b_k)  )
  }
  return(list(l0 = l0, f0 = f0, qg = qg, b = b, a = a))
}

#' @title Empirical Bayes Poisson Matrix Factorization, Background Model (rank 1)
#' @import ebpm
#' @export  rank1_bg
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
    if(is.infinite(fit$log_likelihood)){browser()}
    qf = fit$posterior
    gf = fit$fitted_g
    kl_f = compute_kl_ebpm(y = X_cs, s = s, posterior = qf, ll = fit$log_likelihood)
    #if(is.infinite(kl_f)){browser()}
    rm(fit)
  }
  ## fit for l, and compute kl_l
  if(!fix_option$ql){
    s = sum(f0 * qf$mean) * l0
    fit = do.call(pm_func$l, 
                  c(list(x = X_rs, s = s, g_init = gl, fix_g = fix_option$gl), pm_control))
    if(is.infinite(fit$log_likelihood)){browser()}
    ql = fit$posterior
    gl = fit$fitted_g
    kl_l = compute_kl_ebpm(y = X_rs, s = s, posterior = ql, ll = fit$log_likelihood)
    rm(fit)
  }
  ## list to return
  qg = list(ql = ql, gl = gl, kl_l = kl_l, qf = qf, gf = gf, kl_f = kl_f)
  return(qg)
}




#' @export compute_elbo_bg
compute_elbo_bg <- function(l0, f0, qg, b, a, d, const){
  KL = sum(qg$kl_l) + sum(qg$kl_f)
  elbo = - sum(colSums(l0 * qg$qls_mean) * colSums(f0 * qg$qfs_mean) ) +
            sum(d$x * (log(l0[d$i]) + log(f0[d$j]) + b + a) ) - KL - const
  return(elbo)
}


check_progress_bg <- function(elbo_prev, which_part,k,
                              l0, f0, qg, b, a, d, const){
  ## check progress
  elbo = compute_elbo_bg(l0 = l0, f0 = f0, qg = qg,
                          b = b, a = a, d = d, const = const)
  if(elbo < elbo_prev){
          print(sprintf("k = %d, %s updated, elbo_diff = %f",
                        k, which_part, elbo - elbo_prev))
  }
  return(elbo)
}