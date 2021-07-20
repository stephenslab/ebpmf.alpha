#' @export initialize_wbg3_from_LF
initialize_wbg3_from_LF <- function(L, F ){
  init = ebpmf.alpha::initialize_qgl0f0w_from_LF(L = L, F = F)
  init$l0 <- list(mean = init$l0, mean_log = log(init$l0), g = NULL, kl = 0)
  init$f0 <- list(mean = init$f0, mean_log = log(init$f0), g = NULL, kl = 0)
  return(init)
}

#' @export init_ebpmf_wbg3
init_ebpmf_wbg3 <- function(X, K, init, d, seed = 123){
	n = nrow(X)
  p = ncol(X)
  if(is.null(init)){
 		nnmf_fit = NNLM::nnmf(A = as.matrix(X), k = K,
                        loss = "mkl", method = "scd",
                        max.iter = 50, verbose = FALSE,
                        show.warning = FALSE)
    L = nnmf_fit$W
    F = t(nnmf_fit$H)
		init = ebpmf.alpha::initialize_wbg3_from_LF(L = L, F = F)
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

#' @export compute_elbo_wbg3
compute_elbo_wbg3 <- function(w, l0, f0, qg, b, a, d, const){
  KL = sum(qg$kl_l) + sum(qg$kl_f) + f0$kl + l0$kl
  elbo = - sum(w * colSums(l0$mean * qg$qls_mean) * colSums(f0$mean * qg$qfs_mean) ) +
            sum(d$x * (l0$mean_log[d$i] + f0$mean_log[d$j] + b + a) ) - KL - const
  return(elbo)
}







