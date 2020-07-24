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