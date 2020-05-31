## ================================================================================================================
## Functions updating in ebpmf
## ================================================================================================================

#' @export compute_EZ
compute_EZ <- function(d, b, b_k){
  Ez = sparseMatrix(i = d$i, j = d$j, x = d$x * exp(b_k - b))
  return(list(rs = Matrix::rowSums(Ez), cs = Matrix::colSums(Ez)))
}

#' @export update_qg
update_qg <- function(tmp, qg, k){
  qg$qls_mean[,k] = tmp$ql$mean
  qg$qls_mean_log[,k] = tmp$ql$mean_log
  qg$qfs_mean[,k] = tmp$qf$mean
  qg$qfs_mean_log[,k] = tmp$qf$mean_log
  qg$gls[k] = list(tmp$gl)
  qg$gfs[k] = list(tmp$gf)
	qg$kl_l[k] = tmp$kl_l
	qg$kl_f[k] = tmp$kl_f
  return(qg)
}



## ================================================================================================================
## Functions for initialization
## ================================================================================================================

#' @export bg_prior
bg_prior <- function(){
  aL = c(seq(0.01, 0.10, 0.01), seq(0.2, 0.9, 0.1), seq(1,15,2), 20, 50, 75, 100, 200, 1e3, 1e-8, 1e-16)
  D = length(aL)
  g = ebpm::gammamix(pi = replicate(D, 1/D), shape = aL, scale = 1/aL)
  return(g)
}


#' @export initialize_qg_from_LF
initialize_qg_from_LF <- function(L0,F0){
  K = ncol(L0)
  qls_mean = L0
  qls_mean[qls_mean == 0] = 1e-10
  qfs_mean = F0
  qfs_mean[qfs_mean == 0] = 1e-10
  qls_mean_log = log(L0)
  qfs_mean_log = log(F0)
	kl_l = replicate(K, 0)
	kl_f = replicate(K, 0)
  qg = list(qls_mean = qls_mean, qls_mean_log = qls_mean_log, kl_l = kl_l,
            qfs_mean = qfs_mean, qfs_mean_log = qfs_mean_log,	kl_f = kl_f, 
            gls = replicate(K, list(NULL)),gfs = replicate(K, list(NULL)))
  return(qg)
}

#' @export initialize_qg
initialize_qg <- function(X, K, init_method = "scd", init_iter = 20, seed = 123){
  set.seed(seed)
  nnmf_fit = NNLM::nnmf(A = as.matrix(X), k = K, 
												loss = "mkl", method = init_method,
												max.iter = init_iter, verbose = FALSE,
												show.warning = FALSE)
	qg = initialize_qg_from_LF(L0 = nnmf_fit$W, F0 = t(nnmf_fit$H))
  return(qg)
}


## initialization for `ebpmf`

#' @export init_ebpmf
init_ebpmf <- function(X,K, init, d){
  qg = init$qg
  if(is.null(qg)){
    qg = initialize_qg(X, K, init_method =  init$init_method, init_iter = init$init_iter)
  }
  ## TODO: speedup
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
	return(list(qg = qg, b = b, a = a))
}

##### functions for initializing `ebpmf_wbg`

#' @export initialize_qg_l0f0
initialize_qg_l0f0 <- function(X, K, seed = 123){
	set.seed(seed)
	n = nrow(X)
	p = ncol(X)
	L0 = matrix(exp(runif(n*K, min = -1.5, max = 1)), ncol = K)
  F0 = matrix(exp(runif(p*K, min = -1.5, max = 1)), ncol = K)
  qg = initialize_qg_from_LF(L0, F0)
  l0 = rowSums(X)
  denom <- colSums( t(qg$qfs_mean) * colSums(l0 * qg$qls_mean))
  f0 <- colSums(X)/denom
	## initialize g
	aL = c(seq(0.01, 0.10, 0.01), seq(0.2, 0.9, 0.1), seq(1,15,2), 20, 50, 75, 100, 200, 1e3)
  D = length(aL)
  g = ebpm::gammamix(pi = replicate(D, 1/D), shape = aL, scale = 1/aL)
  qg$gls = replicate(K, list(g))
  qg$gfs = replicate(K, list(g))
	return(list(qg = qg, l0 = l0, f0 = f0))
}

#' @export initialize_qgl0f0_from_LF
initialize_qgl0f0_from_LF <- function(L, F){
  L[L <  1e-8] <- 1e-8
  F[F <  1e-8] <- 1e-8
	l0 = apply(L, 1, mean) ## use median tends to get huge numbers ...
  #l0[l0 == 0] <- 1e-8
  f0 = apply(F, 1, mean)
  #f0[f0 == 0] <- 1e-8
  L = L/l0
  F = F/f0
  qg = ebpmf.alpha::initialize_qg_from_LF(L0 = L, F0 = F)
  ## replace g with mixture of gamma
  K = ncol(L)
  qg$gls = replicate(K, list(bg_prior()))
  qg$gfs = replicate(K, list(bg_prior()))
  return(list(qg = qg, l0 = l0, f0 = f0))
}

## given L_ik, we compute MLE for X_ij ~ Pois(f_j0 sum_k l_ik)
## which gives us `f_j0 = X_.j/l_..`
## then we transform it into background model

#' @export initialize_qgl0f0_from_L
initialize_qgl0f0_from_L <- function(X, L){
	p = ncol(X)
	K = ncol(L)
  L[L <  1e-8] <- 1e-8
  l0 = apply(L, 1, mean) ## use median tends to get huge numbers ...
  f0 = colSums(X)/sum(L)
  #f0[f0 == 0] <- 1e-8
  L = L/l0
  F = matrix(replicate(p*K, 1), ncol = K)
  qg = ebpmf.alpha::initialize_qg_from_LF(L0 = L, F0 = F)
  ## replace g with mixture of gamma
  qg$gls = replicate(K, list(bg_prior()))
  qg$gfs = replicate(K, list(bg_prior()))
  return(list(qg = qg, l0 = l0, f0 = f0))
}

## output: qg, B
## init either NULL, or list(qg, l0, f0)

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

##### functions for initializing `ebpmf_wbg`

#' @export initialize_qgl0f0w_from_LF
initialize_qgl0f0w_from_LF <- function(L, F){
  K = ncol(L)
	w = replicate(K, 1)
  init_tmp = ebpmf.alpha::initialize_qgl0f0_from_LF(L = L, F = F)
	return(list(qg = init_tmp$qg, l0 = init_tmp$l0, f0 = init_tmp$f0, w = w))
}

#' @export initialize_qgl0f0w_from_L
initialize_qgl0f0w_from_L <- function(X, L){
  K = ncol(L)
  w = replicate(K, 1)
  init_tmp = ebpmf.alpha::initialize_qgl0f0_from_L(X = X, L = L)
  return(list(qg = init_tmp$qg, l0 = init_tmp$l0, f0 = init_tmp$f0, w = w))
}

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


## ================================================================================================================
## Functions for computing ELBO
## ================================================================================================================

## compute KL divergence between prior and posterior from `ebpm` outputs

#' @export compute_kl_ebpm
compute_kl_ebpm <- function(y,s, posterior, ll){
  mask <- (y != 0)
  E_loglik = - sum(s * posterior$mean) + sum(y[mask] * log(s[mask])) + sum(y[mask]*posterior$mean_log[mask]) - sum(lgamma(y[mask] + 1))
  return(E_loglik - ll)
}


compute_elbo_bg <- function(l0, f0, qg, b, a, d, const){
  KL = sum(qg$kl_l) + sum(qg$kl_f)
  elbo = - sum(colSums(l0 * qg$qls_mean) * colSums(f0 * qg$qfs_mean) ) +
            sum(d$x * (log(l0[d$i]) + log(f0[d$j]) + b + a) ) - KL - const
  return(elbo)
}


compute_elbo_wbg <- function(w, l0, f0, qg, b, a, d, const){
	KL = sum(qg$kl_l) + sum(qg$kl_f)
	elbo = - sum(w * colSums(l0 * qg$qls_mean) * colSums(f0 * qg$qfs_mean) ) +
            sum(d$x * (log(l0[d$i]) + log(f0[d$j]) + b + a) ) - KL - const
	return(elbo)
}

## ================================================================================================================
## MLE for Poisson Means, with the same format as `ebpm`
## ================================================================================================================

#' @export mle_pm
mle_pm <- function(x, s, g_init, fix_g){
	mask <- (x != 0)
	mle <- replicate(length(x),0)
	mle[mask] <- x/s
	mle[!mask] <- 1e-8 ## set 0s to be some small numbers
	posterior <- list(mean = mle, mean_log = log(mle))
	log_likelihood <- - sum(s * posterior$mean) + sum(x[mask] * log(s[mask])) + sum(x[mask]*posterior$mean_log[mask])- sum(lgamma(x[mask] + 1))
	out = list(fitted_g = list(NULL), posterior = posterior, log_likelihood = log_likelihood)
}


## ================================================================================================================
## Other Functions
## ================================================================================================================

compute_rmse <- function(lam1, lam2){
  return(sqrt(mean((lam1 - lam2)^2)))
}


# Apply operation f to all nonzeros of a sparse matrix.
apply.nonzeros <- function (X, f) {
  d <- summary(X)
  return(sparseMatrix(i = d$i,j = d$j,x = f(d$x),dims = dim(X)))
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



check_progress_wbg <- function(elbo_prev, which_part,k, 
															 w_log, l0, f0, qg, b, a, d, const){
	## check progress
  elbo = compute_elbo_wbg(w = exp(w_log), l0 = l0, f0 = f0, qg = qg,
                          b = b, a = a, d = d, const = const)
  if(elbo < elbo_prev){
					print(sprintf("k = %d, %s updated, elbo_diff = %f", 
												k, which_part, elbo - elbo_prev))
  }
	return(elbo)
}


