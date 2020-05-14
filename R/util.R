## ================================================================================================================
## Functions updating in ebpmf
## ================================================================================================================

compute_EZ <- function(d, B, B_k){
  Ez.val = replicate(length(d$i), 0)
  mask <- (B != 0)
  Ez.val[mask] = d$x[mask] * B_k[mask]/B[mask]
  Ez = sparseMatrix(i = d$i, j = d$j, x = Ez.val)
  return(list(rs = Matrix::rowSums(Ez), cs = Matrix::colSums(Ez)))
}

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
init_ebpmf <- function(X,K, init, d){
  qg = init$qg
  if(is.null(qg)){
    qg = initialize_qg(X, K, init_method =  init$init_method, init_iter = init$init_iter)
  }
  ## TODO: speedup
  B = exp(qg$qls_mean_log[d$i, 1] + qg$qfs_mean_log[d$j, 1])
  for(k in 2:K){
    B <- B + exp(qg$qls_mean_log[d$i, k] + qg$qfs_mean_log[d$j, k])
  }
  return(list(qg = qg, B = B))
}

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
  g = gammamix(pi = replicate(D, 1/D), shape = aL, scale = 1/aL)
  qg$gls = replicate(K, list(g))
  qg$gfs = replicate(K, list(g))
	return(list(qg = qg, l0 = l0, f0 = f0))
}

## output: qg, B
## init either NULL, or list(qg, l0, f0)
init_ebpmf_bg <- function(X,K, init, d, seed = 123){
  n = nrow(X)
	p = ncol(X)
  if(is.null(init)){
		init_tmp = initialize_qg_l0f0(X = X, K = K, seed = seed)
		qg = init_tmp$qg
		l0 = init_tmp$l0
		f0 = init_tmp$f0
	}
  ## TODO: speedup
  B = exp(qg$qls_mean_log[d$i, 1] + qg$qfs_mean_log[d$j, 1])
  for(k in 2:K){
    B <- B + exp(qg$qls_mean_log[d$i, k] + qg$qfs_mean_log[d$j, k])
  }
  return(list(l0 = l0, f0 = f0, qg = qg, B = B))
}

## ================================================================================================================
## Functions for computing ELBO
## ================================================================================================================

## compute KL divergence between prior and posterior from `ebpm` outputs
compute_kl_ebpm <- function(y,s, posterior, ll){
  mask <- (y != 0)
  E_loglik = - sum(s * posterior$mean) + sum(y[mask] * log(s[mask])) + sum(y[mask]*posterior$mean_log[mask]) - sum(lgamma(y[mask] + 1))
  return(E_loglik - ll)
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
