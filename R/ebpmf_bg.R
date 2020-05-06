#' @title Empirical Bayes Poisson Matrix Factorization (Background Model)
#' @import ebpm
#' @import Matrix

#' @param X count matrix (dim(X) = c(n, p)).
#' @param k number of topics
#' @param pm_func function for solving the \code{ebpm} subproblem; can be
#' \code{ebpm_point_gamma, ebpm_two_gamma, ebpm_exponential_mixture, ebpm_gamma_mixture_single_scale}
#' @param pm_control control parameters for pm_func function
#' @param init list(qg, init_method, init_iter)
#' @param fix_g list(l, f) where l, f are either TRUE or FALSE
#' @param maxiter maximum number of iterations
#' @param tol stopping tolerance for ELBO
#' @param verbose print progress if set TRUE
#'
#' @return A list containing elements:
#'     \describe{
#'       \item{\code{l0}}{sample-wise mean}
#'       \item{\code{f0}}{feature-wise mean}
#'       \item{\code{qg}}{approximate posterior for l}
#'       \item{\code{ELBO}}{ELBO objective for this VEB algorithm}
#'      }
#' @examples
#' To add
#' @export  ebpmf_bg

library(ebpm)
library(ebpmf.alpha)
library(Matrix)
source("~/Desktop/git/ebpmf.alpha/r/util.R")
ebpmf_bg <- function(X, K, pm_func = ebpm::ebpm_gamma_mixture,
                   init = list(qg = NULL, init_method = "scd", init_iter = 20), pm_control = NULL,
                   fix_g = list(l = FALSE, f = FALSE), maxiter = 100,
                   tol = 1e-8, verbose = FALSE){
	## TODO: input check, require X_rs, X_cs to be nonzero

  ## transform to sparse matrix 
	X_rs = rowSums(X)
	X_cs = colSums(X)
  X <- as(X, "sparseMatrix") 
  d = summary(X)
  #const = sum(apply.nonzeros(X = X + 1, f = lgamma))
  const = 0 ## TODO: use sum(lgamma(X + 1))

  ## initialization
  init_tmp <- init_ebpmf_bg(X = X, K = K, init = init, d = d)
  qg <- init_tmp$qg
  B <- init_tmp$B
	Lam <- init_tmp$Lam
	l0 <- init_tmp$l0
	f0 <- init_tmp$f0
	l0f0 = f0[d$j] * l0[d$i]
  rm(init_tmp)

  ## update iteratively
  ELBOs <- c()
  for(i in 1:maxiter){
    KL <- 0
    for(k in 1:K){
			print(k)
			if(i == 2 && k == 3){browser()}
			## B_k, Lam_k
 			B_k = exp(qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k]) * l0f0
  		Lam_k = (l0 * qg$qls_mean[,k]) %o% (f0 * qg$qfs_mean[,k])

    	## compute q(Z)
      Ez <- compute_EZ_bg(d = d,B,B_k)

      ## update (qL, gL, qF, gF) 
      init_r1 = list(sf = sum(qg$qls_mean[,k]),
                     gl = qg$gls[[k]], gf = qg$gfs[[k]])
      rank1_tmp <- rank1_bg(d = d, X_rs = Ez$rs, X_cs = Ez$cs,
														l0 = l0, f0 = f0, pm_func = pm_func, pm_control = pm_control,
														init = init_r1, fix_g = fix_g)
      rm(Ez)
      KL = KL + rank1_tmp$kl_l + rank1_tmp$kl_f
      qg = update_qg(rank1_tmp$qg, qg, k)
      rm(rank1_tmp)

			## update B, Lam
			B = B - B_k + exp(qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k]) * l0f0
			Lam = Lam - Lam_k + (l0 * qg$qls_mean[,k]) %o% (f0 * qg$qfs_mean[,k])
    }
		## update l0, f0
		## TODO : fix bug here (X_rs[84] = 1, while denom = 0, because of numerical issue)
    f0 <- (X_cs/colSums(Lam)) * f0
    l0 <- (X_rs/rowSums(Lam)) * l0
    l0f0 = f0[d$j] * l0[d$i]

    ## compute ELBO
    ELBO = -sum( f0 * colSums(l0*Lam) ) + sum(d$x * log(B) ) - KL - const
		ELBOs <- c(ELBOs, ELBO)
    ## verbose
    if(verbose){
      print("iter         ELBO")
      print(sprintf("%d:    %f", i, ELBO))
    }
    ## check convergence
    # diff = ifelse(i > 2, ELBOs[i] - ELBOs[i-1], Inf)
    # if(diff < tol){
    #   if(verbose){print(sprintf("reaches tol %f in %d iterations", tol, i))}
    #   break
    # }
  }
  return(list(l0 = l0, f0 = f0, qg = qg, ELBO = ELBOs))
}




initialize_qg_bg <- function(X, K, init_method = "scd", init_iter = 20, seed = 123){
  set.seed(seed)
  nnmf_fit = NNLM::nnmf(A = as.matrix(X), k = K, loss = "mkl", max.iter = init_iter, verbose = F, method = init_method)
  qls_mean = nnmf_fit$W
  qls_mean[qls_mean == 0] = 1e-10
  qfs_mean = t(nnmf_fit$H)
  qfs_mean[qfs_mean == 0] = 1e-10
  qls_mean_log = log(qls_mean)
  qfs_mean_log = log(qfs_mean)
	#al = c(seq(0.01, 0.10, 0.01), seq(0.2, 0.9, 0.1), seq(1,15,2), 20, 50, 75, 100, 200, 1e3)
	al = c(seq(0.01, 0.10, 0.01), seq(0.2, 0.9, 0.1), seq(1,15,2), 20, 50)
	D = length(al)
	g = gammamix(pi = replicate(D, 1/D), shape = al, scale = 1/al)
  qg = list(qls_mean = qls_mean, qls_mean_log =qls_mean_log,
            qfs_mean = qfs_mean, qfs_mean_log =qfs_mean_log,
            gls = rep(list(g), K),gfs = rep(list(g), K))
  return(qg)
}


## output: qg, B
init_ebpmf_bg <- function(X,K, init, d){
  start = proc.time()
  qg = init$qg
  if(is.null(qg)){
    qg = initialize_qg_bg(X, K, init_method =  init$init_method, init_iter = init$init_iter)
  }
	f0 = apply(qg$qfs_mean ,1, median)
	f0[f0 == 0] = 1e-10
	l0 = apply(qg$qls_mean ,1, median)
	l0[l0 == 0] = 1e-10
	l0f0 = l0[d$i] * f0[d$j]
  ## TODO: speedup
  B = exp(qg$qls_mean_log[d$i, 1] + qg$qfs_mean_log[d$j, 1]) * l0f0
  Lam = (l0 * qg$qls_mean[,1]) %o% (f0 * qg$qfs_mean[,1])
  for(k in 2:K){
    B <- B + exp(qg$qls_mean_log[d$i, k] + qg$qfs_mean_log[d$j, k]) * l0f0
    Lam <- Lam + (l0 * qg$qls_mean[,k]) %o% (f0 * qg$qfs_mean[,k])
  }
  runtime = proc.time() - start
  print(runtime)
  return(list(l0 = l0, f0 = f0, qg = qg, B = B, Lam = Lam))
}

compute_EZ_bg <- function(d, B, B_k){
  Ez.val = replicate(length(d$i), 0)
  mask <- (B != 0)
  Ez.val[mask] = d$x[mask] * B_k[mask]/B[mask]
  Ez = sparseMatrix(i = d$i, j = d$j, x = Ez.val)
  return(list(rs = rowSums(Ez), cs = colSums(Ez)))
}

#' @title Empirical Bayes Poisson Matrix Factorization (rank 1)
#' @import ebpm

#' @param d: summary(X), so it has `i`, `j` and `x` as attributes
#' @param X_rs, rowSums of X,  count matrix (dim(X) = c(n, p)).
#' @param X_cs, colSums of X,  count matrix (dim(X) = c(n, p)).
#' @param pm_func function for solving the \code{ebpm} subproblem
#' @param pm_control control parameters for pm_func function
#' @param init list(sf, gl, gf), where sf is the scale for ebpm problem for F,
#' gl, gf are initialization for gl, gf; all can be set to NULL
#' @param fix_g list(l, f) where l, f are either TRUE or FALSE
#'
#' @return A list containing elements:
#'     \describe{
#'       \item{\code{ql}}{approximate posterior for l}
#'       \item{\code{gl}}{fitted g for l}
#'       \item{\code{kl_l}}{kl divergence between q and g for l}
#'       \item{\code{qf}}{approximate posterior for f}
#'       \item{\code{gf}}{fitted g for f}
#'       \item{\code{kl_f}}{kl divergence between q and g for f}
#'      }
#' @examples
#' To add
#' @export  rank1

## TODO: uupdate KL, Lam
rank1_bg <- function(d, X_rs, X_cs, l0, f0, pm_func,pm_control, init, fix_g){
  p = length(X_cs)
  n = length(X_rs)
  ## initialization (in fact, any non-negative number well do)
  sf = init$sf
  if(is.null(sf)){
    sf = sum(NNLM::nnmf(A = X, k = 1, loss = "mkl", method = "lee", max.iter = 1, verbose = F)$W[,1])
  }
  ## fit for f, and compute kl_f
  fit_f = do.call(pm_func, c(list(x = X_cs, s = sf*f0, g_init = init$gf), pm_control))
  kl_f = compute_kl(X_cs, sf*f0, fit_f)
  ## fit for l, and compute kl_l
  sl = sum(fit_f$posterior$mean)
  fit_l = do.call(pm_func, c(list(x = X_rs, s = sl*l0, g_init = init$gl), pm_control))
  kl_l = compute_kl(X_rs, sl*l0, fit_l)
  ## list to return
  qg = list(ql = fit_l$posterior, gl = fit_l$fitted_g, qf = fit_f$posterior, gf = fit_f$fitted_g)
  out = list(qg = qg, kl_l = kl_l, kl_f = kl_f)
  return(out)
}







