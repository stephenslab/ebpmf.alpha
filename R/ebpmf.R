#' @title Empirical Bayes Poisson Matrix Factorization
#' @import ebpm

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
#'       \item{\code{ql}}{approximate posterior for l}
#'       \item{\code{gl}}{fitted g for l}
#'       \item{\code{kl_l}}{kl divergence between q and g for l}
#'       \item{\code{qf}}{approximate posterior for f}
#'       \item{\code{gf}}{fitted g for f}
#'       \item{\code{kl_f}}{kl divergence between q and g for f}
#'      }
#' @examples
#' To add
#' @export  ebpmf

ebpmf <- function(X, K,pm_func = ebpm::ebpm_point_gamma,
                  init = list(qg = NULL, init_method = "scd", init_iter = 20), pm_control = NULL,
                  fix_g = list(l = FALSE, f = FALSE), maxiter = 100,
                  tol = 1e-8, verbose = FALSE){
  ELBOs = c()
  # initialize: get q, g, Ez
  init_tmp = ebpmf_init(X, K, init)
  qg = init_tmp$qg
  Ez = init_tmp$Ez ## todo: get rowsum, colsum
  rm(init_tmp)
  # update iteratively
  for(i in 1:maxiter){
    #if(i == 30){browser()}
    KL = 0
    for(k in 1:K){
      X_rs = rowSums(Ez[,,k])
      X_cs = colSums(Ez[,,k])
      init_r1 = list(sf = sum(qg$qls_mean[,k]), gl = qg$gls[[k]], gf = qg$gfs[[k]])
      ### rank-1 update for qg
      tmp = ebpmf_rank1(X_rs = X_rs, X_cs = X_cs, pm_func = pm_func, pm_control = pm_control, init = init_r1, fix_g = fix_g)
      qg = update_qg(tmp, qg, k)
      ### get Ez
      Ez = get_Ez(X, qg, K)$Ez
      KL = KL + tmp$kl_l + tmp$kl_f
    }
    ## compute ELBO
    ELBO = compute_ll(X, qg) - KL
    ELBOs <- c(ELBOs, ELBO)
    ## verbose
    if(verbose){
      print("iter         ELBO")
      print(sprintf("%d:    %f", i, ELBO))
    }
    ## check convergence
    diff = ifelse(i > 2, ELBOs[i] - ELBOs[i-1], Inf)
    if(diff < tol){
      if(verbose){print(sprintf("reaches tol %f in %d iterations", tol, i))}
      break
    }
  }
  return(list(qg = qg, ELBO = ELBOs))
}


ebpmf_init <- function(X, K, init){
  qg = init$qg
  if(is.null(qg)){
    qg = initialize_qg(X, K, init_method =  init$init_method, init_iter = init$init_iter)
  }
  Ez = get_Ez(X, qg, K)$Ez
  return(list(qg = qg, Ez = Ez))
}

#' @title Empirical Bayes Poisson Matrix Factorization (rank 1)
#' @import ebpm

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
#' @export  ebpmf_rank1

ebpmf_rank1 <- function(X_rs, X_cs, pm_func,pm_control, init, fix_g){
	p = length(X_cs)
  n = length(X_rs)
  ## initialization (in fact, any non-negative number well do)
  sf = init$sf
	if(is.null(sf)){
		sf = sum(NNLM::nnmf(A = X, k = 1, loss = "mkl", method = "lee", max.iter = 1, verbose = F)$W[,1])
	}
	## fit for f, and compute kl_f
	fit_f = do.call(pm_func, c(list(x = X_cs, s = sf, g_init = init$gf), pm_control))
  kl_f = compute_kl(X_cs, sf, fit_f)
	## fit for l, and compute kl_l
	sl = sum(fit_f$posterior$mean)
	fit_l = do.call(pm_func, c(list(x = X_rs, s = sl, g_init = init$gl), pm_control))
  kl_l = compute_kl(X_rs, sl, fit_l)
  ## list to return
	out = list(ql = fit_l$posterior, gl = fit_l$fitted_g,  kl_l = kl_l, qf = fit_f$posterior, gf = fit_f$fitted_g, kl_f = kl_f)
	return(out)
}

################## testing code ##########################

# ## test code for rank-1
# n = 100; p = 300
# a = 100; b = 10
# l = rgamma(n, shape = a, rate = b)
# f = rgamma(p, shape = a, rate = b)
# lam = outer(l, f, "*")
# X = matrix(rpois(n*p, lam), nrow = n)
#
# library(NNLM)
# fit_nnmf = NNLM::nnmf(X, k = 1, loss = "mkl", method = "lee")
#
# fit_ebpmf_r1 = ebpmf_rank1(X_rs = rowSums(X),X_cs = colSums(X), pm_func = ebpm::ebpm_point_gamma, pm_control = list(NULL), init = list(sf = NULL, gl = NULL, gf = NULL))
# plot(fit_ebpmf_r1$ql$mean, fit_nnmf$W)
# plot(fit_ebpmf_r1$qf$mean, fit_nnmf$H)


## test code for ebpmf
# n = 100; p = 300
# a = 100; b = 10
# k = 2
# l = matrix(rgamma(n*k, shape = a, rate = b), ncol = k)
# f = matrix(rgamma(p*k, shape = a, rate = b), ncol = k)
# lam = l %*% t(f)
# X = matrix(rpois(n*p, lam), nrow = n)

# library(NNLM)
# lf_init = NNLM::nnmf(X, k = k, loss = "mkl", method = "lee", max.iter = 3)
# L0 = lf_init$W
# F0 = t(lf_init$H)
# qg0 = initialize_qg_from_LF(L0 = L0, F0 = F0)


# fit_ebpmf = ebpmf(X = X, K = k)

#
# init_ = list(qg = qg0, init_method = "scd", init_iter = 20)
# start = proc.time()
# fit_ebpmf = ebpmf(X = X, K = k, pm_func = ebpm::ebpm_point_gamma, pm_control = NULL, init = init_, maxiter = 50, tol = 1e-3)
# t1 = proc.time() - start
#
#
# start = proc.time()
# fit_pg = ebpmf.alpha::ebpmf_point_gamma(X = X, K = k, qg = qg0, maxiter.out = 50)
# t2 = proc.time() - start
#
# identical(fit_ebpmf$qg$qls_mean, fit_pg$qg$qls_mean)
#
# t1
# t2

# n = 100; p = 300
# a = 100; b = 10
# k = 2
# l = matrix(rgamma(n*k, shape = a, rate = b), ncol = k)
# f = matrix(rgamma(p*k, shape = a, rate = b), ncol = k)
# lam = l %*% t(f)
# X = matrix(rpois(n*p, lam), nrow = n)
#
# fit_ebpmf = ebpmf(X = X, K = k, pm_func = ebpm::ebpm_gamma_mixture_single_scale)
