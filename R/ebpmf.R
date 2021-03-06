#' @title Empirical Bayes Poisson Matrix Factorization
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
									 fix_option = list(gl = FALSE, ql = FALSE,
                                     gf = FALSE, qf = FALSE),
									 maxiter = 100,
                   tol = 1e-8, verbose = FALSE){
  ## transform to sparse matrix ## TODO: when we prefer to use dense matrix?
  X <- as(X, "sparseMatrix") ## TODO: which format is better?
  d = summary(X)
  const = sum(apply.nonzeros(X = X, f = function(x) lgamma(x + 1)))
	## initialization
  init_tmp <- init_ebpmf(X = X, K = K, init = init, d = d)
  qg <- init_tmp$qg
  b <- init_tmp$b
  a <- init_tmp$a
  rm(init_tmp)
  ## update iteratively
  ELBOs <- c()
  KLs <- c()
  for(i in 1:maxiter){
		b_k_max = replicate(length(d$x), 0) ## max b_k
    for(k in 1:K){
      ## compute Ez 
      ## use q to compute Ez; output list(rs, cs, B_k)
      b_k = qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k] - a
			Ez <- compute_EZ(d = d, b = b, b_k = b_k)
      ## rank1 update 
      ## use Ez to upudate q, g; output (list(B, kl_l, kl_f, qg))
			rank1_qg <- rank1(d = d, X_rs = Ez$rs, X_cs = Ez$cs,
                        pm_func = pm_func, pm_control = pm_control,
												ql = list(mean = qg$qls_mean[,k],
                                  mean_log = qg$qls_mean_log[,k]),
                        qf = list(mean = qg$qfs_mean[,k],
                                  mean_log = qg$qfs_mean_log[,k]),
                        gl = qg$gls[[k]],
                        gf = qg$gfs[[k]],
                        kl_l = qg$kl_l[k],
                        kl_f = qg$kl_f[k],
                        fix_option = fix_option)
      rm(Ez)
      qg = update_qg(rank1_qg, qg, k)
      rm(rank1_qg)
    	b_k0 = b_k
      b_k = qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k] - a
      b = log( exp(b) - exp(b_k0) + exp(b_k)  )
      b_k_max = pmax(b_k, b_k_max)
		}
    ## compute ELBO
    KL = sum(qg$kl_l) + sum(qg$kl_f)
    ELBO = - sum( colSums(qg$qls_mean) * colSums(qg$qfs_mean) ) + sum(d$x * (b + a)) - KL - const
    ELBOs <- c(ELBOs, ELBO)
		KLs <- c(KLs, KL)
		## update a & b
    #a0 = a
    #a = b_k_max + a0
    #b = (b + a0) - a
    b = b - b_k_max
    a = b_k_max + a
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
  return(list(qg = qg, ELBO = ELBOs, KL = KLs))
}


#' @title Empirical Bayes Poisson Matrix Factorization (rank 1)
#' @import ebpm

#' @export  rank1

rank1 <- function(d, X_rs, X_cs, 
									pm_func,pm_control,
									ql, gl, kl_l,
									qf, gf, kl_f,
									fix_option){
  p = length(X_cs)
  n = length(X_rs)

	if(!fix_option$qf){
    s <- sum(ql$mean) 
    fit = do.call(pm_func,
                  c(list(x = X_cs, s = s, g_init = gf, fix_g = fix_option$gf), pm_control))
    qf = fit$posterior
    gf = fit$fitted_g
    kl_f = compute_kl_ebpm(y = X_cs, s = replicate(p, s), posterior = qf, ll = fit$log_likelihood)
    rm(fit)
  }
	## fit for l, and compute kl_l
  if(!fix_option$ql){
    s = sum(qf$mean) 
    fit = do.call(pm_func,
                  c(list(x = X_rs, s = s, g_init = gl, fix_g = fix_option$gl), pm_control))
    ql = fit$posterior
    gl = fit$fitted_g
    kl_l = compute_kl_ebpm(y = X_rs, s = replicate(n, s), posterior = ql, ll = fit$log_likelihood)
    rm(fit)
  }
	qg = list(ql = ql, gl = gl, kl_l = kl_l, qf = qf, gf = gf, kl_f = kl_f)
	return(qg)
}

