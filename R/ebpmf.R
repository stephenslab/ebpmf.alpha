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
                   fix_g = list(l = FALSE, f = FALSE), maxiter = 100,
                   tol = 1e-8, verbose = FALSE){
  ## transform to sparse matrix ## TODO: when we prefer to use dense matrix?
  X <- as(X, "sparseMatrix") ## TODO: which format is better?
  d = summary(X)
  #const = sum(apply.nonzeros(X = X + 1, f = lgamma))
  const = 0 ## TODO: use sum(lgamma(X + 1))

  ## initialization
  init_tmp <- init_ebpmf(X = X, K = K, init = init, d = d)
  qg <- init_tmp$qg
  B <- init_tmp$B
  rm(init_tmp)

  ## update iteratively
  ELBOs <- c()
  for(i in 1:maxiter){
    KL <- 0
    for(k in 1:K){
      ## compute Ez (list(rs, cs, B))
      Ez <- compute_EZ(B = B, d = d,
                       ql_log = qg$qls_mean_log[,k], qf_log = qg$qfs_mean_log[,k])
      ## rank1 update (list(D, kl_l, kl_f, qg))
      init_r1 = list(sf = sum(qg$qls_mean[,k]),
                     gl = qg$gls[[k]], gf = qg$gfs[[k]])
      rank1_tmp <- rank1(d = d, X_rs = Ez$rs, X_cs = Ez$cs,
                         pm_func = pm_func, pm_control = pm_control,
                         init = init_r1, fix_g = fix_g,
                         B = B, B_k = Ez$B_k)
      rm(Ez)
      B = rank1_tmp$B
      KL = KL + rank1_tmp$kl_l + rank1_tmp$kl_f
      qg = update_qg(rank1_tmp$qg, qg, k)
      rm(rank1_tmp)
    }
    ## compute ELBO
    ELBO = - sum( colSums(qg$qls_mean) * colSums(qg$qfs_mean) ) + sum(d$x * log(B)) - KL - const
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
  return(list(qg = qg, ELBO = ELBOs))
}

## output: qg, B
init_ebpmf <- function(X,K, init, d){
  start = proc.time()
  qg = init$qg
  if(is.null(qg)){
    qg = initialize_qg(X, K, init_method =  init$init_method, init_iter = init$init_iter)
  }
  ## TODO: speedup
  B = exp(qg$qls_mean_log[d$i, 1] + qg$qfs_mean_log[d$j, 1])
  for(k in 2:K){
    B <- B + exp(qg$qls_mean_log[d$i, k] + qg$qfs_mean_log[d$j, k])
  }
  runtime = proc.time() - start
  print(runtime)
  return(list(qg = qg, B = B))
}

compute_EZ <- function(d, B, ql_log, qf_log){
  B_k = exp(ql_log[d$i] + qf_log[d$j])
  Ez.val = replicate(length(d$i), 0)
  mask <- (B != 0)
  Ez.val[mask] = d$x[mask] * B_k[mask]/B[mask]
  Ez = sparseMatrix(i = d$i, j = d$j, x = Ez.val)
  return(list(rs = rowSums(Ez), cs = colSums(Ez), B_k = B_k))
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

rank1 <- function(d, X_rs, X_cs, pm_func,pm_control, init, fix_g, B, B_k){
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
  ## update B
  B = B - B_k + exp(fit_l$posterior$mean_log[d$i] + fit_f$posterior$mean_log[d$j])
  ## list to return
  qg = list(ql = fit_l$posterior, gl = fit_l$fitted_g, qf = fit_f$posterior, gf = fit_f$fitted_g)
  out = list(qg = qg, kl_l = kl_l, kl_f = kl_f, B = B)
  return(out)
}

