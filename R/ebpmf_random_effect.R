#' @title Empirical Bayes Poisson Matrix Factorization with Prior Family Point Gamma and with random effect in  mean
#' @description Uses Empirical Bayes to fit the model \deqn{X_{ij}  ~ Poi(\sum_k L_{ik} F_{jk})} with \deqn{L_{.k} ~ g_k()}, where dim(X) = c(n,p); dim(L) = c(n,K); dim(F) = c(p,K)
#' with g_k being either Mixture of Exponential, or Point Gamma
#' @import mixsqp
#' @import ebpm
#' @import gtools
#' @import NNLM

#' @details The model is fit in 2 stages: i) estimate \eqn{g} by maximum likelihood (over pi_k)
#' ii) Compute posterior distributions for \eqn{\lambda_j} given \eqn{x_j,\hat{g}}.
#' @param X count matrix (dim(X) = c(n, p)).
#' @param K number of topics
#' @param  maxiter.out  maximum iterations in the outer  loop
#' @param  maxiter.int  maximum iterations in the inner loop
#' @param verbose T if print ELBO
#' @param init_method used in \code{NNLM::nnmf}. Either `scd` or `lee`
#' @param seed random seed
#'
#' @return A list containing elements (this is freeuqntly updated):
#'     \describe{
#'       \item{\code{qls_mean}}{A n by k matrix: Approximate posterior mean for L}
#'       \item{\code{qls_mean_log}}{A n by k matrix: Approximate posterior log mean for L}
#'       \item{\code{gls}}{A list of K elements, each element is the estimated prior for the kth column of L}
#'       \item{\code{qfs_mean}}{A p by k matrix: Approximate posterior mean for F}
#'       \item{\code{qfs_mean_log}}{A p by k matrix: Approximate posterior log mean for F}
#'       \item{\code{gfs}}{A list of K elements, each element is the estimated prior for the kth column of F}
#'      }
#'
#'
#' @examples
#' To add

#'
#' @export ebpmf_random_effect
#'
ebpmf_random_effect <- function(X, K, qg = NULL, maxiter.out = 10, fix_g = F, verbose = F, init_method = "scd",
                              seed = 123, Lam_true = NULL, pi0_l = "estimate", pi0_f = "estimate", prior_mu = "point"){
  if(identical(pi0_l, "estimate")){pi0_l = replicate(K, "estimate")}
  if(identical(pi0_f, "estimate")){pi0_f = replicate(K, "estimate")}
  set.seed(seed)
  ## init from NNLM::nnmf result
  start = proc.time()
  if(is.null(qg)){
    qg = initialize_qg_random_effect(X, K, init_method)
  }
  runtime_init = (proc.time() - start)[[3]]

  runtime_rank1 = 0
  runtime_ez = 0
  ELBOs = c()
  KLs = c()
  lls = c()
  rmses = c()

  start = proc.time()
  #browser()
  tmp = get_Ez_random_effect(X, qg, K)
  Ez = tmp$Ez
  zeta = tmp$zeta
  rm(tmp)
  runtime_ez = runtime_ez + (proc.time() - start)[[3]]

  for(iter in 1:maxiter.out){
    KL = 0
    for(k in 1:K){
      ## update q, g
      start = proc.time()
      init_l  = list(mean = qg$qls_mean[,k])
      if(is.null(qg$gls[[k]])){
        tmp = ebpmf_rank1_point_gamma_helper(X = Ez[,,k],init = init_l, pi0_l = pi0_l[k], pi0_f = pi0_f[k])
      }else{
        if(fix_g){
          tmp = ebpmf_rank1_point_gamma_helper(X = Ez[,,k],init = init_l,
                                               gl_init = qg$gls[[k]], gf_init = qg$gfs[[k]], fix_gl = T, fix_gf = T, pi0_l = pi0_l[k], pi0_f = pi0_f[k])
        }else{
          tmp = ebpmf_rank1_point_gamma_helper(X = Ez[,,k],init = init_l, gl_init = qg$gls[[k]], gf_init = qg$gfs[[k]], pi0_l = pi0_l[k], pi0_f = pi0_f[k])
        }
      }
      KL = KL + tmp$kl_l + tmp$kl_f
      qg = update_qg(tmp, qg, k)
      ## update Z
      start = proc.time()
      tmp = get_Ez_random_effect(X, qg, K)
      Ez = tmp$Ez
      zeta = tmp$zeta
      rm(tmp)
      runtime_ez = runtime_ez + (proc.time() - start)[[3]]
    }
    ## update qg for mu
    tmp = ebpmf_update_random_effect(Ez[,,K+1], ebpm_fn = prior_mu)
    KL = KL + tmp$kl
    qg = update_qg_random_effect(tmp, qg)

    ll = compute_ll_random_effect(X, qg)
    ELBO = ll - KL
    KLs <- c(KLs, KL)
    ELBOs <- c(ELBOs, ELBO)
    lls <- c(lls, ll)
    if(!is.null(Lam_true)){
      rmse = compute_rmse(Lam_true, qg$qls_mean %*% t(qg$qfs_mean))
      rmses = c(rmses, rmse)
    }
    if(verbose){
      print("iter         ELBO")
      print(sprintf("%d:    %f", iter, ELBO))
      if(!is.null(Lam_true)){
        print(sprintf("%d:    %f", iter, rmse))
      }
    }
  }
  return(list(qg = qg, ELBO = ELBOs, KL = KLs, ll = lls, RMSE = rmses))
}

ebpmf_update_random_effect <- function(X, ebpm_fn = "gamma"){
  if(identical(ebpm_fn, "gamma")){
    fit = ebpm::ebpm_point_gamma(x = as.vector(X), pi0 = 0)
  }
  if(identical(ebpm_fn, "exponential_mix")){
    fit = ebpm::ebpm_exponential_mixture(x = as.vector(X))
  }
  gmu = fit$fitted_g
  qmu_mean = fit$posterior$mean
  qmu_mean_log = fit$posterior$mean_log
  dim(qmu_mean) = dim(X)
  dim(qmu_mean_log) = dim(X)
  kl = compute_kl(x = as.vector(X), s = replicate(length(as.vector(X)), 1), fit)
  return(list(kl = kl, gmu = gmu, qmu_mean = qmu_mean, qmu_mean_log = qmu_mean_log))
}
