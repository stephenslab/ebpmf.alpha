#' @title Empirical Bayes Poisson Matrix Factorization
#' @description Uses Empirical Bayes to fit the model \deqn{X_{ij}  ~ Poi(\sum_k L_{ik} F_{jk})} with \deqn{L_{.k} ~ g_k()}
#' with g_k being either Mixture of Exponential, or Point Gamma
#' @import mixsqp
#' @import ebpm

#' @details The model is fit in 2 stages: i) estimate \eqn{g} by maximum likelihood (over pi_k)
#' ii) Compute posterior distributions for \eqn{\lambda_j} given \eqn{x_j,\hat{g}}.
#' @param X count matrix (dim(X) = c(n, p)).
#' @param K number of topics
#' @param m multiplicative parameter for selecting grid in "ebpm::ebpm_exponential_mixture"

#' @param  maxiter.out  maximum iterations in the outer  loop
#' @param  maxiter.int  maximum iterations in the inner loop
#' @param seed random seed
#'
#'
#' @return A list containing elements:
#'     \describe{
#'       \item{\code{qls_mean}}{A n by k matrix: Approximate posterior mean for L}
#'       \item{\code{qls_mean_log}}{A n by k matrix: Approximate posterior log mean for L}
#'       \item{\code{gls}}{A list of K elements, each element is the estimated prior for the kth column of L}
#'       \item{\code{qfs_mean}}{A p by k matrix: Approximate posterior mean for F}
#'       \item{\code{qfs_mean_log}}{A p by k matrix: Approximate posterior log mean for F}
#'       \item{\code{gfs}}{A list of K elements, each element is the estimated prior for the kth column of F}
#'      }
#' @examples
#' To add

#'
#' @export  ebpmf_exponential_mixture

ebpmf_exponential_mixture <- function(X, K, maxiter.out = 10, fix_g = F, fix_grid = F, verbose = F, m = 2, seed = 123){
  set.seed(seed)
  ## init from NNLM::nnmf result
  start = proc.time()
  qg = initialize_qg(X, K)
  runtime_init = (proc.time() - start)[[3]]

  runtime_rank1 = 0
  runtime_ez = 0
  ELBOs = c()
  KLs = c()
  lls = c()

  start = proc.time()
  tmp = get_Ez(X, qg, K)
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
        tmp = ebpmf_rank1_exponential_helper(X = Ez[,,k],init = init_l,m = m)
      }else{
        if(fix_g){
          tmp = ebpmf_rank1_exponential_helper(X = Ez[,,k],init = init_l,m = m,
                                               gl_init = qg$gls[[k]], gf_init = qg$gfs[[k]], fix_gl = T, fix_gf = T)
        }else{
          if(fix_grid){
            tmp = ebpmf_rank1_exponential_helper(X = Ez[,,k],init = init_l,m = m,
                                                 scale_l = list(a = qg$gls[[k]]$a, b = qg$gls[[k]]$b), scale_f = list(a = qg$gfs[[k]]$a, b = qg$gfs[[k]]$b))
          }else{
            tmp = ebpmf_rank1_exponential_helper(X = Ez[,,k],init = init_l,m = m)
          }
        }
      }
      KL = KL + tmp$kl_l + tmp$kl_f
      qg = update_qg(tmp, qg, k)
      ## update Z
      start = proc.time()
      tmp = get_Ez(X, qg, K)
      Ez = tmp$Ez
      zeta = tmp$zeta
      rm(tmp)
      runtime_ez = runtime_ez + (proc.time() - start)[[3]]
    }
    ll = compute_ll(X, qg)
    ELBO = ll - KL
    KLs <- c(KLs, KL)
    ELBOs <- c(ELBOs, ELBO)
    lls <- c(lls, ll)
    if(verbose){
      print("iter         ELBO")
      print(sprintf("%d:    %f", iter, ELBO))
    }
  }
  # print("summary of  runtime:")
  # print(sprintf("init           : %f", runtime_init))
  # print(sprintf("Ez     per time: %f", runtime_ez/(iter*K)))
  # print(sprintf("rank1  per time: %f", runtime_rank1/(iter*K)))
  return(list(qg = qg, ELBO = ELBOs, KL = KLs, ll = lls))
}

## ================== helper functions ==================================

#' @export ebpmf_rank1_exponential_helper
#'
ebpmf_rank1_exponential_helper <- function(X, init = NULL, m = 2,
                                           scale_l = "estimate", scale_f = "estimate", gl_init = NULL, gf_init = NULL, fix_gl = F, fix_gf = F,
                                           maxiter = 1,verbose = F){
  X_rowsum = rowSums(X)
  X_colsum = colSums(X)
  p = length(X_colsum)
  n = length(X_rowsum)

  if(verbose){
    print(sprintf("%2s   %10s  %10s   %10s   %10s   %10s", "iter", "elbo", "kl_l", "kl_f", "sum_El", "sum_Ef"))
  }

  ## initialization.  It doesn't matter when doing rank-1 case. Maybe useful for rank-k case  when we assess convergence.
  if(is.null(init)){
    nnmf_res = NNLM::nnmf(A = X, k = 1, loss = "mkl", method = "lee", max.iter = 1)
    ql =  list(mean = nnmf_res$W[,1])
  }else{ql = init}

  for(i in 1:maxiter){
    ## update q(f), g(f)
    sum_El = sum(ql$mean)
    tmp_f = ebpm::ebpm_exponential_mixture(x = X_colsum, s = replicate(p,sum_El), m = m,scale = scale_f, g_init = gf_init, fix_g = fix_gl)
    qf = tmp_f$posterior
    gf = tmp_f$fitted_g
    ll_f = tmp_f$log_likelihood
    ### compute KL(q(f) || g(f)) = -Eq log(g(f)/q(f))
    kl_f = compute_kl(X_colsum, replicate(p,sum_El), tmp_f)
    ## update q(l), g(l)
    sum_Ef = sum(qf$mean)
    tmp_l = ebpm::ebpm_exponential_mixture(x = X_rowsum, s = replicate(n,sum_Ef), m = m, scale = scale_l, g_init = gl_init, fix_g = fix_gf)
    ql = tmp_l$posterior
    gl = tmp_l$fitted_g
    ll_l = tmp_l$log_likelihood
    ### compute KL(q(l) || g(l)) = -Eq log(g(l)/q(l))
    kl_l = compute_kl(X_rowsum, replicate(n,sum_Ef), tmp_l)
    qg = list(ql = ql, gl = gl, ll_l = ll_l, kl_l = kl_l, qf = qf, gf = gf, ll_f = ll_f, kl_f = kl_f)

    if(verbose){
      ## compute ELBO (for debugging)
      tmp = X * outer(ql$mean_log, qf$mean_log, "+")
      tmp[X == 0] = 0
      ll = - outer(ql$mean,  qf$mean, "*") + tmp
      ll = sum(ll)
      elbo = ll - kl_l - kl_f
      print(sprintf("%3d   %.10f  %.10f   %.10f   %.10f   %.10f", i, elbo, kl_l, kl_f, sum_El, sum_Ef))
    }
  }
  return(qg)
}


