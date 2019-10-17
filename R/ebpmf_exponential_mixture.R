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
#'
#'
#' @examples
#' To add

#'
#' @export  ebpmf_exponential_mixture

ebpmf_exponential_mixture <- function(X, K, m = 2, maxiter.out = 10, maxiter.int = 1, verbose = F, seed = 123,
                                      fix_g = F, fix_grid = F){
  set.seed(seed)

  ## init from NNLM::nnmf result
  start = proc.time()
  qg = initialize_qg(X, K)
  runtime_init = (proc.time() - start)[[3]]

  runtime_rank1 = 0
  runtime_ez = 0
  ELBOs = c()
  KLs = c()

  for(iter in 1:maxiter.out){
    ELBO = 0
    KL = 0
    ## get <Z_ijk>, zeta_ijk
    start = proc.time()
    tmp = get_Ez(X, qg, K)
    Ez = tmp$Ez
    zeta = tmp$zeta
    rm(tmp)
    runtime_ez = runtime_ez + (proc.time() - start)[[3]]
    #print(sprintf("iter: %d", iter))
    #browser()

    for(k in 1:K){
      ## update q, g
      start = proc.time()
      init_l  = list(mean = qg$qls_mean[,k])
      if(is.null(qg$gls[[k]])){
        tmp = ebpmf_rank1_exponential_helper(X = Ez[,,k],init = init_l,m = m, maxiter = maxiter.int)
      }else{
        if(fix_g){
          tmp = ebpmf_rank1_exponential_helper(X = Ez[,,k],init = init_l,m = m, maxiter = maxiter.int,
                                               fitted_gl = qg$gls[[k]], fitted_gf = qg$gfs[[k]])
        }else{
          if(fix_grid){
            tmp = ebpmf_rank1_exponential_helper(X = Ez[,,k],init = init_l,m = m, maxiter = maxiter.int,
                                                 grid_l = list(a = qg$gls[[k]]$a, b = qg$gls[[k]]$b), grid_f = list(a = qg$gfs[[k]]$a, b = qg$gfs[[k]]$b))
          }else{
            tmp = ebpmf_rank1_exponential_helper(X = Ez[,,k],init = init_l,m = m, maxiter = maxiter.int)
          }
        }
      }
      KL = KL + tmp$kl_l + tmp$kl_f
      runtime_rank1 = runtime_rank1 + (proc.time() - start)[[3]]
      qg = update_qg(tmp, qg, k)
    }
    ELBO = compute_elbo(Ez, zeta, qg, KL)
    KLs <- c(KLs, KL)
    ELBOs <- c(ELBOs, ELBO)
    if(verbose){
      print("iter         ELBO")
      print(sprintf("%d:    %f", iter, ELBO))
    }
  }
  # print("summary of  runtime:")
  # print(sprintf("init           : %f", runtime_init))
  # print(sprintf("Ez     per time: %f", runtime_ez/(iter*K)))
  # print(sprintf("rank1  per time: %f", runtime_rank1/(iter*K)))
  return(list(qg = qg, ELBO = ELBOs, KL = KLs))
}

## ================== helper functions ==================================

#' @export ebpmf_rank1_exponential_helper
ebpmf_rank1_exponential_helper <- function(X, init = NULL, m = 2, grid_l = NULL, grid_f = NULL, fitted_gl = NULL, fitted_gf = NULL, maxiter = 1, verbose = F){
  X_rowsum = rowSums(X)
  X_colsum = colSums(X)
  p = length(X_colsum)
  n = length(X_rowsum)

  ## initialization.  It doesn't matter when doing rank-1 case. Maybe useful for rank-k case  when we assess convergence.
  if(is.null(init)){
    nnmf_res = NNLM::nnmf(A = X, k = 1, loss = "mkl", method = "lee", max.iter = 1)
    ql =  list(mean = nnmf_res$W[,1])
  }else{ql = init}

  for(i in 1:maxiter){
    ## update q(f), g(f)
    sum_El = sum(ql$mean)
    tmp_f = ebpm::ebpm_exponential_mixture(x = X_colsum, s = replicate(p,sum_El), m = m,grid = grid_f, fitted_g = fitted_gf)
    qf = tmp_f$posterior
    gf = tmp_f$fitted_g
    ll_f = tmp_f$log_likelihood
    ### compute KL(q(f) || g(f)) = -Eq log(g(f)/q(f))
    kl_f = compute_kl(X_colsum, replicate(p,sum_El), tmp_f)
    ## update q(l), g(l)
    sum_Ef = sum(qf$mean)
    tmp_l = ebpm::ebpm_exponential_mixture(x = X_rowsum, s = replicate(n,sum_Ef), m = m, grid = grid_l, fitted_g = fitted_gl)
    ql = tmp_l$posterior
    gl = tmp_l$fitted_g
    ll_l = tmp_l$log_likelihood
    ### compute KL(q(l) || g(l)) = -Eq log(g(l)/q(l))
    kl_l = compute_kl(X_rowsum, replicate(p,sum_Ef), tmp_l)
    qg = list(ql = ql, gl = gl, ll_l = ll_l, kl_l = kl_l, qf = qf, gf = gf, ll_f = ll_f, kl_f = kl_f)

    if(verbose){
      ## compute ELBO (for debugging)
      tmp = X * outer(ql$mean_log, qf$mean_log, "+")
      tmp[X == 0] = 0
      ll = - outer(ql$mean,  qf$mean, "*") + tmp
      ll = sum(ll)
      elbo = ll - kl_l - kl_f
      print(sprintf("%3d   %.10f  %.10f   %.10f", i, elbo, kl_l, kl_f))
    }
  }
  return(qg)
}

initialize_qg <- function(X, K, seed = 123){
  nnmf_fit = NNLM::nnmf(A = X, k = K, loss = "mkl", max.iter = 20)
  qls_mean = nnmf_fit$W
  qfs_mean = t(nnmf_fit$H)
  qls_mean_log = log(nnmf_fit$W)
  qfs_mean_log = log(t(nnmf_fit$W))
  qg = list(qls_mean = qls_mean, qls_mean_log =qls_mean_log,
            qfs_mean = qfs_mean, qfs_mean_log =qfs_mean_log,
            gls = replicate(K, list(NULL)),gfs = replicate(K, list(NULL)))
  return(qg)
}

## compute the row & col sum of <Z_ijk> for a given k
get_Ez <- function(X, qg,K){
  n = nrow(X)
  p = ncol(X)
  zeta = array(dim = c(n, p, K))
  ## get <ln l_ik> + <ln f_jk>
  for(d in 1:K){
    zeta[,,d] = outer(qg$qls_mean_log[,d], qg$qfs_mean_log[,d], "+")
  }
  ## do softmax
  zeta = softmax3d(zeta)
  Ez = as.vector(zeta)*as.vector(X)
  dim(Ez) = dim(zeta)
  return(list(Ez = Ez, zeta = zeta))
}

update_qg <- function(tmp, qg, k){
  qg$qls_mean[,k] = tmp$ql$mean
  qg$qls_mean_log[,k] = tmp$ql$mean_log
  qg$qfs_mean[,k] = tmp$qf$mean
  qg$qfs_mean_log[,k] = tmp$qf$mean_log
  qg$gls[[k]] = tmp$gl
  qg$gfs[[k]] = tmp$gf
  return(qg)
}















