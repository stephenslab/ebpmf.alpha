#' @title Empirical Bayes Poisson Matrix Factorization with Prior Family Mixture of Exponential
#' @description Uses Empirical Bayes to fit the model \deqn{X_{ij}  ~ Poi(\sum_k L_{ik} F_{jk})} with \deqn{L_{.k} ~ g_k()}
#' with g_k being either Mixture of Exponential, or Point Gamma
#' @import mixsqp
#' @import ebpm
#' @import gtools


#' @details The model is fit in 2 stages: i) estimate \eqn{g} by maximum likelihood (over pi_k)
#' ii) Compute posterior distributions for \eqn{\lambda_j} given \eqn{x_j,\hat{g}}.
#' @param X count matrix (dim(X) = c(n, p)).
#' @param K number of topics
#' @param  maxiter.out  maximum iterations in the outer  loop
#' @param  maxiter.int  maximum iterations in the inner loop
#' @param verbose T if print ELBO
#' @param seed random seed
#'
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
#' @export
ebpmf_point_gamma <- function(X, K, maxiter.out = 10, maxiter.int = 1, verbose = F, seed = 123){
  set.seed(123)
  start = proc.time()
  qg = initialize_qg(X, K)
  tmp = get_Ez(X, qg, K)
  Ez = tmp$Ez
  zeta = tmp$zeta
  rm(tmp)
  runtime_init = (proc.time() - start)[[3]]
  runtime_rank1 = 0
  runtime_ez = 0

  ELBOs = c()
  KLs = c()
  for(iter in 1:maxiter.out){
    #if(iter == 2){browser()}
    ELBO = 0
    KL = 0
    # ## get <Z_ijk>
    # start = proc.time()
    # tmp = get_Ez(X, qg, K)
    # Ez = tmp$Ez
    # zeta = tmp$zeta
    # rm(tmp)
    # runtime_ez = runtime_ez + (proc.time() - start)[[3]]
    #browser()
    for(k in 1:K){
      ## update q, g
      start = proc.time()
      tmp = ebpmf_rank1_point_gamma_helper(rowSums(Ez[,,k]),colSums(Ez[,,k]),NULL,maxiter.int)
      KL = KL + tmp$kl_l + tmp$kl_f
      runtime_rank1 = runtime_rank1 + (proc.time() - start)[[3]]
      qg = update_qg(tmp, qg, k)
    }
    ## get <Z_ijk>
    start = proc.time()
    tmp = get_Ez(X, qg, K)
    runtime_ez = runtime_ez + (proc.time() - start)[[3]]
    Ez = tmp$Ez
    zeta = tmp$zeta
    rm(tmp)

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
## for each pair of l, f, give them 1/k of the row & col sum
initialize_qg <- function(X, K, seed = 123){
  n = nrow(X)
  p = ncol(X)
  set.seed(seed)
  X_rsum = rowSums(X)
  X_csum = colSums(X)
  prob_r = replicate(n, rdirichlet(1,replicate(K, 1/K)))[1,,] ## K by n
  prob_c = replicate(p, rdirichlet(1,replicate(K, 1/K)))[1,,] ## K  by p
  rsums = matrix(replicate(K*n,0), nrow = K)
  csums = matrix(replicate(K*p,0), nrow = K)
  for(i in  1:n){
    if(X_rsum[i] == 0){rsums[,i] = replicate(K, 0)}
    else{rsums[,i] = rmultinom(1, X_rsum[i],prob_r[,i])}
  }
  for(j in  1:p){
    if(X_csum[j] == 0){csums[,j] = replicate(K, 0)}
    else{csums[,j] = rmultinom(1, X_csum[j],prob_c[,j])}
  }
  qg = list(qls_mean = matrix(replicate(n*K, 0), ncol =  K), qls_mean_log = matrix(replicate(n*K, 0), ncol =  K), gls = replicate(K, list(NaN)),
            qfs_mean = matrix(replicate(p*K, 0), ncol =  K), qfs_mean_log = matrix(replicate(p*K, 0), ncol =  K), gfs = replicate(K, list(NaN))
  )
  for(k in 1:K){
    qg_ = ebpmf_rank1_exponential_helper(rsums[k,], csums[k, ], init = NULL, m = 2, maxiter = 1)
    #qg_ = ebpmf_rank1_point_gamma_helper(rsums[k,], csums[k, ], init = NULL, maxiter = 1) ## THIS IS PROBLEMATIC!!
    qg   = update_qg(qg_, qg, k)
  }
  return(qg)
}

## compute the row & col sum of <Z_ijk> for a given k
get_Ez <- function(X, qg,K){
  n = nrow(X)
  p = ncol(X)
  zeta = array(dim = c(n, p, K))
  ## get <ln l_ik> + <ln f_jk>
  for(d in 1:K){
    zeta[,,d] = outer(qg$qls_mean_log[,d], qg$qfs_mean_log[,d], "+") ##  this can be -Inf
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

ebpmf_rank1_point_gamma_helper <- function(X_rowsum,X_colsum, init = NULL,maxiter = 1){
  p = length(X_colsum)
  n = length(X_rowsum)
  if(is.null(init)){init = list(mean = runif(length(X_rowsum), 0, 1))}
  ql = init
  for(i in 1:maxiter){
    ## update q(f), g(f)
    sum_El = sum(ql$mean)
    tmp_f = ebpm::ebpm_point_gamma(x = X_colsum, s = replicate(p,sum_El))
    qf = tmp_f$posterior
    gf = tmp_f$fitted_g
    ll_f = tmp_f$log_likelihood
    ### compute KL(q(f) || g(f)) = -Eq log(g(f)/q(f))
    kl_f = compute_kl(X_colsum, replicate(p,sum_El), tmp_f)
    ## update q(l), g(l)
    sum_Ef = sum(qf$mean)
    tmp_l = ebpm_point_gamma(x = X_rowsum, s = replicate(n,sum_Ef))
    ql = tmp_l$posterior
    gl = tmp_l$fitted_g
    ll_l = tmp_l$log_likelihood
    ### compute KL(q(l) || g(l)) = -Eq log(g(l)/q(l))
    kl_l = compute_kl(X_rowsum, replicate(p,sum_Ef), tmp_l)
    qg = list(ql = ql, gl = gl, ll_l = ll_l, kl_l = kl_l, qf = qf, gf = gf, ll_f = ll_f, kl_f = kl_f)
  }
  return(qg)
}

#' @export ebpmf_rank1_point_gamma
ebpmf_rank1_point_gamma <- function(X, init = NULL,maxiter = 1, verbose  = F){
  X_rowsum = rowSums(X)
  X_colsum = colSums(X)
  p = length(X_colsum)
  n = length(X_rowsum)
  if(is.null(init)){init = list(mean = 100*runif(length(X_rowsum), 0, 1))}
  ql = init

  elbos = c()
  kls = c()

  if(verbose){
    print(sprintf("%2s   %15s  %15s   %15s  %15s  %15s","iter", "ELBO","KL_L", "KL_F" ,"sum_El", "sum_Ef"))
  }

  for(i in 1:maxiter){
    ## update q(f), g(f)
    sum_El = sum(ql$mean)
    tmp_f = ebpm::ebpm_point_gamma(x = X_colsum, s = replicate(p,sum_El))
    qf = tmp_f$posterior
    gf = tmp_f$fitted_g
    ll_f = tmp_f$log_likelihood
    ### compute KL(q(f) || g(f)) = -Eq log(g(f)/q(f))
    kl_f = compute_kl(X_colsum, replicate(p,sum_El), tmp_f)
    ## update q(l), g(l)
    sum_Ef = sum(qf$mean)
    tmp_l = ebpm::ebpm_point_gamma(x = X_rowsum, s = replicate(n,sum_Ef))
    ql = tmp_l$posterior
    gl = tmp_l$fitted_g
    ll_l = tmp_l$log_likelihood
    ### compute KL(q(l) || g(l)) = -Eq log(g(l)/q(l))
    kl_l = compute_kl(X_rowsum, replicate(p,sum_Ef), tmp_l)

    ## compute loglikelihood
    tmp = X * outer(ql$mean_log, qf$mean_log, "+")
    tmp[X == 0] = 0
    ll = - outer(ql$mean,  qf$mean, "*") + tmp
    ll = sum(ll)
    elbo = ll - kl_l - kl_f
    elbos  = c(elbos, elbo)
    kls = c(kls, kl_l + kl_f)
    qg = list(ql = ql, gl = gl, ll_l = ll_l, kl_l = kl_l, qf = qf, gf = gf, ll_f = ll_f, kl_f = kl_f, elbos = elbos, kls = kls)
    if(verbose){
      print(sprintf("%2d %10f  %.10f  %.10f %.10f %.10f", i, elbo, kl_l, kl_f, sum_El, sum_Ef))
    }
  }
  return(qg)
}












