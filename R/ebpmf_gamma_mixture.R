#' @title Empirical Bayes Poisson Matrix Factorization with Mixture of Gamma as Prior Family
#' @description Uses Empirical Bayes to fit the model \deqn{X_{ij}  ~ Poi(\sum_k L_{ik} F_{jk})} with \deqn{L_{.k} ~ g_k()}
#' with g_k being either Mixture of gamma, or Point Gamma
#' @import mixsqp
#' @import ebpm

#' @details The model is fit in 2 stages: i) estimate \eqn{g} by maximum likelihood (over pi_k)
#' ii) Compute posterior distributions for \eqn{\lambda_j} given \eqn{x_j,\hat{g}}.
#' @param X count matrix (dim(X) = c(n, p)).
#' @param K number of topics
#' @param m multiplicative parameter for selecting grid in "ebpm::ebpm_gamma_mixture"
#' @param maxiter.out  maximum iterations in the outer  loop
#' @param init_method used in \code{NNLM::nnmf}. Either `scd` or `lee`
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


#' @export  ebpmf_gamma_mixture

ebpmf_gamma_mixture <- function(X, K, qg = NULL, maxiter.out = 10, fix_g = F, fix_grid = F, verbose = F, m = 2, init_method = "scd",
                                theta_l = "one", theta_f = "one",
                                threshold =  NULL,uniform_mixture = F,
                                seed = 123, Lam_true = NULL){
  set.seed(seed)
  ## init from NNLM::nnmf result
  start = proc.time()
  if(is.null(qg)){
    qg = initialize_qg(X, K, init_method)
  }
  runtime_init = (proc.time() - start)[[3]]

  runtime_rank1 = 0
  runtime_ez = 0
  ELBOs = c()
  KLs = c()
  lls = c()
  rmses = c()

  start = proc.time()
  tmp = get_Ez(X, qg, K, threshold)
  Ez = tmp$Ez
  zeta = tmp$zeta
  rm(tmp)
  runtime_ez = runtime_ez + (proc.time() - start)[[3]]

  for(iter in 1:maxiter.out){
    #if(iter == 73){browser()}
    KL = 0
    for(k in 1:K){
      ## update q, g
      start = proc.time()
      init_l  = list(mean = qg$qls_mean[,k])
      if(is.null(qg$gls[[k]])){
        tmp = ebpmf_rank1_gamma_helper(X = Ez[,,k],init = init_l,m = m, uniform_mixture = uniform_mixture, theta_l = theta_l, theta_f = theta_f)
      }else{
        if(fix_g){
          tmp = ebpmf_rank1_gamma_helper(X = Ez[,,k],init = init_l,m = m,
                                         gl_init = qg$gls[[k]], gf_init = qg$gfs[[k]], fix_gl = T, fix_gf = T,
                                         uniform_mixture  = uniform_mixture, theta_l = theta_l, theta_f = theta_f)
        }else{
          if(fix_grid){
            tmp = ebpmf_rank1_gamma_helper(X = Ez[,,k],init = init_l,m = m,
                                          scale_l = list(a = qg$gls[[k]]$a, b = qg$gls[[k]]$b),
                                           scale_f = list(a = qg$gfs[[k]]$a, b = qg$gfs[[k]]$b),
                                          uniform_mixture  =  uniform_mixture, theta_l = theta_l, theta_f = theta_f)
          }else{
            tmp = ebpmf_rank1_gamma_helper(X = Ez[,,k],init = init_l,m = m,
                                           uniform_mixture = uniform_mixture, theta_l = theta_l, theta_f = theta_f)
          }
        }
      }
      KL = KL + tmp$kl_l + tmp$kl_f
      qg = update_qg(tmp, qg, k)
      ## update Z
      start = proc.time()
      tmp = get_Ez(X, qg, K, threshold)
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
  # print("summary of  runtime:")
  # print(sprintf("init           : %f", runtime_init))
  # print(sprintf("Ez     per time: %f", runtime_ez/(iter*K)))
  # print(sprintf("rank1  per time: %f", runtime_rank1/(iter*K)))
  return(list(qg = qg, ELBO = ELBOs, KL = KLs, ll = lls, RMSE = rmses))
}

## ================== helper functions ==================================

#' @export ebpmf_rank1_gamma_helper
#'
ebpmf_rank1_gamma_helper <- function(X, init = NULL, m = 2,
                                           scale_l = "estimate", scale_f = "estimate", theta_l = "one", theta_f = "one",
                                     gl_init = NULL, gf_init = NULL, fix_gl = F, fix_gf = F,
                                           uniform_mixture = F,low = NULL,
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
    nnmf_res = NNLM::nnmf(A = X, k = 1, loss = "mkl", method = "lee", max.iter = 1, verbose = F)
    ql =  list(mean = nnmf_res$W[,1])
  }else{ql = init}

  for(i in 1:maxiter){
    ## update q(f), g(f)
    sum_El = sum(ql$mean)

    if(uniform_mixture){
      ## f
      gf_init = get_uniform_mixture(x = X_colsum,s = replicate(p,sum_El),grid_res = gf_init, m = m, low = low)
      tmp_f = ebpm::ebpm_gamma_mixture_single_scale(x = X_colsum, s = replicate(p,sum_El), m = m, g_init = gf_init, fix_g = T, theta = theta_f)
    }else{
      tmp_f = ebpm::ebpm_gamma_mixture_single_scale(x = X_colsum, s = replicate(p,sum_El), m = m,scale = scale_f, g_init = gf_init, fix_g = fix_gl,
                                                    low = low, theta = theta_f)
    }
    qf = tmp_f$posterior
    gf = tmp_f$fitted_g
    ll_f = tmp_f$log_likelihood
    ### compute KL(q(f) || g(f)) = -Eq log(g(f)/q(f))
    kl_f = compute_kl(X_colsum, replicate(p,sum_El), tmp_f)
    ## update q(l), g(l)
    sum_Ef = sum(qf$mean)

    if(uniform_mixture){
      ## f
      gl_init = get_uniform_mixture(x = X_rowsum,s = replicate(n,sum_Ef),grid_res = gl_init, m, low =  low)
      tmp_l = ebpm::ebpm_gamma_mixture_single_scale(x = X_rowsum, s = replicate(n,sum_Ef), m = m, g_init = gl_init, fix_g = T, theta = theta_l)
    }else{
      tmp_l = ebpm::ebpm_gamma_mixture_single_scale(x = X_rowsum, s = replicate(n,sum_Ef), m = m, scale = scale_l, g_init = gl_init, fix_g = fix_gf,
                                       low = low, theta = theta_l)
    }
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


