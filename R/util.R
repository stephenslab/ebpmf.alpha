## ================================================================================================================
## Functions updating in ebpmf
## ================================================================================================================

## compute the row & col sum of <Z_ijk> for a given k
get_Ez <- function(X, qg,K, threshold = NULL){
  n = nrow(X)
  p = ncol(X)
  zeta = array(dim = c(n, p, K))
  ## get <ln l_ik> + <ln f_jk>
  for(d in 1:K){
    zeta[,,d] = outer(qg$qls_mean_log[,d], qg$qfs_mean_log[,d], "+") ##  this can be -Inf
  }
  ## do softmax
  zeta = softmax3d(zeta, threshold)
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

# ## TESTING
# ## only update
# update_qg_fixl <- function(tmp, qg, k, iter){
#   qg$qfs_mean[,k] = tmp$qf$mean
#   qg$qfs_mean_log[,k] = tmp$qf$mean_log
#   if(iter == 1){
#     qg$qls_mean[,k] = tmp$ql$mean
#     qg$qls_mean_log[,k] = tmp$ql$mean_log
#     qg$gls[[k]] = tmp$gl
#   }
#   #qg$gls[[k]] = tmp$gl
#   qg$gfs[[k]] = tmp$gf
#   return(qg)
# }

#' @export initialize_qg_from_LF
initialize_qg_from_LF <- function(L0,F0){
  K = ncol(L0)
  qls_mean = L0
  qfs_mean = F0
  qls_mean_log = log(L0)
  qfs_mean_log = log(F0)
  qg = list(qls_mean = qls_mean, qls_mean_log =qls_mean_log,
            qfs_mean = qfs_mean, qfs_mean_log =qfs_mean_log,
            gls = replicate(K, list(NULL)),gfs = replicate(K, list(NULL)))
  return(qg)
}

initialize_qg <- function(X, K, init_method = "scd", seed = 123){
  nnmf_fit = NNLM::nnmf(A = X, k = K, loss = "mkl", max.iter = 20, verbose = F, method = init_method)
  qls_mean = nnmf_fit$W
  qfs_mean = t(nnmf_fit$H)
  qls_mean_log = log(nnmf_fit$W)
  qfs_mean_log = log(t(nnmf_fit$H))
  qg = list(qls_mean = qls_mean, qls_mean_log =qls_mean_log,
            qfs_mean = qfs_mean, qfs_mean_log =qfs_mean_log,
            gls = replicate(K, list(NULL)),gfs = replicate(K, list(NULL)))
  return(qg)
}

## x is 3d array, and can contain -Inf
softmax3d <- function(x, threshold = NULL){
  score.exp <- exp(x)
  probs <-as.vector(score.exp)/as.vector(rowSums(score.exp,dims=2))
  if(!is.null(threshold)){
    dim(probs) <- dim(x)
    probs[probs < threshold] = 0
    probs <-as.vector(probs)/as.vector(rowSums(probs,dims=2))
  }
  probs[is.na(probs)] = 0 ##  since softmax((0,0,...,0)) = (NA, NA,...,NA), I manually set them to be 0. But it is not safe!!!
  dim(probs) <- dim(x)
  return(probs)
}

compute_rmse <- function(lam1, lam2){
  return(sqrt(mean((lam1 - lam2)^2)))
}

## ================================================================================================================
## Functions for random effect
## ================================================================================================================


initialize_qg_random_effect <- function(X, K, init_method = "scd", seed = 123){
  nnmf_fit = NNLM::nnmf(A = X, k = K+1, loss = "mkl", max.iter = 20, verbose = F, method = init_method)
  qls_mean = nnmf_fit$W[,1:K]
  qfs_mean = t(nnmf_fit$H)[,1:K]
  qmu_mean = outer(nnmf_fit$W[,K+1],nnmf_fit$H[K+1,], "*")
  qls_mean_log = log(qls_mean)
  qfs_mean_log = log(qfs_mean)
  qmu_mean_log = log(qmu_mean)
  qg = list(qls_mean = qls_mean, qls_mean_log = qls_mean_log,
            qfs_mean = qfs_mean, qfs_mean_log = qfs_mean_log,
            qmu_mean = qmu_mean, qmu_mean_log = qmu_mean_log,
            gls = replicate(K, list(NULL)),gfs = replicate(K, list(NULL)),
            gmu = list(NULL))
  return(qg)
}


## compute the row & col sum of <Z_ijk> for a given k
get_Ez_random_effect <- function(X, qg,K){
  n = nrow(X)
  p = ncol(X)
  zeta = array(dim = c(n, p, K+1))
  ## get <ln l_ik> + <ln f_jk>
  for(d in 1:K){
    zeta[,,d] = outer(qg$qls_mean_log[,d], qg$qfs_mean_log[,d], "+") ##  this can be -Inf
  }
  zeta[,,K+1] = qg$qmu_mean_log
  ## do softmax
  zeta = softmax3d(zeta)
  Ez = as.vector(zeta)*as.vector(X)
  dim(Ez) = dim(zeta)
  return(list(Ez = Ez, zeta = zeta))
}

update_qg_random_effect <- function(tmp, qg){
  qg$qmu_mean = tmp$qmu_mean
  qg$qmu_mean_log = tmp$qmu_mean_log
  qg$gmu = tmp$gmu
  return(qg)
}

compute_ll_random_effect <- function(X, qg){
  tmp = X * log(rowSums(exp(operate_lf(qg$qls_mean_log, qg$qfs_mean_log, "+")), dims = 2) + exp(qg$qmu_mean_log))
  tmp[X ==  0] =  0
  out = - sum(operate_lf(qg$qls_mean, qg$qfs_mean, "*")) - sum(qg$qmu_mean) + sum(tmp)
  return(out)
}



## ================================================================================================================
## Functions for computing ELBO
## ================================================================================================================

## compute KL(q || p) = Eq(log(q(lambda)/p(lambda))) from ebpm fitted result
compute_kl <- function(x, s, fit_pm){
  #browser()
  ## -KL = fit_pm$log_likelihood - sum(tmp) + sum(lgamma(x + 1))
  ## tmp =  - s * fit_pm$posterior$mean  +  x * (log(s) + fit_pm$posterior$mean_log)
  tmp = x * (log(s) + fit_pm$posterior$mean_log)
  tmp[x == 0] = 0
  tmp = tmp - s * fit_pm$posterior$mean
  out = fit_pm$log_likelihood - sum(tmp) + sum(lgamma(x + 1))
  return(-out)
}

compute_ll <- function(X, qg){
  tmp = X * log(rowSums(exp(operate_lf(qg$qls_mean_log, qg$qfs_mean_log, "+")), dims = 2))
  tmp[X ==  0] =  0
  out = - sum(operate_lf(qg$qls_mean, qg$qfs_mean, "*")) + sum(tmp)
  return(out)
}


## this is too slow
operate_lf <- function(L, F, operation){
  K = ncol(L)
  out = array(dim = c(nrow(L), nrow(F), K))
  for(k in 1:K){
    out[,,k] = outer(L[,k], F[,k], operation)
  }
  return(out)
}







