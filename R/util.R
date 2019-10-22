## ================================================================================================================
## Functions updating in ebpmf
## ================================================================================================================

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


initialize_qg <- function(X, K, seed = 123){
  nnmf_fit = NNLM::nnmf(A = X, k = K, loss = "mkl", max.iter = 20)
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
softmax3d <- function(x){
  score.exp <- exp(x)
  probs <-as.vector(score.exp)/as.vector(rowSums(score.exp,dims=2))
  probs[is.na(probs)] = 0 ##  since softmax((0,0,...,0)) = (NA, NA,...,NA), I manually set them to be 0. But it is not safe!!!
  dim(probs) <- dim(x)
  return(probs)
}


## ================================================================================================================
## Functions for computing ELBO
## ================================================================================================================

## compute KL(q || p) = Eq(log(q(lambda)/p(lambda))) from ebpm fitted result
compute_kl <- function(x, s, fit_pm){
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







