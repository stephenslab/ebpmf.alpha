##


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


compute_elbo <- function(Ez, zeta, qg, KL){
  ## ELBO = -KL + log lik - Eq log q(Z)
  ELBO = -KL
  tmp = Ez*(operate_lf(qg$qls_mean_log, qg$qfs_mean_log, "+") - log(zeta))
  tmp[Ez == 0] = 0
  ELBO = ELBO + sum(tmp) - sum(operate_lf(qg$qls_mean, qg$qfs_mean, "*"))
}

## add l_ik + f_jk
## this is too slow
operate_lf <- function(L, F, operation){
  K = ncol(L)
  out = array(dim = c(nrow(L), nrow(F), K))
  for(k in 1:K){
    out[,,k] = outer(L[,k], F[,k], operation)
  }
  return(out)
}

## x is 3d array, and can contain -Inf
softmax3d <- function(x){
  score.exp <- exp(x)
  probs <-as.vector(score.exp)/as.vector(rowSums(score.exp,dims=2))
  probs[is.na(probs)] = 0 ##  since softmax((0,0,...,0)) = (NA, NA,...,NA), I  amnually set them to be 0. But it is not safe!!!
  dim(probs) <- dim(x)
  return(probs)
}
