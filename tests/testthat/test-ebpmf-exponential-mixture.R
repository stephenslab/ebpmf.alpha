rm(list  = ls())
context("test ebpmf_exponential_mixture.R")

library(NNLM)
library(gtools)
library(ebpm)

sim_mgamma <- function(dist){
  pi = dist$pi
  a = dist$a
  b = dist$b
  idx = which(rmultinom(1,1,pi) == 1)
  return(rgamma(1, shape = a[idx], rate =  b[idx]))
}

## simulate a poisson mean problem
## to do:
simulate_pm  <-  function(n, p, dl, df, K,scale_b = 10, seed = 1234){
  set.seed(seed)
  ## simulate L
  a = replicate(dl,1)
  b = 10*runif(dl)
  pi <- rdirichlet(1,rep(1/dl, dl))
  gl = list(pi = pi, a = a, b= b)
  L = matrix(replicate(n*K, sim_mgamma(gl)), ncol = K)
  ## simulate F
  a = replicate(df,1)
  b = 10*runif(df)
  pi <- rdirichlet(1,rep(1/df, df))
  gf = list(pi = pi, a = a, b= b)
  F = matrix(replicate(p*K, sim_mgamma(gf)), ncol = K)
  ## simulate X
  lam = L %*% t(F)
  X = matrix(rpois(n*p, lam), nrow = n)
  Y = matrix(rpois(n*p, lam), nrow = n)
  ## prepare output
  g = list(gl = gl, gf = gf)
  out = list(X = X, Y = Y, L = L, F = F, g = g)
  return(out)
}

### simulate data
n = 100
p = 200
K = 2
dl = 10
df = 10
scale_b = 5
sim = simulate_pm(n, p, dl, df, K, scale_b = scale_b)

## ebpmf
#browser()


out_ebpmf = ebpmf::ebpmf_exponential_mixture(sim$X, K, maxiter.out = 100)

## plot ELBOs & KLs
plot(out_ebpmf$ELBO, type = "l")
plot(out_ebpmf$KL, type = "l")


## nnmf
W0 = out_ebpmf$qg$qls_mean
H0 = t(out_ebpmf$qg$qfs_mean)
out_nmf = NNLM::nnmf(sim$X, K,init = list(W0 = W0, H0 = H0), loss = "mkl", method = "lee", max.iter = 100, rel.tol = -1, verbose  =  F)


## testing:
test_that("training loglikelihood beats oracle", {
  lam_pm = out_ebpmf$qg$qls_mean %*% t(out_ebpmf$qg$qfs_mean)
  ll_train_ebpmf = sum(dpois(sim$X, lambda = lam_pm, log = T))
  #ll_val_ebpmf = sum(dpois(sim$Y, lambda = lam_pm, log = T))

  ll_train_oracle = sum(dpois(sim$X, lambda = sim$L %*% t(sim$F), log = T))
  #ll_val_oracle = sum(dpois(sim$Y, lambda = sim$L %*% t(sim$F), log = T))
  expect_gt(ll_train_ebpmf, ll_train_oracle)
})


## testing:
test_that("validation loglikelihood beats nnmf", {
  lam_pm = out_ebpmf$qg$qls_mean %*% t(out_ebpmf$qg$qfs_mean)
  #ll_train_ebpmf = sum(dpois(sim$X, lambda = lam_pm, log = T))
  ll_val_ebpmf = sum(dpois(sim$Y, lambda = lam_pm, log = T))

  lam_nmf = out_nmf$W %*% out_nmf$H
  #ll_train_nnmf = sum(dpois(sim$X, lambda = lam_nmf, log = T))
  ll_val_nnmf = sum(dpois(sim$Y, lambda = lam_nmf, log = T))
  expect_gt(ll_val_ebpmf, ll_val_nnmf)
})


# plot(out_ebpmf$ELBO, ylab = "elbo")
#
# plot(out_ebpmf$ELBO - out_ebpmf$KL, ylab = "log-prob")

# test_that("elbo increases monotonically",{
#   elbos = out_ebpmf$ELBO
#   n = length(elbos)
#   expect_false(any(elbos[1:(n-1)] > elbos[2:n]))
# })


## experiment to see RMSE on Lambda
Lam_true = sim$L %*% t(sim$F)
try_experiment_rmse <- function(iter, Lam_true){
  test = ebpmf::ebpmf_exponential_mixture(sim$X, K, maxiter.out = iter)
  Lam = test$qg$qls_mean %*% t(test$qg$qfs_mean)
  return(sum((Lam - Lam_true)^2))
}

iters = seq(10,50,10)
rmses <- c()
for(iter in iters){
  rmse = try_experiment_rmse(iter, Lam_true)
  rmses = c(rmses, rmse)
}

## testing:
test_that("rmse decreases monotonically", {
  n = length(rmses)
  expect_false(any(rmses[1:(n-1)] < rmses[2:n]))
})

