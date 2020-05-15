rm(list = ls())
library(ebpm)
library(Matrix)
library(ebpmf.alpha)
# test code for ebpmf

set.seed(123)
maxiter = 20
verbose = TRUE
k = 10
n = 500; p = 100
a = 0.3; b = 1
l = matrix(rgamma(n*k, shape = a, rate = b), ncol = k)
f = matrix(rgamma(p*k, shape = a, rate = b), ncol = k)
lam = l %*% t(f)
X = matrix(rpois(n*p, lam), nrow = n)

library(NNLM)
system.time(
  lf_init <- NNLM::nnmf(X, k = k, loss = "mkl", 
												method = "scd", max.iter = maxiter)
)
L0 = lf_init$W 
F0 = t(lf_init$H) 
init = list(L = L0, F = F0)

system.time(
	fit_pmf <- ebpmf.alpha::pmf(X = X, K = k, 
															init = init, maxiter = maxiter, 
															verbose = verbose)
)

