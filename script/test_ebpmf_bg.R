rm(list = ls())
library(ebpm)
library(Matrix)
library(ebpmf.alpha)
# test code for ebpmf
set.seed(123)
maxiter = 100
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
# L0 = lf_init$W + 1e-10
# F0 = t(lf_init$H) + 1e-10
# qg0 = initialize_qg_from_LF(L0 = L0, F0 = F0)

system.time(
  fit_ebpmf <- ebpmf.alpha::ebpmf_bg(X = X, K = k, 
												list(f = ebpm::ebpm_gamma_mixture, l = ebpmf.alpha::mle_pm),
												maxiter = maxiter, verbose = verbose, 
												fix_g = list(l = FALSE, f = FALSE))
)

