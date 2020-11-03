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


system.time(
	fit_pmf <- ebpmf.alpha::pmf_bg(X = X, K = k, 
																 fix_option = list(L = FALSE, F = FALSE, l0 = FALSE, f0 = FALSE),
																 init = NULL, maxiter = maxiter, verbose = verbose)
)

