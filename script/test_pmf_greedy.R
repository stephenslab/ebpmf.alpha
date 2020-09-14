rm(list = ls())
library(ebpm)
library(Matrix)
library(ebpmf.alpha)
# test code for ebpmf

set.seed(123)
maxiter = 150
rate = 0.9
verbose = TRUE
k = 10
n = 500; p = 100
a = 0.3; b = 1
l = matrix(rgamma(n*k, shape = a, rate = b), ncol = k)
f = matrix(rgamma(p*k, shape = a, rate = b), ncol = k)
lam = l %*% t(f)
X = matrix(rpois(n*p, lam), nrow = n)


system.time(
	fit_pmf <- ebpmf.alpha::pmf_greedy(X = X, K = 3 * k, 
															init = NULL, maxiter = maxiter, rate = rate,
															verbose = verbose)
)

