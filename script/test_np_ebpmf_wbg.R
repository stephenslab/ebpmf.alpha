rm(list = ls())
library(ebpm)
library(Matrix)
library(ebpmf.alpha)
# test code for ebpmf
set.seed(123)
maxiter = 100
alpha = 5
verbose = TRUE
k = 10
n = 500; p = 100
a = 0.3; b = 1
l = matrix(rgamma(n*k, shape = a, rate = b), ncol = k)
f = matrix(rgamma(p*k, shape = a, rate = b), ncol = k)
lam = l %*% t(f)
X = matrix(rpois(n*p, lam), nrow = n)

#library(NNLM)
#system.time(
#  lf_init <- NNLM::nnmf(X, k = k, loss = "mkl", 
#												method = "scd", max.iter = maxiter)
#)
#L0 = lf_init$W 
#F0 = t(lf_init$H)
#qg0 = ebpmf.alpha::initialize_qg_bg_from_LF(L0 = L0, F0 = F0)
#

system.time(
  fit_ebpmf <- ebpmf.alpha::np_ebpmf_wbg(X = X, K = 2 *k, alpha = alpha,  
												pm_func = list(f = ebpm::ebpm_gamma_mixture, l = ebpm::ebpm_gamma_mixture),
												init = NULL,
												fix_option = list(l0 = FALSE, ql = FALSE, gl = FALSE,
																					f0 = FALSE, qf = FALSE, gf = FALSE),
												maxiter = maxiter, verbose = verbose) 
)

