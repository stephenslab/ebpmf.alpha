
context("test_ebpmf_wbg3")
library(ebpmf.alpha)
library(fastTopics)
set.seed(123)

n = 100
p = 50
K = 3

## simulate data
L = matrix(100*runif(n*K), ncol = K)
F = matrix(runif(p*K), ncol = K)
X = matrix(rpois(n = n*p, lambda = L %*% t(F)), nrow = n)
## initialization
fit = fastTopics::fit_poisson_nmf(X = X, k = K)
init <- ebpmf.alpha::initialize_wbg3_from_LF(L = fit$L, F = fit$F)
fit = ebpmf.alpha::ebpmf_wbg3(X = X, K = K, init = init, verbose = TRUE)

test_that("ebpmf-wbg3 works", {
  expect_false(any(diff(fit$ELBO) < 0))
})
