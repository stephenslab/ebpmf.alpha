
context("test_ebpmf_wbg2")
library(ebpmf.alpha)
# source("../../R/pmf_subject.R")
# library("Matrix")
# source("../../R/util.R")
# source("../../R/util_subject.R")

library(fastTopics)
set.seed(123)

n = 100
p = 50
K = 3

## simulate data
L = matrix(100*runif(n*K), ncol = K)
F = matrix(runif(p*K), ncol = K)
m = matrix(runif(p*2), nrow = p)
u = c(replicate(n/2, 1), replicate(n/2, 2))
lam = L %*% t(F)
lam[u == 1,] = t(t(lam[u == 1,]) * m[,1])
lam[u == 2,] = t(t(lam[u == 2,]) * m[,2])
X = matrix(rpois(n = n*p, lambda = lam), nrow = n)

## initialization

fit = pmf_subject(X = X, u = u, K = K, verbose = TRUE)

test_that("pmf-subject works", {
  expect_false(any(diff(fit$ELBO) < 0))
})