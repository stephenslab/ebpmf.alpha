
context("test_ebpmf_wbg2")
library(ebpmf.alpha)
# source("../../R/ebpmf_wbg2_subject.R")
# source("../../R/ebpmf_wbg2_subject_util.R")
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
X = matrix(rpois(n = n*p, lambda = L %*% t(F)), nrow = n)
u = c(replicate(n/2, 1), replicate(n/2, 2))
## initialization

fit = ebpmf_wbg2_subject(X = X, u = u, K = K, verbose = TRUE)

test_that("ebpmf-wbg2-subject works", {
  expect_false(any(diff(fit$ELBO) < 0))
})