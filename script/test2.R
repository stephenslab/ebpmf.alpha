sim_mgamma <- function(dist){
  pi = dist$pi
  a = dist$a
  b = dist$b
  idx = which(rmultinom(1,1,pi) == 1)
  return(rgamma(1, shape = a[idx], rate =  b[idx]))
}

## simulate a poisson mean problem
## to do:
simulate_pm  <-  function(n, p, dl, df, K,scale_b = 10, seed = 123){
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

n = 50
p = 100
K = 2
dl = 10
df = 10
scale_b = 5
sim = simulate_pm(n, p, dl, df, K, scale_b = scale_b, seed =12)
print(dim(sim$X))

out_ebpmf_point = ebpmf::ebpmf_point_gamma(sim$X, K, maxiter.out = 100)
plot(out_ebpmf_point$ELBO, type = "l", xlab = "niter", ylab = "ELBO")

m = 2
out_ebpmf_exp = ebpmf::ebpmf_exponential_mixture(sim$X, K, m = m, maxiter.out = 100)
plot(out_ebpmf_exp$ELBO, type = "l", xlab = "niter", ylab = "ELBO")









