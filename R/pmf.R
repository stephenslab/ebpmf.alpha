#' @title Poisson Matrix Factorization

#' @param X count matrix (dim(X) = c(n, p)).
#' @param K number of topics
#' @param init Either NULL or \code{list(L, F)} 
#' @param maxiter
#' @param verbose
#' @param seed only used when init is NULL

#' @return \code{list(L, F, log_liks)} 

#' @export pmf
pmf <- function(X, K, init, maxiter = 100, verbose = FALSE, seed = 123){
	X <- as(X, "sparseMatrix")
	d = summary(X)
	const = sum(apply.nonzeros(X = X, f = function(x) lgamma(x + 1)))
	init_tmp = init_pmf(X = X, K = K, init = init, d= d, seed = seed)
	L = init_tmp$L
	F = init_tmp$F
	b = init_tmp$b
	a = init_tmp$a
	rm(init_tmp)
	log_liks = c()
	for(i in 1:maxiter){
		b_k_max = replicate(length(d$x), 0) ## max b_k
		for(k in 1:K){
			## compute Ez
			b_k = log(L[d$i,k]) + log(F[d$j,k]) - a
			Ez <- compute_EZ(d = d, b = b, b_k = b_k)
			## rank-1 update
			L[,k] = Ez$rs/sum(F[,k])
			F[,k] = Ez$cs/sum(L[,k])
			rm(Ez)
			## update B
			b_k0 = b_k
			b_k = log(L[d$i,k]) + log(F[d$j,k]) - a
      b = log( exp(b) - exp(b_k0) + exp(b_k)  )
      b_k_max = pmax(b_k, b_k_max)
		}
		## compute loglikelihood
		ll = - sum(colSums(L)*colSums(F)) + sum(d$x * (b + a)) - const
		log_liks = c(log_liks, ll)
		## verbose
    if(verbose){
      print("iter         loglik")
      print(sprintf("%d:    %f", i, ll))
    }
	}
	return(list(L = L, F = F, log_liks = log_liks))
}

init_pmf <- function(X, K, init, d, seed = 123){
	n = nrow(X)
	p = ncol(X)
	if(is.null(init)){
		set.seed(seed)
		L <- matrix(runif(n*K), ncol = K)
		F <- matrix(runif(p*K), ncol = K)
	}else{
		L <- init$L
		F <- init$F
	}

	## compute `a`
  a = replicate(length(d$x), 0)
  for(k in 1:K){
    b_k_tmp <- log(L[d$i,k]) + log(F[d$j,k])
    a <- pmax(a, b_k_tmp)
  }
  ## compute b
  b = log(L[d$i,1]) + log(F[d$j,1]) - a
  for(k in 2:K){
    b_k = log(L[d$i,k]) + log(F[d$j,k]) - a
    b <- log( exp(b) + exp(b_k)  )
  }
	return(list(L = L, F = F, b = b, a = a))
}
