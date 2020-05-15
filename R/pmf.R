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
	B = init_tmp$B
	rm(init_tmp)
	log_liks = c()
	for(i in 1:maxiter){
		for(k in 1:K){
			## compute Ez
			B_k = L[d$i,k] * F[d$j,k]
			Ez <- compute_EZ(d = d, B = B, B_k = B_k)
			## rank-1 update
			L[,k] = Ez$rs/sum(F[,k])
			F[,k] = Ez$cs/sum(L[,k])
			rm(Ez)
			## update B
			B = B - B_k + L[d$i,k] * F[d$j,k]
		}
		## compute loglikelihood
		ll = - sum(colSums(L)*colSums(F)) + sum(d$x * log(B)) - const
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
  B = L[d$i,1] * F[d$j,1]
  for(k in 2:K){
    B <- B + L[d$i,k] * F[d$j,k]
  }
	return(list(L = L, F = F, B = B))
}
