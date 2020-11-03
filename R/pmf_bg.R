#' @title Poisson Matrix Factorization with background

#' @param X count matrix (dim(X) = c(n, p)).
#' @param K number of topics
#' @param init Either NULL or \code{list(l0, f0, L, F)} 
#' @param maxiter
#' @param verbose
#' @param seed only used when init is NULL

#' @return \code{list(l0, f0, L, F, log_liks)} 

#' @export pmf_bg
pmf_bg <- function(X, K, init, fix_option,
								maxiter = 100, verbose = FALSE, seed = 123){
	X <- as(X, "sparseMatrix")
	X_rs = Matrix::rowSums(X)
  X_cs = Matrix::colSums(X)
	d = summary(X)
	const = sum(apply.nonzeros(X = X, f = function(x) lgamma(x + 1)))
	init_tmp = init_pmf_bg(X = X, K = K, init = init, d= d, seed = seed)
	l0 = init_tmp$l0
	f0 = init_tmp$f0
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
			if(!fix_option$L){
				L[,k] = Ez$rs/(l0 * sum(f0 * F[,k]))
			}
			if(!fix_option$F){
				F[,k] = Ez$cs/(f0 * sum(l0 * L[,k]))
			}
			rm(Ez)
			## update B
			b_k0 = b_k
			b_k = log(L[d$i,k]) + log(F[d$j,k]) - a
      b = log( exp(b) - exp(b_k0) + exp(b_k)  )
      b_k_max = pmax(b_k, b_k_max)
			if(!fix_option$l0){
        denom <- colSums(t(L) * colSums(f0 * F))
        l0 <- X_rs/denom
      }
			if(!fix_option$f0){
				denom <- colSums(t(F) * colSums(l0 * L))
      	f0 <- X_cs/denom
			}
		}
		## compute loglikelihood
		ll = - sum(colSums(l0 * L)*colSums(f0 * F)) + sum(d$x * (log(l0[d$i]) + log(f0[d$j]) +  b + a)) - const
		log_liks = c(log_liks, ll)
		## verbose
    if(verbose){
      print("iter         loglik")
      print(sprintf("%d:    %f", i, ll))
    }
	}
	return(list(l0 = l0, f0 = f0, L = L, F = F, log_liks = log_liks))
}

init_pmf_bg <- function(X, K, init, d, seed = 123){
	n = nrow(X)
	p = ncol(X)
	if(is.null(init)){
		set.seed(seed)
		l0 = Matrix::rowMeans(X)
  	f0 = Matrix::colMeans(X)
		l0 = l0 * (sum(X)/sum(f0)/K) / sum(l0)
		L <- matrix(2 * runif(n*K), ncol = K)
		F <- matrix(2 * runif(p*K), ncol = K)
	}else{
		l0 <- init$l0
		f0 <- init$f0
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
	return(list(l0 = l0, f0 = f0, L = L, F = F, b = b, a = a))
}
