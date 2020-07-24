#' @title Poisson Matrix Factorization (greedy)

#' @param X count matrix (dim(X) = c(n, p)).
#' @param K number of topics
#' @param init Either NULL or \code{list(L, F)} 
#' @param maxiter for each k
#' @param verbose
#' @param seed only used when init is NULL

#' @return \code{list(L, F, log_liks)} 

#' @export pmf_greedy
pmf_greedy <- function(X, K, init, 
											 maxiter = 100, rate = 0.9, 
											 verbose = FALSE, seed = 123){
	progress = matrix(, nrow = maxiter + 1, ncol = K - 1)
	X <- as(X, "sparseMatrix")
	d = summary(X)
	const = sum(apply.nonzeros(X = X, f = function(x) lgamma(x + 1)))
	init_tmp = init_pmf_greedy(X = X, K = K, init = init, d = d, seed = seed)	
	L = init_tmp$L
	F = init_tmp$F
	Lam = init_tmp$Lam
	zeta = NULL 
	loglik <- - sum(colSums(L)*colSums(F)) + sum(d$x * log(Lam)) - const
	if(verbose){
      print("k         loglik")
      print(sprintf("%d:    %f", 1,  loglik))
	}
	logliks <- c(loglik)
	for(k in 2:K){
		r1_tmp = rank1_pmf_greedy(X = X, Lam = Lam, d = d, k = k, zeta = zeta, 
															maxiter = maxiter, rate = rate)
		if(test_terminate(l = r1_tmp$L, f = r1_tmp$F)){
      print(sprintf("greedy algo terminates after adding %d lf", k))
			break
    }
		zeta = r1_tmp$zeta
		L <- cbind(L, r1_tmp$L)
		F <- cbind(F, r1_tmp$F)
		Lam = r1_tmp$Lam
		loglik <- - sum(colSums(L)*colSums(F)) + sum(d$x * log(Lam)) - const
		if(verbose){
      print(sprintf("%d:    %f", k, loglik))
    }
		logliks <- c(logliks, loglik)
		progress[, k - 1] = r1_tmp$obj
	}
	return(list(L = L, F = F, loglik = logliks, progress = progress))
}

rank1_pmf_greedy <- function(X, Lam, d, k, zeta = NULL, 
														 maxiter = 100, eps = 1e-4, rate = 0.9){
	if(is.null(zeta)){zeta = replicate(length(d$i), 1/k)}
	else{zeta = rate * zeta}
	Ez = sparseMatrix(i = d$i, j = d$j, x = d$x * zeta)
	l = Matrix::rowSums(Ez)
	f = Matrix::colSums(Ez)/sum(l)	
	objs <- c()
	obj = - sum(l)*sum(f) + sum(d$x * log(Lam + l[d$i] * f[d$j]))
	objs <- c(objs, obj)

	for(iter in 1:maxiter){
		zeta_prev = zeta
		zeta = 1 - Lam/(Lam + l[d$i] * f[d$j])
		rel_diff = max(abs(zeta - zeta_prev)/(zeta + zeta_prev))
		Ez = sparseMatrix(i = d$i, j = d$j, x = d$x * zeta)
		l = Matrix::rowSums(Ez)/sum(f)
		f = Matrix::colSums(Ez)/sum(l)

		obj = - sum(l)*sum(f) + sum(d$x * log(Lam + l[d$i] * f[d$j]))
		objs <- c(objs, obj)
	}
	zeta = 1 - Lam/(Lam + l[d$i] * f[d$j])
	Lam = Lam + l[d$i] * f[d$j]
	return(list(L = l, F = f, Lam = Lam, zeta = zeta, obj = objs))
}


test_terminate <- function(l, f, eps = 1e-4){
	(max(l) < eps) && (max(f) < eps)
}

init_pmf_greedy <- function(X, K, init, d, seed = 123){
	n = nrow(X)
	p = ncol(X)
	if(is.null(init)){
		L = Matrix::rowSums(X)
		F = Matrix::colSums(X)/sum(L)
		L = matrix(L, ncol = 1)
		F = matrix(F, ncol = 1)
	}else{
		L <- init$L
		F <- init$F
	}
	Lam = L[d$i,1] * F[d$j,1]
	return(list(L = L, F = F, Lam = Lam))
}
