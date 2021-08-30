#' @title Poisson Matrix Factorization (using subject information)

#' @param X count matrix (dim(X) = c(n, p)).
#' @param u vector encoding subject (length(u) = n; elements need to be from {1,..., length(unique(u))}).
#' @param K number of topics
#' @param init Either NULL or \code{list(L, F)} 
#' @param maxiter
#' @param verbose
#' @param seed only used when init is NULL

#' @return \code{list(L, F, log_liks)} 

#' @export pmf_subject
pmf_subject <- function(X, u, K, init = NULL, maxiter = 100, verbose = FALSE, seed = 123){
	S = get_S_from_u(u)
	E = length(unique(u))
	X <- as(X, "sparseMatrix")
	X_cs_subject = compute_cs_by_subject(X, S) ## E by p matrix
	d = summary(X)
	const = sum(apply.nonzeros(X = X, f = function(x) lgamma(x + 1)))
	init_tmp = init_pmf_subject(X = X, u = u, K = K, init = init, d= d, seed = seed)
	L = init_tmp$L
	F = init_tmp$F
	m = init_tmp$m
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
			s = as.vector((t(m) %*% F[,k])[u])
			L[,k] = Ez$rs/s

			l_sum_subject = compute_sum_by_subject(L[,k], S)
			
			#if(k == 3){browser()}
			s <- as.vector(m %*% l_sum_subject)
			F[,k] = Ez$cs/s
			rm(Ez)
			## update B
			b_k0 = b_k
			b_k = log(L[d$i,k]) + log(F[d$j,k]) - a
      b = log( exp(b) - exp(b_k0) + exp(b_k)  )
      b_k_max = pmax(b_k, b_k_max)
		}
		## update m
		L_cs_subject = compute_cs_by_subject(L, S)
		wlf_subject = L_cs_subject %*% t(F) ## E by p	
		m = t(X_cs_subject / wlf_subject)
		
		## compute loglikelihood
		L_cs_subject = compute_cs_by_subject(L, S) 
		tmp = t(log(m)) * X_cs_subject
		tmp[is.nan(tmp)] <- 0
		ll = - sum((m %*% L_cs_subject) * F) + sum(d$x * (b + a)) + sum(tmp) - const
		log_liks = c(log_liks, ll)
		## verbose
    if(verbose){
      print("iter         loglik")
      print(sprintf("%d:    %f", i, ll))
    }
	}
	return(list(L = L, F = F, m = m, log_liks = log_liks))
}

init_pmf_subject <- function(X, K, u, init, d, seed = 123){
	n = nrow(X)
	p = ncol(X)
	E = length(unique(u))
	if(is.null(init)){
		set.seed(seed)
		L <- matrix(runif(n*K), ncol = K)
		F <- matrix(runif(p*K), ncol = K)
		m <- matrix(runif(p*E), ncol = E)
	}else{
		L <- init$L
		F <- init$F
		m <- init$m
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
	return(list(L = L, F = F, m = m, b = b, a = a))
}
