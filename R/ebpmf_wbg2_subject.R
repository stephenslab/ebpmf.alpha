#' @title Empirical Bayes Poisson Matrix Factorization (Background Model with weights, using subject information)
#' @import ebpm
#' @import Matrix

#' @param X count matrix (dim(X) = c(n, p)).
#' @param u vector encoding subject (length(u) = n; elements need to be from {1,..., length(unique(u))}).
#' @param k number of topics
#' @param pm_func functions for solving the \code{ebpm} subproblem for \code{L} and \code{F};
#' It is a list \code{list(l, f, m)};
#' @param pm_control control parameters for pm_func function
#' @param init Either \code{NULL} or \code{list(qg, m, w)}
#' @param fix_option list(m, gl, ql, gf, qf) where each is either TRUE or FALSE
#' @param maxiter maximum number of iterations
#' @param tol stopping tolerance for ELBO
#' @param seed used when init is NULL
#'
#' @return A list containing elements:
#'     \describe{
#'       \item{\code{m}}{list of "mean" and "mean_log", "g", "kl"; 
#' 					"mean" and "mean_log" are each a p by D background frequency matrix}
#'       \item{\code{qg}}{list(ql, gl,qf, gf)}
#'       \item{\code{ELBO}}{ELBO objective for this VEB algorithm}
#'      }
#' @examples
#' To add
#' @export  ebpmf_wbg2_subject

ebpmf_wbg2_subject <- function(X, u, K,
										 pm_func = list(f = ebpm::ebpm_gamma_mixture,
																		l = ebpm::ebpm_gamma,
																		m = ebpm::ebpm_gamma),
										 init = NULL, pm_control = NULL,
										 fix_option = list(m = FALSE,
										                   gl = FALSE, ql = FALSE,
										                   gf = FALSE, qf = FALSE),
										 maxiter = 100, tol = 1e-8,
										 verbose = FALSE, seed = 123){
	## TODO: input check, require X_rs, X_cs to be nonzero
	S = get_S_from_u(u)
	E = length(unique(u))
  ## transform to sparse matrix
  X <- as(X, "sparseMatrix")
	X_rs = Matrix::rowSums(X)
	X_cs = Matrix::colSums(X)
	X_cs_subject = compute_cs_by_subject(X, S) ## E by p matrix
  d = summary(X)
	const = sum(apply.nonzeros(X = X, f = function(x) lgamma(x + 1)))
  ## initialization
	init_tmp <- init_ebpmf_wbg2_subject(X = X, K = K, u = u, init = init, d = d, seed = seed)
  qg <- init_tmp$qg
  b <- init_tmp$b
  a <- init_tmp$a
	m <- init_tmp$m
	w <- init_tmp$w
	w_log = log(w)
  rm(init_tmp)
  ## update iteratively
  ELBOs <- c()
	KLs <- c()
  for(i in 1:maxiter){
		b_k_max = replicate(length(d$x),0) ## max b_k
		for(k in 1:K){
			## store B_k
 			b_k = w_log[k] + qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k] - a
    	## compute q(Z)
      Ez <- compute_EZ(d = d,b = b,b_k = b_k)
      ## update (qL, gL, qF, gF)
      rank1_qg <- rank1_wbg2_subject(d = d, X_rs = Ez$rs, X_cs = Ez$cs,
      											S = S, u = u, 
														m = m, w_log_k = w_log[k],
														pm_func = pm_func, pm_control = pm_control,
														ql = list(mean = qg$qls_mean[,k],
																			mean_log = qg$qls_mean_log[,k]),
														qf = list(mean = qg$qfs_mean[,k],
																			mean_log = qg$qfs_mean_log[,k]),
													 	gl = qg$gls[[k]],
													 	gf = qg$gfs[[k]],
														kl_l = qg$kl_l[k],
														kl_f = qg$kl_f[k],
														fix_option = fix_option)
			rm(Ez)
      qg = update_qg(tmp = rank1_qg, qg = qg, k = k)
      rm(rank1_qg)
			## update b
      b_k0 = b_k
      b_k = w_log[k] + qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k] - a
			b = log( exp(b) - exp(b_k0) + exp(b_k)  )
    }
		## update w[1:K], and b accordingly
		b_new = b
		L_cs_subject = compute_cs_by_subject(qg$qls_mean, S)
		mlf = colSums((m$mean %*% L_cs_subject) * qg$qfs_mean)
		for(k in 1:K){
			b_k = w_log[k] + qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k] - a
			w_log[k] = log( sum(d$x * exp(b_k - b)) ) - log( mlf[k] ) 
			b_k_new = w_log[k] + qg$qls_mean_log[d$i,k] + qg$qfs_mean_log[d$j, k] - a
			b_new =  log( exp(b_new) - exp(b_k) + exp(b_k_new)  )
			b_k_max = pmax(b_k, b_k_max)
		}
		b = b_new
	  ## update m
	  wlf_subject = L_cs_subject %*% (exp(w_log) * t(qg$qfs_mean)) ## E by p
    if(!fix_option$m){
      for(e in 1:E){
      	s = wlf_subject[e,]
      	m1_tmp = do.call(pm_func$m, c(list(x = X_cs_subject[e,], s = s, g_init = m$g[[e]]), pm_control))
      	m$mean[, e] <- m1_tmp$posterior$mean
      	m$mean_log[, e] <- m1_tmp$posterior$mean_log
      	m$g[[e]] <- m1_tmp$fitted_g
      	m$kl[[e]] <- compute_kl_ebpm(y = X_cs_subject[e,], s = s, 
      		posterior = m1_tmp$posterior, ll = m1_tmp$log_likelihood)
      }
    }
		## compute ELBO
		w = exp(w_log)
		KL = sum(qg$kl_l) + sum(qg$kl_f) + sum(unlist(m$kl))
		ELBO = compute_elbo_wbg2_subject(w = w, m = m, qg = qg,
														b = b, a = a, d = d, X_cs_subject = X_cs_subject, S = S, const = const)
		ELBOs <- c(ELBOs, ELBO)
		KLs <- c(KLs, KL)
		## update a & b
    #a0 = a
    #a = b_k_max + a0
    #b = (b + a0) - a
    b = b - b_k_max
    a = b_k_max + a
		## verbose
    if(verbose){
      print("iter         ELBO")
      print(sprintf("%d:    %f", i, ELBO))
    }
  }
  return(list(w = w, m = m, qg = qg, ELBO = ELBOs, KL = KLs))
}









