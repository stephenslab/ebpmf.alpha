## faster than ebpmf

ebpmf2 <- function(X, K,pm_func = ebpm::ebpm_point_gamma,
                  init = list(qg = NULL, init_method = "scd", init_iter = 20), pm_control = NULL,
                  fix_g = list(l = FALSE, f = FALSE), maxiter = 100,
                  tol = 1e-8, verbose = FALSE){
  ## get nonzero index
  mask_nz <- (X != 0)
  X_nz = X[mask_nz]
  ## initialization
  init_tmp <- init_ebpmf(X = X, K = K, init = init, mask_nz = mask_nz)
  qg <- init_tmp$qg
  D <- init_tmp$D
  rm(init_tmp)

  ## update iteratively
  ELBOs <- c()
  for(i in 1:maxiter){
    KL <- 0
    for(k in 1:K){
      ## compute Ez (list(rs, cs, B))
      Ez <- compute_EZ(X_nz = X_nz, D = D,
                       ql_log = qg$qls_mean_log[,k], qf_log = qg$qfs_mean_log[,k],
                       mask_nz = mask_nz)
      ## rank1 update (list(D, kl_l, kl_f, qg))
      init_r1 = list(sf = sum(qg$qls_mean[,k]),
                     gl = qg$gls[[k]], gf = qg$gfs[[k]])
      rank1_tmp <- rank1(D_minus = D - Ez$B,
                         X_rs = Ez$rs, X_cs = Ez$cs,
                         pm_func = pm_func, pm_control = pm_control,
                         init = init_r1, fix_g = fix_g, mask_nz = mask_nz)
      rm(Ez)
      D = rank1_tmp$D
      KL = KL + rank1_tmp$kl_l + rank1_tmp$kl_f
      qg = update_qg(rank1_tmp$qg, qg, k)
    }
    ## compute ELBO
    ELBO = compute_ll(X, qg) - KL
    ELBOs <- c(ELBOs, ELBO)
    ## verbose
    if(verbose){
      print("iter         ELBO")
      print(sprintf("%d:    %f", i, ELBO))
    }
    ## check convergence
    diff = ifelse(i > 2, ELBOs[i] - ELBOs[i-1], Inf)
    if(diff < tol){
      if(verbose){print(sprintf("reaches tol %f in %d iterations", tol, i))}
      break
    }
  }
  return(list(qg = qg, ELBO = ELBOs))
}

## output: qg, D
init_ebpmf <- function(X,K, init, mask_nz){
  qg = init$qg
  if(is.null(qg)){
    qg = initialize_qg(X, K, init_method =  init$init_method, init_iter = init$init_iter)
  }
  D = (exp(qg$qls_mean_log) %*% t(exp(qg$qfs_mean_log)))[mask_nz]
  return(list(qg = qg, D = D))
}

compute_EZ <- function(X_nz, D, ql_log, qf_log, mask_nz){
  ## todo: speed up by only calculating nonzero ones
  B = exp(outer(ql_log, qf_log, "+"))[mask_nz]
  Ez_sub = X_nz * (B/D)
  Ez = array(0, dim = dim(mask_nz))
  Ez[mask_nz] = Ez_sub
  return(list(rs = rowSums(Ez), cs = colSums(Ez), B = B))
}

rank1 <- function(D_minus, X_rs, X_cs, pm_func,pm_control, init, fix_g, mask_nz){
  p = length(X_cs)
  n = length(X_rs)
  ## initialization (in fact, any non-negative number well do)
  sf = init$sf
  if(is.null(sf)){
    sf = sum(NNLM::nnmf(A = X, k = 1, loss = "mkl", method = "lee", max.iter = 1, verbose = F)$W[,1])
  }
  ## fit for f, and compute kl_f
  fit_f = do.call(pm_func, c(list(x = X_cs, s = sf, g_init = init$gf), pm_control))
  kl_f = compute_kl(X_cs, sf, fit_f)
  ## fit for l, and compute kl_l
  sl = sum(fit_f$posterior$mean)
  fit_l = do.call(pm_func, c(list(x = X_rs, s = sl, g_init = init$gl), pm_control))
  kl_l = compute_kl(X_rs, sl, fit_l)
  ## update D
  D = D_minus + exp(outer(fit_l$posterior$mean_log,fit_f$posterior$mean_log,"+"))[mask_nz]
  ## list to return
  qg = list(ql = fit_l$posterior, gl = fit_l$fitted_g, qf = fit_f$posterior, gf = fit_f$fitted_g)
  out = list(qg = qg, kl_l = kl_l, kl_f = kl_f, D = D)
  return(out)
}



