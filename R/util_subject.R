
get_S_from_u <- function(u){ ## S[[d]] finds all i s.t. u[i] = d
  S = list()
  E = length(unique(u))
  for(e in 1:E){
    S[[e]] = which(u == e)
  }
  return(S)
}


compute_cs_by_subject <- function(L, S){
  E = length(S)
  K = ncol(L)
  out <- matrix(, nrow = E, ncol = K)
  for(e in 1:E){
    out[e,] = Matrix::colSums(L[S[[e]],])
  }
  return(out)
}

compute_sum_by_subject <- function(l, S){
  E = length(S)
  out <- replicate(E, NA)
  for(e in 1:E){
    out[e] = sum(l[S[[e]]])
  }
  return(out)
}
