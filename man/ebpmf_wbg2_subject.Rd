% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ebpmf_wbg2_subject.R
\name{ebpmf_wbg2_subject}
\alias{ebpmf_wbg2_subject}
\title{Empirical Bayes Poisson Matrix Factorization (Background Model with weights, using subject information)}
\usage{
ebpmf_wbg2_subject(
  X,
  u,
  K,
  pm_func = list(f = ebpm::ebpm_gamma_mixture, l = ebpm::ebpm_gamma, m =
    ebpm::ebpm_gamma),
  init = NULL,
  pm_control = NULL,
  fix_option = list(m = FALSE, gl = FALSE, ql = FALSE, gf = FALSE, qf = FALSE),
  maxiter = 100,
  tol = 1e-08,
  verbose = FALSE,
  seed = 123
)
}
\arguments{
\item{X}{count matrix (dim(X) = c(n, p)).}

\item{u}{vector encoding subject (length(u) = n; elements need to be from {1,..., length(unique(u))}).}

\item{pm_func}{functions for solving the \code{ebpm} subproblem for \code{L} and \code{F};
It is a list \code{list(l, f, m)};}

\item{init}{Either \code{NULL} or \code{list(qg, m, w)}}

\item{pm_control}{control parameters for pm_func function}

\item{fix_option}{list(m, gl, ql, gf, qf) where each is either TRUE or FALSE}

\item{maxiter}{maximum number of iterations}

\item{tol}{stopping tolerance for ELBO}

\item{seed}{used when init is NULL}

\item{k}{number of topics}
}
\value{
A list containing elements:
    \describe{
      \item{\code{m}}{list of "mean" and "mean_log", "g", "kl"; 
					"mean" and "mean_log" are each a p by D background frequency matrix}
      \item{\code{qg}}{list(ql, gl,qf, gf)}
      \item{\code{ELBO}}{ELBO objective for this VEB algorithm}
     }
}
\description{
Empirical Bayes Poisson Matrix Factorization (Background Model with weights, using subject information)
}
\examples{
To add
}
