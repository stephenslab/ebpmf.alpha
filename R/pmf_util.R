#' @export mle_pm
mle_pm <- function(x, s, g_init, fix_g){
	mask <- (x != 0)
	mle <- replicate(length(x),0)
	mle[mask] <- x/s
	mle[!mask] <- 1e-8 ## set 0s to be some small numbers
	posterior <- list(mean = mle, mean_log = log(mle))
	log_likelihood <- - sum(s * posterior$mean) + sum(x[mask] * log(s[mask])) + sum(x[mask]*posterior$mean_log[mask])- sum(lgamma(x[mask] + 1))
	out = list(fitted_g = list(NULL), posterior = posterior, log_likelihood = log_likelihood)
}