boxcox <- function(D, A, d = 2, X1 = NULL, cmds_start = TRUE, lambda = 1, mu = 1, nu = 0, tau = 1, niter = 500, sample_rate = 100, bfgs = FALSE, adaptive = TRUE, scale_out = TRUE, verbose = FALSE)
{
	random_start <- FALSE
	if(is.null(X1))
	{
		X1 <- matrix(0, nrow(D), d)
		if(cmds_start == FALSE)
		{
			random_start <- TRUE;
		}
	}
	tmp <- (INTERNAL_BOXCOX(D, A, X1, cmds_start, random_start, d, lambda, mu, nu, tau, niter, sample_rate, bfgs, adaptive, scale_out, verbose))
	tmp$description <- paste("LMDS, lambda = ", lambda, ", mu = ", mu ,", nu = ", nu ,", tau = ", tau, sep = "")
	return(tmp)
}