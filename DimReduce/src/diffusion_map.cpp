#include <RcppArmadillo.h>
#include "fastPdist.h"
#include "neighbor_graph.h"


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

SEXP diffusion_map(arma::mat X, unsigned int d = 2, double t = 1.0, double sigma = -1.0)
{
	unsigned int n = X.n_rows;
	arma::mat D = fastPdist(X, X);
	if (sigma == -1)
	{
		sigma = arma::max(arma::max(D, 1));
	}

	D = arma::exp(((-D % D) / (2 * (sigma * sigma))));
	for (int i = 0; i < n; ++i)
	{
		D.row(i) = D.row(i) / arma::sum(D.row(i));
	}
	D = D * D;
	arma::cx_vec cx_eigval;
	arma::cx_mat cx_eigvec;
	arma::eig_gen(cx_eigval, cx_eigvec, D);
	arma::vec eigval = arma::conv_to<arma::vec>::from(cx_eigval);
	arma::mat eigvec = arma::conv_to<arma::mat>::from(cx_eigvec);
	
	
	arma::uvec sorted = arma::sort_index(eigval, 1);
	sorted = sorted.subvec(1, d);

	D =(eigvec.cols(sorted));
	return Rcpp::List::create(Rcpp::Named("Y") = D);
}