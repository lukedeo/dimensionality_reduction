#include <RcppArmadillo.h>
#include "fast_scale.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

SEXP cmds(arma::mat D, double d, bool scale = true, bool verbose = false) 
{
	if (verbose)
	{
		std::cout << "Embedding...";
	}
	int n = D.n_rows;

	arma::mat J(n, n); 
	J.ones();
	J = J / (-n);
	J.diag() += 1;

	arma::vec eigval;
	arma::mat eigvec;


	arma::eig_sym(eigval, eigvec, ((-1.0/2) * J * pow(D, 2) * J), "dc");
	std::string desc = "Classical MDS";
	arma::mat Y = arma::real(eigvec.cols((n - d), (n - 1)));
	if (scale)
	{
		fast_scale(Y);
	}
	if (verbose)
	{
		std::cout << "Done." << std::endl;
	}
    return Rcpp::List::create(
    	Rcpp::Named("description") = desc,
        Rcpp::Named("Y") = Y);
}

