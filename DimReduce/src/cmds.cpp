#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

SEXP cmds(arma::mat D, double d) 
{

	int n = D.n_rows;

	arma::mat J(n, n); 
	J.ones();
	J = J / (-n);
	J.diag() += 1;

	arma::vec eigval;
	arma::mat eigvec;


	arma::eig_sym(eigval, eigvec, ((-1.0/2) * J * pow(D, 2) * J), "dc");

    return Rcpp::List::create(
    	Rcpp::Named("vectors") = eigvec,
        Rcpp::Named("values") = eigval,
        Rcpp::Named("embedding") = arma::mat(real(eigvec.cols((n - d), (n - 1))))
    );
}

