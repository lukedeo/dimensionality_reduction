#include <RcppArmadillo.h>
#include "fastPdist.h"
#include "neighbor_graph.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]



SEXP local_linear_embedding(arma::mat X, int k = 6, int d = 2, bool verbose = false)
{
	unsigned int n = X.n_rows;
	unsigned int p = X.n_cols;
	if (verbose)
	{
		std::cout << "Constructing k-NN Graph...";
	}
	arma::mat W = fastPdist(X, X);
	arma::mat _adj = sort(W);
	arma::rowvec _proxy = (_adj.row(k) + _adj.row(k+1)) / 2;
	arma::umat Adj(n, n);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			if (W(i, j) > _proxy(i))
			{
				Adj(i, j) = 0;
			}
			else
			{
				Adj(i, j) = 1;
			}
		}
	}
	Adj.diag().zeros();
	if (verbose)
	{
		std::cout << "Done." << std::endl;
		std::cout << "Forming Local Reconstructions...";
	}
	W.zeros();

	arma::mat U, V;

	arma::vec ones_vec(k), s;
	arma::uvec i_vec(1);
	ones_vec.ones();

	for (int i = 0; i < n; ++i)
	{
		arma::mat Z = X.rows(arma::find(Adj.row(i) > 0)).t();
		Z.each_col() -= X.row(i).t();
		arma::mat C = Z.t() * Z;
		C.diag() += 0.001 * arma::trace(C);
		arma::vec w = arma::solve(C, ones_vec);
		w = w / arma::sum(w);
		i_vec(0) = i;
		W(i_vec, arma::find(Adj.row(i) > 0)) = -w.t();
	}

	W.diag() += 1;
	if (verbose)
	{
		std::cout << "Done." << std::endl;
		std::cout << "Embedding...";
	}

	arma::svd(U, s, V, W, "dc");

	W = (V.cols((n - d - 1), (n - 2)));
	if (verbose)
	{
		std::cout << "Done." << std::endl;
	}
	std::stringstream ss;
	ss << "Locally Linear Embedding, k = " << k;

	std::string desc = ss.str();

    return Rcpp::List::create(Rcpp::Named("Y") = W,
					    	  Rcpp::Named("description") = desc);

}