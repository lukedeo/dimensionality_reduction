#include <RcppArmadillo.h>
#include "fastPdist.h"
#include "neighbor_graph.h"





// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


SEXP laplacian_eigenmap(arma::mat X, int d, int k, double heat = 2.0)
{
	unsigned int n = X.n_rows;
	std::cout << "Constructing heat kernel...";
	arma::mat D = fastPdist(X, X);
	arma::mat graph = neighbor_graph(D, k, true);
	graph.diag().zeros();
	if (heat != 0)
	{
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < i; ++j)
			{
				if (graph(i, j) > 0)
				{
					graph(i, j) = exp((-D(i, j)) / heat);
				}
			}
			for (int j = (i + 1); j < n; ++j)
			{
				if (graph(i, j) > 0)
				{
					graph(i, j) = exp((-D(i, j)) / heat);
				}
			}
		}
	}
	std::cout << "Done.\nConstructing laplacian...";
	arma::vec weight = (arma::sum(graph, 1));
	arma::mat laplacian = -graph;
	laplacian.diag() = laplacian.diag() + weight;
	for (int i = 0; i < n; ++i)
	{
		laplacian.row(i) = laplacian.row(i) / weight(i);
	}
	std::cout << "Done. \nCreating embedding...";
	
	arma::cx_vec cx_eigval;
	arma::cx_mat cx_eigvec;
	arma::eig_gen(cx_eigval, cx_eigvec, laplacian);

	arma::vec eigval = arma::conv_to<arma::vec>::from(cx_eigval);
	arma::mat eigvec = arma::conv_to<arma::mat>::from(cx_eigvec);
	
	
	arma::uvec sorted = arma::sort_index(eigval, 0);
	sorted = sorted.subvec(1, d);

	arma::mat Y = (eigvec.cols(sorted));
	std::cout << "Done." << std::endl;
	return Rcpp::List::create(Rcpp::Named("Y") = Y);


}