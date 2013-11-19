#include <RcppArmadillo.h>
#include "fastPdist.h"
#include "neighbor_graph.h"





// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


SEXP laplacian_eigenmap(arma::mat X, int d, int k, double heat = 2.0)
{
	unsigned int n = X.n_rows;
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
					graph(i, j) = -exp((-D(i, j)) / heat);
				}
			}
		}
	}
	arma::vec weight = (arma::sum(graph, 1));
	arma::mat laplacian = graph;
	laplacian.diag() = laplacian.diag() + weight;
	for (int i = 0; i < n; ++i)
	{
		laplacian.row(i) = laplacian.row(i) / weight(i);
	}

	arma::cx_vec eigval;
	arma::cx_mat eigvec;
	arma::eig_gen(eigval, eigvec, laplacian);
	// laplacian = arma::real(eigvec.cols((n - d - 1), n - 2));
	laplacian = arma::real(eigvec);
	return Rcpp::List::create(Rcpp::Named("Y") = laplacian);


}





// Standard Diffusion Map a la Coifman et al.

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
	arma::cx_vec eigval;
	arma::cx_mat eigvec;
	arma::eig_gen(eigval, eigvec, D);
	D = arma::real(eigvec.cols(1, d));
	return Rcpp::List::create(Rcpp::Named("Y") = D);
}







SEXP local_linear_embedding(arma::mat X, int k, int d = 2, bool verbose = false)
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

    return Rcpp::List::create(Rcpp::Named("Y") = W);

}