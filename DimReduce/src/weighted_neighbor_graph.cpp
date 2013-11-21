#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


arma::mat weighted_neighbor_graph(arma::mat D, int k = 4)
{
	int n = D.n_rows;
	arma::mat _adj = sort(D);
	arma::rowvec _proxy = (_adj.row(k) + _adj.row(k+1)) / 2;
	_adj.zeros();

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			if ((D(i, j) <= _proxy(i)))
			{
				_adj(i, j) = D(i, j);
				_adj(j, i) = D(j, i);
			}
		}
	}
	return (_adj);
}