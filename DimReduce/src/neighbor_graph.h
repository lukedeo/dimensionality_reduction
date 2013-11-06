#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


arma::mat neighbor_graph(arma::mat D, int k = 4) //not symmetric
{
	int n = D.n_rows;
	arma::mat _adj = sort(D);
	arma::rowvec _proxy = _adj.row(k);
	_adj.fill(0.0);

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			if (D(i, j) < _proxy(j))
			{
				_adj(i, j) = 1;
			}
		}
	}
	return (_adj);
}



