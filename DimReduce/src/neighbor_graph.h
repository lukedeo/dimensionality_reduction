#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


inline arma::mat neighbor_graph(arma::mat D, int k = 4, bool sym = false) //not symmetric
{
	int n = D.n_rows;
	arma::mat _adj = sort(D);
	arma::rowvec _proxy = (_adj.row(k) + _adj.row(k+1)) / 2;
	_adj.zeros();

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			if ((D(i, j) <= _proxy(i)))
			{
				_adj(i, j) = 1;
				if (sym)
				{
					_adj(j, i) = 1;
				}
			}
		}
	}
	return (_adj);
}



