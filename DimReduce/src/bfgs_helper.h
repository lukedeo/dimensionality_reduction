#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec to_vec(arma::mat X, unsigned int rows, unsigned int columns)
{
	arma::vec out(rows * columns);

	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < columns; ++j)
		{
			out(i * columns + j) = X(i, j);
		}
	}
	return out;
}

arma::mat to_mat(arma::vec X, unsigned int rows, unsigned int columns)
{
	arma::mat out(rows, columns);

	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < columns; ++j)
		{
			out(i, j) = X(i * columns + j);
		}
	}
	return out;
}