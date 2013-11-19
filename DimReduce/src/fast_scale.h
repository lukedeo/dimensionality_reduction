#ifndef FAST_SCALE_H
#define FAST_SCALE_H 
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

inline void fast_scale(arma::mat &X)
{
	unsigned int n = X.n_rows;
	unsigned int d = X.n_cols;
	arma::mat _sproxy(n, d), _mproxy(n, d);
	_sproxy.each_row() = arma::stddev(X);
	_mproxy.each_row() = arma::mean(X);
	X = (X - _mproxy) / _sproxy;
}
#endif