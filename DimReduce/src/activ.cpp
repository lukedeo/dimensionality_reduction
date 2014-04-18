#include <RcppArmadillo.h>
#include "activation_functions.h"
#include "layer.h"
#include <stdexcept>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat pred(arma::mat A, int num_outs = 2)
{
	layer L(A.n_cols, num_outs, sigmoid);
	return L.predict(A);	
}