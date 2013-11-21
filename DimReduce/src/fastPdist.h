#ifndef FASTPDIST_H
#define FASTPDIST_H 

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat fastPdist(arma::mat A, arma::mat B);
#endif