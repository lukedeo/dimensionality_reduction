#ifndef FASTPDIST_H
#define FASTPDIST_H 

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

inline arma::mat fastPdist(arma::mat A, arma::mat B) 
{
 
    arma::colvec An =  sum(square(A),1);
    arma::colvec Bn =  sum(square(B),1);
 
    arma::mat C = -2 * (A * B.t());
    C.each_col() += An;
    C.each_row() += Bn.t();
 
    return (sqrt(C)); 
}
#endif