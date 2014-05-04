#define ARMA_64BIT_WORD
#ifndef __ACTIV_FUNC__
#define __ACTIV_FUNC__ 
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


//----------------------------------------------------------------------------
arma::mat _sigmoid(arma::mat A)
{
    A = arma::exp(-A);
    for (arma::mat::iterator elem = A.begin(); elem != A.end(); ++elem)
    {   
        *elem = 1 / (1 + *elem);
    }
    return A;
}
//----------------------------------------------------------------------------
template <typename T>
T _sigmoid_derivative(T A)
{
    return A % (1 - A);
}
//----------------------------------------------------------------------------
arma::mat _softmax(arma::mat A, unsigned int axis = 1)
{
    A = arma::exp(A);
    if (axis == 1)
    {
        arma::vec expsum = arma::sum(A, axis);
        for (int i = 0; i < A.n_rows; ++i)
        {
            A.row(i) /= expsum(i);
        }
        return A;
    } 
    else if (axis == 0)
    {
        arma::rowvec expsum = arma::sum(A, axis);   
        for (int i = 0; i < A.n_cols; ++i)
        {
            A.col(i) /= expsum(i);
        }
        return A;
    }
    throw std::runtime_error("axis for softmax can only be 1 or 0");
}
//----------------------------------------------------------------------------
template <class T>
T _identity(T &O)
{
    return O;
}
//----------------------------------------------------------------------------
template <class T>
double _identity_derivative(const T &O)
{
    return 1.0;
}
//----------------------------------------------------------------------------






#endif