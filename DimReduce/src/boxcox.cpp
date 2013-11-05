#include <RcppArmadillo.h>
#include "fastPdist.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

SEXP BOXCOX(arma::mat D, 
            arma::mat A, 
            arma::mat X1,
            int cmds_start = 1,
            int random_start = 0, 
            int d = 3, 
            double lambda = 1, 
            double mu = 1, 
            double nu = 0, 
            double tau = 1, 
            int niter = 1000) 
{

    unsigned int n = D.n_rows;

    std::vector<double> median_holder;
    arma::vec neighbors = sum(A, 1);

    double k_avg = (sum(neighbors) - n) / n;
    double Mka = 0, temp;

    arma::mat A_sum(n, n), D_nu(n, n), D_nu_lambda(n, n);
    A_sum.each_col() = neighbors;

    arma::mat A_sym = A;

    lambda = 1 / lambda;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            if ((A_sym(i, j) > 0) || (A_sym(j, i) > 0))
            {
                A_sym(i, j) = 1;
                A_sym(j, i) = 1;
                D_nu(i, j) = pow(D(i, j), nu);
                D_nu(j, i) = pow(D(j, i), nu);
                D_nu_lambda(i, j) = pow(D(i, j), (nu + lambda));
                D_nu_lambda(j, i) = pow(D(j, i), (nu + lambda));
                median_holder.push_back(D_nu_lambda(i, j));
                median_holder.push_back(D_nu_lambda(j, i));
            }
        }
    }
    
    
    arma::vec zeros(n); zeros.fill(0.0);
    D_nu.diag() = zeros;
    D_nu_lambda.diag() = zeros;

    double median = arma::median(arma::vec(median_holder));
    double cc = (accu(A_sym) - n); cc /= n; cc /= (n * median);
    double t = tau * cc;

    arma::mat Gradient(n, d);

    if (random_start == 1)
    {
        X1.randn(n, d);
    }

    if (cmds_start && (random_start == 0))
    {
        //Do CDMS here 
        // (don't use a function because this reduces overhead and eval time)
        //----------------------------------------------------------------------------
        arma::vec values;
        arma::mat eigvec;
        arma::mat J(n, n); 
        J.ones();
        J = J / (-n);
        J.diag() += 1;
        if (n > 500)
        {
            arma::eig_sym(values, eigvec, ((-1.0/2) * J * pow(D, 2) * J), "dc");
        }
        else
        {
            arma::eig_sym(values, eigvec, ((-1.0/2) * J * pow(D, 2) * J), "standard");
        }
        X1 = arma::mat(real(eigvec.cols((n - d), (n - 1))));
        //----------------------------------------------------------------------------
        // Done

        double adj_norm = (arma::norm(D, 2) / n / n * 0.01);
        for (int i = 0; i < d; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                X1(j, i) *= values(i);
                X1(j, i) += adj_norm *  Rcpp::as<double>(Rcpp::rnorm(1));
            }
        }
    }

    arma::mat D1 = fastPdist(X1, X1);

    X1 = (X1 * arma::norm(dist, 2)) / arma::norm(D1, 2);

    arma::mat X0 = X1;
    arma::mat normgrad = X1;

    double s1 = 999999999.0, 
           s0 = 2, 
           stepsize = 0.1;

    while ((stepsize > 1e-5) && (i < niter))
    {
        if ((s1 >= s0) && (i > 1))
        {
            stepsize *= 0.5;


            X1 <- X0 - stepsize*normgrad
        }








    }



    X1 <- (X1 * norm(dist)) / norm(D1)
    s1 <- Inf
    s0 <- 2
    stepsize <- 0.1
    i <- 0










    return Rcpp::List::create(
        Rcpp::Named("A") = A
    );

}