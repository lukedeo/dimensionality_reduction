#include <RcppArmadillo.h>

// [[Rcpp::export]]

SEXP BOXCOX(SEXP dist, SEXP adjacency, SEXP X1 = false, int random_start = 1, int d = 3, double lambda = 1, double mu = 1, double nu = 0, double tau = 1, int niter = 1000) 
{
    Rcpp::NumericMatrix _D(dist);
    Rcpp::NumericMatrix _X(X1);
    Rcpp::NumericMatrix _A(adjacency);

    int n = _D.nrow();

    std::vector<double> median_holder;

    arma::mat D(_D.begin(), n, n, false); 
    arma::mat A(_A.begin(), n, n, false); 
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

    if (
    {
        std::cout << "hi" << std::endl;
    }




    return Rcpp::List::create(
        Rcpp::Named("A") = A
    );

}