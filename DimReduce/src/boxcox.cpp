#include <RcppArmadillo.h>
#include "fastPdist.h"
#include "neighbor_graph.h"
#include "fast_scale.h"
#include "bfgs_helper.h"
#include "boxcox_helper.h"
#include <cmath>

// get MNIST

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


SEXP boxcox_mds(arma::mat D, arma::umat A, arma::mat X1,int cmds_start = 1,
            int random_start = 0, int d = 3, double lambda = 1, 
            double mu = 1, double nu = 0, double tau = 1, int niter = 1000, 
            int sample_rate = 100, bool bfgs = 0, bool adaptive = 1, bool scale_out = 1)
{


    unsigned int n = D.n_rows;
    std::vector<double> median_holder;

    arma::uvec neighbors = sum(A, 1);
    arma::vec column_sum;

    int k_avg = ((int)(sum(neighbors) - n)) / n;

    arma::mat A_reduced, A_pure = neighbor_graph(D, k_avg);

    double Mka = 0, Mka_best = -1;

    arma::mat D_nu(n, n), 
              D_nu_lambda(n, n), 
              D_mu_emb(n, n),
              D_mu_lambda_emb(n, n), 
              M(n, n), 
              not_A(n, n), 
              X_best(n, d),
              repulsion(n, n);


    arma::umat A_sym = A;

    lambda = 1 / lambda;

    // encapsulate a bunch of stuff in a for loop so we can do 
    // it quickly without returning to old indices.

    for (int i = 0; i < n; ++i)
    {
        // We only need to visit below the main diagonal (n^2 - n) / 2 instead of n^2
        for (int j = 0; j < i; ++j)
        {
            if ((A_sym(i, j) > 0) || (A_sym(j, i) > 0))
            {
                A_sym(i, j) = 1;
                A_sym(j, i) = 1;
                not_A(i, j) = 0;
                not_A(j, i) = 0;
                D_nu(i, j) = pow(D(i, j), nu);
                D_nu(j, i) = pow(D(j, i), nu);
                D_nu_lambda(i, j) = pow(D(i, j), (nu + lambda));
                D_nu_lambda(j, i) = pow(D(j, i), (nu + lambda));
                median_holder.push_back(D_nu_lambda(i, j));
                median_holder.push_back(D_nu_lambda(j, i));
            }
            // only perform ops if we have to
            else 
            {
                A_sym(i, j) = 0;
                A_sym(j, i) = 0;
                D_nu(i, j) = 0;
                D_nu(j, i) = 0;
                not_A(i, j) = 1;
                not_A(j, i) = 1;
                D_nu_lambda(i, j) = 0;
                D_nu_lambda(j, i) = 0;
            }
        }

        // We make the diagonals zero.
        D_nu(i, i) = 0;
        D_nu_lambda(i, i) = 0;
        not_A(i, i) = 0;
    }

    
    D_nu.diag().fill(0.0);
    D_nu_lambda.diag().fill(0.0);

    double median = arma::median(arma::vec(median_holder));
    double cc = (accu(A_sym) - n); cc /= n; cc /= (n / median);
    double t = tau * cc;

    for (int i = 0; i < n; ++i)
    {
        // We can exploit symmetry again!
        for (int j = 0; j < i; ++j)
        {
            if (A_sym(i, j) == 1)
            {
                repulsion(i, j) = 0.0;
                repulsion(j, i) = 0.0;
            }
            else
            {
                repulsion(i, j) = t;
                repulsion(j, i) = t;
            }
        }
        repulsion(i, i) = 0.0;
    }

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

    X1 = (X1 * arma::norm(D, 2)) / arma::norm(D1, 2);

    arma::mat X0 = X1;

    double s1 = 999999999.0, 
           s0 = 2, 
           stepsize = 0.1;
    int i = 0;


    if (bfgs)
    {
        std::cout << "Using BFGS optimization." << std::endl;
        arma::vec s(n*d), y(n*d), search_direction(n*d), gradient_vector(n*d), gradient_vector_old(n*d);
        arma::mat B_inv(n*d, n*d);
        B_inv.eye();
        while ((stepsize > 1e-4) && (i < (niter - 1)))
        {
            X0 = X1;
            gradient_dump(X0, 
                          D1, 
                          D_nu, 
                          D_nu_lambda, 
                          repulsion, 
                          D_mu_emb, 
                          D_mu_lambda_emb, 
                          column_sum, 
                          M, 
                          Gradient, 
                          mu, 
                          lambda);
            
            gradient_vector = (to_vec(Gradient, n, d));

            if (i > 0)
            {
                y = gradient_vector - gradient_vector_old;

                double sTy = arma::dot(s, y);

                arma::mat ysT = y * s.t();

                B_inv = B_inv + ((sTy + arma::dot(y.t() * B_inv, y)) * (s * s.t())) / (sTy * sTy) - (B_inv * ysT + ysT.t() * B_inv) / sTy;
            }
            // Line search
            //----------------------------------------------------------------------------
            double base_stress = stress(D_mu_emb, D1, D_mu_lambda_emb, not_A, D_nu, D_nu_lambda, t, mu, lambda);
            double gamma = 0.78, c = 0.01;
            
            search_direction = -B_inv * gradient_vector;
            double grad_dot_dir = c * arma::dot(gradient_vector.t(), search_direction);

            stepsize = 1;
            s = search_direction;
            X1 = X0 + to_mat(s, n, d);
            D1 = fastPdist(X1,X1);
            double new_stress = stress(D_mu_emb, D1, D_mu_lambda_emb, not_A, D_nu, D_nu_lambda, t, mu, lambda);
            int itermax = 30, iter = 0;
            while((new_stress > (base_stress + stepsize * grad_dot_dir)) && (iter < itermax))
            {
                stepsize *= gamma;
                s = stepsize * search_direction;
                X1 = X0 + to_mat(s, n, d);
                D1 = fastPdist(X1,X1);
                new_stress = stress(D_mu_emb, D1, D_mu_lambda_emb, not_A, D_nu, D_nu_lambda, t, mu, lambda);
                ++iter;
            }
            gradient_vector_old = gradient_vector;

            ++i;

            s0 = s1;
            D1 = fastPdist(X1,X1);
            s1 = new_stress;
            if (((i + 1) % sample_rate) == 0)
            {
                A_reduced = neighbor_graph(D1, k_avg);
                A_reduced = A_reduced % A_pure;
                arma::vec Nk = (arma::sum(A_reduced, 1) - 1);
                Nk = Nk / (k_avg - 1);
                Mka = (double)arma::sum(Nk) / n - k_avg/n;
                if ((Mka >= Mka_best) && (i > 1))
                {
                    X_best = X1;
                    Mka_best = Mka;
                }

                std::cout << "Iteration Number: " << i + 1 << ", Stress: " << s1 << ", Mk_adj: " << Mka << std::endl;
            }
        }
    }
    else
    {
        std::cout << "Using standard gradient descent optimization";
        if (adaptive)
        {
            std::cout << " with adaptive stepsize." << std::endl;
        }
        else
        {
            std::cout << " with backbracking line search stepsize." << std::endl;
        }
        while ((stepsize > 1e-5) && (i < niter))
        {
            if (adaptive)
            {
                if ((s1 >= s0) && (i > 1))
                {
                    stepsize *= 0.5;
                    X1 = X0 - stepsize * Gradient;
                }
                else 
                {
                    stepsize *= 1.05;
                    X0 = X1;

                    D_mu_emb = arma::pow(D1, (mu - 2));
                    D_mu_emb.diag().fill(0.0);

                    D_mu_lambda_emb = arma::pow(D1, (mu + lambda - 2));
                    D_mu_lambda_emb.diag().fill(0.0);

                    M = D_nu % D_mu_lambda_emb - D_mu_emb % (D_nu_lambda + repulsion);
                    
                    column_sum = sum(M, 1);
                    Gradient.each_col() = column_sum;
                    Gradient = X0 % Gradient - M * X0;

                    Gradient = (norm(X0, 2) / norm(Gradient, 2)) * Gradient;

                    X1 = X0 - stepsize * Gradient;
                }
                D1 = fastPdist(X1,X1);
                s0 = s1;
                s1 = stress(D_mu_emb, D1, D_mu_lambda_emb, not_A, D_nu, D_nu_lambda, t, mu, lambda);
            }
            else
            {
                arma::vec s(n*d), search_direction(n*d), gradient_vector(n*d);
                X0 = X1;
                gradient_dump(X0, D1, D_nu, D_nu_lambda, repulsion, D_mu_emb, D_mu_lambda_emb, column_sum, M, Gradient, mu, lambda);
                gradient_vector = (to_vec(Gradient, n, d));

                // Line search
                //----------------------------------------------------------------------------
                double base_stress = stress(D_mu_emb, D1, D_mu_lambda_emb, not_A, D_nu, D_nu_lambda, t, mu, lambda);
                double gamma = 0.9, c = 0.002;
                
                search_direction = -1 * gradient_vector;
                double grad_dot_dir = c * arma::dot(gradient_vector.t(), search_direction);

                stepsize = 1;
                s = search_direction;
                X1 = X0 + to_mat(s, n, d);
                D1 = fastPdist(X1,X1);
                double new_stress = stress(D_mu_emb, D1, D_mu_lambda_emb, not_A, D_nu, D_nu_lambda, t, mu, lambda);
                int itermax = 30, iter = 0;
                while((new_stress > (base_stress + stepsize * grad_dot_dir)) && (iter < itermax))
                {
                    stepsize *= gamma;
                    s = stepsize * search_direction;
                    X1 = X0 + to_mat(s, n, d);
                    D1 = fastPdist(X1,X1);
                    new_stress = stress(D_mu_emb, D1, D_mu_lambda_emb, not_A, D_nu, D_nu_lambda, t, mu, lambda);
                    ++iter;
                }
                s0 = s1;
                s1 = new_stress;
            }
            ++i;
            if (((i + 1) % sample_rate) == 0)
            {
                A_reduced = neighbor_graph(D1, k_avg);
                A_reduced = A_reduced % A_pure;
                arma::vec Nk = (arma::sum(A_reduced, 1) - 1);
                Nk = Nk / (k_avg - 1);
                Mka = (double)arma::sum(Nk) / n - k_avg/n;
                if ((Mka >= Mka_best) && (i > 1))
                {
                    X_best = X1;
                    Mka_best = Mka;
                }
                std::cout << "Iteration Number: " << i + 1 << ", Stress: " << s1 << ", Mk_adj: " << Mka << std::endl;
            }
        }
    }
    if (scale_out)
    {
        fast_scale(X1);
        fast_scale(X_best);
    }
    return Rcpp::List::create(
        Rcpp::Named("embedding") = (X1),
        Rcpp::Named("best") = (X_best)
    );

}