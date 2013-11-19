#include <RcppArmadillo.h>


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

void gradient_dump(arma::mat &X0, 
				   arma::mat &D1, 
				   arma::mat &D_nu, 
				   arma::mat &D_nu_lambda, 
				   arma::mat &repulsion, 
				   arma::mat &D_mu_emb, 
				   arma::mat &D_mu_lambda_emb, 
				   arma::mat &column_sum,
				   arma::mat &M, 
				   arma::mat &Gradient,
				   double mu,
				   double lambda)
{
	D_mu_emb = arma::pow(D1, (mu - 2));
	D_mu_emb.diag().fill(0.0);
	D_mu_lambda_emb = arma::pow(D1, (mu + lambda - 2));
	D_mu_lambda_emb.diag().fill(0.0);
	M = D_nu % D_mu_lambda_emb - D_mu_emb % (D_nu_lambda + repulsion);
	column_sum = sum(M, 1);
	Gradient.each_col() = column_sum;
	Gradient = X0 % Gradient - M * X0;
	Gradient = (norm(X0, 2) / norm(Gradient, 2)) * Gradient;
}




double stress(arma::mat &D_mu_emb,
			  arma::mat &D1,
			  arma::mat &D_mu_lambda_emb,
			  arma::mat &not_A,
			  arma::mat &D_nu,
			  arma::mat &D_nu_lambda,
			  double t,
			  double mu,
			  double lambda)
{
	D_mu_emb = arma::pow(D1, (mu));
    D_mu_emb.diag().fill(0.0);

    D_mu_lambda_emb = arma::pow(D1, (mu + lambda));
    D_mu_lambda_emb.diag().fill(0.0);

    if ((mu + lambda) == 0)
    {
        D1.diag().fill(1.0);
        D_mu_emb = D_mu_emb - 1;
        return arma::accu(D_nu % arma::log(D1)) - arma::accu((D_mu_emb) % D_nu_lambda) / mu - t * arma::accu((D_mu_emb) % (not_A)) / mu;
    }
    if (mu == 0)
    {
        D1.diag().fill(1.0);
        return arma::accu(D_nu % (D_mu_lambda_emb - 1)) / (mu + lambda) - arma::accu(arma::log(D1) % D_nu_lambda) - t * arma::accu(arma::log(D1) % (not_A));
    }
    else
    {
        D_mu_emb = D_mu_emb - 1;
        return arma::accu(D_nu % (D_mu_lambda_emb - 1)) / (mu + lambda) - arma::accu((D_mu_emb) % D_nu_lambda)/mu-t*arma::accu((D_mu_emb) % (not_A)) / mu;
    }
}




