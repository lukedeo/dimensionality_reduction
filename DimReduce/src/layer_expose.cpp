#include <RcppArmadillo.h>
#include "autoencoder.h"


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP learn_autoencoder(arma::mat X, int n_hidden, std::string activation_type, int epochs = 10, int batch = 2, double learning = 0.02, double momentum = 0.9, double regularization = 0.001)
{
	activ_func type;
	if (activation_type == "linear") type = linear;
	else if (activation_type == "softmax") type = softmax;
	else if (activation_type == "sigmoid") type = sigmoid;
	else
	{
		throw std::runtime_error("unregonized activation function \'" + activation_type + "\'.");
	}

	autoencoder L(X.n_cols, n_hidden, sigmoid, type);
	L.set_learning(learning).set_momentum(momentum).set_regularization(regularization);
	int n_passes = X.n_rows * epochs;
	int ctr = 0;
	std::cout << std::endl;

	for (int i = 0; i < epochs; ++i)
	{
		for (int row = batch - 1; row < X.n_rows; ++row)
		{
			L.sgd_update(X.rows(row-batch+1, row));
			if (row % 2 == 0)
			{
				std::cout << "\r" << ((double)ctr / n_passes) * 100.0 << "% complete." << std::flush;
			}
			++ctr;
		}
	}
	std::cout << std::endl;
	return L.to_list();
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP learn_bfgs_autoencoder(arma::mat X, int n_hidden, std::string activation_type, int epochs = 1, int batch = 10, double learning = 0.02, double momentum = 0.9, double regularization = 0.001)
{
	activ_func type;
	if (activation_type == "linear") type = linear;
	else if (activation_type == "softmax") type = softmax;
	else if (activation_type == "sigmoid") type = sigmoid;
	else
	{
		throw std::runtime_error("unregonized activation function \'" + activation_type + "\'.");
	}

	autoencoder L(X.n_cols, n_hidden, sigmoid, type);
	int n_passes = X.n_rows * epochs;
	int ctr = 0;
	std::cout << std::endl;
	L.set_learning(learning).set_momentum(momentum).set_regularization(regularization);
	for (int i = 0; i < epochs; ++i)
	{
		for (int row = batch - 1; row < X.n_rows; ++row)
		{
			L.lbfgs_update(X.rows(row-batch+1, row));
			if (row % 2 == 0)
			{
				std::cout << "\r" << ((double)ctr / n_passes) * 100.0 << "% complete." << std::flush;
			}
			++ctr;
		}
	}
	std::cout << std::endl;
	return L.to_list();
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP learn_denoising_autoencoder(arma::mat X, int n_hidden, std::string activation_type, int epochs = 10, int batch = 2, double learning = 0.02, double momentum = 0.9, double regularization = 0.001, double noise = 0.02)
{
	activ_func type;
	if (activation_type == "linear") type = linear;
	else if (activation_type == "softmax") type = softmax;
	else if (activation_type == "sigmoid") type = sigmoid;
	else
	{
		throw std::runtime_error("unregonized activation function \'" + activation_type + "\'.");
	}

	autoencoder L(X.n_cols, n_hidden, sigmoid, type);
	int n_passes = X.n_rows * epochs;
	int ctr = 0;
	std::cout << std::endl;
	L.set_learning(learning).set_momentum(momentum).set_regularization(regularization).set_noise(noise);
	for (int i = 0; i < epochs; ++i)
	{
		for (int row = batch - 1; row < X.n_rows; ++row)
		{
			L.sgd_update(X.rows(row - batch + 1, row), true);
			if (row % 2 == 0)
			{
				std::cout << "\r" << ((double)ctr / n_passes) * 100.0 << "% complete." << std::flush;
			}
			++ctr;
		}
	}
	std::cout << std::endl;
	return L.to_list();
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat reconstruct_autoencoder(Rcpp::List autoenc, arma::mat X)
{
	autoencoder L;
	L.from_list(autoenc);
	return L.reconstruct(X);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat encode_autoencoder(Rcpp::List autoenc, arma::mat X)
{
	autoencoder L;
	L.from_list(autoenc);
	return L.encode(X);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat decode_autoencoder(Rcpp::List autoenc, arma::mat X)
{
	autoencoder L;
	L.from_list(autoenc);
	return L.decode(X);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP continue_learn_autoencoder(Rcpp::List autoenc, arma::mat X, int epochs = 10, int batch = 2, double learning = 0.02, double momentum = 0.9, double regularization = 0.001)
{
	autoencoder L;
	L.from_list(autoenc);
	L.set_learning(learning).set_momentum(momentum).set_regularization(regularization);
	int n_passes = X.n_rows * epochs;
	int ctr = 0;
	std::cout << std::endl;
	for (int i = 0; i < epochs; ++i)
	{
		for (int row = batch - 1; row < X.n_rows; ++row)
		{
			L.sgd_update(X.rows(row-batch+1, row));
			if (row % 2 == 0)
			{
				std::cout << "\r" << ((double)ctr / n_passes) * 100.0 << "% complete." << std::flush;
			}
			++ctr;
		}
	}
	std::cout << std::endl;
	return L.to_list();
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP continue_learn_bfgs_autoencoder(Rcpp::List autoenc, arma::mat X, int epochs = 1, int batch = 10, double learning = 0.02, double momentum = 0.9, double regularization = 0.001)
{
	autoencoder L;
	L.from_list(autoenc);
	L.set_learning(learning).set_momentum(momentum).set_regularization(regularization);
	int n_passes = X.n_rows * epochs;
	int ctr = 0;
	std::cout << std::endl;
	for (int i = 0; i < epochs; ++i)
	{
		for (int row = batch - 1; row < X.n_rows; ++row)
		{
			L.lbfgs_update(X.rows(row-batch+1, row));
			if (row % 2 == 0)
			{
				std::cout << "\r" << ((double)ctr / n_passes) * 100.0 << "% complete." << std::flush;
			}
			++ctr;
		}
	}
	std::cout << std::endl;
	return L.to_list();
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP continue_learn_denoising_autoencoder(Rcpp::List autoenc, arma::mat X, int epochs = 10, int batch = 2, double learning = 0.02, double momentum = 0.9, double regularization = 0.001, double noise = 0.02)
{

	autoencoder L;
	L.from_list(autoenc);
	L.set_learning(learning).set_momentum(momentum).set_regularization(regularization).set_noise(noise);
	int n_passes = X.n_rows * epochs;
	int ctr = 0;
	std::cout << std::endl;
	for (int i = 0; i < epochs; ++i)
	{
		for (int row = batch - 1; row < X.n_rows; ++row)
		{
			L.sgd_update(X.rows(row - batch + 1, row), true);
			if (row % 2 == 0)
			{
				std::cout << "\r" << ((double)ctr / n_passes) * 100.0 << "% complete." << std::flush;
			}
			++ctr;
		}
	}
	std::cout << std::endl;
	return L.to_list();
}