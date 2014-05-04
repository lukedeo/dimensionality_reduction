#ifndef __AUTOENCODER__HH__
#define __AUTOENCODER__HH__ 

#include "activation_functions.h"
#include "layer.h"
#include <deque>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

class autoencoder
{
public:
	autoencoder(int n_visible = 1, int n_hidden = 1, activ_func enc_func = sigmoid, activ_func dec_func = linear);
	~autoencoder() = default;
	arma::mat encode(const arma::mat &X);
	arma::mat decode(const arma::mat &R);
	arma::mat reconstruct(const arma::mat &X);

	void sgd_update(const arma::mat &X, bool denoising = false);
	void bfgs_update(const arma::mat &X);
	void lbfgs_update(const arma::mat &X);

	template <typename T>
	double gradient_magnitude(const arma::mat &X, const T &type);
	double reconstruction_error(const arma::mat &X);

	autoencoder& set_mem(const int &v)
	{
		mem = v;
		return *this;
	}
	autoencoder& set_learning(const double &v)
	{
		learning = v;
		return *this;
	}
	autoencoder& set_tau(const double &v)
	{
		tau = v;
		return *this;
	}
	autoencoder& set_lambda(const double &v)
	{
		lambda = v;
		return *this;
	}
	autoencoder& set_c(const double &v)
	{
		c = v;
		return *this;
	}

	Rcpp::List to_list()
	{
		return Rcpp::List::create(Rcpp::Named("encoder") = encoder.to_list(),
								  Rcpp::Named("decoder") = decoder.to_list());
	}

	void from_list(const Rcpp::List &list)
	{
		encoder.from_list(list["encoder"]);
		decoder.from_list(list["decoder"]);
	}

private:
	layer encoder, decoder;
	arma::mat B;
	std::deque<arma::vec> S_store, Y_store;
	std::deque<double> alphas;

	int t;
	int n_above_thresh;
	int m_max_iters;
	int mem;
	double learning;
	double tau;
	double lambda;
	double c;
	double thresh;
	int k;

	arma::vec s;
	arma::vec y;

};
//----------------------------------------------------------------------------
autoencoder::autoencoder(int n_visible, int n_hidden, activ_func enc_func, activ_func dec_func)
: encoder(n_visible, n_hidden, enc_func), decoder(n_hidden, n_visible, dec_func)
{
	int n_parm = 2 * n_visible * n_hidden + n_hidden + n_visible;
	if (n_parm < ARMA_MAX_UWORD)
	{
			// B.eye(n_parm, n_parm);
			// B.diag() *= 1e-10;
	}


	t = 0;
	n_above_thresh = 0;
	m_max_iters = 400;
	mem = 15;
	learning = 0.1;
	tau = 10;
	lambda = 0.1;
	c = 0.1;
	thresh = 1e-4;
	k = 5;
}
//----------------------------------------------------------------------------
arma::mat autoencoder::encode(const arma::mat &X)
{
	return encoder.predict(X);
}
//----------------------------------------------------------------------------
arma::mat autoencoder::decode(const arma::mat &R)
{
	return decoder.predict(R);
}
//----------------------------------------------------------------------------
arma::mat autoencoder::reconstruct(const arma::mat &X)
{
	return decoder.predict(encoder.predict(X));
}
//----------------------------------------------------------------------------
void autoencoder::sgd_update(const arma::mat &X, bool denoising)
{
	arma::mat X_tilde;
	if (!denoising)
	{
		X_tilde = reconstruct(X);
	}
	else
	{
		arma::mat X_plus_noise = X;
		X_plus_noise += arma::randn<arma::mat>(X.n_rows, X.n_cols) * 0.02;
		X_tilde = reconstruct(X_plus_noise);
	}
	 
	arma::mat error = X_tilde - X;

	decoder.backpropagate(error);
	encoder.backpropagate(decoder.pass_down().t());
}
//----------------------------------------------------------------------------
template <typename T>
double autoencoder::gradient_magnitude(const arma::mat &X, const T &type)
{
	arma::mat X_tilde = reconstruct(X);
	arma::mat error = X_tilde - X;

	decoder.calculate_derivatives(error);
	encoder.calculate_derivatives(decoder.pass_down().t());
	arma::vec gradient = arma::join_vert(encoder.grad_params, decoder.grad_params);
	return arma::norm(gradient, type);
}
//----------------------------------------------------------------------------
double autoencoder::reconstruction_error(const arma::mat &X)
{
	arma::mat X_tilde = reconstruct(X);
	arma::mat error = X_tilde - X;
	return arma::norm(error, "fro") / X.n_rows;
}
//----------------------------------------------------------------------------
void autoencoder::bfgs_update(const arma::mat &X)
{
	arma::mat X_tilde = reconstruct(X);
	arma::mat error = X_tilde - X;

	decoder.calculate_derivatives(error);
	encoder.calculate_derivatives(decoder.pass_down().t());
	arma::vec gradient = arma::join_vert(encoder.grad_params, decoder.grad_params);

	// Get search direction
	arma::vec p = -B * gradient;
	arma::vec s = (learning / c) * (tau /(tau + t + 1)) * p;


	// step in search direction
	encoder.params += s.subvec(0, encoder.grad_params.size() - 1);
	decoder.params += s.subvec(encoder.grad_params.size(), gradient.size() - 1);

	// find new error
	X_tilde = reconstruct(X);
	error = X_tilde - X;
	decoder.calculate_derivatives(error);
	encoder.calculate_derivatives(decoder.pass_down().t());

	// get gradient difference
	arma::vec y = -gradient + lambda * s;
	gradient = arma::join_vert(encoder.grad_params, decoder.grad_params);
	y += gradient;

	// check if our gradient difference is small
	if (arma::norm(y, 2) < thresh)
	{
		n_above_thresh++;
	}
	else
	{
		n_above_thresh = 0;
	}

	if (t == 0)
	{
		B.eye();
		B *= (arma::dot(s, y) / arma::dot(y, y));
	}

	double xi = 1 / arma::dot(s, y);

	arma::mat temp = -xi * s * y.t();
	temp.diag() += 1;

	B = temp * B;

	temp = -xi * y * s.t();
	temp.diag() += 1;

	B = B * temp + c * xi * s * s.t();

	++t;
}
//----------------------------------------------------------------------------
void autoencoder::lbfgs_update(const arma::mat &X)
{
	arma::mat X_tilde = reconstruct(X);
	arma::mat error = X_tilde - X;

	decoder.calculate_derivatives(error);
	encoder.calculate_derivatives(decoder.pass_down().t());
	arma::vec gradient = arma::join_vert(encoder.grad_params, decoder.grad_params);

	int n_above_thresh = 0;

	// Get search direction
	arma::vec p = -gradient;

	alphas.clear();
	for (int i = 0; i < std::min(t, mem); ++i)
	{
		double alpha = arma::dot(S_store.at(i), p) / arma::dot(S_store.at(i), Y_store.at(i));
		alphas.push_back(alpha);
		p -= alpha * Y_store.at(i);
	}

	if (t > 0)
	{
		double tmp = 0;
		for (int i = 0; i < std::min(t, mem); ++i)
		{
			tmp += arma::dot(S_store.at(i), Y_store.at(i)) / arma::dot(Y_store.at(i), Y_store.at(i));
		}
		p *= (tmp / std::min(t, mem));
	}
	else
	{
		p *= 1e-10;
	}

	for (int i = 0; i < std::min(t, mem); ++i)
	{
		p += S_store.at(i) * (alphas.at(i) - arma::dot(Y_store.at(i), p) / arma::dot(S_store.at(i), Y_store.at(i)));
	}


	s = (learning / c) * (tau /(tau + t + 1)) * p;

	// step in search direction
	encoder.params += s.subvec(0, encoder.grad_params.size() - 1);
	decoder.params += s.subvec(encoder.grad_params.size(), gradient.size() - 1);

	// find new error
	X_tilde = reconstruct(X);
	error = X_tilde - X;
	decoder.calculate_derivatives(error);
	encoder.calculate_derivatives(decoder.pass_down().t());

	// get gradient difference
	y = -gradient + lambda * s;
	gradient = arma::join_vert(encoder.grad_params, decoder.grad_params);
	y += gradient;

	// update memory bank

	S_store.push_back(s);
	Y_store.push_back(y);

	if (S_store.size() > mem)
	{
		S_store.pop_front();
		Y_store.pop_front();
	}

	// check if our gradient difference is small
	if (arma::norm(y, 2) < thresh)
	{
		n_above_thresh++;
	}
	else
	{
		n_above_thresh = 0;
	}
	++t;

}



#endif