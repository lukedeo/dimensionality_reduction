#ifndef DEEP__AUTO__HH__
#define DEEP__AUTO__HH__ 

#include "autoencoder.h"
#include <vector>
#include <iomanip>

class deep_autoencoder
{
public:
	deep_autoencoder(std::vector<int> s, std::vector<std::string> types);
	deep_autoencoder();
	
	~deep_autoencoder();

	arma::mat encode(const arma::mat &X);
	arma::mat decode(const arma::mat &R);

	arma::mat reconstruct(const arma::mat &X);

	void learn(const arma::mat &X, int epochs = 2, int batch = 1, bool denoising = false);

	Rcpp::List to_list();
	void from_list(const Rcpp::List &spec);

	deep_autoencoder& set_learning(double v)
	{
		for (int i = 0; i < (stack.size()); ++i)
		{
			stack[i]->set_learning(v);
		}
		return *this;
	}
	deep_autoencoder& set_momentum(double v)
	{
		for (int i = 0; i < (stack.size()); ++i)
		{
			stack[i]->set_momentum(v);
		}
		return *this;
	}
	deep_autoencoder& set_regularization(double v)
	{
		for (int i = 0; i < (stack.size()); ++i)
		{
			stack[i]->set_regularization(v);
		}
		return *this;
	}
	deep_autoencoder& set_noise(double v)
	{
		for (int i = 0; i < (stack.size()); ++i)
		{
			stack[i]->set_noise(v);
		}
		return *this;
	}




private:

	std::vector<autoencoder*> stack;
};



deep_autoencoder::deep_autoencoder(std::vector<int> s, std::vector<std::string> types)
{
	for (int i = 0; i < (s.size() - 1); ++i)
	{
		activ_func type;
		if (types[i] == "linear") type = linear;
		else if (types[i] == "softmax") type = softmax;
		else if (types[i] == "sigmoid") type = sigmoid;
		else
		{
			throw std::runtime_error("unregonized activation function \'" + types[i] + "\'.");
		}
		stack.push_back(new autoencoder(s[i], s[i + 1], sigmoid, type));
	}
}
deep_autoencoder::deep_autoencoder()
{
}
deep_autoencoder::~deep_autoencoder()
{
	for (int i = 0; i < (stack.size()); ++i)
	{
		delete stack[i];
	}
}

inline arma::mat deep_autoencoder::encode(const arma::mat &X)
{
	arma::mat R(X);
	for (int i = 0; i < (stack.size()); ++i)
	{
		R = stack[i]->encode(R);
	}
	return R;
}

inline arma::mat deep_autoencoder::reconstruct(const arma::mat &X)
{
	return decode(encode(X));
}

inline arma::mat deep_autoencoder::decode(const arma::mat &X)
{
	arma::mat R(X);
	for (int i = (stack.size() - 1); i >= 0; --i)
	{
		R = stack[i]->decode(R);
	}
	return R;
}

inline void deep_autoencoder::learn(const arma::mat &X, int epochs, int batch, bool denoising)
{
	int n_passes = X.n_rows * epochs * stack.size();
	int ctr = 0;
	arma::mat R(X);
	for (int layer = 0; layer < (stack.size()); ++layer)
	{
		for (int i = 0; i < epochs; ++i)
		{
			for (int row = batch - 1; row < X.n_rows; ++row)
			{
				stack[layer]->sgd_update(R.rows(row - batch + 1, row), denoising);
				if (row % 2 == 0)
				{
					// std::cout << "\r" << std::setprecision(3) << ((double)ctr / n_passes) * 100.0 << "% complete." << std::flush;
					printf ("\r%.2f percent complete.", ((double)ctr / n_passes) * 100.0);
					std::cout << std::flush;
				}
				++ctr;
			}
		}
		R = stack[layer]->encode(R);
	}
}

inline Rcpp::List deep_autoencoder::to_list()
{
	Rcpp::List DAE;
	for (int i = 0; i < (stack.size()); ++i)
	{
		DAE.push_back(stack[i]->to_list());
	}
	return DAE;
}

inline void deep_autoencoder::from_list(const Rcpp::List &spec)
{
	for (int i = 0; i < (stack.size()); ++i)
	{
		delete stack[i];
	}
	stack.clear();

	for (int i = 0; i < spec.size(); ++i)
	{
		stack.push_back(new autoencoder());
		stack[i]->from_list(spec[i]);
	}
}







#endif