#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

namespace nnet
{

//-----------------------------------------------------------------------------
//	Helper functions.
//-----------------------------------------------------------------------------
double sigmoid(double x) 
{
	return 1 / (1 + exp(-x));
}
double dsigmoid(double x) 
{
	return x * (1 - x);
}

arma::mat sigmoid(arma::mat A)
{
	int m = A.n_rows;
	int n = A.n_cols;

	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			A(i, j) = sigmoid(A(i, j));
		}
	}
	return A;
}

arma::mat dsigmoid(arma::mat A)
{
	int m = A.n_rows;
	int n = A.n_cols;

	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			A(i, j) = dsigmoid(A(i, j));
		}
	}
	return A;
}

arma::vec sigmoid(arma::vec A)
{
	int n = A.n_elem;

	for (int i = 0; i < n; ++i)
	{
		A(i) = sigmoid(A(i));
	}
	return A;
}

arma::vec dsigmoid(arma::vec A)
{
	int n = A.n_elem;

	for (int i = 0; i < n; ++i)
	{
		A(i) = dsigmoid(A(i));
	}
	return A;
}

//-----------------------------------------------------------------------------
//	CLass definintion
//-----------------------------------------------------------------------------

class layer
{
public:
	layer(int inputs, int outputs);

	void set_learning(double x);
	void set_momentum(double x);

	void pass(arma::vec v);
	arma::vec get_output(bool transformed = false);
	arma::vec get_delta();
	arma::mat W();



	void update_weights();
	void update_gradient(arma::mat W_grad, arma::vec b_grad);
	void set_delta(arma::vec v);


	~layer();
private:
	arma::mat W_gradient, delta_W, W;
	arma::vec b_gradient, delta_b, b, output, delta;
	double learning, momentum, regularizer;
	int ctr, batch;
};


//-----------------------------------------------------------------------------
//	Class Implementation
//-----------------------------------------------------------------------------



layer::layer(int inputs, int outputs) : W ( inputs, outputs ), W_gradient(inputs, outputs), 
										delta_W(inputs, outputs), b(outputs), 
										b_gradient(outputs), delta_b(outputs), output(outputs),
										learning(0.01), momentum(0.5), ctr(0), batch(1), regularizer(0.01)
{
	W = 2 * W.randu() - 1;
	b = 2 * b.randu() - 1;
	W_gradient.zeros();
	b_gradient.zeros();
	delta_b.zeros();
	delta_W.zeros();
}

void layer::set_learning(double x)
{
	learning = x;
}
void layer::set_momentum(double x)
{
	momentum = x;
}

void layer::pass(arma::vec v)
{
	output = W * v + b;
}
arma::vec layer::get_output(bool transformed)
{
	if (!transformed)
	{
		return output;
	}
	else
	{
		return sigmoid(output);
	}
}
void update_gradient(arma::mat W_grad, arma::vec b_grad)
{
	W_gradient = W_gradient + W_grad;
	b_gradient = b_gradient + b_grad;
	++ctr;
	if (ctr == batch)
	{
		ctr = 0;
		W_gradient = W_gradient / batch;
		b_gradient = b_gradient / batch;
		update_weights();
		W_gradient.zeros();
		b_gradient.zeros();
	}
}
void layer::update_weights()
{
	delta_W = (1 - momentum) * W_gradient + momentum * delta_W;
	delta_b = (1 - momentum) * b_gradient + momentum * delta_b;

	W = W - learning * (delta_W + regularizer * W);
	b = b - learning * delta_b;
}
void layer::set_delta(arma::vec v)
{
	delta = 
}
arma::vec get_delta()
arma::mat W()
layer::~layer()
{
}


arma::mat W_gradient, delta_W, W;
arma::vec b_gradient, delta_b, b, output;
double learning, momentum;







}






