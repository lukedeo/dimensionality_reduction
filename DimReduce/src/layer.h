#ifndef __LAYER__H__
#define __LAYER__H__ 

#include "activation_functions.h"
#include <stdexcept>

enum activ_func { linear, sigmoid, softmax };

//----------------------------------------------------------------------------
class layer
{
public:
	explicit layer(int inputs, int outputs, activ_func type, 
		double learning = 0.01, double momentum = 0.7, 
		double regularizer = 0.00001);

	~layer() = default;

	template <class T>
	T predict(const T &M);

private:
	arma::mat W, W_old, W_change;
    arma::colvec b, b_old, b_change, m_out, m_in, delta, m_dump_below;
    int m_batch_size, m_inputs, m_outputs;
    double learning, momentum, regularization;
    activ_func layer_type;
};

inline layer::layer(int inputs, int outputs, activ_func type, 
	double learning, double momentum, double regularizer):
W_old(outputs, inputs), W_change(outputs, inputs), b_old(outputs), b_change(outputs), m_out(outputs), 
m_in(inputs), delta(outputs), learning(learning), momentum(momentum), 
regularization(regularizer), layer_type(type)
{
	W = arma::randn<arma::mat>(outputs, inputs) * 0.1;
	b = arma::randn<arma::colvec>(outputs) * 0.1;
	b_old(;
		b_change;
	m_out;
	m_in;
	delta;

}

//----------------------------------------------------------------------------
inline template <class T>
T layer::predict(const T &M)
{
	throw std::invalid_argument("unrecognized input type.");
}
//----------------------------------------------------------------------------
template <>
arma::mat layer::predict(const arma::mat &M)
{
	arma::mat tmp = W * M.t();
	tmp.each_col() += b;
	switch(layer_type)
	{
		case linear: return tmp.t();
		case sigmoid: return _sigmoid(tmp.t());
		case softmax: return _softmax(tmp.t());
		default: throw std::logic_error("invalid layer type.");
	}
}
//----------------------------------------------------------------------------
template <>
arma::vec layer::predict(const arma::vec &v)
{
	switch(layer_type)
	{
		case linear: return W * v + b;
		case sigmoid: return _sigmoid(W * v + b);
		case softmax: return _softmax(W * v + b);
		default: throw std::logic_error("invalid layer type.");
	}
}

#endif