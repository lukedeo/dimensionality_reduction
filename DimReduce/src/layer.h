#ifndef __LAYER__H__
#define __LAYER__H__ 

#include <RcppArmadillo.h>
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

    arma::mat predict(const arma::mat &M);

    // arma::mat correct(const arma::mat &M);

    arma::mat backpropagate(const arma::mat &V);

    virtual learn_encoding(const arma::mat &X)
    {
        throw std::logic_error("standar layer cannot learn encoding.");
    }

    virtual encode(const arma::mat &X)
    {
        throw std::logic_error("standar layer cannot encode.");
    }

    arma::mat get_W(){return W;}
    arma::mat get_b(){return b;}

    Rcpp:List to_list();



private:
    arma::mat W, W_old, W_change, m_out, m_in;
    arma::colvec b, b_old, b_change;
    int m_batch_size, m_inputs, m_outputs;
    double learning, momentum, regularization;
    activ_func layer_type;
};

inline layer::layer(int inputs, int outputs, activ_func type, 
    double learning, double momentum, double regularizer):
W_old(outputs, inputs), W_change(outputs, inputs), b_old(outputs), b_change(outputs), 
learning(learning), momentum(momentum), regularization(regularizer), layer_type(type)
{
    W = arma::randn<arma::mat>(outputs, inputs) * 0.1;
    b = arma::randn<arma::colvec>(outputs) * 0.1;
    b_old.fill(0.0);
    b_change.fill(0.0);
    W_old.fill(0.0);
    W_change.fill(0.0);

}
//---------------------------------------------------------------------------
arma::mat layer::predict(const arma::mat &M)
{
    m_in = M.t();
    m_out = W * M.t();
    m_out.each_col() += b;
    switch(layer_type)
    {
        case linear: return m_out.t();
        case sigmoid: return _sigmoid(m_out.t());
        case softmax: return _softmax(m_out.t());
        default: throw std::logic_error("invalid layer type.");
    }
}
//----------------------------------------------------------------------------

arma::mat layer::backpropagate(const arma::mat &V)
{
    arma::mat delta = V.t();
    if (layer_type == sigmoid)
    {
        delta = delta % _sigmoid_derivative(_sigmoid(m_out));
    }
    W_change = delta * m_in.t() / V.n_rows;
    b_change = arma::sum(delta, 1) / V.n_rows;
    W_old = momentum * W_old - learning * (W_change + regularization * W);
    b_old = momentum * b_old - learning * b_change;
    W += W_old;
    b += b_old;

    return W.t() * delta;
}

//----------------------------------------------------------------------------
Rcpp:List layer::to_list()
{
    return Rcpp::List::create(Rcpp::Named("W") = L.get_W(), 
                              Rcpp::Named("b") = L.get_b());
}





#endif