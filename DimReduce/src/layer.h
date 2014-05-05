#ifndef __LAYER__H__
#define __LAYER__H__ 

#include "activation_functions.h"
#include <stdexcept>

enum activ_func { linear, sigmoid, softmax };

class autoencoder;

//----------------------------------------------------------------------------
class layer
{
public:
    explicit layer(int inputs = 1, int outputs = 1, activ_func type = linear, 
        double learning = 0.02, double momentum = 0.9, 
        double regularizer = 0.001);

    ~layer( ) { }

    arma::mat predict(const arma::mat &M);

    void backpropagate(const arma::mat &V);

    void calculate_derivatives(const arma::mat &V);

    arma::mat pass_down();

    arma::mat get_W()
    {
        return arma::reshape(params.subvec(0, m_divide), m_outputs, m_inputs);
    }
    arma::mat get_b()
    {
        return params.subvec(m_divide + 1, m_total);
    }

    Rcpp::List to_list();

    void from_list(const Rcpp::List &list);

    void set_learning(double x)
    {
        learning = x;
    }
    void set_momentum(double x)
    {
        momentum = x;
    }
    void set_regularization(double x)
    {
        regularization = x;
    }

private:
    friend class autoencoder;

    arma::vec params, grad_params, old_params;

    arma::mat m_out, m_in, B;
    int m_batch_size, m_inputs, m_outputs, m_divide, m_total;
    double learning, momentum, regularization;
    activ_func layer_type;
    arma::mat to_pass_down;
};

//----------------------------------------------------------------------------
layer::layer(int inputs, int outputs, activ_func type, 
    double learning, double momentum, double regularizer):
grad_params(inputs * outputs + outputs), 
old_params(inputs * outputs + outputs), m_inputs(inputs), m_outputs(outputs), 
m_divide(inputs * outputs - 1), m_total(inputs * outputs + outputs - 1)
,learning(learning), momentum(momentum), regularization(regularizer), layer_type(type)
{
    params = arma::randn<arma::vec>(inputs * outputs + outputs) * 0.1;
    old_params.fill(0.0);
    grad_params.fill(0.0);
}
//---------------------------------------------------------------------------
arma::mat layer::predict(const arma::mat &M)
{
    m_in = M.t();
    // m_out = W * M.t();
    // arma::mat W = arma::reshape(params.subvec(0, m_divide), m_outputs, m_inputs);
    // arma::colvec b = params.subvec(m_divide + 1, m_total);
    m_out = arma::reshape(params.subvec(0, m_divide), m_outputs, m_inputs) * M.t();
    m_out.each_col() += params.subvec(m_divide + 1, m_total);

    switch(layer_type)
    {
        case linear: return m_out.t();
        case sigmoid: return _sigmoid(m_out.t());
        case softmax: return _softmax(m_out.t());
        default: throw std::logic_error("invalid layer type.");
    }
}
// //----------------------------------------------------------------------------

void layer::backpropagate(const arma::mat &V)
{
    arma::mat delta = V.t();
    if (layer_type == sigmoid)
    {
        delta = delta % _sigmoid_derivative(_sigmoid(m_out));
    }
    // std::cout << "here" << std::endl;
    arma::mat W = arma::reshape(params.subvec(0, m_divide), m_outputs, m_inputs);
    grad_params.subvec(0, m_divide) = arma::vectorise(delta * m_in.t() / V.n_rows + regularization * W);
    grad_params.subvec(m_divide + 1, m_total) = arma::sum(delta, 1) / V.n_rows;

    old_params *= momentum;
    old_params -= learning * grad_params;
    old_params.subvec(0, m_divide) -= learning * regularization * params.subvec(0, m_divide);

    params += old_params;

    to_pass_down = W.t() * delta;

}


void layer::calculate_derivatives(const arma::mat &V)
{
    arma::mat delta = V.t();
    if (layer_type == sigmoid)
    {
        delta = delta % _sigmoid_derivative(_sigmoid(m_out));
    }
    arma::mat W = arma::reshape(params.subvec(0, m_divide), m_outputs, m_inputs);
    grad_params.subvec(0, m_divide) = arma::vectorise(delta * m_in.t() / V.n_rows + regularization * W);
    grad_params.subvec(m_divide + 1, m_total) = arma::sum(delta, 1) / V.n_rows;
    to_pass_down = W.t() * delta;
}

arma::mat layer::pass_down()
{
    return to_pass_down;
}


Rcpp::List layer::to_list()
{
    std::string type;

    if (layer_type == linear) type = "linear"; 
    else if (layer_type == sigmoid) type = "sigmoid"; 
    else if (layer_type == softmax) type = "softmax"; 

    return Rcpp::List::create(Rcpp::Named("W") = get_W(), 
                              Rcpp::Named("b") = get_b(),
                              Rcpp::Named("activation") = type,
                              Rcpp::Named("learning") = learning,
                              Rcpp::Named("momentum") = momentum,
                              Rcpp::Named("regularization") = regularization);
}
//----------------------------------------------------------------------------

void layer::from_list(const Rcpp::List &list)
{
    std::string type = list["activation"];

    if (type == "linear") layer_type = linear;
    else if (type == "sigmoid") layer_type = sigmoid;
    else if (type == "softmax") layer_type = sigmoid;


    arma::mat W = list["W"];
    arma::mat b = list["b"];

    int inputs = W.n_cols;
    int outputs = W.n_rows;

    grad_params.set_size(inputs * outputs + outputs);
    params.set_size(inputs * outputs + outputs);
    old_params.set_size(inputs * outputs + outputs);

    m_inputs = (inputs);
    m_outputs = (outputs);
    m_divide = (inputs * outputs - 1);
    m_total = (inputs * outputs + outputs - 1);
    learning = (list["learning"]);
    momentum = (list["momentum"]);
    regularization = (list["regularization"]);

    params.subvec(0, m_divide) = arma::vectorise(W);

    params.subvec(m_divide + 1, m_total) = b;

    old_params.fill(0.0);
    grad_params.fill(0.0);


    // std::cout << "the type is " << type << std::endl;
}

#endif