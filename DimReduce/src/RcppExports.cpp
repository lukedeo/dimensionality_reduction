// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// INTERNAL_BOXCOX
SEXP INTERNAL_BOXCOX(arma::mat D, arma::umat A, arma::mat X1, int cmds_start = 1, int random_start = 0, int d = 3, double lambda = 1, double mu = 1, double nu = 0, double tau = 1, int niter = 1000, int sample_rate = 100, bool bfgs = 0, bool adaptive = 1, bool scale_out = 1, bool verbose = false);
RcppExport SEXP DimReduce_INTERNAL_BOXCOX(SEXP DSEXP, SEXP ASEXP, SEXP X1SEXP, SEXP cmds_startSEXP, SEXP random_startSEXP, SEXP dSEXP, SEXP lambdaSEXP, SEXP muSEXP, SEXP nuSEXP, SEXP tauSEXP, SEXP niterSEXP, SEXP sample_rateSEXP, SEXP bfgsSEXP, SEXP adaptiveSEXP, SEXP scale_outSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP );
        Rcpp::traits::input_parameter< arma::umat >::type A(ASEXP );
        Rcpp::traits::input_parameter< arma::mat >::type X1(X1SEXP );
        Rcpp::traits::input_parameter< int >::type cmds_start(cmds_startSEXP );
        Rcpp::traits::input_parameter< int >::type random_start(random_startSEXP );
        Rcpp::traits::input_parameter< int >::type d(dSEXP );
        Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP );
        Rcpp::traits::input_parameter< double >::type mu(muSEXP );
        Rcpp::traits::input_parameter< double >::type nu(nuSEXP );
        Rcpp::traits::input_parameter< double >::type tau(tauSEXP );
        Rcpp::traits::input_parameter< int >::type niter(niterSEXP );
        Rcpp::traits::input_parameter< int >::type sample_rate(sample_rateSEXP );
        Rcpp::traits::input_parameter< bool >::type bfgs(bfgsSEXP );
        Rcpp::traits::input_parameter< bool >::type adaptive(adaptiveSEXP );
        Rcpp::traits::input_parameter< bool >::type scale_out(scale_outSEXP );
        Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP );
        SEXP __result = INTERNAL_BOXCOX(D, A, X1, cmds_start, random_start, d, lambda, mu, nu, tau, niter, sample_rate, bfgs, adaptive, scale_out, verbose);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cmds
SEXP cmds(arma::mat D, double d, bool scale = true, bool verbose = false);
RcppExport SEXP DimReduce_cmds(SEXP DSEXP, SEXP dSEXP, SEXP scaleSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP );
        Rcpp::traits::input_parameter< double >::type d(dSEXP );
        Rcpp::traits::input_parameter< bool >::type scale(scaleSEXP );
        Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP );
        SEXP __result = cmds(D, d, scale, verbose);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// diffusion_map
SEXP diffusion_map(arma::mat X, unsigned int d = 2, double t = 1.0, double sigma = -1.0, bool verbose = false);
RcppExport SEXP DimReduce_diffusion_map(SEXP XSEXP, SEXP dSEXP, SEXP tSEXP, SEXP sigmaSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        Rcpp::traits::input_parameter< unsigned int >::type d(dSEXP );
        Rcpp::traits::input_parameter< double >::type t(tSEXP );
        Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP );
        Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP );
        SEXP __result = diffusion_map(X, d, t, sigma, verbose);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// fastPdist
arma::mat fastPdist(arma::mat A, arma::mat B);
RcppExport SEXP DimReduce_fastPdist(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP );
        Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP );
        arma::mat __result = fastPdist(A, B);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// laplacian_eigenmap
SEXP laplacian_eigenmap(arma::mat X, int d, int k, double heat = 2.0, bool verbose = false);
RcppExport SEXP DimReduce_laplacian_eigenmap(SEXP XSEXP, SEXP dSEXP, SEXP kSEXP, SEXP heatSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        Rcpp::traits::input_parameter< int >::type d(dSEXP );
        Rcpp::traits::input_parameter< int >::type k(kSEXP );
        Rcpp::traits::input_parameter< double >::type heat(heatSEXP );
        Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP );
        SEXP __result = laplacian_eigenmap(X, d, k, heat, verbose);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// learn_layer
SEXP learn_layer(arma::mat X, arma::mat Y);
RcppExport SEXP DimReduce_learn_layer(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP );
        SEXP __result = learn_layer(X, Y);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// layer_predict
arma::mat layer_predict(Rcpp::List list, arma::mat X);
RcppExport SEXP DimReduce_layer_predict(SEXP listSEXP, SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP );
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        arma::mat __result = layer_predict(list, X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// local_linear_embedding
SEXP local_linear_embedding(arma::mat X, int k = 6, int d = 2, bool verbose = false);
RcppExport SEXP DimReduce_local_linear_embedding(SEXP XSEXP, SEXP kSEXP, SEXP dSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        Rcpp::traits::input_parameter< int >::type k(kSEXP );
        Rcpp::traits::input_parameter< int >::type d(dSEXP );
        Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP );
        SEXP __result = local_linear_embedding(X, k, d, verbose);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// neighbor_graph
arma::mat neighbor_graph(arma::mat D, int k, bool sym);
RcppExport SEXP DimReduce_neighbor_graph(SEXP DSEXP, SEXP kSEXP, SEXP symSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP );
        Rcpp::traits::input_parameter< int >::type k(kSEXP );
        Rcpp::traits::input_parameter< bool >::type sym(symSEXP );
        arma::mat __result = neighbor_graph(D, k, sym);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// weighted_neighbor_graph
arma::mat weighted_neighbor_graph(arma::mat D, int k = 4);
RcppExport SEXP DimReduce_weighted_neighbor_graph(SEXP DSEXP, SEXP kSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP );
        Rcpp::traits::input_parameter< int >::type k(kSEXP );
        arma::mat __result = weighted_neighbor_graph(D, k);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
