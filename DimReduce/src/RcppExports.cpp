// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// BOXCOX
SEXP BOXCOX(arma::mat D, arma::umat A, arma::mat X1, int cmds_start = 1, int random_start = 0, int d = 3, double lambda = 1, double mu = 1, double nu = 0, double tau = 1, int niter = 1000, int sample_rate = 100, bool bfgs = 0, bool adaptive = 1, bool scale_out = 1);
RcppExport SEXP DimReduce_BOXCOX(SEXP DSEXP, SEXP ASEXP, SEXP X1SEXP, SEXP cmds_startSEXP, SEXP random_startSEXP, SEXP dSEXP, SEXP lambdaSEXP, SEXP muSEXP, SEXP nuSEXP, SEXP tauSEXP, SEXP niterSEXP, SEXP sample_rateSEXP, SEXP bfgsSEXP, SEXP adaptiveSEXP, SEXP scale_outSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    arma::mat D = Rcpp::as<arma::mat >(DSEXP);
    arma::umat A = Rcpp::as<arma::umat >(ASEXP);
    arma::mat X1 = Rcpp::as<arma::mat >(X1SEXP);
    int cmds_start = Rcpp::as<int >(cmds_startSEXP);
    int random_start = Rcpp::as<int >(random_startSEXP);
    int d = Rcpp::as<int >(dSEXP);
    double lambda = Rcpp::as<double >(lambdaSEXP);
    double mu = Rcpp::as<double >(muSEXP);
    double nu = Rcpp::as<double >(nuSEXP);
    double tau = Rcpp::as<double >(tauSEXP);
    int niter = Rcpp::as<int >(niterSEXP);
    int sample_rate = Rcpp::as<int >(sample_rateSEXP);
    bool bfgs = Rcpp::as<bool >(bfgsSEXP);
    bool adaptive = Rcpp::as<bool >(adaptiveSEXP);
    bool scale_out = Rcpp::as<bool >(scale_outSEXP);
    SEXP __result = BOXCOX(D, A, X1, cmds_start, random_start, d, lambda, mu, nu, tau, niter, sample_rate, bfgs, adaptive, scale_out);
    return Rcpp::wrap(__result);
END_RCPP
}
// cmds
SEXP cmds(arma::mat D, double d);
RcppExport SEXP DimReduce_cmds(SEXP DSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    arma::mat D = Rcpp::as<arma::mat >(DSEXP);
    double d = Rcpp::as<double >(dSEXP);
    SEXP __result = cmds(D, d);
    return Rcpp::wrap(__result);
END_RCPP
}
// laplacian_eigenmap
SEXP laplacian_eigenmap(arma::mat X, int d, int k, double heat = 2.0);
RcppExport SEXP DimReduce_laplacian_eigenmap(SEXP XSEXP, SEXP dSEXP, SEXP kSEXP, SEXP heatSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    arma::mat X = Rcpp::as<arma::mat >(XSEXP);
    int d = Rcpp::as<int >(dSEXP);
    int k = Rcpp::as<int >(kSEXP);
    double heat = Rcpp::as<double >(heatSEXP);
    SEXP __result = laplacian_eigenmap(X, d, k, heat);
    return Rcpp::wrap(__result);
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP DimReduce_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    List __result = rcpp_hello_world();
    return Rcpp::wrap(__result);
END_RCPP
}
