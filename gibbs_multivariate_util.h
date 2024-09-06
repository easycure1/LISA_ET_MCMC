#ifndef __GIBBS_MULTIVARIATE_UTIL_INCLUDED__
#define __GIBBS_MULTIVARIATE_UTIL_INCLUDED__

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace std;
using namespace Rcpp;

// I/O: Only use within Rcpp
// [[Rcpp::export]]
arma::cx_cube cx_cube_from_ComplexVector(ComplexVector x);

#endif