// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cholesky_xFromPhi(NumericVector phi) {
  // lenghy of phi: d^2 - 1
  // length of x:   d^2
#include <math.h> // sin, cos
  const unsigned N = phi.size();
  double phiProd = 1.0;
  NumericVector x(N+1);
  for (unsigned j=0; j<N; ++j) {
    x(j) = std::cos(phi(j)) * phiProd;
    phiProd *= std::sin(phi(j));
  }
  x(N) = phiProd;
  return x;
}

// [[Rcpp::export]]
arma::cx_mat cholseky_LFromx(arma::vec x) {
#include <math.h> // sqrt
  const unsigned d = std::sqrt(x.size());
  arma::cx_mat res(d, d, arma::fill::zeros);
  unsigned k=0;
  for (int i=0; i<d; ++i) {
    for (int j=0; j<i; ++j) {
      res(i,j) = arma::cx_double(x(k), -x(k+1));
      k += 2;
    }
    res(i,i) = x(k);
    ++k;
  }
  return res;
}

// [[Rcpp::export]]
double cholesky_jacobianLogDeterminant(NumericVector phi) {
#include <math.h> // abs, cos, sine, sqrt, log
  const int N = phi.size();
  const unsigned n = std::sqrt(N+1); // actually d
  double res = 0.0;
  int i = 1;
  for (int l=1; l <= N; ++l) { // !! take care to index phi by l-1 !!
    if (l == i*i) {
      // cosine part
      const int p_l = 2*(n-i)+1;
      res += (double)p_l * std::log(std::abs(std::cos(phi[l-1])));
      ++i;
    }
    // sine part
    const int i_2 = i-1; // !! take care not to use i (but i_2) for sine part !!
    const int m = l - (i_2)*(i_2);
    const int kappa_l = n-i_2-1;
    const int lambda_l = (i_2-1)*n + 1 + m;
    const int q_l = n*n + kappa_l*n - lambda_l;
    res += (double)q_l * std::log(std::abs(std::sin(phi[l-1])));
  }
  return res;
}

// [[Rcpp::export]]
NumericVector cholesky_pVec(unsigned d) { // (67) in Mittelbach
  const unsigned N=d*d-1;
  NumericVector res(N);
  int i=1;
  for (unsigned l=1; l<=N; ++l) {
    if (l==i*i) {
      res[l-1] = 2*(d-i)+1;
      ++i;
    } else {
      res[l-1] = 0;
    }
  }
  return res;
}

// [[Rcpp::export]]
NumericVector cholesky_qVec(unsigned d) { // (68) in Mittelbach
  const unsigned N=d*d-1;
  NumericVector res(N);
  int i=1;
  for (unsigned l=1; l<=N; ++l) {
    if (l == i*i) {
      ++i;
    }
    const int i_2 = i-1; // !! take care not to use i (but i_2) for sine part !!
    const int m = l - (i_2)*(i_2);
    const int kappa_l = d-i_2-1;
    const int lambda_l = (i_2-1)*d + 1 + m;
    res[l-1] = d*d + kappa_l*d - lambda_l;
  }
  return res;
}
