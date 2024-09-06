// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "gibbs_multivariate_util.h" // cx_cube_from_ComplexVector
using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
arma::cx_cube expand_matrix_process(NumericVector x,
                                    NumericVector r, 
                                    ComplexVector U_,
                                    NumericVector omega) {
  // Note: omega denotes nodes in [0,1] here (not the AGamma parameter)
  // assert x is sorted
  // assert the indices of x, r and U are consistent
  const arma::cx_cube U = cx_cube_from_ComplexVector(U_);
  const unsigned d = U.n_rows;
  const unsigned N = omega.length();
  arma::cx_cube res(d,d,N);
  unsigned x_ind = 0;
  arma::cx_mat currentSum(d,d,arma::fill::zeros);
  for (unsigned j=0; j<N; ++j) {
    if (x_ind < x.length() && x[x_ind] <= omega[j]) {
      currentSum += r[x_ind] * U.slice(x_ind);
      ++x_ind;
    }
    res.slice(j) = currentSum;
  }
  return res;
}

// [[Rcpp::export]]
NumericVector beta_fun_AGamma_process_cube(ComplexVector U_, 
                                           ComplexVector Sigma_) {
  
  const arma::cx_cube U = cx_cube_from_ComplexVector(U_);
  const arma::cx_cube Sigma = cx_cube_from_ComplexVector(Sigma_);
  
  const unsigned d = U.n_rows;
  const unsigned N = U.n_slices;
  NumericVector res(N);
  for (unsigned j=0; j<N; ++j) {
    arma::cx_mat Sigma_inv;
    
    // BEGIN debuggy
    const arma::cx_mat eye_mat = arma::eye<arma::cx_mat>(d,d);
    // if (!arma::approx_equal(Sigma.slice(j), eye_mat, "absdiff", 1e-15)) {
    //   //std::stringstream ss;
    //   Rcout << "Unexpected entry at at j=" << j << ": ";
    //   Sigma.slice(j).print(Rcout);
    //   Rcout << std::endl;
    //   //throw(std::runtime_error(""));
    // }
    try {
      Sigma_inv = Sigma.slice(j).i(); //arma::inv_sympd(Sigma.slice(j));
    } catch (const std::exception &e) {
      for (unsigned jj=0; jj<U.n_slices; ++jj) {
        Rcout << Sigma.slice(jj) << std::endl;
      }
      Rcout << j << std::endl;
      throw(e);
    }
    // END debuggy
    
    res[j] = arma::trace(Sigma_inv * U.slice(j)).real();
  }
  return res;
}
