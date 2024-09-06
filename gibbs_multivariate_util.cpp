#include <cmath> // acos
#include "gibbs_multivariate_util.h" // cx_cube_from_ComplexVector

// Vector extension of ar_forecast (see gibbs_util.R)
// [[Rcpp::export]]
arma::mat var_forecast(arma::mat x,
                           arma::mat ar,
                           arma::mat noise) {
  const int n = x.n_rows;
  const int d = x.n_cols;
  const int p = ar.n_cols / d;
  const int k = noise.n_rows; // forecast length
  arma::mat y(k+p, d);
  // starting values for auto-regression
  for (int j=0; j<p; ++j) {
    y.row(j) = x.row(n-p+j);
  }
  // fill up the rest
  for (int i=0; i<k; ++i) {
    y.row(p+i) = noise.row(i);
    for (int j=0; j<p; ++j) {
      y.row(p+i) += (ar.submat(0,j*d,d-1,(j+1)*d-1) * y.row(p+i-(j+1)).t()).t();
    }
  }
  return y.submat(p,0,p+k-1,d-1);
}

//
// I/O: Only use *within* Rcpp
//
// [[Rcpp::export]]
arma::cx_cube cx_cube_from_ComplexVector(ComplexVector x) {
  const IntegerVector dim_x = x.attr("dim");
  arma::cx_vec x_vec(x);
  return arma::cx_cube(x_vec.begin(), dim_x[0], dim_x[1], 
                       dim_x[2], true); // re-allocate memory
}
arma::cube cube_from_NumericVector(NumericVector x) {
  const IntegerVector dim_x = x.attr("dim");
  arma::vec x_vec(x);
  return arma::cube(x_vec.begin(), dim_x[0], dim_x[1], 
                    dim_x[2], true); // re-allocate memory
}
//
//
//

// [[Rcpp::export]]
arma::cx_double tr(arma::cx_mat A) {
  arma::cx_double res(0.0, 0.0);
  for (unsigned j=0; j < A.n_cols; ++j) {
    res += A(j,j);
  }
  return res;
}

// [[Rcpp::export]]
double matCond(arma::cx_mat A) { // inverse condition number
  return arma::rcond(A);
}

// [[Rcpp::export]]
bool numericalUnstable(ComplexVector f_, bool excludeBoundary, double TOL=1e-12) {
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  for (unsigned j=excludeBoundary; j<f.n_slices-excludeBoundary; ++j) {
    const double current_cond = arma::rcond(f.slice(j));
    if (current_cond < TOL) {
      return true;
    }
  }
  return false;
}

// [[Rcpp::export]]
arma::cx_cube varma_transfer2psd(ComplexVector transfer_ar_, ComplexVector transfer_ma_, arma::cx_mat sigma) {
  const arma::cx_cube transfer_ar = cx_cube_from_ComplexVector(transfer_ar_);
  const arma::cx_cube transfer_ma = cx_cube_from_ComplexVector(transfer_ma_);
  const unsigned d = transfer_ar.n_rows;
  const unsigned N = transfer_ar.n_slices;
  const double pi = std::acos(-1.0);
  arma::cx_cube res(d,d,N);
  for (unsigned j=0; j<N; ++j) {
    const arma::cx_mat transfer_ar_inv = arma::inv(transfer_ar.slice(j));
    res.slice(j) = transfer_ar_inv * transfer_ma.slice(j) * sigma * 
      transfer_ma.slice(j).t() * transfer_ar_inv.t() / 2.0 / pi;
  }
  return res;
}

// [[Rcpp::export]]
arma::cx_cube transfer_polynomial(NumericVector lambda, arma::mat coef) {
  const unsigned d = coef.n_rows;
  const unsigned p = coef.n_cols / d;
  const unsigned N = lambda.size();
  const arma::cx_mat eye(d,d,arma::fill::eye);
  arma::cx_cube res(d,d,N);
  for (unsigned l=0; l < N; ++l) {
    res.slice(l) = eye;
    for (unsigned j=0; j < p; ++j) {
      res.slice(l) += coef.submat(0,j*d,d-1,(j+1)*d-1) * std::polar<double>(1.0, -lambda[l]*(double)(j+1));
    }
  }
  return res;
}

// [[Rcpp::export]]
arma::cube realValuedPsd(ComplexVector f_) { 
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  // convert to {f11, Re(f12), Im(f12), f22}-representation:
  arma::cube res(f.n_rows, f.n_cols, f.n_slices, arma::fill::zeros);
  for (unsigned j=0; j < f.n_slices; ++j) {
    for (unsigned r=0; r < f.n_rows; ++r) {
      for (unsigned s=0; s < f.n_cols; ++s) {
        if (r <= s) {
          res(r,s,j) = f(r,s,j).real();
        } else {
          res(r,s,j) = f(r,s,j).imag();
        }
      }
    }
  }
  return res;
}


// log Whittle likelihood (multivariate), unnormalized
// Note: omega=0 and omega=pi excluded!
// [[Rcpp::export]]
double llike_whittle(const arma::cx_mat& FZ, const arma::cx_cube& f) {
  const int N = FZ.n_rows;
  double res(0.0);
  for (unsigned j=1; j < N-1; ++j) {
    const arma::cx_mat Sigma(2.0 * M_PI * f.slice(j));
    const arma::cx_vec z(FZ.row(j).st()); // transpose without hermitian
    std::complex<double> log_det_val;
    double log_det_sign;
    arma::log_det(log_det_val,log_det_sign,Sigma);
    const double log_determ = log_det_val.real();
    const arma::cx_mat zSz = arma::trans(z) * arma::inv(Sigma) * z;
    res += (-log_determ - zSz(0,0).real());
  }
  return res;
}

// Sum of the segmented log Whittle likelihood over all frequencies
// Based on the complex Wishart distribution
//[[Rcpp::export]]
double llike_whittle_sum(const arma::cx_cube& FZ, const arma::cx_cube& f) {
  const int d = FZ.n_cols;
  const int N = FZ.n_rows;
  const int K = FZ.n_slices;
  double res(0.0);
  for (int j=1; j<N-1; ++j) {
    arma::cx_mat mpg_sum(d, d, arma::fill::zeros);
    for (int k=0; k<K; ++k) {
      mpg_sum += arma::trans(FZ.slice(k).row(j)) * FZ.slice(k).row(j);
    }
    const arma::cx_mat f_new(2 * M_PI / K * f.slice(j));
    std::complex<double> tr = arma::trace(arma::inv(f_new) * mpg_sum / K);
    res += K * arma::log_det(f_new).real() + tr.real();
  }
  
  return -res;
}

// [[Rcpp::export]]
double sldcmvnorm_t(arma::cx_mat z_t, ComplexVector Sigma_t_) {
  const arma::cx_cube Sigma_t = cx_cube_from_ComplexVector(Sigma_t_);
  // sum of complex multivariate normal log densities with mean 0 and inhomogeneous sigma_t
  double res(0.0);
  for (unsigned j=0; j < z_t.n_rows; ++j) {
    const arma::cx_mat Sigma(Sigma_t.slice(j));
    const arma::cx_vec z(z_t.row(j).st()); // transpose without hermitian
    std::complex<double> log_det_val;
    double log_det_sign;
    arma::log_det(log_det_val,log_det_sign,Sigma);
    const double log_determ = log_det_val.real();
    const arma::cx_mat zSz = arma::trans(z) * arma::inv(Sigma) * z; // consider arma::inv_sympd() for larger dimensions ?
    res += (-log_determ - zSz(0,0).real());
  }
  return(res);
}

// [[Rcpp::export]]
double sldmvnorm_t(arma::mat z_t, NumericVector Sigma_t_) {
  const arma::cube Sigma_t = cube_from_NumericVector(Sigma_t_);
  // sum of multivariate normal log densities with mean 0 and inhomogeneous sigma_t
  double res(0.0);
  for (unsigned j=0; j < z_t.n_rows; ++j) {
    const arma::mat Sigma(Sigma_t.slice(j));
    const arma::vec z(z_t.row(j).t()); // transpose of real row vector
    double log_det_val;
    double log_det_sign;
    arma::log_det(log_det_val,log_det_sign,Sigma);
    //const double log_determ = log_det_val.real();
    const arma::mat zSz = arma::trans(z) * arma::inv(Sigma) * z; // consider arma::inv_sympd() for larger dimensions ?
    res += (-log_det_val - zSz(0,0));
  }
  return(res);
}

// [[Rcpp::export]]
arma::mat epsilon_var(arma::mat zt, arma::mat ar) {
  const unsigned d = zt.n_cols;
  const unsigned n = zt.n_rows;
  const unsigned p = ar.n_cols / d;
  arma::mat res(n-p, d, arma::fill::zeros);
  for (unsigned t=p; t < n; ++t) {
    res.row(t-p) = zt.row(t);
    for (unsigned tt=1; tt<=p; ++tt) {
      res.row(t-p) -= zt.row(t-tt) * ar.submat(0,(tt-1)*d,d-1,tt*d-1).t();
    }
  }
  return res;
}

// [[Rcpp::export]]
double sldmvnorm(arma::mat z_t, arma::mat Sigma) { 
  // sum of multivariate normal log densities with mean 0 and homogeneous sigma (!! log 2 pi stuff missing !!)
  double res(0.0);
  double log_det_val;
  double log_det_sign;
  arma::log_det(log_det_val,log_det_sign,Sigma);
  const arma::mat Sigma_inv = arma::inv(Sigma);
  for (unsigned j=0; j < z_t.n_rows; ++j) {
    const arma::vec z(z_t.row(j).st());
    arma::mat zSz = trans(z) * Sigma_inv * z; // consider arma::inv_sympd() for larger dimensions ?
    res += 0.5 * (-log_det_val - zSz(0,0));
  }
  return(res);
}

// [[Rcpp::export]]
arma::mat acvToeplitz(arma::mat acv) {
  const unsigned d = acv.n_rows;
  const unsigned p = acv.n_cols / d;
  arma::mat m(p*d, p*d);
  for (int i=0; i < p; ++i) {
    for (int j=0; j < p; ++j) {
      unsigned index = abs(i-j);
      m.submat(i*d,j*d,(i+1)*d-1,(j+1)*d-1) = acv.submat(0, index*d, d-1, (index+1)*d-1);
    }
  }
  return(m);
}

// [[Rcpp::export]]
arma::cx_mat get_CFZ(arma::cx_mat FZ, ComplexVector f_half_inv_, 
                     ComplexVector f_param_half_, bool excludeBoundary) {
  const arma::cx_cube f_half_inv = cx_cube_from_ComplexVector(f_half_inv_);
  const arma::cx_cube f_param_half = cx_cube_from_ComplexVector(f_param_half_);
  arma::cx_mat res(FZ.n_rows, FZ.n_cols);
  if (excludeBoundary) {
    res.row(0) = arma::cx_rowvec(FZ.n_cols, arma::fill::zeros);
    res.row(FZ.n_rows-1) = arma::cx_rowvec(FZ.n_cols, arma::fill::zeros);
  }
  for (unsigned j=excludeBoundary; j<FZ.n_rows-excludeBoundary; ++j) {
    res.row(j) = (f_param_half.slice(j) * f_half_inv.slice(j) * FZ.row(j).st()).st();
  }
  return res;
}

// [[Rcpp::export]]
arma::cx_mat get_CFZ_q(arma::cx_mat FZ, ComplexVector q_,
                       ComplexVector f_param_half_, bool excludeBoundary) {
  const arma::cx_cube q = cx_cube_from_ComplexVector(q_);
  const arma::cx_cube f_param_half = cx_cube_from_ComplexVector(f_param_half_);
  arma::cx_mat res(FZ.n_rows, FZ.n_cols);
  if (excludeBoundary) {
    res.row(0) = arma::cx_rowvec(FZ.n_cols, arma::fill::zeros);
    res.row(FZ.n_rows-1) = arma::cx_rowvec(FZ.n_cols, arma::fill::zeros);
  }
  for (unsigned j=excludeBoundary; j<FZ.n_rows-excludeBoundary; ++j) {
    res.row(j) = (
      f_param_half.slice(j) *
        arma::inv(f_param_half.slice(j) * trans(arma::chol(q.slice(j)))) *
        FZ.row(j).st()).st();
  }
  return res;
}

// 
// Some cube algebra
//
// [[Rcpp::export]]
arma::cx_cube cubeTimesMatrix(ComplexVector f_, arma::cx_mat A) {
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  arma::cx_cube res(f.n_rows, f.n_cols, f.n_slices); // Carful: No fill
  for (unsigned j=0; j<f.n_slices; ++j) {
    res.slice(j) = f.slice(j) * A;
  }
  return res;
}
// [[Rcpp::export]]
arma::cx_cube chol_cube(ComplexVector f_, bool excludeBoundary) { // ok
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  arma::cx_cube f_half(f.n_rows, f.n_cols, f.n_slices); // Carful: No fill
  if (excludeBoundary) {
    f_half.slice(0) = arma::cx_mat(f.n_rows, f.n_cols, arma::fill::zeros);
    f_half.slice(f.n_slices-1) = arma::cx_mat(f.n_rows, f.n_cols, arma::fill::zeros);
  }
  for (unsigned j=excludeBoundary; j < f.n_slices-excludeBoundary; ++j) {
    f_half.slice(j) = trans(arma::chol(f.slice(j)));
  }
  return f_half;
}
// [[Rcpp::export]]
arma::cx_cube inv_cube(ComplexVector f_, bool excludeBoundary) { // ok
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  arma::cx_cube f_inv(f.n_rows, f.n_cols, f.n_slices); // Carful: No fill
  if (excludeBoundary) {
    f_inv.slice(0) = arma::cx_mat(f.n_rows, f.n_cols, arma::fill::zeros);
    f_inv.slice(f.n_slices-1) = arma::cx_mat(f.n_rows, f.n_cols, arma::fill::zeros);
  }
  for (unsigned j=excludeBoundary; j < f.n_slices-excludeBoundary; ++j) {
    f_inv.slice(j) = arma::inv(f.slice(j));
  }
  return f_inv;
}
// [[Rcpp::export]]
arma::cx_cube mult_cube(ComplexVector a_, ComplexVector b_) { // ok
  const arma::cx_cube a = cx_cube_from_ComplexVector(a_);
  const arma::cx_cube b = cx_cube_from_ComplexVector(b_);
  arma::cx_cube c(a.n_rows, a.n_cols, a.n_slices);
  for (unsigned j=0; j<a.n_slices; ++j) {
    c.slice(j) = a.slice(j) * b.slice(j);
  }
  return c;
}
// [[Rcpp::export]]
NumericVector logdet_cube(ComplexVector f_, bool excludeBoundary) { // ok
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  const unsigned N = f.n_slices;
  NumericVector res(N);
  if (excludeBoundary) {
    res(0) = 0;
    res(N-1) = 0;
  }
  for (unsigned j=excludeBoundary; j<N-excludeBoundary; ++j) {
    std::complex<double> log_det_val;
    double log_det_sign;
    arma::log_det(log_det_val,log_det_sign,f.slice(j));
    res(j) = log_det_val.real();
  }
  return res;
}
// [[Rcpp::export]]
arma::cx_cube const_cube(arma::cx_mat sigma, unsigned N) { // ok
  arma::cx_cube res(sigma.n_rows, sigma.n_cols, N);
  for (unsigned j=0; j<N; ++j) {
    res.slice(j) = sigma;
  }
  return res;
}
// [[Rcpp::export]]
arma::cx_cube trans_cube(ComplexVector f_) { // ok
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  arma::cx_cube res(f.n_rows, f.n_cols, f.n_slices);
  for (unsigned j=0; j<f.n_slices; ++j) {
    res.slice(j) = arma::trans(f.slice(j)); // Hermitian conjugate
  }
  return res;
}
// [[Rcpp::export]]
arma::cx_cube rev_cube(ComplexVector f_) {
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  const unsigned N = f.n_slices;
  arma::cx_cube res(f.n_rows, f.n_cols, N);
  for (unsigned j=0; j<N; ++j) {
    res.slice(j) = f.slice(N-1-j);
  }
  return res;
}
// [[Rcpp::export]]
arma::cx_cube c_cube(ComplexVector f_, ComplexVector g_) {
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  const arma::cx_cube g = cx_cube_from_ComplexVector(g_);
  // TODO Check that dim(f)[i]==dim(g)[i] for i=1,2
  arma::cx_cube res(f.n_rows, f.n_cols, f.n_slices+g.n_slices);
  for (unsigned j=0; j<f.n_slices; ++j) {
    res.slice(j) = f.slice(j);
  }
  for (unsigned k=0; k<g.n_slices; ++k) {
    res.slice(f.n_slices+k) = g.slice(k);
  }
  return res;
}

// [[Rcpp::export]]
arma::cx_mat rcWishart(unsigned nu, arma::cx_mat Sigma_half) {
  const unsigned d = Sigma_half.n_rows;
  arma::cx_mat res(d, d, arma::fill::zeros);
  const arma::cx_double i(0,1);
  for (unsigned j=0; j < nu; ++j) {
    arma::cx_vec innov(d);
    for (unsigned l=0; l < d; ++l) {
      const arma::vec real_innov = arma::randn(2);
      innov(l) = (arma::cx_double(real_innov[0], real_innov[1])) / std::sqrt(2.0);
    }
    const arma::cx_vec X(Sigma_half * innov);
    res += X * arma::trans(X);
  }
  return res;
}

// [[Rcpp::export]]
arma::cx_mat chol_cpp(arma::cx_mat A) {
  // Carful: Cholesky defined as lower triangular L here! (LL*=A)
  return trans(arma::chol(A));
}

// [[Rcpp::export]]
bool hasEigenValueSmallerZero(arma::cx_mat A, double TOL=0.0) {
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  arma::eig_gen(eigval, eigvec, A);
  bool smallerZero = false;
  for (unsigned j=0; j < eigval.size(); ++j) {
    if (eigval(j).real() < TOL) {
      smallerZero = true;
    }
  }
  return smallerZero;
}
