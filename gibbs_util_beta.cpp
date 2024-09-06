
//[[Rcpp::export]]
double llike_whittle_sum(const arma::cx_cube& FZ, const arma::cx_cube& f) {
  const int d = FZ.n_cols;
  const int N = FZ.n_rows;
  const int K = FZ.n_slices;
  double res(0.0);
  for (j=0; j<N; ++j) {
    arma::cx_mat mpg_sum(d, d, fill::zeros);
    for (k=0; k<K; ++k) {
      mpg_mean += arma::trans(FZ.slice(k).row(j)) * FZ.slice(k).row(j);
    }
    double tr = arma::trace(arma::inv(f.slice(j)) * mpg_mean / K);
    res += K * log_det(f.slice(j)).real() + tr;
  }
  
  return -res;
}
