// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
// AR Forecast -- noise parsed by R
NumericVector ar_forecast(NumericVector x, 
                          NumericVector ar, 
                          NumericVector noise) {
  const int k = noise.size(); // forecast length
  const int p = ar.size();
  NumericVector res(k);
  for (int i=0; i<k; ++i) {
    res[i] = noise[i];
    for (int j=0; j<p; ++j) {
      res[i] += ar[j]*x[x.length()-1-j];
    }
    x.push_back(res[i]);
  }
  return res;
}

// [[Rcpp::export]]
double pacf2ARacv(NumericVector pacf, double sigma2) {
  const unsigned p = pacf.size();
  double acv = sigma2;
  for (unsigned j = 0; j < p; ++j) {
    acv *= (1 / (1 - pacf[j]*pacf[j]));
  }
  return(acv);
}

// [[Rcpp::export]]
NumericMatrix pacf2AR(NumericVector pacf) {
  //#include <cassert>
  unsigned p = pacf.size();
  NumericMatrix arCoef(p, p);
  arCoef(p-1,p-1) = pacf[p-1];
  if (p==1) return(arCoef);
  NumericVector pacfTmp(p-1);
  for(unsigned j=0; j<p-1; ++j) pacfTmp[j] = pacf[j];
  NumericMatrix coefTmp = pacf2AR(pacfTmp);
  //assert(coefTmp.nrow() == p-1 && coefTmp.ncol == p-1);
  for(unsigned i=0; i < p-1; ++i) {
    for(unsigned j=0; j < p-1; ++j) {
      arCoef(i,j)=coefTmp(i,j);
    }
  }
  if (p==2) {
    arCoef(p-1,0) = pacf[0] * (1 - pacf[1]);
  } 
  if (p > 2) {
    for (int j=p-1; j >= 1; --j) {
      arCoef(p-1,j-1) = pacf[j-1];
      for (int r=1; r <= p-j; ++r) {
        arCoef(p-1,j-1) = arCoef(p-1,j-1) - pacf[j+r-1]*arCoef(j+r-2,r-1);
      }
    }
  }
  return(arCoef);
}

// [[Rcpp::export]]
NumericVector genEpsMAC(NumericVector zt, NumericVector ma) {
  int m = zt.size(), q = ma.size();
  NumericVector epsilon_t(m + q), masum(m);
  for (int ii = 0; ii < q; ++ii) {
    epsilon_t[ii] = 0;
  }
  for (int tt = 0; tt < m; ++tt) {
    for (int jj = 0; jj < q; ++jj) {
      masum[tt] += ma[jj] * epsilon_t[tt + q - jj - 1];
    }
    epsilon_t[tt + q] = zt[tt] - masum[tt];
  }          
  NumericVector::const_iterator first = epsilon_t.begin() + q;
  NumericVector::const_iterator last = epsilon_t.begin() + m + q;
  NumericVector epsilon_s(first, last);         
  return epsilon_s;
}

// [[Rcpp::export]]
NumericVector genEpsARC(NumericVector zt, NumericVector ar) {
  int m = zt.size(), p = ar.size();
  NumericVector epsilon_s(m - p), arsum(m - p);
  for (int tt = 0; tt < m - p; ++tt) {
    for (int jj = 0; jj < p; ++jj) {
      arsum[tt] += ar[jj] * zt[tt + p - jj - 1];
    }
    epsilon_s[tt] = zt[tt + p] - arsum[tt];
  }            
  return epsilon_s;
}

// [[Rcpp::export]]
NumericVector genEpsARMAC(NumericVector zt, NumericVector ar, NumericVector ma) {
  int m = zt.size(), p = ar.size(), q = ma.size();
  NumericVector epsilon_t(m + q - p), arsum(m - p), masum(m - q);
  for (int ii = 0; ii < q - p; ++ii) {
    epsilon_t[ii] = 0;
  }
  for (int tt = 0; tt < m - p; ++tt) {
    for (int jj = 0; jj < p; ++jj) {
      arsum[tt] += ar[jj] * zt[tt + p - jj - 1];
    }
    for (int kk = 0; kk < q; ++kk) {
      masum[tt] += ma[kk] * epsilon_t[tt + q - kk - 1];
    }
    epsilon_t[tt + q] = zt[tt + p] - arsum[tt] - masum[tt];
  }            
  NumericVector::const_iterator first = epsilon_t.begin() + q;
  NumericVector::const_iterator last = epsilon_t.begin() + m + q - p;
  NumericVector epsilon_s(first, last);         
  return epsilon_s;
}

// [[Rcpp::export]]
NumericVector psd_arma(NumericVector freq, NumericVector ar, NumericVector ma, double sigma2 = 1.0) {
  const unsigned n = freq.size();
  const unsigned p = ar.size();
  const unsigned q = ma.size();
  const double constant = sigma2 / (2*M_PI);
  NumericVector psd(n);
  for (unsigned j = 0; j < n; ++j) {
    const double lambda = freq[j];
    // compute numerator (MA part)
    std::complex<double> numerator_c(1.0, 0.0);
    for (unsigned i = 0; i < q; ++i) {
      numerator_c += ma[i]*std::polar<double>(1.0, -lambda*(double)(i+1));
    }
    // compute denominator (AR part)
    std::complex<double> denominator_c(1.0, 0.0);
    for (unsigned i = 0; i < p; ++i) {
      denominator_c -= ar[i]*std::polar<double>(1.0, -lambda*(double)(i+1));
    }
    psd[j] = constant * std::norm(numerator_c) / std::norm(denominator_c);
  }
  return(psd);
}

// [[Rcpp::export]]
NumericMatrix acvMatrix(NumericVector acv) {
  unsigned p = acv.size();
  NumericMatrix m(p, p);
  for (int i=0; i < p; ++i) {
    for (int j=0; j < p; ++j) {
      unsigned index = abs(i-j);
      m(i,j) = acv[index];
    }
  }
  return(m);
}

// [[Rcpp::export]]
double acceptanceRate(NumericVector trace) {
  unsigned rejections = 0;
  for (unsigned i=1; i < trace.length(); ++i) {
    rejections += (trace[i]==trace[i-1]);
  }
  double rejectionRate = (double)rejections / (double)trace.length();
  return 1 - rejectionRate;
}