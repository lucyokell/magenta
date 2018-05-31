#include <Rcpp.h>
#include "util.h"

SEXP bitset_to_sexp(boost::dynamic_bitset<> x) {
  
  SEXP ret = PROTECT(Rf_allocVector(RAWSXP, x.size()));
  unsigned char *data = RAW(ret);
  for (size_t i = 0; i < x.size(); ++i) {
    data[i] = x[i];
  }
  UNPROTECT(1);
  return ret;
}

boost::dynamic_bitset<> sexp_to_bitset(SEXP x, unsigned int n) {
  const unsigned char * data = RAW(x);
  boost::dynamic_bitset<> ret(n);
  for (size_t i = 0; i < n; ++i) {
    ret[i] = data[i];
  }
  return ret;
}

// [[Rcpp::export]]
SEXP test_bitset_serialisation(SEXP x, unsigned int n) {
  boost::dynamic_bitset<> y = sexp_to_bitset(x, n);
  for (size_t i = 0; i < n; ++i) {
    Rcpp::Rcout << y[i];
  }
  Rcpp::Rcout << std::endl;
  return bitset_to_sexp(y);
}
