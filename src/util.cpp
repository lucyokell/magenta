#include <Rcpp.h>
#include "util.h"

SEXP bitset_to_sexp(std::bitset<barcode_length> x) {
  SEXP ret = PROTECT(Rf_allocVector(RAWSXP, barcode_length));
  unsigned char *data = RAW(ret);
  for (size_t i = 0; i < barcode_length; ++i) {
    data[i] = x[i];
  }
  UNPROTECT(1);
  return ret;
}

std::bitset<barcode_length> sexp_to_bitset(SEXP x) {
  const unsigned char * data = RAW(x);
  std::bitset<barcode_length> ret;
  for (size_t i = 0; i < barcode_length; ++i) {
    ret[i] = data[i];
  }
  return ret;
}

// [[Rcpp::export]]
SEXP test_bitset_serialisation(SEXP x) {
  std::bitset<barcode_length> y = sexp_to_bitset(x);
  for (size_t i = 0; i < barcode_length; ++i) {
    Rcpp::Rcout << y[i];
  }
  Rcpp::Rcout << std::endl;
  return bitset_to_sexp(y);
}
