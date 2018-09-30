#include "util.h"

SEXP bitset_to_sexp(boost::dynamic_bitset<> x, unsigned int n) {
  
  SEXP ret = PROTECT(Rf_allocVector(RAWSXP, n));
  unsigned char *data = RAW(ret);
  for (size_t i = 0; i < n; ++i) {
    data[i] = x[i];
  }
  UNPROTECT(1);
  return ret;
}

SEXP bitset_vector_to_sexp(std::vector<boost::dynamic_bitset<> > x, unsigned int n) {
  
  SEXP ret = PROTECT(Rf_allocVector(RAWSXP, n*x.size()));
  unsigned char *data = RAW(ret);
  unsigned int count = 0;
  for(unsigned int j = 0; j < x.size() ; j++){
  for (size_t i = 0; i < n; ++i, ++count) {
    data[count] = x[j][i];
  }
  }
  UNPROTECT(1);
  return ret;
}

Rcpp::RawVector bitset_vector_to_rawvector(std::vector<boost::dynamic_bitset<> > x, unsigned int n) {
  
  Rcpp::RawVector ret = Rcpp::no_init(n*x.size());
  unsigned int count = 0;
  for(unsigned int j = 0; j < x.size() ; j++){
    for (size_t i = 0; i < n; ++i, ++count) {
      ret[count] = x[j][i];
    }
  }
  return ret;
}

Rcpp::RawMatrix bitset_vector_to_raw_matrix(std::vector<boost::dynamic_bitset<> > x, unsigned int n) {
  
  Rcpp::RawMatrix ret = Rcpp::no_init(n, x.size());
  unsigned int count = 0;
  for(unsigned int j = 0; j < x.size() ; j++){
    for (unsigned int i = 0; i < n; ++i, ++count) {
      ret(i,j) = x[j][i];
    }
  }
  return ret;
}

Rcpp::RawMatrix vector_of_bitset_vectors_to_raw_matrix(std::vector<std::vector<boost::dynamic_bitset<> > > x, unsigned int n) {
  
  // size of matrix
  unsigned int rows = 0;
  for(unsigned int m = 0; m< x.size(); m++){
    rows += x[m].size();    
  }
  
  Rcpp::RawMatrix ret = Rcpp::no_init(rows, n);
  unsigned int count = 0;
  for(unsigned int i = 0; i < x.size() ; i++){
    for(unsigned int j = 0; j < x[i].size() ; j++, ++count){
      for (unsigned int k = 0; k < n; ++k) {
        ret(count,k) = x[i][j][k];
      }
    }
  }
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

void rcpp_out(bool quiet, std::string message){
  
  if(!quiet){
    Rcpp::Rcout << message;
  }
  
}

std::string enum_spatial_convert(Parameters::g_spatial_type_enum v) {
  
  switch (v)
  {
  case Parameters::NON:   return std::string("NON");
  case Parameters::ISLAND:   return std::string("ISLAND");
  case Parameters::METAPOPULATION: return std::string("METAPOPULATION");
  case Parameters::NUMBER_OF_SPATIAL_TYPE_OPTIONS: return std::string("NUMBER_OF_SPATIAL_TYPE_OPTIONS");
  default: return("Out of spatial types: serious error");
  }
  
}

// Tests
// ------------------------------
// [[Rcpp::export]]
SEXP test_bitset_serialisation(SEXP x, unsigned int n) {
  boost::dynamic_bitset<> y = sexp_to_bitset(x, n);
  for (size_t i = 0; i < n; ++i) {
    Rcpp::Rcout << y[i];
  }
  Rcpp::Rcout << std::endl;
  return bitset_to_sexp(y, n);
}
