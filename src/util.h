// -*-c++-*-
#ifndef UTIL_H
#define UTIL_H

#include <Rcpp.h>
#include <bitset>
#include <R.h>
#include <Rinternals.h>
#include <boost/dynamic_bitset.hpp>
#include "parameters.h"

SEXP bitset_to_sexp(boost::dynamic_bitset<> x);

boost::dynamic_bitset<> sexp_to_bitset(SEXP x, unsigned int n);

void rcpp_out(bool quiet, std::string messsage);

std::string enum_spatial_convert(Parameters::g_spatial_type_enum v);

#endif
