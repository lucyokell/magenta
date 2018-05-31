// -*-c++-*-
#ifndef UTIL_H
#define UTIL_H

#include <bitset>
#include <R.h>
#include <Rinternals.h>
#include <boost/dynamic_bitset.hpp>

SEXP bitset_to_sexp(boost::dynamic_bitset<> x);
boost::dynamic_bitset<> sexp_to_bitset(SEXP x);

#endif
