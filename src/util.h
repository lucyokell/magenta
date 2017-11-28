// -*-c++-*-
#ifndef UTIL_H
#define UTIL_H

#include <bitset>
#include <R.h>
#include <Rinternals.h>
#include "strain.h"

SEXP bitset_to_sexp(std::bitset<barcode_length> x);
std::bitset<barcode_length> sexp_to_bitset(SEXP x);

#endif
