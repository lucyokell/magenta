// -*-c++-*-
#ifndef UNIVERSE_UTILS_H
#define UNIVERSE_UTILS_H

#include <Rcpp.h>
#include "person.h"

// Returns the population's parasite genetics for ibd style summarised by pibd for given sample size and state
Rcpp::List population_get_genetics_ibd_df_n(Rcpp::List paramList);

// Returns the population's parasite genetics summarised by coi for given sample size and state
Rcpp::List population_get_genetics_df_n(Rcpp::List paramList);


#endif
