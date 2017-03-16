//
//  MAGENTA
//  probability.h
//
//  Created: OJ on 12/01/2017
//	Most recent edits: OJ on 10/03/2017
//
//  Distributed under the MIT software licence
//
//  Collection of required probability functions.
//
// ---------------------------------------------------------------------------

#ifndef __MAGENTA__PROBABILITY__
#define __MAGENTA__PROBABILITY__

#include <vector>

//------------------------------------------------
// draw from bernoulli distribution for cD
bool rbernoulli_cD();

//------------------------------------------------
// draw from bernoulli distribution for cT
bool rbernoulli_cT();

//------------------------------------------------
// draw from bernoulli distribution for cU
bool rbernoulli_cU();

//------------------------------------------------
// draw from exponenetial distribution for mosquito death
int rexpint_mu0_1();

//------------------------------------------------
// sample from uniform(a,b) distribution
double runif1(double a = 0, double b = 1.0);

//------------------------------------------------
// sample from uniform integer(a,b) distribution
int runiform_int_1(int a = 0, int b = 3);

//------------------------------------------------
// sample from given probability vector (that sums to pSum)
int sample1(std::vector<double> &p, double pSum = 1);

//------------------------------------------------
// sample from gamma(alpha,beta) distribution
double rgamma1(double shape, double rate);

//------------------------------------------------
// draw from beta(alpha,beta) distribution
double rbeta1(double alpha, double beta);

//------------------------------------------------
// draw from univariate normal distribution
double rnorm1(double mean, double sd);

//------------------------------------------------
// draw from lognormal distribution
double rlognorm1(double mean, double sd);

//------------------------------------------------
// draw from poisson distribution
int rpoisson1(double mean);

//------------------------------------------------
// draw from binomial distribution
int rbinomial1(int trials, double p);

//------------------------------------------------
// draw from bernoulli distribution
bool rbernoulli1(double p);

//------------------------------------------------
// draw from exponenetial distribution
double rexpdouble1(double lambda);

//------------------------------------------------
// draw from exponenetial distribution
int rexpint1(double lambda);

//------------------------------------------------
// sample without replacement from vector
std::vector<int> sample_without_replacement(std::vector<int> &v, int c);

//------------------------------------------------
// sample with replacement from vector
std::vector<int> sample_with_replacement(std::vector<int> &v, int c);

//------------------------------------------------
// sample with replacement from max integer
std::vector<int> sample_with_replacement_from_max(int max, int c);

//------------------------------------------------
// shuffle integer vector
void shuffle_integer_vector(std::vector<int> &vec);

#endif