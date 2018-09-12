//
//  MAGENTA
//  probability.h
//
//  Created: OJ on 12/01/2017
//	Most recent edits: OJ on 29/03/2017
//
//  Distributed under the MIT software licence
//
//  Collection of required probability functions.
//
// ---------------------------------------------------------------------------

#ifndef __MAGENTA__PROBABILITY__
#define __MAGENTA__PROBABILITY__

#include <vector>
#include <Rcpp.h>
//#include <boost/math/special_functions.hpp>
#include <random>
#include <algorithm>


//------------------------------------------------
// draw from bernoulli distribution for cD
bool rbernoulli_cD();

//------------------------------------------------
// draw from binomial distribution for cD
bool rbernoulli_cD_rcpp();

//------------------------------------------------
// draw from bernoulli distribution for cT
bool rbernoulli_cT();

//------------------------------------------------
// draw from bernoulli distribution for cU
bool rbernoulli_cU();

//------------------------------------------------
// sample from uniform(0,1) distribution
double runif0_1();

//------------------------------------------------
// sample from uniform(a,b) distribution
double runif1(double a = 0, double b = 1.0);

//------------------------------------------------
// sample from uniform integer(a,b) distribution
int runiform_int_1(int a, int b);
    
//------------------------------------------------
// sample from given probability vector (that sums to pSum)
int sample1(const std::vector<double> &p, double pSum = 1);

//------------------------------------------------
// sample given sorted random vector, cumsum vector, and n to sammple
void samplerandoms(std::vector<double> &r, std::vector<double> &p, int n, std::vector<int> &bite_storage_queue); 

//------------------------------------------------
// sample from gamma(alpha,beta) distribution
double rgamma1(double shape, double rate);

//------------------------------------------------
// sample from nbinom(size,prob) distribution
double rnbinom_mu(double size, double mu);

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

// Simple implementation of the rmultinom.c that emplaces the draws to the output vector provided
void rmultinomN(int n, std::vector<double> &prob, double p_tot, int K, std::vector<int> &output);

//------------------------------------------------
// draw from smarter binomial distribution
int rbinomialsmarter1(unsigned int trials, double p);

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