
#include "probability.h"

// --------------------------------------------------------------------------------------------------------------------------------
// Predefined distibutions with fixed means and one trial only
// --------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------
// draw from binomial distribution for cD
bool rbernoulli_cD() {
    return(R::unif_rand()<0.0676909);
} 

//------------------------------------------------
// draw from binomial distribution for cT
bool rbernoulli_cT() {
    return(R::unif_rand()<0.02179647);
}

//------------------------------------------------
// draw from binomial distribution for cU
std::bernoulli_distribution bernoulli_dist_cU(0.0062);
bool rbernoulli_cU() {
    return(R::unif_rand()<0.0062);
}

// --------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------
// sample from uniform(a,b) distribution
double runif1(double a, double b) {
    return(R::runif(a, b));
}

//------------------------------------------------
// sample from uniform integer(a,b) distribution
int runiform_int_1(int a, int b) {
    return(static_cast<int>(R::runif(a,b) + 0.5));
}


//------------------------------------------------
// sample from given probability vector (that sums to pSum)
int sample1(std::vector<double> &p, double pSum) {
    double rand = pSum*R::unif_rand();
    double z = 0;
    for (int i = 0; i<int(p.size()); i++) {
        z += p[i];
        if (rand < z)
            return i;
    }
    return(0);
}

//------------------------------------------------
// sample from gamma(alpha,beta) distribution
double rgamma1(double shape, double rate) {

    double x = R::rgamma(shape, 1.0 / rate);
    // check for zero or infinite values (corrects for bug in some compilers)
    while (x == 0 || (1.0 / x) == 0)
    {
        x = R::rgamma(shape, 1.0 / rate);
    }
    
    return(x);
}

//------------------------------------------------
// draw from beta(alpha,beta) distribution
double rbeta1(double alpha, double beta) {
    return(R::rbeta(alpha, beta));
}

//------------------------------------------------
// draw from univariate normal distribution
double rnorm1(double mean, double sd) {
    return(R::rnorm(mean,sd));
}

//------------------------------------------------
// draw from lognormal distribution
double rlognorm1(double mean, double sd) {
    return(R::rlnorm(mean,sd));
}

//------------------------------------------------
// draw from poisson distribution
int rpoisson1(double mean) {
    return(R::rpois(mean));
}

//------------------------------------------------
// draw from binomial distribution
int rbinomial1(int trials, double p) {
    return(R::rbinom(trials,p));
}

//------------------------------------------------
// draw from bernoulli distribution
bool rbernoulli1(double p) {
    return(R::unif_rand()<p);
}

//------------------------------------------------
// draw from exponenetial distribution
double rexpdouble1(double lambda) {
    return(R::rexp(lambda));
}

//------------------------------------------------
// draw from exponenetial distribution
int rexpint1(double lambda) {
    return(static_cast<int>(R::rexp(lambda)));
}

//------------------------------------------------
// sample without replacement from vector
std::vector<int> sample_without_replacement(std::vector<int> &v, int c) {
    
    // Create return vector
    std::vector<int> samples( c );
    
    // Based on Knuth sample algorithm
    int n = c;
    int N = v.size();
    
    int t = 0; // total input records dealt with
    int m = 0; // number of items selected so far
    double u;
    
    while (m < n)
    {
        u = R::unif_rand(); // call a uniform(0,1) random number generator
        
        if ((N - t)*u >= n - m)
        {
            t++;
        }
        else
        {
            samples[m] = t;
            t++; m++;
        }
    }
    return(samples);
}

//------------------------------------------------
// sample with replacement from vector
std::vector<int> sample_with_replacement(std::vector<int> &v, int c) {
    
    // Create return vector
    std::vector<int> samples( c );
    int N = v.size();
    
    for(int i=0 ; i<c ; i++)
    { 
        samples[i] = v[runiform_int_1(0, N)]; // call a uniform(0,1) random integer
    }
    
    return(samples);
}

//------------------------------------------------
// sample with replacement from max integer
std::vector<int> sample_with_replacement_from_max(int max, int c) {
    
    // Create return vector
    std::vector<int> samples(c);
    
    for (int i=0; i<c; i++)
    {
        samples[i] = runiform_int_1(0, max); // call a uniform(0,N) random intenger
    }
    
    return(samples);
}

//------------------------------------------------
// shuffle integer vector
inline int randWrapper(const int n) { return std::floor(R::unif_rand()*n); }

void shuffle_integer_vector(std::vector<int> &vec) {

    std::random_shuffle(vec.begin(), vec.end(), randWrapper);

}
