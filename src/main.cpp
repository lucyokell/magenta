//
//  MAGENTA
//  main.cpp
//
//  Created: Bob Verity on 06/12/2015
//
//  Distributed under the MIT software licence
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include <Rcpp.h>

using namespace std;

// [[Rcpp::export]]
Rcpp::List dummy1_cpp(Rcpp::List paramList) {
    
    // prove that C++ code is being run
    cout << "Rcpp function is working!\n";
    
    // convert inputs to native C++ classes
    int foo = Rcpp::as<int>(paramList["foo"]);
    int bar = Rcpp::as<int>(paramList["bar"]);
    
    // prove that they are here
    cout << "foo=" << foo << ", bar=" << bar << "\n";
    
    // calculate some output
    int goo = foo * bar;
    double gar = foo/double(bar);
    
    // output Rcpp list
    return Rcpp::List::create(Rcpp::Named("goo")=goo, Rcpp::Named("gar")=gar);
    
}
