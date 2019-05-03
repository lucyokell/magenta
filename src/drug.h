//
//  magenta
//  drug.h
//
//  Created: OJ on 13/01/2017
//	Most recent edits: OJ on 10/03/2017
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Class describing a parasite strain.
//
// ---------------------------------------------------------------------------


#ifndef DRUG_H
#define DRUG_H
#include <vector>
#include <queue>
#include <bitset>
#include <random>
#include <boost/dynamic_bitset.hpp>
#include "Rcpp.h"

class Drug {
  
public:
  
  
private:
  
  std::vector<double> m_lpf; // prob of lpf for each bitset combination related to this drug moda
  std::vector<unsigned int> m_barcode_positions; // which positions in the barcode correspond to this drug
  double m_dur_P;					// duration of prophylaxis
  double m_dur_SPC;					// duration of slow parasite clearance
  
public:
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // CONSTRUCTORS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Class constructor using required parameters
  Drug(std::vector<double> lpf, std::vector<unsigned int> barcode_positions, double dur_P = 25, double m_dur_SPC = 5) : 
  
  m_lpf(lpf), 
  m_barcode_positions(barcode_positions),
  m_dur_P(dur_P),
  m_dur_SPC(m_dur_SPC)
  
  {};

  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // GETTERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Get lpf
  std::vector<double> get_m_lpf() const { return(m_lpf); }
  
  // Get barcode positions
  std::vector<unsigned int> get_m_barcode_positions() const { return(m_barcode_positions); }
  
  // Get drug duration of prophylaxis
  double get_m_dur_P() const { return(m_dur_P); }		
  
  // Get  duration of slow parasite clearance
  double get_m_dur_SPC() const { return(m_dur_SPC); }		
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Setters
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Set lpf
  void set_m_lpf(std::vector<double> x) { m_lpf = x; }
  
  // Set barcode positions
  void set_m_barcode_positions(std::vector<unsigned int> x) { m_barcode_positions = x; }
  
  // Set drug duration of prophylaxis
  void set_m_dur_P(double x) { m_dur_P = x; }
  
  // Set  duration of slow parasite clearance
  void set_m_dur_SPC(double x) { m_dur_SPC = x; }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Extra
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Get prob of lpf
  double get_prob_of_lpf_x(unsigned int x) const { return(m_lpf[x]); }
  
  // Get prob of lpf
  double get_prob_of_lpf_barcode(boost::dynamic_bitset<> &x) const { 
    
    unsigned long mask = 1;
    unsigned long result = 0;
    for (unsigned int start_bit = 0 ; start_bit < m_barcode_positions.size(); start_bit++) {
      if (x[m_barcode_positions[start_bit]]) {
        result |= mask;
      }
      mask <<= 1;
    }
    return m_lpf[result];
    
    }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // RcppList Conversion
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Class converter
  Rcpp::List drug_to_rcpp_list() {
    return(
      Rcpp::List::create(
        Rcpp::Named("m_lpf")=m_lpf, 
        Rcpp::Named("m_barcode_positions")=m_barcode_positions,
        Rcpp::Named("m_dur_P")=m_dur_P,
        Rcpp::Named("m_dur_SPC")=m_dur_SPC
      )
    );
  }
  
  // Class constructor from RcppList
  Drug(Rcpp::List list) : 
    
    m_lpf(Rcpp::as<std::vector<double> >(list["m_lpf"])), 
    m_barcode_positions(Rcpp::as<std::vector<unsigned int> >(list["m_barcode_positions"])), 
    m_dur_P(Rcpp::as<double>(list["dur_P"])),
    m_dur_SPC(Rcpp::as<double>(list["dur_SPC"]))
    
  {};
  
  
  
};

#endif
