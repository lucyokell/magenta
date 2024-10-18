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
#include "probability.h"
#include "Rcpp.h"

class Drug {
  
public:
  
  
private:
  
  std::vector<double> m_lpf; // prob of lpf for each bitset combination related to this drug moda
  std::vector<unsigned int> m_barcode_positions; // which positions in the barcode correspond to this drug
  std::vector<unsigned int> m_prophylactic_positions; // which positions in the barcode correspond to the prophylactic
  double m_dur_P;					// duration of prophylaxis
  double m_dur_SPC;					// duration of slow parasite clearance
  double m_hill_n; // Hill function parameter n for curve detailing pretective efficacy
  double m_hill_kA; // Hill function parameter kA for curve detailing pretective efficacy
  double m_hill_res_n; // Hill function parameter n for curve detailing pretective efficacy when challenged by resistant parasite
  double m_hill_res_kA; // Hill function parameter kA for curve detailing pretective efficacy when challenged by resistant parasite
  std::vector<double> m_prophylactic_probability; // prob of lpf for each bitset combination related to this drug moda
  std::vector<double> m_prophylactic_resistant_probability; // prob of lpf for each bitset combination related to this drug moda
  int m_drug_clearance_max_time;
  
public:
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // CONSTRUCTORS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Class constructor using required parameters
  Drug(std::vector<double> lpf, 
       std::vector<unsigned int> barcode_positions, 
       std::vector<unsigned int> prophylactic_positions, 
       double dur_P, 
       double dur_SPC,
       double hill_n,
       double hill_kA,
       double hill_res_n,
       double hill_res_kA,
       std::vector<double> prophylactic_probability, 
       std::vector<double> prophylactic_resistant_probability,
       int drug_clearance_max_time) : 
  
  m_lpf(lpf), 
  m_barcode_positions(barcode_positions),
  m_prophylactic_positions(prophylactic_positions),
  m_dur_P(dur_P),
  m_dur_SPC(dur_SPC),
  m_hill_n(hill_n),
  m_hill_kA(hill_kA),
  m_hill_res_n(hill_res_n),
  m_hill_res_kA(hill_res_kA),
  m_prophylactic_probability(prophylactic_probability),
  m_prophylactic_resistant_probability(prophylactic_resistant_probability),
  m_drug_clearance_max_time(drug_clearance_max_time)
  {};

  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // GETTERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Get lpf
  std::vector<double> get_m_lpf() const { return(m_lpf); }
  
  // Get barcode positions
  std::vector<unsigned int> get_m_barcode_positions() const { return(m_barcode_positions); }
  
  // Get prophylactic positions
  std::vector<unsigned int> get_m_prophylactic_positions() const { return(m_prophylactic_positions); }
  
  // Get drug duration of prophylaxis
  double get_m_dur_P() const { return(m_dur_P); }		
  
  // Get max duration of drug prophylaxis
  double get_m_drug_clearance_max_time() const { return(m_drug_clearance_max_time); }		
  
  // Getduration of slow parasite clearance
  double get_m_dur_SPC() const { return(m_dur_SPC); }		
  
  // Get Hill function parameter n for curve detailing pretective efficacy
  double get_m_hill_n() const { return(m_hill_n); }	
  
  // Get Hill function parameter kA for curve detailing pretective efficacy
  double get_m_hill_kA() const { return(m_hill_kA); }	
  
  // Get Hill function parameter n for curve detailing pretective efficacy when challenged by resistant parasite
  double get_m_hill_res_n() const { return(m_hill_res_n); }	
  
  // Get Hill function parameter kA for curve detailing pretective efficacy when challenged by resistant parasite
  double get_m_hill_res_kA() const { return(m_hill_res_kA); }	
  
  // Get prophylactic_probability
  std::vector<double> get_m_prophylactic_probability() const { return(m_prophylactic_probability); }
  
  // Get m_prophylactic_resistant_probability
  std::vector<double> get_m_prophylactic_resistant_probability() const { return(m_prophylactic_resistant_probability); }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Setters
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Set lpf
  void set_m_lpf(std::vector<double> x) { m_lpf = x; }
  
  // Set barcode positions
  void set_m_barcode_positions(std::vector<unsigned int> x) { m_barcode_positions = x; }
  
  // Set prophylactic positions
  void set_m_prophylactic_positions(std::vector<unsigned int> x) { m_prophylactic_positions = x; }
  
  // Set drug duration of prophylaxis
  void set_m_dur_P(double x) { m_dur_P = x; }
  
  // Set  duration of slow parasite clearance
  void set_m_dur_SPC(double x) { m_dur_SPC = x; }
  
  // Set Hill function parameter n for curve detailing pretective efficacy
  void set_m_hill_n(double x) { m_hill_n = x; }	
  
  // Set Hill function parameter kA for curve detailing pretective efficacy
  void set_m_hill_kA(double x) { m_hill_kA = x; }	
  
  // Set Hill function parameter n for curve detailing pretective efficacy when challenged by resistant parasite
  void set_m_hill_res_n(double x) { m_hill_res_n = x; }	
  
  // Set Hill function parameter kA for curve detailing pretective efficacy when challenged by resistant parasite
  void set_m_hill_res_kA(double x) { m_hill_res_kA = x; }	
  
  // Set prophylactic_probability
  void set_m_prophylactic_probability(std::vector<double> x) { m_prophylactic_probability = x; }
  
  // Set prophylactic_resistant_probability
  void set_m_prophylactic_resistant_probability(std::vector<double> x) { m_prophylactic_resistant_probability = x; }
  
  // Set max drug prophylaxis time
  void set_m_drug_clearance_max_time(int x) { m_drug_clearance_max_time = x; }	
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Extra
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Get prob of lpf
  double get_prob_of_lpf_x(unsigned int x) const { return(m_lpf[x]); }
  
  // Get prob of lpf at position
  double get_prob_of_lpf_barcode(boost::dynamic_bitset<> &x) const { 
    
    unsigned long mask = 1;
    unsigned long result = 0;
    for (unsigned int start_bit = 0 ; start_bit < m_barcode_positions.size(); start_bit++) {
      if (x[m_barcode_positions[start_bit]]) {
        result |= mask;
      }
      mask <<= 1;
    }
    return (m_lpf[result]);
    
    }
  
  // Get prob of prophylaxis at time x
  double get_prob_of_prophylaxis_x(unsigned int x) const { return(m_prophylactic_probability[x]); }
  
  // Get prob of prophylaxis against resistant strain at time x
  double get_prob_of_prophylaxis_resistant_x(unsigned int x) const { return(m_prophylactic_resistant_probability[x]); }
  
  // Does the provided barcode have resistant loci
  bool resistance_to_drug(boost::dynamic_bitset<> &x) const { 
    
    // start as false and then and assign
    bool resistant = false;
    
    // what is the strain's lpf and def(ault)
    double def = get_prob_of_lpf_x(0);
    double lpf = get_prob_of_lpf_barcode(x);
    
    // if it is lower then check if this is because it has prophylactic resistance positions
    if (lpf < def) {
    
    //The next loop is checking if any of the positions in m_prophylactic_positions correspond to a true value in x. If so, resistant will be set to true. This is a way of determining if "resistance" exists at any of the specified positions in x.
    // p loops through all the values within m_prophylactic_positions
    for (unsigned int p : m_prophylactic_positions) {
      resistant = resistant || x[p];
    }
    
    }
    
    return(resistant);
    
  }
  
  // Get early reinfection
  bool early_reinfection(boost::dynamic_bitset<> &x,
                         unsigned int current_time,
                         unsigned int day_of_change,
                         unsigned int day_treated) const { 
    
    // is the infecting strain resistant at prophylactic positions
    bool resistant = resistance_to_drug(x);
    
    // what is their drug concentration
    // double drug_conc = R::dexp(current_time - day_treated,
    //                            1.0/(day_of_change - day_treated),
    //                            false);
    // changing this as confusingly R::dexp uses
    double drug_conc = R::dexp(current_time - day_treated,
                               (day_of_change - day_treated),
                               false);

    
    // normalise so that drug conc is 1 at t = 0
    //drug_conc = drug_conc / R::dexp(0.0, 1.0/(day_of_change - day_treated), false);
    drug_conc = drug_conc / R::dexp(0.0, (day_of_change - day_treated), false);
    
    
    std::cout << "drug_conc after normalising = " << drug_conc << "\n";
    std::cout << "about to calculate early reinfection yes/no \n ";
    std::cout << "current_time - day_treated = " << current_time - day_treated << "\n";
    
    // if it was resistant then check for early reinfection using resistant hill parameters
    if (resistant) {
      std::cout << "hill func res = " << hill_function(drug_conc, m_hill_res_n, m_hill_res_kA) << "\n";
      
      return(rbernoulli1(hill_function(drug_conc, m_hill_res_n, m_hill_res_kA)));
      
    } else {
      //print("no resistance in early reinfection!\n")
      std::cout << "hill func wt = " << hill_function(drug_conc, m_hill_n, m_hill_kA) << "\n";
      
      return(rbernoulli1(hill_function(drug_conc, m_hill_n, m_hill_kA)));
      
    }
    
  }
  
  
   //// LO make new simpler early_reinfection function which does not use drug conc, just user-entered curve.
   // Get early reinfection
   bool early_reinfection_prophylactic_probability(boost::dynamic_bitset<> &x,
                          unsigned int current_time,
                          unsigned int day_treated) const { 
     
     // is the infecting strain resistant at prophylactic positions
     bool resistant = resistance_to_drug(x);
     
     // if it was resistant then check for early reinfection using resistant hill parameters
     
     if (resistant) {
       
       return(rbernoulli1(get_prob_of_prophylaxis_x(current_time - day_treated)));
       
     } else {
       
       return(rbernoulli1(get_prob_of_prophylaxis_resistant_x(current_time - day_treated)));
       
     }
     
   }
  
  
  
  //
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // RcppList Conversion
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Class converter
  Rcpp::List drug_to_rcpp_list() {
    return(
      Rcpp::List::create(
        Rcpp::Named("m_lpf")=m_lpf, 
        Rcpp::Named("m_barcode_positions")=m_barcode_positions,
        Rcpp::Named("m_prophylactic_positions")=m_prophylactic_positions,
        Rcpp::Named("m_dur_P")=m_dur_P,
        Rcpp::Named("m_dur_SPC")=m_dur_SPC,
        Rcpp::Named("m_hill_n")=m_hill_n,
        Rcpp::Named("m_hill_kA")=m_hill_kA,
        Rcpp::Named("m_hill_res_n")=m_hill_res_n,
        Rcpp::Named("m_hill_res_kA")=m_hill_res_kA,
        Rcpp::Named("m_prophylactic_probability")=m_prophylactic_probability,
        Rcpp::Named("m_prophylactic_resistant_probability")=m_prophylactic_resistant_probability,
        Rcpp::Named("m_drug_clearance_max_time")=m_drug_clearance_max_time
      )
    );
  }
  
  // Class constructor from RcppList
  Drug(Rcpp::List list) : 
    
    m_lpf(Rcpp::as<std::vector<double> >(list["m_lpf"])), 
    m_barcode_positions(Rcpp::as<std::vector<unsigned int> >(list["m_barcode_positions"])), 
    m_prophylactic_positions(Rcpp::as<std::vector<unsigned int> >(list["m_prophylactic_positions"])), 
    m_dur_P(Rcpp::as<double>(list["dur_P"])),
    m_dur_SPC(Rcpp::as<double>(list["dur_SPC"])),
    m_hill_n(Rcpp::as<double>(list["m_hill_n"])),
    m_hill_kA(Rcpp::as<double>(list["m_hill_kA"])),
    m_hill_res_n(Rcpp::as<double>(list["m_hill_res_n"])),
    m_hill_res_kA(Rcpp::as<double>(list["m_hill_res_kA"])),
    m_prophylactic_probability(Rcpp::as<std::vector<double> >(list["m_prophylactic_probability"])), 
    m_prophylactic_resistant_probability(Rcpp::as<std::vector<double> >(list["m_prophylactic_resistant_probability"])), 
    m_drug_clearance_max_time(Rcpp::as<int>(list["m_drug_clearance_max_time"]))
    
  {};
  
  
  
};

#endif
