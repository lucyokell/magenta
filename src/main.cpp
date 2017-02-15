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

#include <RcppArmadillo.h>
#include "stdafx.h"
#include <iostream>
#include "parameters.h"
#include "probability.h"
#include <cassert> // for error checking
#include "strain.h"
#include "person.h"
#include <chrono>
#include <functional>
#include <numeric>  
#include <algorithm>

using namespace std;
using namespace Rcpp;
//#define NDEBUG; // This causes all asserts to not be used - good for use after debugging and development has finished

// Static initialisation
int Person::s_person_ID_generator = 0;
int Strain::s_strain_ID_generator = 0;
double Person::s_psi_sum = 0;

// Maternal statics
double Person::s_sum_maternal_immunity = 0;
int Person::s_total_mums = 0;
double Person::s_mean_maternal_immunity = 70; // TODO: chang this to be initialised by R read ins

// Create vector of all barcode sequences
std::vector <barcode_t*> g_barcodes{ 1000000 };


// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// START: MAIN
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List Simulation_cpp(Rcpp::List paramList)
{
  
  // prove that C++ code is being run
  Rcpp::Rcout << "Rcpp function is working!\n";
  
  // start timer
  chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
  
  // Initialise file scope parameters object
  // This will need to be declared at the beginning of each header/cpp file where parameters are needed
  extern Parameters parameters;
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  Rcpp::List eqSS = paramList["eqSS"];
  Rcpp::NumericMatrix Smat(Rcpp::as<NumericMatrix>(eqSS["Smat"]));
  Rcpp::NumericMatrix Dmat(Rcpp::as<NumericMatrix>(eqSS["Dmat"]));
  Rcpp::NumericMatrix Amat(Rcpp::as<NumericMatrix>(eqSS["Amat"]));
  Rcpp::NumericMatrix Umat(Rcpp::as<NumericMatrix>(eqSS["Umat"]));
  Rcpp::NumericMatrix Tmat(Rcpp::as<NumericMatrix>(eqSS["Tmat"]));
  Rcpp::NumericMatrix Pmat(Rcpp::as<NumericMatrix>(eqSS["Pmat"]));
  Rcpp::NumericMatrix IBmat(Rcpp::as<NumericMatrix>(eqSS["IBmat"]));
  Rcpp::NumericMatrix ICAmat(Rcpp::as<NumericMatrix>(eqSS["ICAmat"]));
  Rcpp::NumericMatrix ICMmat(Rcpp::as<NumericMatrix>(eqSS["ICMmat"]));
  Rcpp::NumericMatrix IDmat(Rcpp::as<NumericMatrix>(eqSS["IDmat"]));
  vector<double> age_brackets = Rcpp::as<vector<double> >(eqSS["age_brackets"]);
  vector<double> het_brackets = Rcpp::as<vector<double> >(eqSS["het_brackets"]);
  double Iv = Rcpp::as<double>(eqSS["Iv"]);
  Person::s_mean_maternal_immunity = Rcpp::as<double>(eqSS["MaternalImmunity"]);
  parameters.g_N = Rcpp::as<unsigned int>(paramList["N"]);
  parameters.g_years = Rcpp::as<int>(paramList["years"]);

  
  Rcpp::Rcout << "Matrix unpacking working!\n";
  // could do this //
  // int nr = Smat.nrow(), nc = Smat.ncol() ;
  // std::vector< std::vector<double> > vec( nc ) ;
  // for( int i=0; i<nc; i++){
  //   NumericMatrix::Column col = Smat(_,i) ;
  //   vec[i].assign( col.begin() , col.end() ) ;
  // }
  // Rcpp::Rcout << Smat(1,1) << "\n";
  // Rcpp::Rcout << Smat(1,4) << "\n";
  // Rcpp::Rcout << Smat(2,3) << "\n";
  // Rcpp::Rcout << Smat(3,2) << "\n";
  // Rcpp::Rcout << Smat(4,2) << "\n";
  // Rcpp::Rcout << Smat(21,2) << "\n";
  // Rcpp::Rcout << "Number of rows = " << Smat.nrow() <<"\n";
  // Rcpp::Rcout << "Number of cols = " << Smat.ncol() <<"\n";
  // Rcpp::Rcout << vec[0][1] << " " << vec[4][1] << "\n";
  
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: INITIALISATION FROM EQUILIBRIUM
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Initialise static vectors for keeping biting related variables
  static std::vector<double> s_psi_vector(parameters.g_N);
  static std::vector<double> s_zeta_vector(parameters.g_N);
  static std::vector<double> s_pi_vector(parameters.g_N);
  static std::vector<double> s_cum_pi_vector(parameters.g_N);
  static double s_mean_psi;
  static std::queue<int> s_bite_storage{};
  
  // Initialise humans
  std::vector<Person> Population{ parameters.g_N };
  
  // Initialise strains and mosquitos
  // TODO: Same as above for initialising humans
  
  // Use read in R equilibrium state to then allocate each individual accordingly given tehir age and biting heterogeneity
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // std::vector<double> status_eq{ 0.112426, 0.00640168, 0.7717, 0.0884877, 0.00426778, 0.0169667 };
  // std::vector<double> ID_Eq{ 0.005253793 ,41.30406 };
  // std::vector<double> ICA_Eq{ 0.01841544 ,39.89687 };
  // std::vector<double> IB_Eq{ 0.01059006 ,80.31266 };
  
  int num_age_brackets = age_brackets.size();
  int num_het_brackets = het_brackets.size();
  int age_bracket_in = 0;
  int het_bracket_in = 0;
  std::vector<double> infection_state_probability(6);
  //std::vector<std::vector<std::vector<double> > > infection_state_probability (6,vector<vector<double> >(num_age_brackets,vector <double>(3)));
  Rcpp::Rcout << "Pre human-initialisation working!\n";
  
  for (auto &element : Population) {
    
    // first find out what age bracket they are in
    // are they in the last - do this first for those people with ages that happen to be above the max age bracket
    if(element.get_m_person_age() >= (age_brackets[num_age_brackets-1])){
      age_bracket_in = num_age_brackets;
    } 
    // if not then loop up to last-1
    for(int age_i = 0 ; age_i < (num_age_brackets-1) ; age_i++)
    {
      if(element.get_m_person_age() >= age_brackets[age_i] && element.get_m_person_age() < age_brackets[age_i+1]){
        age_bracket_in = age_i;
      }
    }
    
    // second find out what heterogeneity bracket they are in
    // are they in the last - do this first for those people with heterogeneities that happen to be above the max 
    if(element.m_individual_biting_rate >= het_brackets[num_het_brackets-1]) {
      het_bracket_in = num_het_brackets;
    }
    if(element.m_individual_biting_rate < het_brackets[1]){
      het_bracket_in = 0;
    }
    for(int het_i = 1 ; het_i < (num_het_brackets-1) ; het_i++)
    {
      if(element.m_individual_biting_rate < het_brackets[het_i+1]){
        het_bracket_in = het_i;
      }
    }
    
  
    infection_state_probability[0] = Smat(age_bracket_in, het_bracket_in);
    infection_state_probability[1] = Dmat(age_bracket_in, het_bracket_in);
    infection_state_probability[2] = Amat(age_bracket_in, het_bracket_in);
    infection_state_probability[3] = Umat(age_bracket_in, het_bracket_in);
    infection_state_probability[4] = Tmat(age_bracket_in, het_bracket_in);
    infection_state_probability[5] = Pmat(age_bracket_in, het_bracket_in);
    
    // Sample their infection state given the probabilities of being in any state given their age and het compartment
    element.set_m_infection_state(static_cast<Person::InfectionStatus>(sample1(infection_state_probability, std::accumulate(infection_state_probability.begin(), infection_state_probability.end(),0))));
    element.set_m_IB(IBmat(age_bracket_in, het_bracket_in));
    element.set_m_ICA(ICAmat(age_bracket_in, het_bracket_in));
    element.set_m_ICM(ICMmat(age_bracket_in, het_bracket_in));
    element.set_m_ID(IDmat(age_bracket_in, het_bracket_in));
    
    
    // Using pseudo initialisation from above then assign their infection state and immunity
    // element.set_m_infection_state(static_cast<Person::InfectionStatus>(sample1(status_eq, 1)));
    // element.set_m_ID(1.4 * (element.m_individual_biting_rate * ID_Eq[0] * element.get_m_person_age()) + ID_Eq[1]);
    // element.set_m_ICA(1.3 * (element.m_individual_biting_rate * ICA_Eq[0] * element.get_m_person_age()) + ICA_Eq[1]);
    // element.set_m_IB(0.9 * (element.m_individual_biting_rate * IB_Eq[0] * element.get_m_person_age()) + IB_Eq[1]);
    // element.set_m_ICM(parameters.g_PM * 70 * exp(-element.get_m_person_age() / parameters.g_dCM)); // 70 known mean mother ICA
    // 
    
    // Schedule change for those who are not susceptible
    if (element.get_m_infection_state() != Person::SUSCEPTIBLE) {
      element.schedule_m_day_of_InfectionStatus_change(parameters);
    }
    
    // Check if mother and if so increase the maternal immunity sum and total number
    if (element.get_m_person_age() > 20 * 365 && element.get_m_person_age() < 21 * 365) {
      Person::s_sum_maternal_immunity += element.get_m_ICA();
      Person::s_total_mums++;
    }
    
    // If they are infected, i.e. not S or P, then assign their strains and next strain clearance date
    if (element.get_m_infection_state() != Person::SUSCEPTIBLE && element.get_m_infection_state() != Person::PROPHYLAXIS) {
      // TODO: Think about how we can correctly initialise MOI for a given EIR. Presumably there is a rarefaction of MOI vs EIR, and the MOI is lognormal*age_dependency
      element.set_m_number_of_strains(runiform_int_1(1, 10));
      element.schedule_m_day_of_strain_clearance(parameters);
    }
    
    // Set the next event day
    element.set_m_day_of_next_event();
    
    // Add to the static zeta vector which is required for calculating the overall probability of being bitten, pi
    s_zeta_vector[element.m_person_ID] = element.m_individual_biting_rate;
  }
  
  Rcpp::Rcout << "Human initilisation working\n";
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: INITIALISATION FROM EQUILIBRIUM
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: TIMERS AND PRE LOOP
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
  std::cout << "Time elapsed in initialisation: " << duration << " seconds" << std::endl;
  
  // Start from day 2
  parameters.g_current_time++;
  
  // End of simulation time
  int g_end_time = parameters.g_years * 365;
  
  // Preallocations;
  int num_bites = 0;
  int increasing_bites = 0;
  int individual_binomial_bite_draw = 0;
  double pi_cum_sum = 0;
  double pi_sum = 0;
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START SIMULATION LOOP
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  for (parameters.g_current_time ; parameters.g_current_time < g_end_time ; parameters.g_current_time++)
  {
    
    // Counter print
    if (parameters.g_current_time % 100 == 0) { 
      std::cout << parameters.g_current_time << " days" << std::endl; 
    }
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // DAILY UPDATING AND EVENT HANDLING
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Calculate yesterday's mean maternal immunity
    Person::s_mean_maternal_immunity = Person::s_sum_maternal_immunity / Person::s_total_mums;
    
    // Clear static counters for maternal immunity which will be reset during the update loop
    Person::s_sum_maternal_immunity = 0;
    Person::s_total_mums = 0;
    
    // Clear static sum of psi
    Person::s_psi_sum = 0;
    
    // Loop through each person and mosquito and update
    // PARALLEL_TODO: This loop could easily be parallelised as each person will not require any shared memory (except for parameters)
    for (unsigned int n = 0; n < parameters.g_N; n++) {
      s_psi_vector[n] = Population[n].update(parameters);
    }
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // BITE HANDLING
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // First calculate mean age dependent biting heterogeneity (psi)
    s_mean_psi = Person::s_psi_sum / parameters.g_N;
    
    // Create normalised psi by dividing by the mean age dependent biting rate
    std::transform(s_psi_vector.begin(), s_psi_vector.end(), s_psi_vector.begin(),
                   std::bind1st(std::multiplies<double>(), 1 / s_mean_psi));
    
    // Create overall relative biting rate, pi, i.e. the product of individual biting heterogeneity and age dependent heterogeneity
    std::transform(s_psi_vector.begin(), s_psi_vector.end(),
                   s_zeta_vector.begin(), s_pi_vector.begin(),
                   std::multiplies<double>());
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // BITE ALLOCATIONS
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Multinomial step //
    /////////////////////////////////////////////////////////////////////////
    
    pi_sum = std::accumulate(s_pi_vector.begin(), s_pi_vector.end(), 0.0);
    num_bites = rpoisson1(pi_sum * Iv * 0.30677);
    
    // multinomial procedure for everyone and breaking when n_bites is reached
    for (unsigned int n = 0; n < parameters.g_N - 1; n++) {
      
      individual_binomial_bite_draw = rbinomial1(num_bites - increasing_bites, s_pi_vector[n] / (pi_sum - pi_cum_sum));
      
      for (int element = 0; element < individual_binomial_bite_draw; element++) {
        s_bite_storage.push(n);
      }
      
      pi_cum_sum += s_pi_vector[n];
      increasing_bites += individual_binomial_bite_draw;
      if (increasing_bites >= num_bites) break;
      
    }
    
    // catch rounding errors so just place this here outside loop
    if (increasing_bites != num_bites) {
      individual_binomial_bite_draw = num_bites - increasing_bites;
      for (int element = 0; element < individual_binomial_bite_draw; element++) {
        s_bite_storage.push(parameters.g_N - 1);
      }
    }
    
    // Reset bite and sum of biting rates
    increasing_bites = 0;
    pi_cum_sum = 0;
    
    
    /*
     // OLD BITE SAMPLING// 		
     // Create cumulative pi
     std::partial_sum(s_pi_vector.begin(), s_pi_vector.end(), s_cum_pi_vector.begin());
     
     // Draw total number of bites from adult mosquito population
     // final numbers are simply what EIR of 120 would look like for Iv and 3 day biting
     int num_bites = rpoisson1(s_cum_pi_vector.back() * 0.808 * 0.30677);
     
     // Loop through bites TODO: FIX THIS FOR SPEED
     // PARALLEL_TODO: This could also be very eaily parallelised, except need to assess the storage for bites into a thread safe object
     for (int bite = 0; bite < num_bites; bite++) {
     
     // Locate bites
     s_bite_storage.push(int(std::lower_bound(s_cum_pi_vector.begin(), s_cum_pi_vector.end(), runif1(0.0, s_cum_pi_vector.back())) - s_cum_pi_vector.begin()));
     
     // THE ABOVE DOES WHAT IS COMMENTED BELOW ESSENTIALLY
     // double test = runif1(0.0, s_cum_pi_vector.back());
     // auto it = std::lower_bound(s_cum_pi_vector.begin(), s_cum_pi_vector.end(), test);
     // s_bite_storage.push(std::lower_bound(s_cum_pi_vector.begin(), s_cum_pi_vector.end(), test));
     // int location(it - s_cum_pi_vector.begin());	
     // s_bite_storage.push(location);
     }
     */
    
    // Loop through each bite and allocate accordingly
    // PARALLEL_TODO: Don't know how this could be parallelised yet - come back to with mosquitos in.
    for (int n = 0; n < num_bites; n++) {
      Population[s_bite_storage.front()].allocate_bite(parameters);
      s_bite_storage.pop();
    }
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // LOGGERS
    // TODO: Looping idea
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
  };
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END TIMERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  t1 = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
  std::cout << "Time elapsed total: " << duration << " seconds" << std::endl;
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // SUMMARY LOGGING
  // TODO: As above for including Rcpp call out here
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Final infection states
  // TODO: This needs to be moved into R side when state comes in...
  std::vector<double> status_eq(6);
  std::vector<int> Infection_States(parameters.g_N);
  std::vector<double> Ages(parameters.g_N);
  int total_incidence = 0;
  int total_incidence_05 = 0;
  int daily_incidence_return = 0;
  
  
  for (int element = 0; element < parameters.g_N ; element++) {
    // Match infection state and schedule associated next state change
    switch (Population[element].get_m_infection_state())
    {
    case Person::SUSCEPTIBLE:
      status_eq[0]++; 
      Infection_States[element] = 0;
      break;
    case Person::DISEASED:
      status_eq[1]++;
      Infection_States[element] = 1;
      break;
    case Person::ASYMPTOMATIC:
      status_eq[2]++;
      Infection_States[element] = 2;
      break;
    case Person::SUBPATENT:
      status_eq[3]++;
      Infection_States[element] = 3;
      break;
    case Person::TREATED:
      status_eq[4]++;
      Infection_States[element] = 4;
      break;
    case Person::PROPHYLAXIS:
      status_eq[5]++;
      Infection_States[element] = 5;
      break;
    default:
      assert(NULL && "Infection state not recognised");
    break;
    }
    
    // log incidence
    daily_incidence_return = Population[element].log_daily_incidence();
    if(daily_incidence_return == 2){
      total_incidence++;
      total_incidence_05++;
    } else if (daily_incidence_return == 1) {
      total_incidence++;
    }
    
    // Ages
    Ages[element] = Population[element].get_m_person_age();

  }
  
  // divide by population size
  for (int element = 0; element < 6; element++) {
    status_eq[element] /= parameters.g_N;
    std::cout << status_eq[element] << " | " << std::endl;
  }
  
  // output Rcpp list
  return Rcpp::List::create(Rcpp::Named("S")=status_eq[0], 
                            Rcpp::Named("D")=status_eq[1],
                            Rcpp::Named("A")=status_eq[2],
                            Rcpp::Named("U")=status_eq[3],
                            Rcpp::Named("T")=status_eq[4],
                            Rcpp::Named("P")=status_eq[5],
                            Rcpp::Named("Incidence")=total_incidence,
                            Rcpp::Named("Incidence_05")=total_incidence_05,
                            Rcpp::Named("InfectionStates")=Infection_States,
                            Rcpp::Named("Ages")=Ages);
  
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // fini
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  return 0;
}


