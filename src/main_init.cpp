//
//  MAGENTA
//  main_init.cpp
//
//  Created: OJ Watson on 06/12/2015
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

// Create universe structure for all important variables
// TODO: Move this perhaps to own file...
struct Universe {
  // Human storage
  std::vector<Person> Population;
  std::vector<double> psi_vector;
  std::vector<double> zeta_vector;
  std::vector<double> pi_vector;
  // Mosquito storage
  double Iv;
  // Parameter storage
  Parameters parameters;
};

// Static initialisation - only used for ID generation, mostly irrelevant apart from with for_each statements as this
// then allows the vector position tpo be determined
int Person::s_person_ID_generator = 0;
int Strain::s_strain_ID_generator = 0;

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// START: MAIN
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List Simulation_Init_cpp(Rcpp::List paramList)
{
  
  // prove that C++ code is being run
  Rcpp::Rcout << "Rcpp function is working!\n";
  
  // start timer
  chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
  
  // Initialise file scope parameters object
  // This will need to be declared at the beginning of each header/cpp file where parameters are needed
  extern Parameters parameters;
  
  // Update these first so that they are used in the universe initialisation
  parameters.g_N = Rcpp::as<unsigned int>(paramList["N"]);
  parameters.g_years = Rcpp::as<double>(paramList["years"]);
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: INITIALISATION
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Add the human population
  std::vector<Person> Population{ parameters.g_N };
  
  // Initialise vectors for keeping biting related variables
  std::vector<double> psi_vector(parameters.g_N);
  std::vector<double> zeta_vector(parameters.g_N);
  std::vector<double> pi_vector(parameters.g_N);
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: INITIALISATION
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
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
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  Rcpp::Rcout << "Matrix unpacking working!\n";

  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: INITIALISATION FROM EQUILIBRIUM
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Use read in R equilibrium state to then allocate each individual accordingly given tehir age and biting heterogeneity
  
  int num_age_brackets = age_brackets.size();
  int num_het_brackets = het_brackets.size();
  int age_bracket_in = 0;
  int het_bracket_in = 0;
  std::vector<double> infection_state_probability(6);
  Rcpp::Rcout << "Pre human-initialisation working!\n";
  
  for (auto &element : Population) 
  {
    
    // first find out what age bracket they are in
    // are they in the last - do this first for those people with ages that happen to be above the max age bracket
    if(element.get_m_person_age() >= (age_brackets[num_age_brackets-1]))
    {
      age_bracket_in = num_age_brackets;
    } 
    // if not then loop up to last-1
    for(int age_i = 0 ; age_i < (num_age_brackets-1) ; age_i++)
    {
      if(element.get_m_person_age() >= age_brackets[age_i] && element.get_m_person_age() < age_brackets[age_i+1])
      {
        age_bracket_in = age_i;
      }
    }
    
    // second find out what heterogeneity bracket they are in
    // are they in the last - do this first for those people with heterogeneities that happen to be above the max 
    if(element.m_individual_biting_rate >= het_brackets[num_het_brackets-1]) 
    {
      het_bracket_in = num_het_brackets;
    }
    if(element.m_individual_biting_rate < het_brackets[1])
    {
      het_bracket_in = 0;
    }
    for(int het_i = 1 ; het_i < (num_het_brackets-1) ; het_i++)
    {
      if(element.m_individual_biting_rate < het_brackets[het_i+1])
      {
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
    
    // Schedule change for those who are not susceptible
    if (element.get_m_infection_state() != Person::SUSCEPTIBLE)
    {
      element.schedule_m_day_of_InfectionStatus_change(parameters);
    }
    
    // Check if mother and if so increase the maternal immunity sum and total number
    if (element.get_m_person_age() > 20 * 365 && element.get_m_person_age() < 21 * 365)
    {
      parameters.g_sum_maternal_immunity += element.get_m_ICA();
      parameters.g_total_mums++;
    }
    
    // If they are infected, i.e. not S or P, then assign their strains and next strain clearance date
    if (element.get_m_infection_state() != Person::SUSCEPTIBLE && element.get_m_infection_state() != Person::PROPHYLAXIS) 
    {
      // TODO: Think about how we can correctly initialise MOI for a given EIR. Presumably there is a rarefaction of MOI vs EIR, and the MOI is lognormal*age_dependency
      element.set_m_number_of_strains(runiform_int_1(1, 10));
      element.schedule_m_day_of_strain_clearance(parameters);
    }
    
    // Set the next event day
    element.set_m_day_of_next_event();
    
    // Add to the static zeta vector which is required for calculating the overall probability of being bitten, pi
    zeta_vector[element.m_person_ID] = element.m_individual_biting_rate;
  }
  
  Rcpp::Rcout << "Human initilisation working\n";
  
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: TIMERS AND PRE LOOP
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
  Rcpp::Rcout << "Time elapsed in initialisation: " << duration << " seconds\n";
  
  // Start from day 2
  parameters.g_current_time++;
  
  // End of simulation time
  int g_end_time = parameters.g_current_time + (parameters.g_years * 365);
  
  // Preallocations;
  int num_bites = 0;
  int increasing_bites = 0;
  int individual_binomial_bite_draw = 0;
  double psi_sum = 0;
  // Maternal 
  double mean_psi = 0;
  double pi_cum_sum = 0;
  double pi_sum = 0;
  // Bites
  std::queue<int> bite_storage_queue{};
  
  
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START SIMULATION LOOP
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  for ( ; parameters.g_current_time < g_end_time ; parameters.g_current_time++)
  {
    
    // Counter print
    if (parameters.g_current_time % 100 == 0) 
    { 
      Rcpp::Rcout << parameters.g_current_time << " days" << "\n"; 
    }
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // DAILY UPDATING AND EVENT HANDLING
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // First calculate the mean maternal immunity from yesterday
    // --------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Calculate yesterday's mean maternal immunity
    parameters.g_mean_maternal_immunity = parameters.g_sum_maternal_immunity / parameters.g_total_mums;
    
    // Reset maternal immunity sums
    parameters.g_sum_maternal_immunity = 0;
    parameters.g_total_mums = 0;
    
    // Reset age dependent biting rate sum
    psi_sum = 0;
    
    // Loop through each person and mosquito and update
    // --------------------------------------------------------------------------------------------------------------------------------------------------
    // PARALLEL_TODO: This loop could easily be parallelised as each person will not require any shared memory (except for parameters)
    
    for (unsigned int n = 0; n < parameters.g_N; n++) 
    {
      psi_sum += psi_vector[n] = Population[n].update(parameters);
    }
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // BITE HANDLING
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // First calculate mean age dependent biting heterogeneity (psi)
    // --------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Calculate mean age dependent biting rate
    mean_psi = psi_sum / parameters.g_N;
    
    // Create normalised psi by dividing by the mean age dependent biting rate
    std::transform(psi_vector.begin(), psi_vector.end(), psi_vector.begin(),
                   std::bind1st(std::multiplies<double>(), 1 / mean_psi));
    
    // Create overall relative biting rate, pi, i.e. the product of individual biting heterogeneity and age dependent heterogeneity
    std::transform(psi_vector.begin(), psi_vector.end(),
                   zeta_vector.begin(), pi_vector.begin(),
                   std::multiplies<double>());
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // BITE ALLOCATIONS
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Multinomial step //
    // --------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Caluclate total probability of being bitten within population
    pi_sum = std::accumulate(pi_vector.begin(), pi_vector.end(), 0.0);
    num_bites = rpoisson1(pi_sum * Iv * 0.30677);
    
    // multinomial procedure for everyone and breaking when n_bites is reached
    for (unsigned int n = 0; n < parameters.g_N - 1; n++)
    {
      
      individual_binomial_bite_draw = rbinomial1(num_bites - increasing_bites, pi_vector[n] / (pi_sum - pi_cum_sum));
      
      for (int element = 0; element < individual_binomial_bite_draw; element++)
      {
        bite_storage_queue.push(n);
      }
      
      pi_cum_sum += pi_vector[n];
      increasing_bites += individual_binomial_bite_draw;
      if (increasing_bites >= num_bites) break;
      
    }
    
    // catch rounding errors so just place this here outside loop
    if (increasing_bites != num_bites) 
    {
      individual_binomial_bite_draw = num_bites - increasing_bites;
      for (int element = 0; element < individual_binomial_bite_draw; element++) 
      {
        bite_storage_queue.push(parameters.g_N - 1);
      }
    }
    
    // Reset bite and sum of biting rates
    increasing_bites = 0;
    pi_cum_sum = 0;
    
    // ALLOCATE BITES
    // --------------------------------------------------------------------------------------------------------------------------------------------------
    // PARALLEL_TODO: Don't know how this could be parallelised yet - come back to with mosquitos in.
    
    for (int n = 0; n < num_bites; n++)
    {
      Population[bite_storage_queue.front()].allocate_bite(parameters);
      bite_storage_queue.pop();
    }
    
  };
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END TIMERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  t1 = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
  Rcpp::Rcout << "Time elapsed total: " << duration << " seconds\n";
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // SUMMARY LOGGING
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Final infection states
  std::vector<double> status_eq(6);
  std::vector<int> Infection_States(parameters.g_N);
  std::vector<double> Ages(parameters.g_N);
  std::vector<double> IB(parameters.g_N);
  std::vector<double> ICA(parameters.g_N);
  std::vector<double> ICM(parameters.g_N);
  std::vector<double> ID(parameters.g_N);
  int total_incidence = 0;
  int total_incidence_05 = 0;
  int daily_incidence_return = 0;
  
  for (unsigned int element = 0; element < parameters.g_N ; element++) 
  {
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
    if(daily_incidence_return == 2)
    {
      total_incidence++;
      total_incidence_05++;
    } 
    if (daily_incidence_return == 1) 
    {
      total_incidence++;
    }
    
    // Ages and immunity 
    // TODO: Figure out the best way of standardising this logging 
    // Something like passing in a function name within the paramList which is the 
    // name for a logger written else where which then returns the Loggers obeject below
    Ages[element] = Population[element].get_m_person_age();
    IB[element] = Population[element].get_m_IB();
    ICA[element] = Population[element].get_m_ICA();
    ICM[element] = Population[element].get_m_ICM();
    ID[element] = Population[element].get_m_ID();
    
  }
  
  // divide by population size
  Rcpp::Rcout << "S | D | A | U | T | P:\n" ;
  
  for (int element = 0; element < 6; element++) 
  {
    status_eq[element] /= parameters.g_N;
    Rcpp::Rcout << status_eq[element] << " | ";
  }
  
  // Create Rcpp loggers list
  Rcpp::List Loggers = Rcpp::List::create(Rcpp::Named("S")=status_eq[0],Rcpp::Named("D")=status_eq[1],Rcpp::Named("A")=status_eq[2],
                                          Rcpp::Named("U")=status_eq[3],Rcpp::Named("T")=status_eq[4],Rcpp::Named("P")=status_eq[5],Rcpp::Named("Incidence")=total_incidence,
                                          Rcpp::Named("Incidence_05")=total_incidence_05,Rcpp::Named("InfectionStates")=Infection_States,Rcpp::Named("Ages")=Ages,
                                                      Rcpp::Named("IB")=IB,Rcpp::Named("ICA")=ICA,Rcpp::Named("ICM")=ICM,Rcpp::Named("ID")=ID);
  
  
  // Create universe ptr for memory-continuiation
  Rcpp::XPtr<Universe> universe_ptr(new Universe{ Population, psi_vector, zeta_vector, pi_vector, Iv, parameters},
                                    true);
  
  // Return Named List with pointer and loggers
  return Rcpp::List::create(Rcpp::Named("Ptr") = universe_ptr, Rcpp::Named("Loggers")=Loggers);
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // fini
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
}


