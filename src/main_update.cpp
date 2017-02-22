//
//  MAGENTA
//  main_update.cpp
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

// Create universe structure for all important variables
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


// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// START: MAIN
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List Simulation_Update_cpp(Rcpp::List paramList)
{
  
  // prove that C++ code is being run
  Rcpp::Rcout << "Rcpp function is working!\n";
  
  // start timer
  chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Create universe pointer from paramList statePtr
  Rcpp::XPtr<Universe> universe_ptr = Rcpp::as<Rcpp::XPtr<Universe> > (paramList["statePtr"]);
  // Initialise all the universal variables from the statePtr provided
  //std::vector<Person> Population = universe_ptr->Population;
  //double Iv = universe_ptr->Iv;
   extern Parameters parameters;
   parameters = universe_ptr->parameters;
  parameters.g_years = Rcpp::as<double>(paramList["years"]);
  // std::vector<double> pi_vector = universe_ptr->pi_vector;
  // std::vector<double> psi_vector = universe_ptr->psi_vector;
  // std::vector<double> zeta_vector = universe_ptr->zeta_vector;
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  Rcpp::Rcout << "Pointer unpacking working!\n";
  double Scount = 0;
  for(auto &element : universe_ptr->Population){
    if(element.get_m_infection_state()==Person::SUSCEPTIBLE){
      Scount++;
    }
  }
  Rcpp::Rcout << "Starting susceptible population: " << Scount/parameters.g_N << "\n";
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: TIMERS AND PRE LOOP
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
  Rcpp:Rcout << "Time elapsed in initialisation: " << duration << " seconds\n";
  
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
      psi_sum += universe_ptr->psi_vector[n] = universe_ptr->Population[n].update(parameters);
    }
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // BITE HANDLING
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // First calculate mean age dependent biting heterogeneity (psi)
    // --------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Calculate mean age dependent biting rate
    mean_psi = psi_sum / parameters.g_N;
    
    
    
    // Create normalised psi by dividing by the mean age dependent biting rate
    std::transform(universe_ptr->psi_vector.begin(), universe_ptr->psi_vector.end(), universe_ptr->psi_vector.begin(),
                   std::bind1st(std::multiplies<double>(), 1 / mean_psi));
    
    // Create overall relative biting rate, pi, i.e. the product of individual biting heterogeneity and age dependent heterogeneity
    std::transform(universe_ptr->psi_vector.begin(), universe_ptr->psi_vector.end(),
                   universe_ptr->zeta_vector.begin(), universe_ptr->pi_vector.begin(),
                   std::multiplies<double>());
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // BITE ALLOCATIONS
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Multinomial step //
    /////////////////////////////////////////////////////////////////////////
    
    pi_sum = std::accumulate(universe_ptr->pi_vector.begin(), universe_ptr->pi_vector.end(), 0.0);
    num_bites = rpoisson1(pi_sum * universe_ptr->Iv * 0.30677);
    
    
    // multinomial procedure for everyone and breaking when n_bites is reached
    for (unsigned int n = 0; n < parameters.g_N - 1; n++)
    {
      
      individual_binomial_bite_draw = rbinomial1(num_bites - increasing_bites, universe_ptr->pi_vector[n] / (pi_sum - pi_cum_sum));
      
      for (int element = 0; element < individual_binomial_bite_draw; element++)
      {
        bite_storage_queue.push(n);
      }
      
      pi_cum_sum += universe_ptr->pi_vector[n];
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
    
    // PARALLEL_TODO: Don't know how this could be parallelised yet - come back to with mosquitos in.
    for (int n = 0; n < num_bites; n++)
    {
      universe_ptr->Population[bite_storage_queue.front()].allocate_bite(parameters);
      bite_storage_queue.pop();
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
  Rcpp::Rcout << "Time elapsed total: " << duration << " seconds\n";
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // SUMMARY LOGGING
  // TODO: As above for including Rcpp call out here
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Final infection states
  // TODO: This needs to be moved into R side when state comes in...
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
    switch (universe_ptr->Population[element].get_m_infection_state())
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
    daily_incidence_return = universe_ptr->Population[element].log_daily_incidence();
    if(daily_incidence_return == 2)
    {
      total_incidence++;
      total_incidence_05++;
    } 
    if (daily_incidence_return == 1) 
    {
      total_incidence++;
    }
    
    // Ages
    Ages[element] = universe_ptr->Population[element].get_m_person_age();
    IB[element] = universe_ptr->Population[element].get_m_IB();
    ICA[element] = universe_ptr->Population[element].get_m_ICA();
    ICM[element] = universe_ptr->Population[element].get_m_ICM();
    ID[element] = universe_ptr->Population[element].get_m_ID();
    
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
  
  // Update the universe parameters
  // TODO: Figure out why I can't do parameters the same way as I do the population
  universe_ptr->parameters = parameters;
  
  // Return Named List with pointer and loggers
  return Rcpp::List::create(Rcpp::Named("Ptr") = universe_ptr, Rcpp::Named("Loggers")=Loggers);
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // fini
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
}


