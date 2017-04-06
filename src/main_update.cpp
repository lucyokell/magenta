//
//  MAGENTA
//  main_update.cpp
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
  std::vector<Person> population;
  std::vector<double> psi_vector;
  std::vector<double> zeta_vector;
  std::vector<double> pi_vector;
  // Mosquito storage
  std::vector<Mosquito> scourge;
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
  universe_ptr->parameters.g_years = Rcpp::as<double>(paramList["years"]);
  universe_ptr->parameters.g_ft = Rcpp::as<double>(paramList["ft"]);
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  Rcpp::Rcout << "Pointer unpacking working!\n";
  double Scount = 0;
  for(auto &element : universe_ptr->population){
    if(element.get_m_infection_state()==Person::SUSCEPTIBLE){
      Scount++;
    }
  }
  Rcpp::Rcout << "Starting susceptible population: " << Scount/universe_ptr->parameters.g_N << "\n";
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: TIMERS AND PRE LOOP
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
  Rcpp::Rcout << "Time elapsed in initialisation: " << duration << " seconds\n";
  
  // End of simulation time
  int g_end_time = universe_ptr->parameters.g_current_time + (universe_ptr->parameters.g_years * 365);
  
  // Preallocations;
  unsigned int num_bites = 0;
  unsigned int increasing_bites = 0;
  unsigned int individual_binomial_bite_draw = 0;
  double psi_sum = 0;
  unsigned int scourge_size = universe_ptr->scourge.size();
  
  // status eq for logging and other logging variables
  std::vector<double> status_eq = { 0,0,0,0,0,0 };
  unsigned int log_counter = 0;
  double total_incidence = 0;
  double total_incidence_05 = 0;
  int daily_incidence_return = 0;
  int daily_bite_counters = 0;
  
  // For loop preallocations
  unsigned int human_update_i = 0;
  unsigned int mosquito_update_i = 0;
  unsigned int bite_sampling_i = 0;
  unsigned int bite_sampling_internal_i = 0;
  unsigned int num_bites_i = 0;
  
  // Maternal 
  double mean_psi = 0;
  double pi_cum_sum = 0;
  double pi_sum = 0;
  
  // Bites
  std::vector<int> mosquito_biting_queue;
  mosquito_biting_queue.reserve(scourge_size);
  std::vector<int> bite_storage_queue;
  bite_storage_queue.reserve(scourge_size);
  
  
  // resetart timer
  t0 = std::chrono::high_resolution_clock::now();
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START SIMULATION LOOP
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  for ( ; universe_ptr->parameters.g_current_time < g_end_time ; universe_ptr->parameters.g_current_time++)
  {
    
    // Counter print
    if (universe_ptr->parameters.g_current_time % 100 == 0) 
    { 
      Rcpp::Rcout << universe_ptr->parameters.g_current_time << " days" << "\n"; 
    }
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // DAILY UPDATING AND EVENT HANDLING
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // First calculate the mean maternal immunity from yesterday
    // --------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Calculate yesterday's mean maternal immunity
    universe_ptr->parameters.g_mean_maternal_immunity = universe_ptr->parameters.g_sum_maternal_immunity / universe_ptr->parameters.g_total_mums;
    
    // Reset maternal immunity sums
    universe_ptr->parameters.g_sum_maternal_immunity = 0;
    universe_ptr->parameters.g_total_mums = 0;
    
    // Reset age dependent biting rate sum
    psi_sum = 0;

    // Loop through each person and mosquito and update
    // --------------------------------------------------------------------------------------------------------------------------------------------------
    // PARALLEL_TODO: This loop could easily be parallelised as each person will not require any shared memory (except for universe_ptr->parameters)
    
    // Human update loop
    for (human_update_i = 0; human_update_i < universe_ptr->parameters.g_N; human_update_i++)
    {
      psi_sum += universe_ptr->psi_vector[human_update_i] = universe_ptr->population[human_update_i].update(universe_ptr->parameters);
    }

    // Reset number of bites for each day and update each mosquito. Update returns whether the mosquito is biting today
    num_bites = 0;
    
    // Mosquito update loop
    for (mosquito_update_i = 0; mosquito_update_i < scourge_size; mosquito_update_i++)
    {
      if (universe_ptr->scourge[mosquito_update_i].update(universe_ptr->parameters))
      {
        mosquito_biting_queue.emplace_back(mosquito_update_i);
        num_bites++;
      }
    }
    
    // shuffle the bite queue otherwise you will introduce stepping-stone-esque genetic structuring
    shuffle_integer_vector(mosquito_biting_queue);
    
    // Adjust the number of bites to account for anthrophagy
    num_bites *= universe_ptr->parameters.g_Q0;

    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // BITE HANDLING
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // First calculate mean age dependent biting heterogeneity (psi)
    // --------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Calculate mean age dependent biting rate
    mean_psi = psi_sum / universe_ptr->parameters.g_N;
    
    // Create normalised psi by dividing by the mean age dependent biting rate
    std::transform(universe_ptr->psi_vector.begin(), universe_ptr->psi_vector.end(), universe_ptr->psi_vector.begin(),
                   std::bind1st(std::multiplies<double>(), 1 / mean_psi));
    
    // Create overall relative biting rate, pi, i.e. the product of individual biting heterogeneity and age dependent heterogeneity
    std::transform(universe_ptr->psi_vector.begin(), universe_ptr->psi_vector.end(),
                   universe_ptr->zeta_vector.begin(), universe_ptr->pi_vector.begin(),
                   std::multiplies<double>());
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // START: BITE ALLOCATIONS
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Multinomial step //
    // --------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Caluclate total probability of being bitten within population
    pi_sum = std::accumulate(universe_ptr->pi_vector.begin(), universe_ptr->pi_vector.end(), 0.0);
    
    // multinomial procedure for everyone and breaking when n_bites is reached
    for (bite_sampling_i = 0; bite_sampling_i < universe_ptr->parameters.g_N - 1; bite_sampling_i++)
    {
      
      individual_binomial_bite_draw = rbinomial1(num_bites - increasing_bites, universe_ptr->pi_vector[bite_sampling_i] / (pi_sum - pi_cum_sum));
      
      for (bite_sampling_internal_i = 0; bite_sampling_internal_i < individual_binomial_bite_draw; bite_sampling_internal_i++)
      {
        /*	bite_storage_queue.push(n);*/
        bite_storage_queue.emplace_back(bite_sampling_i);
      }
      
      pi_cum_sum += universe_ptr->pi_vector[bite_sampling_i];
      increasing_bites += individual_binomial_bite_draw;
      if (increasing_bites >= num_bites) break;
      
    }
    
    // catch rounding errors so just place this here outside loop
    if (increasing_bites < num_bites)
    {
      individual_binomial_bite_draw = num_bites - increasing_bites;
      for (bite_sampling_internal_i = 0; bite_sampling_internal_i < individual_binomial_bite_draw; bite_sampling_internal_i++)
      {
        //bite_storage_queue.push(parameters.g_N - 1);
        bite_storage_queue.emplace_back(universe_ptr->parameters.g_N - 1);
      }
    }
    
    // Reset bite and sum of biting rates
    increasing_bites = 0;
    pi_cum_sum = 0;
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // END: BITE ALLOCATION SAMPLING
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // --------------------------------------------------------------------------------------------------------------------------------------------------
    // ALLOCATE BITES
    // --------------------------------------------------------------------------------------------------------------------------------------------------
    
    // PARALLEL_TODO: Don't know how this could be parallelised yet - come back to with mosquitos in.
    for (num_bites_i = 0; num_bites_i < num_bites; num_bites_i++)
    {
      // allocate bite to human if mosquito is infected
      if (universe_ptr->scourge[mosquito_biting_queue[num_bites_i]].m_mosquito_infected) {
        universe_ptr->population[bite_storage_queue[num_bites_i]].allocate_bite(universe_ptr->parameters, universe_ptr->scourge[mosquito_biting_queue[num_bites_i]]);
        if ( universe_ptr->parameters.g_current_time > g_end_time - 7)
        {
          daily_bite_counters++;
        }
      }
      
      // if human would cause infection to mosquito then allocate gametocytes
      if (universe_ptr->population[bite_storage_queue[num_bites_i]].reciprocal_infection_boolean(universe_ptr->parameters)) {
        universe_ptr->scourge[mosquito_biting_queue[num_bites_i]].allocate_gametocytes(universe_ptr->parameters, universe_ptr->population[bite_storage_queue[num_bites_i]].sample_two_barcodes(universe_ptr->parameters));
      }
      
    }

    // clear biting storage vectors
    bite_storage_queue.clear();
    mosquito_biting_queue.clear();
    
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // START: SUMMARY LOGGING
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Log the last week //
    
    if (universe_ptr->parameters.g_current_time > g_end_time - 7)
    {
      
      log_counter++;
      
      for (auto &element : universe_ptr->population)
      {
        // Match infection state and schedule associated next state change
        switch (element.get_m_infection_state())
        {
        case Person::SUSCEPTIBLE:
          status_eq[0]++;
          break;
        case Person::DISEASED:
          status_eq[1]++;
          break;
        case Person::ASYMPTOMATIC:
          status_eq[2]++;
          break;
        case Person::SUBPATENT:
          status_eq[3]++;
          break;
        case Person::TREATED:
          status_eq[4]++;
          break;
        case Person::PROPHYLAXIS:
          status_eq[5]++;
          break;
        default:
          assert(NULL && "Schedule Infection Status Change Error - person's infection status not S, D, A, U, T or P");
        break;
        }
        
        // log incidence
        daily_incidence_return = element.log_daily_incidence(universe_ptr->parameters);
        if(daily_incidence_return == 2)
        {
          total_incidence++;
          total_incidence_05++;
        } 
        if (daily_incidence_return == 1) 
        {
          total_incidence++;
        }
        
      }
    }
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // END: SUMMARY LOGGING
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
  };
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END TIMERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  t1 = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
  Rcpp::Rcout << "Time elapsed total: " << duration << " seconds\n";
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // SUMMARY LOGGING AVERAGING AND VARIABLE RETURN
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // divide by population size and log counter and print to give overview
  Rcpp::Rcout << "S | D | A | U | T | P:\n" ;
  
  for (int element = 0; element < 6; element++) 
  {
    status_eq[element] /= (universe_ptr->parameters.g_N * log_counter);
    Rcpp::Rcout << status_eq[element] << " | ";
  }
  
  
  // Final infection states
  std::vector<int> Infection_States(universe_ptr->parameters.g_N);
  std::vector<double> Ages(universe_ptr->parameters.g_N);
  std::vector<double> IB(universe_ptr->parameters.g_N);
  std::vector<double> ICA(universe_ptr->parameters.g_N);
  std::vector<double> ICM(universe_ptr->parameters.g_N);
  std::vector<double> ID(universe_ptr->parameters.g_N);
  
  // loop through population and grabages and immunities for error checking
  for (unsigned int element = 0; element < universe_ptr->parameters.g_N ; element++) 
  {
    
    // Ages and immunity 
    // TODO: Figure out the best way of standardising this logging 
    // Something like passing in a function name within the paramList which is the 
    // name for a logger written else where which then returns the Loggers obeject below
    Infection_States[element] = static_cast<int>(universe_ptr->population[element].get_m_infection_state());
    universe_ptr->population[element].update_immunities_to_today(universe_ptr->parameters);
    Ages[element] = universe_ptr->population[element].get_m_person_age();
    IB[element] = universe_ptr->population[element].get_m_IB();
    ICA[element] = universe_ptr->population[element].get_m_ICA();
    ICM[element] = universe_ptr->population[element].get_m_ICM();
    ID[element] = universe_ptr->population[element].get_m_ID();
    
  }
  
  // Create Rcpp loggers list
  Rcpp::List Loggers = Rcpp::List::create(Rcpp::Named("S")=status_eq[0],Rcpp::Named("D")=status_eq[1],Rcpp::Named("A")=status_eq[2],
                                          Rcpp::Named("U")=status_eq[3],Rcpp::Named("T")=status_eq[4],Rcpp::Named("P")=status_eq[5],Rcpp::Named("Incidence")=total_incidence/log_counter,
                                          Rcpp::Named("Incidence_05")=total_incidence_05/log_counter,Rcpp::Named("InfectionStates")=Infection_States,Rcpp::Named("Ages")=Ages,
                                                      Rcpp::Named("IB")=IB,Rcpp::Named("ICA")=ICA,Rcpp::Named("ICM")=ICM,Rcpp::Named("ID")=ID,
                                                        Rcpp::Named("Daily_Bites")=daily_bite_counters/log_counter);
  
  
  // Return Named List with pointer and loggers
  return Rcpp::List::create(Rcpp::Named("Ptr") = universe_ptr, Rcpp::Named("Loggers")=Loggers);
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // fini
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
}


