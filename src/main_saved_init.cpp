//
//  MAGENTA
//  main_saved_init.cpp
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
struct Universe {
  // Human storage
  std::vector<Person> population;
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
Rcpp::List Simulation_Saved_Init_cpp(Rcpp::List paramList)
{
  
  // prove that C++ code is being run
  Rcpp::Rcout << "Rcpp function is working!\n";
  
  // start timer
  chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
  
  // Initialise parameters
  Parameters parameters;
  parameters.g_years = Rcpp::as<double>(paramList["years"]);
  
  // Unpack R List to Rcpp Lists
  Rcpp::List savedState = paramList["savedState"];
  Rcpp::List population_List = savedState["population_List"];
  Rcpp::List parameters_List = savedState["parameters_List"];
  double Iv = Rcpp::as<double>(savedState["Iv"]);
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: R -> C++ CONVERSIONS: parameters
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Update this first as needed for initialising the population
  parameters.g_N = Rcpp::as<unsigned int>(parameters_List["g_N"]);
  parameters.g_current_time = Rcpp::as<int>(parameters_List["g_current_time"]);
  parameters.g_mean_maternal_immunity = Rcpp::as<double>(parameters_List["g_mean_maternal_immunity"]);
  parameters.g_sum_maternal_immunity = Rcpp::as<double>(parameters_List["g_sum_maternal_immunity"]);
  parameters.g_total_mums = Rcpp::as<int>(parameters_List["g_total_mums"]);
  parameters.g_years = Rcpp::as<double>(paramList["years"]);
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS: parameters
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: INITIALISATION
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Add the human population
  std::vector<Person> population;
  population.reserve(parameters.g_N);
  
  // Initialise vectors for keeping biting related variables
  std::vector<double> psi_vector(parameters.g_N);
  std::vector<double> zeta_vector(parameters.g_N);
  std::vector<double> pi_vector(parameters.g_N);
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: R -> C++ CONVERSIONS: population
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  std::vector<double> Zetas = Rcpp::as<vector<double> >(population_List["Zetas"]);
  std::vector<int> Infection_States_v = Rcpp::as<vector<int> >(population_List["Infection_States"]);
  std::vector<int> Ages_v = Rcpp::as<vector<int> >(population_List["Ages"]);
  std::vector<double> IB_v = Rcpp::as<vector<double> >(population_List["IB"]);
  std::vector<double> ICA_v = Rcpp::as<vector<double> >(population_List["ICA"]);
  std::vector<double> ICM_v = Rcpp::as<vector<double> >(population_List["ICM"]);
  std::vector<double> ID_v = Rcpp::as<vector<double> >(population_List["ID"]);
  std::vector<int> IB_last_boost_time = Rcpp::as<vector<int> >(population_List["IB_last_boost_time"]);
  std::vector<int> ICA_last_boost_time = Rcpp::as<vector<int> >(population_List["ICA_last_boost_time"]);
  std::vector<int> ID_last_boost_time = Rcpp::as<vector<int> >(population_List["ID_last_boost_time"]);
  std::vector<int> IB_last_calculated_time = Rcpp::as<vector<int> >(population_List["IB_last_calculated_time"]);
  std::vector<int> I_C_D_CM_last_calculated_time = Rcpp::as<vector<int> >(population_List["I_C_D_CM_last_calculated_time"]);
  std::vector<double> Immunity_boost_float = Rcpp::as<vector<double> >(population_List["Immunity_boost_float"]);
  std::vector<int> Day_of_InfectionStatus_change = Rcpp::as<vector<int> >(population_List["Day_of_InfectionStatus_change"]);
  std::vector<int> Day_of_strain_clearance = Rcpp::as<vector<int> >(population_List["Day_of_strain_clearance"]);
  std::vector<int> Day_of_death = Rcpp::as<vector<int> >(population_List["Day_of_death"]);
  std::vector<int> Number_of_Strains = Rcpp::as<vector<int> >(population_List["Number_of_Strains"]);
  
  Rcpp::List Infection_time_realisation_queues_List = population_List["Infection_time_realisation_queues"];
  Rcpp::List Infection_state_realisation_queues_List = population_List["Infection_state_realisation_queues"];
  
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS: population
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  Rcpp::Rcout << "Matrix unpacking working!\n";
  
  // Use read in R equilibrium state to then allocate each individual accordingly given tehir age and biting heterogeneity
  
  int id_counter = 0;
  Rcpp::Rcout << "Pre human-initialisation working!\n";
  
  for (int n=0; n < parameters.g_N; n++) 
  {
    
    // Set their id
    population.emplace_back(Person(parameters));
    population[n].set_m_person_ID(id_counter++);
    
    // Set infection state
    population[n].set_m_infection_state(static_cast<Person::InfectionStatus>(Infection_States_v[n]));
    
    // Set Age
    population[n].set_m_person_age(Ages_v[n]);
    
    // Set Immunities
    population[n].set_m_IB(IB_v[n]);
    population[n].set_m_ICA(ICA_v[n]);
    population[n].set_m_ICM(ICM_v[n]);
    population[n].set_m_ID(ID_v[n]);
    
    // Set Boost times
    population[n].set_m_IB_last_boost_time(IB_last_boost_time[n]);
    population[n].set_m_ICA_last_boost_time(ICA_last_boost_time[n]);
    population[n].set_m_ID_last_boost_time(ID_last_boost_time[n]);
    
    // Set Calc times
    population[n].set_m_IB_last_calculated_time(IB_last_calculated_time[n]);
    population[n].set_m_I_C_D_CM_last_calculated_time(I_C_D_CM_last_calculated_time[n]);
    
    // Set immunity boost float
    population[n].set_m_immunity_boost_float(Immunity_boost_float[n]);
    
    // Set Day Changes
    population[n].set_m_day_of_InfectionStatus_change(Day_of_InfectionStatus_change[n]);
    population[n].set_m_day_of_strain_clearance(Day_of_strain_clearance[n]);
    population[n].set_m_day_of_death(Day_of_death[n]);
    
    // Strain Numbers
    population[n].set_m_number_of_strains(Number_of_Strains[n]);
    
    // Set queues
    population[n].set_m_infection_time_realisation_queue_from_vector(Rcpp::as<std::vector<int> >(Infection_time_realisation_queues_List[n]));
    population[n].set_m_infection_state_realisation_queue_from_vector(Rcpp::as<std::vector<int> >(Infection_state_realisation_queues_List[n]));
    
    // Set the next event day
    population[n].set_m_day_of_next_event();
    
    // Asign zeta and add to the zeta vector which is required for calculating the overall probability of being bitten, pi
    population[n].set_m_individual_biting_rate(Zetas[n]);
    zeta_vector[n] = population[n].get_m_individual_biting_rate();
    
  }
  
  Rcpp::Rcout << "Human initilisation working\n";
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: INITIALISATION
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  Rcpp::Rcout << "Pointer unpacking working!\n";
  double Scount = 0;
  for(auto &element : population){
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
  Rcpp::Rcout << "Time elapsed in initialisation: " << duration << " seconds\n";
  
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
      psi_sum += psi_vector[n] = population[n].update(parameters);
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
      population[bite_storage_queue.front()].allocate_bite(parameters);
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
    switch (population[element].get_m_infection_state())
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
    daily_incidence_return = population[element].log_daily_incidence(parameters);
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
    Ages[element] = population[element].get_m_person_age();
    IB[element] = population[element].get_m_IB();
    ICA[element] = population[element].get_m_ICA();
    ICM[element] = population[element].get_m_ICM();
    ID[element] = population[element].get_m_ID();
    
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
  Rcpp::XPtr<Universe> universe_ptr(new Universe{ population, psi_vector, zeta_vector, pi_vector, Iv, parameters},
                                    true);
  
  // Return Named List with pointer and loggers
  return Rcpp::List::create(Rcpp::Named("Ptr") = universe_ptr, Rcpp::Named("Loggers")=Loggers);
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // fini
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
}


