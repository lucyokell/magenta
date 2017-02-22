//
//  MAGENTA
//  main_get.cpp
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
Rcpp::List Simulation_Get_cpp(Rcpp::List paramList)
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
  
  // Initialise variables from the statePtr provided that are needing to be saved
  std::vector<double> Zeta = universe_ptr->zeta_vector;
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  Rcpp::Rcout << "Pointer unpacking working!\n";
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: MODEL STATE GET
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  std::vector<int> Infection_States(universe_ptr->parameters.g_N);
  std::vector<double> Ages(universe_ptr->parameters.g_N);
  std::vector<double> IB(universe_ptr->parameters.g_N);
  std::vector<double> ICA(universe_ptr->parameters.g_N);
  std::vector<double> ICM(universe_ptr->parameters.g_N);
  std::vector<double> ID(universe_ptr->parameters.g_N);
  std::vector<int> IB_last_boost_time(universe_ptr->parameters.g_N);
  std::vector<int> ICA_last_boost_time(universe_ptr->parameters.g_N);
  std::vector<int> ID_last_boost_time(universe_ptr->parameters.g_N);
  std::vector<int> IB_last_calculated_time(universe_ptr->parameters.g_N);
  std::vector<int> I_C_D_CM_last_calculated_time(universe_ptr->parameters.g_N);
  std::vector<int> Immunity_boost_float(universe_ptr->parameters.g_N);
  std::vector<int>  Day_of_InfectionStatus_change(universe_ptr->parameters.g_N);	
  std::vector<int>  Day_of_strain_clearance(universe_ptr->parameters.g_N);		
  std::vector<int>  Day_of_death(universe_ptr->parameters.g_N);					
  
  std::vector<std::vector<int> > Infection_time_realisation_queues(universe_ptr->parameters.g_N);
  std::vector<std::vector<int> > Infection_state_realisation_queues(universe_ptr->parameters.g_N);
  
  // Temporary necessities for unpacking queues
  std::queue<int> temp_infection_time_realisation_queue{};
  std::queue<Person::InfectionStatus> temp_infection_state_realisation_queue{};
  
  std::vector<int> temp_infection_time_realisation_vector{};
  std::vector<int> temp_infection_state_realisation_vector{};
  // TODO: Strains and mosquito saving
  //std::vector<std::queue<int> > Infection_strain_realisation_queues(universe_ptr->parameters.g_N);
  
  for (unsigned int element = 0; element < universe_ptr->parameters.g_N ; element++) 
  {
    
    // Infection States
    Infection_States[element] = static_cast<int>(universe_ptr->Population[element].get_m_infection_state());
    
    // Ages
    Ages[element] = universe_ptr->Population[element].get_m_person_age();
    
    // Immunities
    IB[element] = universe_ptr->Population[element].get_m_IB();
    ICA[element] = universe_ptr->Population[element].get_m_ICA();
    ICM[element] = universe_ptr->Population[element].get_m_ICM();
    ID[element] = universe_ptr->Population[element].get_m_ID();
    
    // Boost times
    IB_last_boost_time[element] = universe_ptr->Population[element].get_m_IB_last_boost_time();
    ICA_last_boost_time[element] = universe_ptr->Population[element].get_m_ICA_last_boost_time();
    ID_last_boost_time[element] = universe_ptr->Population[element].get_m_ID_last_boost_time();
    
    // Calc times
    IB_last_calculated_time[element] = universe_ptr->Population[element].get_m_IB_last_calculated_time();
    I_C_D_CM_last_calculated_time[element] = universe_ptr->Population[element].get_m_I_C_D_CM_last_calculated_time();
    
    // Immunity float
    Immunity_boost_float[element] = universe_ptr->Population[element].get_m_immunity_boost_float();
    
    // Day Changes
    Day_of_InfectionStatus_change[element] = universe_ptr->Population[element].get_m_day_of_InfectionStatus_change();
    Day_of_strain_clearance[element] = universe_ptr->Population[element].get_m_day_of_strain_clearance();
    Day_of_death[element]	 = universe_ptr->Population[element].get_m_day_of_death();
    
    // Queue temp creates
    temp_infection_time_realisation_queue = universe_ptr->Population[element].get_m_infection_time_realisation_queue();
    temp_infection_state_realisation_queue = universe_ptr->Population[element].get_m_infection_state_realisation_queue();
    
    // Push back to the temp vector all the queue elements and then assign before clearing the temp vector
    while(!temp_infection_time_realisation_queue.empty()){
      temp_infection_time_realisation_vector.push_back(temp_infection_time_realisation_queue.front());
      temp_infection_time_realisation_queue.pop();
    }
    Infection_time_realisation_queues[element] = temp_infection_time_realisation_vector;
    temp_infection_time_realisation_vector.clear();
    
    // Push back to the temp vector all the queue elements and then assign before clearing the temp vector
    while(!temp_infection_state_realisation_queue.empty()){
      temp_infection_state_realisation_vector.push_back(temp_infection_state_realisation_queue.front());
      temp_infection_state_realisation_queue.pop();
    }
    Infection_state_realisation_queues[element] = temp_infection_state_realisation_vector;
    temp_infection_state_realisation_vector.clear();
    
  }
  
  // Create Rcpp Population list
  Rcpp::List Population_List = Rcpp::List::create(Rcpp::Named("Infection_States")=Infection_States,
                                                  Rcpp::Named("Zetas")=Zeta,
                                                  Rcpp::Named("Ages")=Ages,
                                                  Rcpp::Named("IB")=IB,
                                                  Rcpp::Named("ICA")=ICA,
                                                  Rcpp::Named("ICM")=ICM,
                                                  Rcpp::Named("ID")=ID,
                                                  Rcpp::Named("IB_last_boost_time")=IB_last_boost_time,
                                                  Rcpp::Named("ICA_last_boost_time")=ICA_last_boost_time,
                                                  Rcpp::Named("ID_last_boost_time")=ID_last_boost_time,
                                                  Rcpp::Named("IB_last_calculated_time")=IB_last_calculated_time,
                                                  Rcpp::Named("I_C_D_CM_last_calculated_time")=I_C_D_CM_last_calculated_time,
                                                  Rcpp::Named("Immunity_boost_float")=Immunity_boost_float,
                                                  Rcpp::Named("Day_of_InfectionStatus_change")=Day_of_InfectionStatus_change,
                                                  Rcpp::Named("Day_of_strain_clearance")=Day_of_strain_clearance,
                                                  Rcpp::Named("Day_of_death")=Day_of_death,
                                                  Rcpp::Named("Infection_time_realisation_queues")=Infection_time_realisation_queues,
                                                  Rcpp::Named("Infection_state_realisation_queues")=Infection_state_realisation_queues
  );
  
  // Create Rcpp Parameters list
  // TODO: Chat to Rich about why this fails on long lists and how to circumvent
  Rcpp::List Parameters_List = Rcpp::List::create(Rcpp::Named("g_current_time")=universe_ptr->parameters.g_current_time,
                                                  Rcpp::Named("g_mean_maternal_immunity")=universe_ptr->parameters.g_mean_maternal_immunity,
                                                  Rcpp::Named("g_sum_maternal_immunity")=universe_ptr->parameters.g_sum_maternal_immunity,
                                                  Rcpp::Named("g_total_mums")=universe_ptr->parameters.g_total_mums
  );
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: MODEL STATE GET
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
  std::cout << "Time elapsed in fetching state: " << duration << " seconds" << std::endl;
  
  // Return Named List with Population and parameters
  return Rcpp::List::create(Rcpp::Named("Population") = Population_List, Rcpp::Named("Parameters")=Parameters_List);
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // fini
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
}


