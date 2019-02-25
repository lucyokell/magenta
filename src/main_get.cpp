//
//  magenta
//  main_get.cpp
//
//  Created: OJ Watson on 06/12/2015
//
//  Distributed under the MIT software licence
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

//#include <RcppArmadillo.h>
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
//' Returns whole model to R in series of nested lists
//'
//' @param param_list parameter list generated with \code{Param_List_Simulation_Get_Create}
//' @return list of 4 lists with the entire model state
//' @export
// [[Rcpp::export]]
Rcpp::List Simulation_Get_cpp(Rcpp::List param_list)
{
  

  // start timer
  chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Create universe pointer from param_list statePtr
  Rcpp::XPtr<Universe> universe_ptr = Rcpp::as<Rcpp::XPtr<Universe> > (param_list["statePtr"]);
  
    // prove that C++ code is being run
  rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Rcpp function is working!\n");
  
  // Initialise variables from the statePtr provided that are needing to be saved
  std::vector<double> Zeta = universe_ptr->zeta_vector;
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Pointer unpacking working!\n");
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: MODEL STATE GET
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: HUMAN POPULATION
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Create model state storage
  std::vector<int> Infection_States(universe_ptr->parameters.g_N);
  std::vector<int> Ages(universe_ptr->parameters.g_N);
  std::vector<double> IB(universe_ptr->parameters.g_N);
  std::vector<double> ICA(universe_ptr->parameters.g_N);
  std::vector<double> ICM(universe_ptr->parameters.g_N);
  std::vector<double> ID(universe_ptr->parameters.g_N);
  std::vector<double> IB_last_boost_time(universe_ptr->parameters.g_N);
  std::vector<double> ICA_last_boost_time(universe_ptr->parameters.g_N);
  std::vector<double> ID_last_boost_time(universe_ptr->parameters.g_N);
  std::vector<int> IB_last_calculated_time(universe_ptr->parameters.g_N);
  std::vector<int> I_C_D_CM_last_calculated_time(universe_ptr->parameters.g_N);
  std::vector<double> Immunity_boost_float(universe_ptr->parameters.g_N);
  std::vector<int>  Day_of_InfectionStatus_change(universe_ptr->parameters.g_N);	
  std::vector<int>  Day_of_strain_clearance(universe_ptr->parameters.g_N);		
  std::vector<int>  Day_of_death(universe_ptr->parameters.g_N);					
  std::vector<int>  Number_of_Strains(universe_ptr->parameters.g_N);	
  std::vector<int> Number_of_Realised_Infections(universe_ptr->parameters.g_N);
  std::vector<int> Day_of_next_strain_state_change(universe_ptr->parameters.g_N);
  std::vector<int> Day_of_next_event(universe_ptr->parameters.g_N);
  std::vector<int> Day_of_last_treatment(universe_ptr->parameters.g_N);
  std::vector<std::vector<int> > Infection_time_realisation_vectors(universe_ptr->parameters.g_N);
  std::vector<std::vector<int> > Infection_state_realisation_vectors(universe_ptr->parameters.g_N);
  std::vector<std::vector<std::vector<bool> > > Infection_barcode_realisation_vectors(universe_ptr->parameters.g_N);
  std::vector<std::vector<std::vector<bool> > > Strain_barcode_vectors(universe_ptr->parameters.g_N);
  std::vector<std::vector<int> > Strain_infection_state_vectors(universe_ptr->parameters.g_N);
  std::vector<std::vector<int> > Strain_day_of_infection_state_change_vectors(universe_ptr->parameters.g_N);
  std::vector<std::vector<int> > Strain_day_of_acquisition_vectors(universe_ptr->parameters.g_N);
  std::vector<std::vector<bool> > Strain_cotransmission(universe_ptr->parameters.g_N);
  
  std::vector<std::vector<unsigned int> > recent_barcode_integers(universe_ptr->parameters.g_N);
  
  rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Vector initialisation working!\n");
  
  // Temporary necessities for casting vectors for pending states
  std::vector<Person::InfectionStatus> temp_infection_state_realisation_vector{};
  std::vector<bool> temp_cotransmission_vector{};
  unsigned int temp_status_iterator = 0;
  
  // Temporary necessities for pending barcodes
  std::vector<boost::dynamic_bitset<> > temp_infection_barcode_realisation_vector = {};
  std::vector<bool> temp_barcode_bool_vector(universe_ptr->parameters.g_barcode_length,false);
  unsigned int temp_barcode_iterator = 0;
  
  // Temporary necessities for strains
  Strain temp_strain;
  boost::dynamic_bitset<> temp_barcode;
  int temp_strain_iterator = 0;

  rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Preloop working!\n");

  
  for (unsigned int element = 0; element < universe_ptr->parameters.g_N ; element++) 
  {
    
    // Infection States
    Infection_States[element] = static_cast<int>(universe_ptr->population[element].get_m_infection_state());
    
    // Ages
    Ages[element] = universe_ptr->population[element].get_m_person_age();
    
    // Immunities
    IB[element] = universe_ptr->population[element].get_m_IB();
    ICA[element] = universe_ptr->population[element].get_m_ICA();
    ICM[element] = universe_ptr->population[element].get_m_ICM();
    ID[element] = universe_ptr->population[element].get_m_ID();
    
    // Boost times
    IB_last_boost_time[element] = universe_ptr->population[element].get_m_IB_last_boost_time();
    ICA_last_boost_time[element] = universe_ptr->population[element].get_m_ICA_last_boost_time();
    ID_last_boost_time[element] = universe_ptr->population[element].get_m_ID_last_boost_time();
    
    // Calc times
    IB_last_calculated_time[element] = universe_ptr->population[element].get_m_IB_last_calculated_time();
    I_C_D_CM_last_calculated_time[element] = universe_ptr->population[element].get_m_I_C_D_CM_last_calculated_time();
    
    // Immunity float
    Immunity_boost_float[element] = universe_ptr->population[element].get_m_immunity_boost_float();
    
    // Day Changes
    Day_of_InfectionStatus_change[element] = universe_ptr->population[element].get_m_day_of_InfectionStatus_change();
    Day_of_strain_clearance[element] = universe_ptr->population[element].get_m_day_of_strain_clearance();
    Day_of_death[element]	 = universe_ptr->population[element].get_m_day_of_death();
    Day_of_last_treatment[element] = universe_ptr->population[element].get_m_day_last_treated(); 
    
    // Strain Numbers
    Number_of_Strains[element] = universe_ptr->population[element].get_m_number_of_strains();
    
    // Temp strain to next change state
    Day_of_next_strain_state_change[element] = universe_ptr->population[element].get_m_day_of_next_strain_state_change();
    Day_of_next_event[element] = universe_ptr->population[element].get_m_day_of_next_event();
    
    // Realised Infections
    Number_of_Realised_Infections[element] = universe_ptr->population[element].get_m_number_of_realised_infections();
    
    // Pending Infection time vector
    Infection_time_realisation_vectors[element] = universe_ptr->population[element].get_m_infection_time_realisation_vector();
    
    // rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Prepending working!" + std::to_string(element) "\n");
    
    // Pending Infection barcode and state vector
    // ---------------------------------------
    temp_infection_barcode_realisation_vector = universe_ptr->population[element].get_m_infection_barcode_realisation_vector();
    Infection_barcode_realisation_vectors[element].reserve(temp_infection_barcode_realisation_vector.size());
    
    temp_infection_state_realisation_vector = universe_ptr->population[element].get_m_infection_state_realisation_vector();
    temp_cotransmission_vector = universe_ptr->population[element].get_m_cotransmission_realisation_vector();

    Infection_state_realisation_vectors[element].reserve(temp_infection_state_realisation_vector.size());

    for(temp_status_iterator = 0 ; temp_status_iterator < static_cast<unsigned int>(temp_infection_state_realisation_vector.size()) ; )
    {
      Infection_state_realisation_vectors[element].emplace_back(temp_infection_state_realisation_vector[temp_status_iterator]);
      Strain_cotransmission[element].emplace_back(temp_cotransmission_vector[temp_status_iterator]);

      for(temp_barcode_iterator = 0; temp_barcode_iterator < universe_ptr->parameters.g_barcode_length ; temp_barcode_iterator++ )
      {
        temp_barcode_bool_vector[temp_barcode_iterator] = static_cast<bool>(temp_infection_barcode_realisation_vector[temp_status_iterator][temp_barcode_iterator]);
      }

      Infection_barcode_realisation_vectors[element].push_back(temp_barcode_bool_vector);
      temp_status_iterator++;
    }
    
    temp_infection_barcode_realisation_vector.clear();
    temp_infection_state_realisation_vector.clear();

    // Active strains
    // ---------------------------------------
    // Reserve space for active strains
    Strain_barcode_vectors[element].reserve(Number_of_Strains[element]);
    Strain_infection_state_vectors[element].reserve(Number_of_Strains[element]);
    Strain_day_of_infection_state_change_vectors[element].reserve(Number_of_Strains[element]);
    
    std::vector<Strain> temp_strain_vector;
    temp_strain_vector = universe_ptr->population[element].get_m_active_strains();
    
    // rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Prestrains loop working!\n");
    
    // Loop through each strain converting into barcodes, states and days of infection state changes
    for(temp_strain_iterator = 0 ; temp_strain_iterator < Number_of_Strains[element] ; temp_strain_iterator++)
    {
      
      temp_strain = universe_ptr->population[element].get_m_person_strain_x(temp_strain_iterator);
      
      Strain_infection_state_vectors[element].push_back(temp_strain.get_m_strain_infection_status());
      Strain_day_of_infection_state_change_vectors[element].push_back(temp_strain.get_m_day_of_strain_infection_status_change());
      Strain_day_of_acquisition_vectors[element].push_back(temp_strain.get_m_day_of_strain_acquisition());
      // fetch barcode and turn into vector<bool>
      temp_barcode = temp_strain.get_m_barcode();
      for(temp_barcode_iterator = 0; temp_barcode_iterator < universe_ptr->parameters.g_barcode_length ; temp_barcode_iterator++ )
      {
        temp_barcode_bool_vector[temp_barcode_iterator] = temp_barcode[temp_barcode_iterator];
      }
      Strain_barcode_vectors[element].push_back(temp_barcode_bool_vector);
      
      if(temp_strain_iterator == (Number_of_Strains[element]-1)){
        if(universe_ptr->parameters.g_barcode_type == Parameters::IBD){
          recent_barcode_integers[element].reserve(1);
          recent_barcode_integers[element] = Strain::ibd_barcode_to_integer_vector(temp_barcode);
        }
      }
    }
    
    // rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Poststrains loop working!\n");
    
  }
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // END: HUMAN POPULATION
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // START: MOSQUITO POPULATION
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Create model state storage
    unsigned int scourge_size = universe_ptr->scourge.size();
  
    std::vector<int> Mosquito_Infection_States(scourge_size);
    std::vector<int> Mosquito_Day_of_next_blood_meal(scourge_size);	
    std::vector<bool> Mosquito_Off_Season(scourge_size);	
    std::vector<int> Mosquito_Day_of_death(scourge_size);
    std::vector<int> Mosquito_Number_of_ruptured_oocysts (scourge_size);
    std::vector<std::vector<int> > Mosquito_Oocyst_rupture_time_vectors(scourge_size);
    std::vector<std::vector<std::vector<bool> > > Mosquito_Oocyst_barcode_male_vectors(scourge_size);
    std::vector<std::vector<std::vector<bool> > > Mosquito_Oocyst_barcode_female_vectors(scourge_size);
    
    // Temporary necessities for pending barcodes
    std::vector<boost::dynamic_bitset<>> temp_male_barcode_realisation_vector = {};
    std::vector<boost::dynamic_bitset<>> temp_female_barcode_realisation_vector = {};

    std::vector<bool> temp_barcode_male_bool_vector{};
    std::vector<bool> temp_barcode_female_bool_vector{};
    temp_barcode_female_bool_vector.reserve(universe_ptr->parameters.g_barcode_length);
    temp_barcode_male_bool_vector.reserve(universe_ptr->parameters.g_barcode_length);
    
    rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Premosquito loop working!\n");

    for (unsigned int element = 0; element < scourge_size ; element++) 
    {
      
      // Infection States
      Mosquito_Infection_States[element] = static_cast<int>(universe_ptr->scourge[element].get_m_mosquito_infection_state());

      // Day Changes
      Mosquito_Day_of_death[element] = universe_ptr->scourge[element].get_m_day_of_death();

      Mosquito_Day_of_next_blood_meal[element] = universe_ptr->scourge[element].get_m_day_of_next_blood_meal();

      Mosquito_Off_Season[element] = universe_ptr->scourge[element].get_m_mosquito_off_season();
      
      // Ruptured oocyst Numbers 
      Mosquito_Number_of_ruptured_oocysts[element] = static_cast<int>(universe_ptr->scourge[element].get_m_ruptured_oocyst_count());

      // Pending Infection time vector
      Mosquito_Oocyst_rupture_time_vectors[element] = universe_ptr->scourge[element].get_m_oocyst_rupture_time_vector();

      // Pending oocyst barcode vector
      // ---------------------------------------
      temp_male_barcode_realisation_vector = universe_ptr->scourge[element].get_m_oocyst_barcode_male_vector();

      temp_female_barcode_realisation_vector = universe_ptr->scourge[element].get_m_oocyst_barcode_female_vector();

      Mosquito_Oocyst_barcode_male_vectors[element].reserve(temp_male_barcode_realisation_vector.size());

      Mosquito_Oocyst_barcode_female_vectors[element].reserve(temp_female_barcode_realisation_vector.size());

      temp_status_iterator = 0;

      for(temp_status_iterator = 0; temp_status_iterator < temp_male_barcode_realisation_vector.size() ; temp_status_iterator++)
      {

        // fetch barcode and turn into vector<bool>
        for(temp_barcode_iterator = 0; temp_barcode_iterator < universe_ptr->parameters.g_barcode_length ; temp_barcode_iterator++ )
        {
          temp_barcode_male_bool_vector.push_back(temp_male_barcode_realisation_vector[temp_status_iterator][temp_barcode_iterator]);
          temp_barcode_female_bool_vector.push_back(temp_female_barcode_realisation_vector[temp_status_iterator][temp_barcode_iterator]);
        }
        
        Mosquito_Oocyst_barcode_male_vectors[element].push_back(temp_barcode_male_bool_vector);
        Mosquito_Oocyst_barcode_female_vectors[element].push_back(temp_barcode_female_bool_vector);
        temp_barcode_male_bool_vector.clear();
        temp_barcode_female_bool_vector.clear();

      }

      temp_male_barcode_realisation_vector.clear();
      temp_female_barcode_realisation_vector.clear();

      // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      // END: MOSQUITO POPULATION
      // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    }
    
  
  
  rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Postloop working!\n");
    
  
  // Create Rcpp population list
  Rcpp::List population_List = Rcpp::List::create(
    Rcpp::Named("Infection_States")=Infection_States,
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
    Rcpp::Named("Immunity_boost_float")=Immunity_boost_float
  );
  
  // Create Rcpp population event_and_strain list
  Rcpp::List populations_event_and_strains_List = Rcpp::List::create(
    Rcpp::Named("Day_of_InfectionStatus_change")=Day_of_InfectionStatus_change,
    Rcpp::Named("Day_of_strain_clearance")=Day_of_strain_clearance,
    Rcpp::Named("Day_of_death")=Day_of_death,
    Rcpp::Named("Number_of_Strains")=Number_of_Strains,
    Rcpp::Named("Day_of_next_strain_state_change")=Day_of_next_strain_state_change,
    Rcpp::Named("Day_of_next_event")=Day_of_next_event,
    Rcpp::Named("Day_of_last_treatment")=Day_of_last_treatment,
    Rcpp::Named("Number_of_Realised_Infections")=Number_of_Realised_Infections,
    Rcpp::Named("Infection_time_realisation_vectors")=Infection_time_realisation_vectors,
    Rcpp::Named("Infection_state_realisation_vectors")=Infection_state_realisation_vectors,
    Rcpp::Named("Infection_barcode_realisation_vectors")=Infection_barcode_realisation_vectors,
    Rcpp::Named("Strain_infection_state_vectors")=Strain_infection_state_vectors,
    Rcpp::Named("Strain_day_of_infection_state_change_vectors")=Strain_day_of_infection_state_change_vectors,
    Rcpp::Named("Strain_day_of_acquisition_vectors")=Strain_day_of_acquisition_vectors,
    Rcpp::Named("Strain_barcode_vectors")=Strain_barcode_vectors,
    Rcpp::Named("Strain_cotransmission")=Strain_cotransmission,
    Rcpp::Named("Recent_identity_vectors")=recent_barcode_integers
  );
  
  // Create Mosquito population list
  Rcpp::List scourge_List = Rcpp::List::create(
    Rcpp::Named("Mosquito_Infection_States")=Mosquito_Infection_States,
    Rcpp::Named("Mosquito_Day_of_next_blood_meal")=Mosquito_Day_of_next_blood_meal,
    Rcpp::Named("Mosquito_Off_Season")=Mosquito_Off_Season,
    Rcpp::Named("Mosquito_Day_of_death")=Mosquito_Day_of_death,
    Rcpp::Named("Mosquito_Number_of_ruptured_oocysts")=Mosquito_Number_of_ruptured_oocysts,
    Rcpp::Named("Mosquito_Oocyst_rupture_time_vectors")=Mosquito_Oocyst_rupture_time_vectors,
    Rcpp::Named("Mosquito_Oocyst_barcode_male_vectors")=Mosquito_Oocyst_barcode_male_vectors,
    Rcpp::Named("Mosquito_Oocyst_barcode_female_vectors")=Mosquito_Oocyst_barcode_female_vectors,
    Rcpp::Named("Scourge_size")=universe_ptr->scourge.size()
  );
  
  // Create Rcpp Parameters list
  Rcpp::List parameters_List = Rcpp::List::create(
    Rcpp::Named("g_current_time")=universe_ptr->parameters.g_current_time,
    Rcpp::Named("g_mean_maternal_immunity")=universe_ptr->parameters.g_mean_maternal_immunity,
    Rcpp::Named("g_sum_maternal_immunity")=universe_ptr->parameters.g_sum_maternal_immunity,
    Rcpp::Named("g_total_mums")=universe_ptr->parameters.g_total_mums,
    Rcpp::Named("g_N")=universe_ptr->parameters.g_N,
    Rcpp::Named("g_theta")=universe_ptr->parameters.g_theta,
    Rcpp::Named("g_calendar_day")=universe_ptr->parameters.g_calendar_day,
    Rcpp::Named("g_mosquito_deficit")=universe_ptr->parameters.g_mosquito_deficit,
    Rcpp::Named("g_scourge_today")=universe_ptr->parameters.g_scourge_today,
    Rcpp::Named("g_mean_mv")=universe_ptr->parameters.g_mean_mv,
    // genetics
    Rcpp::Named("g_identity_id")=universe_ptr->parameters.g_identity_id,
    Rcpp::Named("g_num_loci")=universe_ptr->parameters.g_num_loci,
    Rcpp::Named("g_ibd_length")=universe_ptr->parameters.g_ibd_length,
    Rcpp::Named("g_barcode_length")=universe_ptr->parameters.g_barcode_length,
    Rcpp::Named("g_plaf")=universe_ptr->parameters.g_plaf,
    Rcpp::Named("g_prob_crossover")=universe_ptr->parameters.g_prob_crossover,
    Rcpp::Named("g_barcode_type")=static_cast<unsigned int>(universe_ptr->parameters.g_barcode_type),
    Rcpp::Named("g_spatial_type")=static_cast<unsigned int>(universe_ptr->parameters.g_spatial_type),
    // housekeeping
    Rcpp::Named("g_h_quiet_print")=universe_ptr->parameters.g_h_quiet_print
  );
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: MODEL STATE GET
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
   rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Time elapsed in fetching state: " + std::to_string(duration) + " seconds\n");
  
  // Return Named List with population and parameters
  return Rcpp::List::create(
    Rcpp::Named("population_List") = population_List, 
    Rcpp::Named("populations_event_and_strains_List") = populations_event_and_strains_List,
    Rcpp::Named("scourge_List") = scourge_List,
    Rcpp::Named("parameters_List")=parameters_List);
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // fini
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
}


