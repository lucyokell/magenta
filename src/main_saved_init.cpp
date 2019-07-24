//
//  magenta
//  main_saved_init.cpp
//
//  Created: OJ Watson on 06/12/2015
//
//  NOT ACCURATE 
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
//' Creates initial model simulation using a saved model state
//'
//' @param param_list parameter list generated with \code{Param_List_Simulation_Get_Create}
//' @return list with ptr to model state and loggers describing the current model state
//' @export
//' 
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::List Simulation_Saved_Init_cpp(Rcpp::List param_list)
{
  
  // Initialise parameters
  Parameters parameters;
  
  // Unpack R List to Rcpp Lists
  Rcpp::List savedState = param_list["savedState"];
  Rcpp::List population_List = savedState["population_List"];
  Rcpp::List populations_event_and_strains_List = savedState["populations_event_and_strains_List"];
  Rcpp::List scourge_List = savedState["scourge_List"];
  Rcpp::List parameters_List = savedState["parameters_List"];
  
  
  // prove that C++ code is being run
  parameters.g_h_quiet_print = Rcpp::as<bool>(parameters_List["g_h_quiet_print"]);
  rcpp_out(parameters.g_h_quiet_print, "Rcpp function is working!\n");
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: R -> C++ CONVERSIONS: parameters
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Update this first as needed for initialising the population
  parameters.g_N = Rcpp::as<unsigned int>(parameters_List["g_N"]);
  parameters.g_current_time = Rcpp::as<int>(parameters_List["g_current_time"]);
  parameters.g_mean_maternal_immunity = Rcpp::as<double>(parameters_List["g_mean_maternal_immunity"]);
  parameters.g_sum_maternal_immunity = Rcpp::as<double>(parameters_List["g_sum_maternal_immunity"]);
  parameters.g_total_mums = Rcpp::as<int>(parameters_List["g_total_mums"]);
  parameters.g_theta = Rcpp::as<vector<double> >(parameters_List["g_theta"]);
  parameters.g_calendar_day = Rcpp::as<int>(parameters_List["g_calendar_day"]);
  parameters.g_mosquito_deficit = Rcpp::as<int>(parameters_List["g_mosquito_deficit"]);
  parameters.g_scourge_today = Rcpp::as<int>(parameters_List["g_scourge_today"]);
  parameters.g_mean_mv = Rcpp::as<int>(parameters_List["g_mean_mv"]);
  
  parameters.g_identity_id = Rcpp::as<int>(parameters_List["g_identity_id"]);
  parameters.g_num_loci = Rcpp::as<int>(parameters_List["g_num_loci"]);
  parameters.g_ibd_length = Rcpp::as<int>(parameters_List["g_ibd_length"]);
  parameters.g_barcode_length = Rcpp::as<int>(parameters_List["g_barcode_length"]);
  parameters.g_plaf = Rcpp::as<std::vector<double> >(parameters_List["g_plaf"]);
  parameters.g_prob_crossover = Rcpp::as<std::vector<double> >(parameters_List["g_prob_crossover"]);
  parameters.g_barcode_type = static_cast<Parameters::g_barcode_type_enum>(
    Rcpp::as<unsigned int>(parameters_List["g_barcode_type"])
  );
  parameters.g_spatial_type = static_cast<Parameters::g_spatial_type_enum>(
    Rcpp::as<unsigned int>(parameters_List["g_spatial_type"])
  );
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS: parameters
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  rcpp_out(parameters.g_h_quiet_print, "Paremeter list conversion is working!\n");
  
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
  
  // Add the mosquito population, i.e. scourge
  std::vector<Mosquito> scourge;
  unsigned int scourge_size = Rcpp::as<unsigned int>(scourge_List["Scourge_size"]);
  scourge.reserve(static_cast<int>(scourge_size));
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: R -> C++ CONVERSIONS: humans
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
  
  std::vector<int> Day_of_InfectionStatus_change = Rcpp::as<vector<int> >(populations_event_and_strains_List["Day_of_InfectionStatus_change"]);
  std::vector<int> Day_of_strain_clearance = Rcpp::as<vector<int> >(populations_event_and_strains_List["Day_of_strain_clearance"]);
  std::vector<int> Day_of_death = Rcpp::as<vector<int> >(populations_event_and_strains_List["Day_of_death"]);
  std::vector<int> Day_of_last_treatment = Rcpp::as<vector<int> >(populations_event_and_strains_List["Day_of_last_treatment"]);
  std::vector<int> Number_of_Strains = Rcpp::as<vector<int> >(populations_event_and_strains_List["Number_of_Strains"]);
  std::vector<int> Number_of_Realised_Infections = Rcpp::as<vector<int> >(populations_event_and_strains_List["Number_of_Realised_Infections"]);
  
  Rcpp::List Infection_time_realisation_vectors = populations_event_and_strains_List["Infection_time_realisation_vectors"];
  Rcpp::List Infection_state_realisation_vectors = populations_event_and_strains_List["Infection_state_realisation_vectors"];
  Rcpp::List Infection_barcode_realisation_vectors = populations_event_and_strains_List["Infection_barcode_realisation_vectors"];
  
  Rcpp::List Strain_infection_state_vectors = populations_event_and_strains_List["Strain_infection_state_vectors"];
  Rcpp::List Strain_day_of_infection_state_change_vectors = populations_event_and_strains_List["Strain_day_of_infection_state_change_vectors"];
  Rcpp::List Strain_day_of_acquisition_vectors = populations_event_and_strains_List["Strain_day_of_acquisition_vectors"];
  Rcpp::List Strain_barcode_vectors = populations_event_and_strains_List["Strain_barcode_vectors"];
  
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS: humans
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  rcpp_out(parameters.g_h_quiet_print, "Human list conversions are working!\n");
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: R -> C++ CONVERSIONS: mosquitos
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  std::vector<int> Mosquito_Infection_States = Rcpp::as<vector<int> >(scourge_List["Mosquito_Infection_States"]);
  std::vector<bool> Mosquito_Off_Season = Rcpp::as<vector<bool> >(scourge_List["Mosquito_Off_Season"]);
  std::vector<unsigned int> Mosquito_Day_of_next_blood_meal = Rcpp::as<vector<unsigned int> >(scourge_List["Mosquito_Day_of_next_blood_meal"]);
  std::vector<unsigned int> Mosquito_Day_of_death = Rcpp::as<vector<unsigned int> >(scourge_List["Mosquito_Day_of_death"]);
  std::vector<unsigned int> Mosquito_Number_of_ruptured_oocysts = Rcpp::as<vector<unsigned int> >(scourge_List["Mosquito_Number_of_ruptured_oocysts"]);
  
  Rcpp::List Mosquito_Oocyst_rupture_time_vectors = scourge_List["Mosquito_Oocyst_rupture_time_vectors"];
  Rcpp::List Mosquito_Oocyst_remaining_spz_vectors = scourge_List["Mosquito_Oocyst_remaining_spz_vectors"];
  Rcpp::List Mosquito_Oocyst_barcode_male_vectors = scourge_List["Mosquito_Oocyst_barcode_male_vectors"];
  Rcpp::List Mosquito_Oocyst_barcode_female_vectors = scourge_List["Mosquito_Oocyst_barcode_female_vectors"];

  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS: mosquitos
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  rcpp_out(parameters.g_h_quiet_print, "Mosquito list unpacking working!\n");
  
  // Use read in R equilibrium state to then allocate each individual accordingly given tehir age and biting heterogeneity
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: HUMAN FETCHING
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Necessary temporary veriables for unpacking strains within humans
  Strain temp_strain;
  std::vector<std::vector<bool> > temp_strain_barcode_vector;
  std::vector<int>  temp_strain_state_vector;
  std::vector<int>  temp_strain_state_change_time_vector;
  std::vector<int>  temp_strain_day_of_acquisition;
  unsigned int temp_barcode_iterator = 0;
  
  for (unsigned int n=0; n < parameters.g_N; n++) 
  {
    
    // Set their id
    population.emplace_back(parameters);
    population[n].set_m_person_ID(n);
    
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
    population[n].set_m_day_last_treated(Day_of_last_treatment[n]);
    
    // Strain Numbers
    population[n].set_m_number_of_strains(Number_of_Strains[n]);
    
    // Strain Numbers
    population[n].set_m_number_of_realised_infections(Number_of_Realised_Infections[n]);
    
    // Set vectors
    population[n].set_m_infection_time_realisation_vector_from_vector(Rcpp::as<std::vector<int> >(Infection_time_realisation_vectors[n]));
    population[n].set_m_infection_state_realisation_vector_from_vector(Rcpp::as<std::vector<int> >(Infection_state_realisation_vectors[n]));
    population[n].set_m_infection_barcode_realisation_vector_from_vector_of_vector_bool(Rcpp::as<std::vector<std::vector<bool> > >(Infection_barcode_realisation_vectors[n]));
    
    // Set strains
    // -----------------------------------------
    
    // First grab the vector of times, states and barcode vectors for the indivdual
    temp_strain_state_vector = Rcpp::as<std::vector<int> >(Strain_infection_state_vectors[n]);
    temp_strain_state_change_time_vector = Rcpp::as<std::vector<int> >(Strain_day_of_infection_state_change_vectors[n]);
    temp_strain_day_of_acquisition = Rcpp::as<std::vector<int> >(Strain_day_of_acquisition_vectors[n]);
    temp_strain_barcode_vector = Rcpp::as<std::vector<std::vector<bool> > >(Strain_barcode_vectors[n]);
    
    // Loop over each vector and making a temporary strain which is to be pushed onto the human
    for(int s = 0; s < Number_of_Strains[n]; s++)
    {
    // Set the temp strains infection state and day of strain state change
    temp_strain.set_m_strain_infection_status(static_cast<Strain::InfectionStatus>(temp_strain_state_vector[s]));
    temp_strain.set_m_day_of_strain_infection_status_change(temp_strain_state_change_time_vector[s]);
    temp_strain.set_m_day_of_strain_acquisition(temp_strain_day_of_acquisition[s]);
    
    // fetch vector<bool> and turn into barcode and then add to strain
    for(temp_barcode_iterator = 0; temp_barcode_iterator < Parameters::g_barcode_length ; temp_barcode_iterator++ )
    {
      Strain::temp_barcode[temp_barcode_iterator] = temp_strain_barcode_vector[s][temp_barcode_iterator];
    }
    temp_strain.set_m_barcode(Strain::temp_barcode);
    
    // Push the created temp strain onto the person
    population[n].allocate_strain_with_push(temp_strain);
    }
    
    // -----------------------------------------
    
    // Set the next strain change day and next event day
    population[n].set_m_day_of_next_strain_state_change();
    population[n].set_m_day_of_next_event();
    
    // Asign zeta and add to the zeta vector which is required for calculating the overall probability of being bitten, pi
    population[n].set_m_individual_biting_rate(Zetas[n]);
    zeta_vector[n] = population[n].get_m_individual_biting_rate();
    
  }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: HUMAN FETCHING
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  rcpp_out(parameters.g_h_quiet_print, "Human fetching working\n");
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: MOSQUITO FETCHING
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  rcpp_out(parameters.g_h_quiet_print, "Pre mosquito-fetching working!\n");
  int test = 0;
  for (unsigned int n=0; n < scourge_size; n++) 
  {
    
    // Set their id
    scourge.emplace_back(parameters);
    scourge[n].set_m_mosquito_ID(n);
    
    // Set infection state
    scourge[n].set_m_mosquito_infection_state(static_cast<Mosquito::InfectionStatus>(Mosquito_Infection_States[n]));
    if(scourge[n].get_m_mosquito_infection_state() == Mosquito::INFECTED){
      scourge[n].m_mosquito_infected = true;
    }
    
    // Set number of ruptured oocysts
    scourge[n].set_m_ruptured_oocyst_count(Mosquito_Number_of_ruptured_oocysts[n]);
    
    // Set Day Changes
    scourge[n].set_m_day_of_next_blood_meal(Mosquito_Day_of_next_blood_meal[n]);
    scourge[n].set_m_day_of_death(Mosquito_Day_of_death[n]);
    scourge[n].set_m_mosquito_off_season(Mosquito_Off_Season[n]);
    
    if(test ==0){
      rcpp_out(parameters.g_h_quiet_print, "Pre-mosquito_vectors working\n");
      test=1;
    }
    
    // Set vectors;
    scourge[n].set_m_oocyst_rupture_time_vector(Rcpp::as<std::vector<int> >(Mosquito_Oocyst_rupture_time_vectors[n]));
    scourge[n].set_m_oocyst_remaining_spz_count(Rcpp::as<std::vector<int> >(Mosquito_Oocyst_remaining_spz_vectors[n]));
    scourge[n].set_m_oocyst_barcode_male_vector_from_vector_of_vector_bool(Rcpp::as<std::vector<std::vector<bool> > >(Mosquito_Oocyst_barcode_male_vectors[n]));
    scourge[n].set_m_oocyst_barcode_female_vector_from_vector_of_vector_bool(Rcpp::as<std::vector<std::vector<bool> > >(Mosquito_Oocyst_barcode_female_vectors[n]));
    
    // Set the next event day
    scourge[n].schedule_m_day_of_next_event();

    
  }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: MOSQUITO FETCHING
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  rcpp_out(parameters.g_h_quiet_print, "Mosquito fetching working\n");
  
  /*
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: INITIALISATION
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  */
  
  // Final infection states
  std::vector<int> Infection_States(parameters.g_N);
  std::vector<double> Ages(parameters.g_N);
  std::vector<double> IB(parameters.g_N);
  std::vector<double> ICA(parameters.g_N);
  std::vector<double> ICM(parameters.g_N);
  std::vector<double> ID(parameters.g_N);
  
  // status eq for logging and other logging variables
  std::vector<double> status_eq = { 0,0,0,0,0,0 };
  
  // loop through population and grab ages and immunities for error checking
  for (unsigned int element = 0; element < parameters.g_N ; element++) 
  {
    
    // Match infection state and schedule associated next state change
    switch (population[element].get_m_infection_state())
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
    
    // Ages and immunity 
    // TODO: Figure out the best way of standardising this logging 
    // Something like passing in a function name within the param_list which is the 
    // name for a logger written else where which then returns the Loggers obeject below
    Infection_States[element] = static_cast<int>(population[element].get_m_infection_state());
    population[element].update_immunities_to_today(parameters);
    Ages[element] = population[element].get_m_person_age();
    IB[element] = population[element].get_m_IB();
    ICA[element] = population[element].get_m_ICA();
    ICM[element] = population[element].get_m_ICM();
    ID[element] = population[element].get_m_ID();
    
  }
  
  // divide by population size and log counter and print to give overview
  rcpp_out(parameters.g_h_quiet_print, "S | D | A | U | T | P:\n");
  
  for (int element = 0; element < 6; element++) 
  {
    status_eq[element] /= (parameters.g_N);
    rcpp_out(parameters.g_h_quiet_print, std::to_string(status_eq[element]) + " | ");
  }
  
  // Create Rcpp loggers list
  Rcpp::List Loggers = Rcpp::List::create(Rcpp::Named("S")=status_eq[0], 
                                          Rcpp::Named("D")=status_eq[1],Rcpp::Named("A")=status_eq[2],
                                                                                                  Rcpp::Named("U")=status_eq[3],Rcpp::Named("T")=status_eq[4],Rcpp::Named("P")=status_eq[5],
                                                                                                                                                                                        Rcpp::Named("InfectionStates")=Infection_States, Rcpp::Named("Ages")=Ages, 
                                                                                                                                                                                        Rcpp::Named("IB")=IB,Rcpp::Named("ICA")=ICA,Rcpp::Named("ICM")=ICM,Rcpp::Named("ID")=ID
                                            );
  
  
  // Create universe ptr for memory-continuiation
  Rcpp::XPtr<Universe> universe_ptr(new Universe{ population, psi_vector, zeta_vector, pi_vector, scourge, parameters},
                                    true);
  
  // Return Named List with pointer and loggers
  return Rcpp::List::create(Rcpp::Named("Ptr") = universe_ptr, Rcpp::Named("Loggers")=Loggers);
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // fini
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
}


