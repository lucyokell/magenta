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
  Rcpp::List populations_event_and_strains_List = savedState["populations_event_and_strains_List"];
  Rcpp::List scourge_List = savedState["scourge_List"];
  Rcpp::List parameters_List = savedState["parameters_List"];
  
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
  
  Rcpp::Rcout << "Paremeter list conversion is working!\n";
  
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
  std::vector<int> Number_of_Strains = Rcpp::as<vector<int> >(populations_event_and_strains_List["Number_of_Strains"]);
  std::vector<int> Number_of_Realised_Infections = Rcpp::as<vector<int> >(populations_event_and_strains_List["Number_of_Realised_Infections"]);
  
  Rcpp::List Infection_time_realisation_vectors = populations_event_and_strains_List["Infection_time_realisation_vectors"];
  Rcpp::List Infection_state_realisation_vectors = populations_event_and_strains_List["Infection_state_realisation_vectors"];
  Rcpp::List Infection_barcode_realisation_vectors = populations_event_and_strains_List["Infection_barcode_realisation_vectors"];
  
  Rcpp::List Strain_infection_state_vectors = populations_event_and_strains_List["Strain_infection_state_vectors"];
  Rcpp::List Strain_day_of_infection_state_change_vectors = populations_event_and_strains_List["Strain_day_of_infection_state_change_vectors"];
  Rcpp::List Strain_barcode_vectors = populations_event_and_strains_List["Strain_barcode_vectors"];
  
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS: humans
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  Rcpp::Rcout << "Human list conversions are working!\n";
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: R -> C++ CONVERSIONS: mosquitos
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  std::vector<int> Mosquito_Infection_States = Rcpp::as<vector<int> >(scourge_List["Mosquito_Infection_States"]);
  std::vector<unsigned short int> Mosquito_Day_of_next_blood_meal = Rcpp::as<vector<unsigned short int> >(scourge_List["Mosquito_Day_of_next_blood_meal"]);
  std::vector<unsigned short int> Mosquito_Day_of_death = Rcpp::as<vector<unsigned short int> >(scourge_List["Mosquito_Day_of_death"]);
  std::vector<unsigned short int> Mosquito_Number_of_ruptured_oocysts = Rcpp::as<vector<unsigned short int> >(scourge_List["Mosquito_Number_of_ruptured_oocysts"]);
  
  Rcpp::List Mosquito_Oocyst_rupture_time_vectors = scourge_List["Mosquito_Oocyst_rupture_time_vectors"];
  Rcpp::List Mosquito_Oocyst_barcode_male_vectors = scourge_List["Mosquito_Oocyst_barcode_male_vectors"];
  Rcpp::List Mosquito_Oocyst_barcode_female_vectors = scourge_List["Mosquito_Oocyst_barcode_female_vectors"];

  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS: mosquitos
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  Rcpp::Rcout << "Mosquito list unpacking working!\n";
  
  // Use read in R equilibrium state to then allocate each individual accordingly given tehir age and biting heterogeneity
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: HUMAN FETCHING
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Necessary temporary veriables for unpacking strains within humans
  Strain temp_strain;
  barcode_t temp_barcode = Strain::generate_random_barcode();
  std::vector<std::vector<bool> > temp_strain_barcode_vector;
  std::vector<int>  temp_strain_state_vector;
  std::vector<int>  temp_strain_state_change_time_vector;
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
    temp_strain_barcode_vector = Rcpp::as<std::vector<std::vector<bool> > >(Strain_barcode_vectors[n]);
    
    // Loop over each vector and making a temporary strain which is to be pushed onto the human
    for(int s = 0; s < Number_of_Strains[n]; s++)
    {
    // Set the temp strains infection state and day of strain state change
    temp_strain.set_m_strain_infection_status(static_cast<Strain::InfectionStatus>(temp_strain_state_vector[s]));
    temp_strain.set_m_day_of_strain_infection_status_change(temp_strain_state_change_time_vector[s]);
    
    // fetch vector<bool> and turn into barcode and then add to strain
    for(temp_barcode_iterator = 0; temp_barcode_iterator < barcode_length ; temp_barcode_iterator++ )
    {
      temp_barcode[temp_barcode_iterator] = temp_strain_barcode_vector[s][temp_barcode_iterator];
    }
    temp_strain.set_m_barcode(temp_barcode);
    
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
  
  Rcpp::Rcout << "Human fetching working\n";
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: MOSQUITO FETCHING
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  Rcpp::Rcout << "Pre mosquito-fetching working!\n";
  int test = 0;
  for (unsigned int n=0; n < scourge_size; n++) 
  {
    
    // Set their id
    scourge.emplace_back(parameters);
    scourge[n].set_m_mosquito_ID(n);
    
    // Set infection state
    scourge[n].set_m_mosquito_infection_state(static_cast<Mosquito::InfectionStatus>(Mosquito_Infection_States[n]));
    
    // Set number of ruptured oocysts
    scourge[n].set_m_ruptured_oocyst_count(Mosquito_Number_of_ruptured_oocysts[n]);
    
    // Set Day Changes
    scourge[n].set_m_day_of_next_blood_meal(Mosquito_Day_of_next_blood_meal[n]);
    scourge[n].set_m_day_of_death(Mosquito_Day_of_death[n]);
    
    if(test ==0){
      Rcpp::Rcout << "Pre-mosquito_vectors working\n";
      test=1;
    }
    
    // Set vectors;
    scourge[n].set_m_oocyst_rupture_time_vector(Rcpp::as<std::vector<unsigned short int> >(Mosquito_Oocyst_rupture_time_vectors[n]));
    scourge[n].set_m_oocyst_barcode_male_vector_from_vector_of_vector_bool(Rcpp::as<std::vector<std::vector<bool> > >(Mosquito_Oocyst_barcode_male_vectors[n]));
    scourge[n].set_m_oocyst_barcode_female_vector_from_vector_of_vector_bool(Rcpp::as<std::vector<std::vector<bool> > >(Mosquito_Oocyst_barcode_female_vectors[n]));
    
    // Set the next event day
    scourge[n].schedule_m_day_of_next_event();

    
  }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: MOSQUITO FETCHING
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  Rcpp::Rcout << "Mosquito fetching working\n";
  
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
  int g_end_time = static_cast<int>(parameters.g_current_time + (parameters.g_years * 365));
  
  // Preallocations;
  unsigned int num_bites = 0;
  unsigned int increasing_bites = 0;
  unsigned int individual_binomial_bite_draw = 0;
  double psi_sum = 0;
  
  // status eq for logging and other logging variables
  std::vector<double> status_eq = { 0,0,0,0,0,0 };
  unsigned int log_counter = 0;
  double total_incidence = 0;
  double total_incidence_05 = 0;
  int daily_incidence_return = 0;
  
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
  mosquito_biting_queue.reserve(scourge.size());
  std::vector<int> bite_storage_queue;
  bite_storage_queue.reserve(scourge.size());
  
  
  // resetart timer
  t0 = std::chrono::high_resolution_clock::now();
  
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
    
    // Human update loop
    for (human_update_i = 0; human_update_i < parameters.g_N; human_update_i++)
    {
      psi_sum += psi_vector[human_update_i] = population[human_update_i].update(parameters);
    }
    
    // Reset number of bites for each day and update each mosquito. Update returns whether the mosquito is biting today
    num_bites = 0;
    
    // Mosquito update loop
    for (mosquito_update_i = 0; mosquito_update_i < scourge_size; mosquito_update_i++)
    {
      if (scourge[mosquito_update_i].update(parameters))
      {
        mosquito_biting_queue.emplace_back(mosquito_update_i);
        num_bites++;
      }
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
    // START: BITE ALLOCATIONS
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Multinomial step //
    // --------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Caluclate total probability of being bitten within population
    pi_sum = std::accumulate(pi_vector.begin(), pi_vector.end(), 0.0);
    //num_bites = rpoisson1(pi_sum * Iv * 0.30677);
    
    // multinomial procedure for everyone and breaking when n_bites is reached
    for (bite_sampling_i = 0; bite_sampling_i < parameters.g_N - 1; bite_sampling_i++)
    {
      
      individual_binomial_bite_draw = rbinomial1(num_bites - increasing_bites, pi_vector[bite_sampling_i] / (pi_sum - pi_cum_sum));
      
      for (bite_sampling_internal_i = 0; bite_sampling_internal_i < individual_binomial_bite_draw; bite_sampling_internal_i++)
      {
        /*	bite_storage_queue.push(n);*/
        bite_storage_queue.emplace_back(bite_sampling_i);
      }
      
      pi_cum_sum += pi_vector[bite_sampling_i];
      increasing_bites += individual_binomial_bite_draw;
      if (increasing_bites >= num_bites) break;
      
    }
    
    // catch rounding errors so just place this here outside loop
    if (increasing_bites != num_bites)
    {
      individual_binomial_bite_draw = num_bites - increasing_bites;
      for (bite_sampling_internal_i = 0; bite_sampling_internal_i < individual_binomial_bite_draw; bite_sampling_internal_i++)
      {
        //bite_storage_queue.push(parameters.g_N - 1);
        bite_storage_queue.emplace_back(parameters.g_N - 1);
      }
    }
    
    // Reset bite and sum of biting rates
    increasing_bites = 0;
    pi_cum_sum = 0;
    
    // shuffle the bite queue otherwise you will introduce stepping-stone-esque genetic structuring
    // shuffle_integer_vector(bite_storage_queue);
    
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
      if (scourge[mosquito_biting_queue[num_bites_i]].m_mosquito_infected) {
        population[bite_storage_queue[num_bites_i]].allocate_bite(parameters, scourge[mosquito_biting_queue[num_bites_i]]);
      }
      
      // if human would cause infection to mosquito then allocate gametocytes
      if (population[bite_storage_queue[num_bites_i]].reciprocal_infection_boolean(parameters)) {
        scourge[mosquito_biting_queue[num_bites_i]].allocate_gametocytes(parameters, population[bite_storage_queue[num_bites_i]].sample_two_barcodes(parameters));
      }
      
    }
    
    // clear biting storage vectors
    bite_storage_queue.clear();
    mosquito_biting_queue.clear();
    
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // START: SUMMARY LOGGING
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Log the last week //
    
    if (parameters.g_current_time > g_end_time - 7)
    {
      
      log_counter++;
      
      for (auto &element : population)
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
        daily_incidence_return = element.log_daily_incidence(parameters);
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
    status_eq[element] /= (parameters.g_N * log_counter);
    Rcpp::Rcout << status_eq[element] << " | ";
  }
  
  
  // Final infection states
  std::vector<int> Infection_States(parameters.g_N);
  std::vector<double> Ages(parameters.g_N);
  std::vector<double> IB(parameters.g_N);
  std::vector<double> ICA(parameters.g_N);
  std::vector<double> ICM(parameters.g_N);
  std::vector<double> ID(parameters.g_N);
  
  // loop through population and grabages and immunities for error checking
  for (unsigned int element = 0; element < parameters.g_N ; element++) 
  {
    
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
  
  // Create Rcpp loggers list
  Rcpp::List Loggers = Rcpp::List::create(Rcpp::Named("S")=status_eq[0],Rcpp::Named("D")=status_eq[1],Rcpp::Named("A")=status_eq[2],
                                          Rcpp::Named("U")=status_eq[3],Rcpp::Named("T")=status_eq[4],Rcpp::Named("P")=status_eq[5],Rcpp::Named("Incidence")=total_incidence/log_counter,
                                          Rcpp::Named("Incidence_05")=total_incidence_05/log_counter,Rcpp::Named("InfectionStates")=Infection_States,Rcpp::Named("Ages")=Ages,
                                                      Rcpp::Named("IB")=IB,Rcpp::Named("ICA")=ICA,Rcpp::Named("ICM")=ICM,Rcpp::Named("ID")=ID);
  
  
  // Create universe ptr for memory-continuiation
  Rcpp::XPtr<Universe> universe_ptr(new Universe{ population, psi_vector, zeta_vector, pi_vector, scourge, parameters},
                                    true);
  
  // Return Named List with pointer and loggers
  return Rcpp::List::create(Rcpp::Named("Ptr") = universe_ptr, Rcpp::Named("Loggers")=Loggers);
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // fini
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
}


