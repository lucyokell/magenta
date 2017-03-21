//
//  MAGENTA
//  main_init.cpp
//
//  Created: OJ Watson on 06/12/2015
//	Most recent edits: OJ on 10/03/2017
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
Rcpp::List Simulation_Init_cpp(Rcpp::List paramList)
{
  
  // prove that C++ code is being run
  Rcpp::Rcout << "Rcpp function is working!\n";
  
  // start timer
  chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
  
  // Initialise parameters object
  Parameters parameters;
  
  // Grab these R parameters first so that they are used in the initialisation
  parameters.g_N = Rcpp::as<unsigned int>(paramList["N"]);
  parameters.g_years = Rcpp::as<double>(paramList["years"]);
  
  // Unpack the R equilibirum state parameter list object
  Rcpp::List eqSS = paramList["eqSS"];
  
  // Mosquito steady state values at a population level
  double Sv = Rcpp::as<double>(eqSS["Sv"]);
  double Ev = Rcpp::as<double>(eqSS["Ev"]);
  double Iv = Rcpp::as<double>(eqSS["Iv"]);
  
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
  scourge.reserve(static_cast<int>(parameters.g_N * (Sv + Ev + Iv)));
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: INITIALISATION
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
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
  double ICM_Init = Rcpp::as<double>(eqSS["MaternalImmunity"]);
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  Rcpp::Rcout << "Matrix unpacking working!\n";
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: HUMAN INITIALISATION FROM EQUILIBRIUM
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Use read in R equilibrium state to then allocate each individual accordingly given tehir age and biting heterogeneity
  
  // Initialisation temp variables required
  std::vector<int> infected_human_count;
  infected_human_count.reserve(parameters.g_N);
  int id_counter = 0;
  int num_age_brackets = age_brackets.size();
  int num_het_brackets = het_brackets.size();
  int age_bracket_in = 0;
  int het_bracket_in = 0;
  std::vector<double> infection_state_probability(6);
  
  // Temporary strain vector needed for initialisation
  std::vector<Strain> temp_strains;
  temp_strains.reserve(100);
  
  for (unsigned int n=0; n < parameters.g_N; n++) 
  {
    
    // Set their id
    population.push_back(Person(parameters));
    population[n].set_m_person_ID(id_counter++);
    
    // first find out what age bracket they are in
    // are they in the last - do this first for those people with ages that happen to be above the max age bracket
    if(population[n].get_m_person_age() >= (age_brackets[num_age_brackets-1]))
    {
      age_bracket_in = num_age_brackets-1;
    } 
    // if not then loop up to last-1
    for(int age_i = 0 ; age_i < (num_age_brackets-1) ; age_i++)
    {
      if(population[n].get_m_person_age() >= age_brackets[age_i] && population[n].get_m_person_age() < age_brackets[age_i+1])
      {
        age_bracket_in = age_i;
      }
    }
    
    // second find out what heterogeneity bracket they are in
    // are they in the last - do this first for those people with heterogeneities that happen to be above the max 
    if(population[n].get_m_individual_biting_rate() >= het_brackets[num_het_brackets-1]) 
    {
      het_bracket_in = num_het_brackets-1;
    }
    if(population[n].get_m_individual_biting_rate() < het_brackets[1])
    {
      het_bracket_in = 0;
    }
    for(int het_i = 1 ; het_i < (num_het_brackets-1) ; het_i++)
    {
      if(population[n].get_m_individual_biting_rate() < het_brackets[het_i+1])
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
    population[n].set_m_infection_state(static_cast<Person::InfectionStatus>(sample1(infection_state_probability, std::accumulate(infection_state_probability.begin(), infection_state_probability.end(),0.0))));
    population[n].set_m_IB(IBmat(age_bracket_in, het_bracket_in));
    population[n].set_m_ICA(ICAmat(age_bracket_in, het_bracket_in));
    population[n].set_m_ICM_init(ICM_Init);
    population[n].set_m_ICM(ICM_Init * exp(-population[n].get_m_person_age() / parameters.g_dCM));
    population[n].set_m_ID(IDmat(age_bracket_in, het_bracket_in));
    
    // Schedule change for those who are not susceptible
    if (population[n].get_m_infection_state() != Person::SUSCEPTIBLE)
    {
      population[n].schedule_m_day_of_InfectionStatus_change(parameters);
    }
    
    // Check if mother and if so increase the maternal immunity sum and total number
    if (population[n].get_m_person_age() > 20 * 365 && population[n].get_m_person_age() < 21 * 365)
    {
      parameters.g_sum_maternal_immunity += population[n].get_m_ICA();
      parameters.g_total_mums++;
    }
    
    // If they are infected, i.e. not S or P, then assign their strains and next strain clearance date
    if (population[n].get_m_infection_state() != Person::SUSCEPTIBLE && population[n].get_m_infection_state() != Person::PROPHYLAXIS) 
    {
      // TODO: Think about how we can correctly initialise MOI for a given EIR. Presumably there is a rarefaction of MOI vs EIR, and the MOI is lognormal*age_dependency
      population[n].set_m_number_of_strains(static_cast<int>(population[n].get_m_individual_biting_rate() * 2) + 1);
      population[n].schedule_m_day_of_strain_clearance(parameters);
      
      // If they are subpatent or treated the strains we allocate should have no state change
      if (population[n].get_m_infection_state() == Person::SUBPATENT  || population[n].get_m_infection_state() == Person::TREATED)
      {
        for (int s = 0; s < population[n].get_m_number_of_strains(); s++)
        {
          temp_strains.push_back( 
            Strain(
              Strain::generate_random_barcode(),
              Strain::m_transition_vector[static_cast<int>(population[n].get_m_infection_state())],
                                         0
            )
          );
        }
      }
      
      // If they are not then they will have state changes assumed equal to the human
      else
      {
        for (int s = 0; s < population[n].get_m_number_of_strains(); s++)
        {
          temp_strains.push_back(
            Strain(
              Strain::generate_random_barcode(),
              Strain::m_transition_vector[static_cast<int>(population[n].get_m_infection_state())],
                                         population[n].get_m_day_of_InfectionStatus_change()
            )
          );
        }
      }
      
      // Set human strains and associated times of strain state changing and clearing
      population[n].set_m_active_strains(temp_strains);
      population[n].set_m_day_of_next_strain_state_change();
      temp_strains.clear();
      infected_human_count.push_back(n);
      population[n].set_m_day_of_next_strain_state_change();
      
    }
    
    // Set the next event day
    population[n].set_m_day_of_next_event();
    
    // Add to the static zeta vector which is required for calculating the overall probability of being bitten, pi
    zeta_vector[n] = population[n].get_m_individual_biting_rate();
  }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: HUMAN INITIALISATION FROM EQUILIBRIUM
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  Rcpp::Rcout << "Human initilisation working\n";
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: MOSQUITO INITIALISATION FROM EQUILIBRIUM
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Mosquito initialisation preallocation
  std::vector<double> mosquito_status_eq{ Sv, Ev, Iv };
  double mosquito_status_eq_sum = mosquito_status_eq[0] + mosquito_status_eq[1] + mosquito_status_eq[2];
  int parent_source = 0;
  std::vector<unsigned short int> pending_oocyst_time{ 0 };
  std::vector<barcode_t> pending_oocyst_barcode_male{ Strain::generate_random_barcode() };
  std::vector<barcode_t> pending_oocyst_barcode_female{ Strain::generate_random_barcode() };
  
  Rcpp::Rcout << "Mosquito preallocation initilisation working\n";
  Rcpp::Rcout << scourge.capacity() << "\n";
  Rcpp::Rcout << infected_human_count.size() << "\n";
  
  // mosquito initialisation
  for (unsigned int n = 0; n < scourge.capacity(); n++)
  {
    // Set their id
    scourge.push_back(Mosquito(parameters));
    scourge[n].set_m_mosquito_ID(n);
    
    // Set infection status
    scourge[n].set_m_mosquito_infection_state(static_cast<Mosquito::InfectionStatus>(sample1(mosquito_status_eq, mosquito_status_eq_sum)));
    
    // Deal with if the mosquito is infected
    // Here we assume they have no pending oocysts but one already burst
    if (scourge[n].get_m_mosquito_infection_state() == Mosquito::INFECTED)
    {
      // Pick a random human for the source of the strain
      parent_source = runiform_int_1(0, infected_human_count.size() - 1);
      
      // If human has only one strain then assign the first strain
      if (population[infected_human_count[parent_source]].get_m_number_of_strains() == 1)
      {
        scourge[n].set_m_oocyst_barcode_male_vector(population[infected_human_count[parent_source]].get_m_person_strain_x(0).get_m_barcode());
        scourge[n].set_m_oocyst_barcode_female_vector(population[infected_human_count[parent_source]].get_m_person_strain_x(0).get_m_barcode());
      }
      // If human has more than one strain then pick random strain for both male and female mosquito barcode
      else
      {
        scourge[n].set_m_oocyst_barcode_male_vector(population[infected_human_count[parent_source]].get_m_person_strain_x(runiform_int_1(0, population[infected_human_count[parent_source]].get_m_number_of_strains() - 1)).get_m_barcode());
        scourge[n].set_m_oocyst_barcode_female_vector(population[infected_human_count[parent_source]].get_m_person_strain_x(runiform_int_1(0, population[infected_human_count[parent_source]].get_m_number_of_strains() - 1)).get_m_barcode());
      }

      // Set the oocyst rupture count to 1
      scourge[n].set_m_ruptured_oocyst_count(1);
      scourge[n].m_mosquito_infected = true;

    }
    
    // Deal with if the mosquito is exposed
    // Here we will assume there is one pending oocyst with a random realisation time of up to delay.mos days for the purpose of initialisation
    if (scourge[n].get_m_mosquito_infection_state() == Mosquito::EXPOSED)
    {
      // Pick a random human for the source of the strain
      parent_source = runiform_int_1(0, infected_human_count.size() - 1);

      // If human has only one strain then assign the first strain
      if (population[infected_human_count[parent_source]].get_m_number_of_strains() == 1)
      {
        pending_oocyst_time[0] = static_cast<unsigned short int>(runiform_int_1(parameters.g_current_time, static_cast<int>(parameters.g_delay_mos)) + 1);
        pending_oocyst_barcode_male[0] = population[infected_human_count[parent_source]].get_m_person_strain_x(0).get_m_barcode();
        pending_oocyst_barcode_female[0] = population[infected_human_count[parent_source]].get_m_person_strain_x(0).get_m_barcode();
        scourge[n].set_m_oocyst_rupture_time_vector(pending_oocyst_time);
        scourge[n].set_m_oocyst_barcode_male_vector(pending_oocyst_barcode_male);
        scourge[n].set_m_oocyst_barcode_female_vector(pending_oocyst_barcode_female);

      }
      // If human has more than one strain then pick random strain for both male and female mosquito barcode
      else
      {

        pending_oocyst_time[0] = runiform_int_1(parameters.g_current_time, static_cast<int>(parameters.g_delay_mos)) + 1;
        pending_oocyst_barcode_male[0] = population[infected_human_count[parent_source]].get_m_person_strain_x(runiform_int_1(0, population[infected_human_count[parent_source]].get_m_number_of_strains() - 1)).get_m_barcode();
        pending_oocyst_barcode_female[0] = population[infected_human_count[parent_source]].get_m_person_strain_x(runiform_int_1(0, population[infected_human_count[parent_source]].get_m_number_of_strains() - 1)).get_m_barcode();
        scourge[n].set_m_oocyst_rupture_time_vector(pending_oocyst_time);
        scourge[n].set_m_oocyst_barcode_male_vector(pending_oocyst_barcode_male);
        scourge[n].set_m_oocyst_barcode_female_vector(pending_oocyst_barcode_female);
      }
      
    }

    // schedule mosquito's next day of event
    scourge[n].schedule_m_day_of_next_event();
    
  }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: MOSQUITO INITIALISATION FROM EQUILIBRIUM
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  Rcpp::Rcout << "Mosquito initilisation working\n";
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: TIMERS AND PRE LOOP
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
  Rcpp::Rcout << "Time elapsed in initialisation: " << duration << " seconds\n";
  
  // Start from day 2
  parameters.g_current_time++;
  
  // End of simulation time
  int g_end_time = static_cast<int>(parameters.g_current_time + (parameters.g_years * 365));
  
  // Preallocations;
  unsigned int num_bites = 0;
  unsigned int increasing_bites = 0;
  unsigned int individual_binomial_bite_draw = 0;
  double psi_sum = 0;
  unsigned int scourge_size = scourge.size();
  
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
  int negative_immunity_check = 0;
  
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
    
    // shuffle the bite queue otherwise you will introduce stepping-stone-esque genetic structuring
    shuffle_integer_vector(mosquito_biting_queue);
    
    // Adjust the number of bites to account for anthrophagy
    num_bites *= parameters.g_Q0;
    
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
    
    // multinomial procedure for everyone and breaking when n_bites is reached
    for (bite_sampling_i = 0; bite_sampling_i < parameters.g_N - 1; bite_sampling_i++)
    {
      
      individual_binomial_bite_draw = rbinomial1(num_bites - increasing_bites, pi_vector[bite_sampling_i] / (pi_sum - pi_cum_sum));
      
      for (bite_sampling_internal_i = 0; bite_sampling_internal_i < individual_binomial_bite_draw; bite_sampling_internal_i++)
      {
        bite_storage_queue.push_back(bite_sampling_i);
      }
      
      pi_cum_sum += pi_vector[bite_sampling_i];
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
        bite_storage_queue.push_back(parameters.g_N - 1);
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
      if (scourge[mosquito_biting_queue[num_bites_i]].m_mosquito_infected) {
        population[bite_storage_queue[num_bites_i]].allocate_bite(parameters, scourge[mosquito_biting_queue[num_bites_i]]);
        if (parameters.g_current_time > g_end_time - 7)
        {
        daily_bite_counters++;
        }
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
        
        if(element.get_m_ICM() < 0 || element.get_m_ICA() < 0 || element.get_m_IB() < 0 || element.get_m_ID() < 0 ){
          negative_immunity_check = 1;
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
  
  Rcpp::Rcout << "Neg imm = " << negative_immunity_check << "\n";
  
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
    Infection_States[element] = static_cast<int>(population[element].get_m_infection_state());
    population[element].update_immunities_to_today(parameters);
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
                                                      Rcpp::Named("IB")=IB,Rcpp::Named("ICA")=ICA,Rcpp::Named("ICM")=ICM,Rcpp::Named("ID")=ID,
                                                        Rcpp::Named("Daily_Bites")=daily_bite_counters/log_counter);
  
  
  // Create universe ptr for memory-continuiation
  Rcpp::XPtr<Universe> universe_ptr(new Universe{ population, psi_vector, zeta_vector, pi_vector, scourge, parameters},
                                    true);
  
  // Return Named List with pointer and loggers
  return Rcpp::List::create(Rcpp::Named("Ptr") = universe_ptr, Rcpp::Named("Loggers")=Loggers);
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // fini
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
}


