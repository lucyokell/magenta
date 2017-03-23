#include "stdafx.h"
#include "person.h"

// Only Constructor
Person::Person(const Parameters &parameters) :
  
  m_individual_biting_rate{ set_initial_m_individual_biting_rate(parameters.g_zeta_meanlog, parameters.g_zeta_sdlog) },
  m_person_age{ set_initial_m_person_age(parameters.g_average_age) },
  m_age_dependent_biting_rate{ set_m_age_dependent_biting_rate(parameters.g_rho, parameters.g_a0) }
  {
    // Allocate death day
    set_initial_m_day_of_death(parameters);
    
    // Reserve storage space for vectors - trialing suggests the below should fit about right in terms of minimising vector resizing
    m_active_strains.reserve(static_cast<int>(100*m_individual_biting_rate));
    m_active_strain_contribution.reserve(static_cast<int>(100*m_individual_biting_rate));
    m_infection_time_realisation_vector.reserve(static_cast<int>(100*m_individual_biting_rate));			// First pending infection time in position 0 to handle multiple infections times that have not been realised yet
    m_infection_state_realisation_vector.reserve(static_cast<int>(100*m_individual_biting_rate));			// First pending infection state in position 0 to handle multiple infections states that have not been realised yet
    m_infection_barcode_realisation_vector.reserve(static_cast<int>(100*m_individual_biting_rate));		// First pending infection barcode in position 0 to handle multiple infections states that have not been realised yet
  }

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// STATIC INITIALIZERS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Initialise static class const vectors
const std::vector<Person::InfectionStatus> Person::m_transition_vector{ DISEASED, TREATED, ASYMPTOMATIC };

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// SEMI-GETTERS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Work out if reciprocal infection happened. 1 if yes, 0 if no
bool Person::reciprocal_infection_boolean(const Parameters &parameters)
{
  
  if (m_cA_counter)
  {
    // double fD = (1 - ((1 - parameters.g_fD0) / (1 + (pow((m_person_age / parameters.g_aD), parameters.g_gD)))));
    // double q = (parameters.g_d1 + ((1 - parameters.g_d1) / (1 + fD *(pow((m_ID / parameters.g_ID0), parameters.g_kD)))));
    m_cA = (parameters.g_cU + (parameters.g_cD - parameters.g_cU) * (pow((parameters.g_d1 + ((1 - parameters.g_d1) / (1 + (1 - ((1 - parameters.g_fD0) / (1 + (pow((m_person_age / parameters.g_aD), parameters.g_gD))))) * (pow((m_ID / parameters.g_ID0), parameters.g_kD))))), parameters.g_gamma1)));
    m_cA_counter = false;
  }
  
  // Match infection state and asses whether onward infection would have happened
  switch (m_infection_state)
  {
  case SUSCEPTIBLE:
    return(false);
  case DISEASED:
    return(rbernoulli_cD());
  case ASYMPTOMATIC:
    /*return(rbinomial1(1,m_cA));*/
    return(rbernoulli1(m_cA));
  case SUBPATENT:
    return(rbernoulli_cU());
  case TREATED:
    return(rbernoulli_cT());
  case PROPHYLAXIS:
    return(false);
  default:
    assert("Schedule Infection Status Change Error - person's infection status not S, D, A, U, T or P");
  return(false);
  }
  
}

// Work out which strain barcodes are passed on as gametocytes
std::vector<barcode_t> Person::sample_two_barcodes(const Parameters &parameters)
{
  // Work out what the contribution for each strain is 
  if (m_contribution_counter == 0) {
    
    for (int n = 0; n < m_number_of_strains; n++)
    {
      // Match infection state and schedule associated next state change
      switch (m_active_strains[n].get_m_strain_infection_status())
      {
      case Strain::DISEASED:
        m_active_strain_contribution.emplace_back(parameters.g_cD);
        break;
      case Strain::ASYMPTOMATIC:
        // We might need to recalculate the contribution from an asymptomatic for the first time here if the person is clinically diseased but is also conifecte with asymptomatic strains
        if (m_cA_counter)
        {
          // double fD = (1 - ((1 - parameters.g_fD0) / (1 + (pow((m_person_age / parameters.g_aD), parameters.g_gD)))));
          // double q = (parameters.g_d1 + ((1 - parameters.g_d1) / (1 + fD *(pow((m_ID / parameters.g_ID0), parameters.g_kD)))));
          m_cA = (parameters.g_cU + (parameters.g_cD - parameters.g_cU) * (pow((parameters.g_d1 + ((1 - parameters.g_d1) / (1 + (1 - ((1 - parameters.g_fD0) / (1 + (pow((m_person_age / parameters.g_aD), parameters.g_gD))))) * (pow((m_ID / parameters.g_ID0), parameters.g_kD))))), parameters.g_gamma1)));
          m_cA_counter = false;
        }
        m_active_strain_contribution.emplace_back(m_cA);
        break;
      case Strain::SUBPATENT:
        m_active_strain_contribution.emplace_back(parameters.g_cU);
        break;
      case Strain::TREATED:
        m_active_strain_contribution.emplace_back(parameters.g_cT);
        break;
      default:
        assert(NULL && "Strain infection status not D, A, U or T");
      break;
      }
    }
    m_contribution_counter = 1;
    m_contribution_sum = std::accumulate(m_active_strain_contribution.begin(), m_active_strain_contribution.end(), 1.0);
  }
  
  // If individual only has one strain then catch for this as we won't need to draw a random sample for what strain is drawn
  if (m_number_of_strains == 1)
  {
    return(std::vector<barcode_t> { m_active_strains[0].get_m_barcode(), m_active_strains[0].get_m_barcode() });
  }
  else
  {
    return(std::vector<barcode_t> { m_active_strains[sample1(m_active_strain_contribution, m_contribution_sum)].get_m_barcode(), m_active_strains[sample1(m_active_strain_contribution, m_contribution_sum)].get_m_barcode() });
  }
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// SETTERS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Set person's individual biting rate
double Person::set_initial_m_individual_biting_rate(double zeta_meanlog, double zeta_sdlog) {
  
  return(rlognorm1(zeta_meanlog, zeta_sdlog));
  
}

// Set person's age
int Person::set_initial_m_person_age(double average_age) {
  
  return(rexpint1(1.0 / average_age));
  
}

// Set person's initial death day
void Person::set_initial_m_day_of_death(const Parameters &parameters)
{
  // Exppnential waiting time plus current day and 1 so not the same day
  m_day_of_death = rexpint1(1.0 / parameters.g_average_age) + parameters.g_current_time + 1;
  
}

// Set person's relative chance of biting (psi)
double Person::set_m_age_dependent_biting_rate(double rho, double a0) {
  
  return(1 - (rho*exp(-m_person_age / a0)));
  
}

// Set day of next strain state change
void Person::set_m_day_of_next_strain_state_change()
{
  
  // First let's reset the more than one strain to change boolean
  m_more_than_one_strain_to_change_today_bool = false;
  
  // Reset the next day to max
  m_day_of_next_strain_state_change = std::numeric_limits<int>::max();
  
  // loop through and keep minimum non zero
  for (int s = 0; s < m_number_of_strains; s++) {
    if (m_active_strains[s].get_m_strain_infection_status() != Strain::TREATED) {
      m_temp_day_of_next_strain_state_change = m_active_strains[s].get_m_day_of_strain_infection_status_change();
      if (m_temp_day_of_next_strain_state_change <= m_day_of_next_strain_state_change && m_temp_day_of_next_strain_state_change != 0)
      {
        
        // if the next strain state change day is the same as the temp flag the bool
        if (m_day_of_next_strain_state_change == m_temp_day_of_next_strain_state_change)
        {
          m_more_than_one_strain_to_change_today_bool = true;
        }
        
        // keep earliest time and strain position 
        m_day_of_next_strain_state_change = m_temp_day_of_next_strain_state_change;
        m_temp_strain_to_next_change = s;
        
      }
    }
  }
}


// Set day of next event
void Person::set_m_day_of_next_event() {
  
  // First if catches if person is completely susceptible, i.e. not with a pending infection as this is the initialisation criteria
  if (m_infection_time_realisation_vector.empty() && m_day_of_strain_clearance == 0 && m_day_of_InfectionStatus_change == 0) {
    m_day_of_next_event = m_day_of_death;
  }
  
  else // find minimum that is not zero
  {
    // Start next event as strain state change
    m_day_of_next_event = m_day_of_death;
    
    // Compare against strain state change
    if (m_day_of_next_event > m_day_of_next_strain_state_change && m_day_of_next_strain_state_change != 0) {
      m_day_of_next_event = m_day_of_next_strain_state_change;
    }
    
    // Catch if there are no pending infections first, i.e. if the number of realized infections is the same as the length of the infection realisation vector
    if (m_infection_time_realisation_vector.size() > m_number_of_realised_infections)
    {
      if (m_day_of_next_event > m_infection_time_realisation_vector[m_number_of_realised_infections] && m_infection_time_realisation_vector[m_number_of_realised_infections] != 0) {
        m_day_of_next_event = m_infection_time_realisation_vector[m_number_of_realised_infections];
      }
    }
    
    // Compare against strain clearance day
    if (m_day_of_next_event > m_day_of_strain_clearance && m_day_of_strain_clearance != 0) {
      m_day_of_next_event = m_day_of_strain_clearance;
    }
    
    // Compare against infection status change day
    if (m_day_of_next_event > m_day_of_InfectionStatus_change && m_day_of_InfectionStatus_change != 0) {
      m_day_of_next_event = m_day_of_InfectionStatus_change;
    }
    
    
  }
  
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// ALLOCATIONS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Allocate bite to person
void Person::allocate_bite(const Parameters &parameters, Mosquito &mosquito)
{
  // If we have already allocated a bite to this indivdual in this time step then we know we won't need to assess the immunity boosting again
  if (!m_number_of_bites)
  {
    
    // If the last bite was less than g_uB (time in which boosting of IB cannot happen) then increse IB
    if (m_IB_last_boost_time < parameters.g_current_time - parameters.g_uB) {
      
      // Increase IB and update IB last boost time
      m_IB++;
      m_IB_last_boost_time = parameters.g_current_time + modf(m_IB_last_boost_time, &m_IB_last_boost_time);
      
    }
    
    // First calculate what the current IB should be given when it was last calculated
    m_IB *= exp((m_IB_last_calculated_time - parameters.g_current_time) / parameters.g_dB);
    
    // Then update the biting success rate 
    m_biting_success_rate = parameters.g_b0 * (parameters.g_b1 + ((1 - parameters.g_b1) / (1 + (pow((m_IB / parameters.g_IB0), parameters.g_kB)))));
    
    // Update last calculated time
    m_IB_last_calculated_time = parameters.g_current_time;
  }
  
  // Work out if the bite has led to an infection
  // Firstly if the human is treated or in prophylaxis then cannot be infected
  if (m_infection_state != TREATED || m_infection_state != PROPHYLAXIS)
  {
    
    // Random draw to see if the bite led to an infection
    if (rbernoulli1(m_biting_success_rate)) {
      
      // Allocate infection
      allocate_infection(parameters, mosquito);
    }
  }

  // Increase number of bites 
  m_number_of_bites++;
  
  }

// Allocate an infection to person, i.e. individuals who return >0 from allocate_force_of_infection()
void Person::allocate_infection(const Parameters &parameters, Mosquito &mosquito)
{
  // Only need to calculate this once per day step
  if (!m_number_of_succesful_bites)
  {
    
    // First calculate what the current immunities should be given when it was last calculated
    m_ID *= exp((m_I_C_D_CM_last_calculated_time - parameters.g_current_time) / parameters.g_dB);
    m_ICA *= exp((m_I_C_D_CM_last_calculated_time - parameters.g_current_time) / parameters.g_dCA);
    m_ICM = m_ICM_init * exp(-m_person_age / parameters.g_dCM);
    
    // Update last calculated time
    m_I_C_D_CM_last_calculated_time = parameters.g_current_time;
    
    // If the last bite was less than g_uD (time in which boosting of ID cannot happen) then increse ID
    if (m_ID_last_boost_time < parameters.g_current_time - parameters.g_uD) {
      
      // Increase ID and update ID last boost time
      m_ID++;
      m_ID_last_boost_time = parameters.g_current_time + modf(m_ID_last_boost_time, &m_ID_last_boost_time);
      
    }
    
    // If the last bite was less than g_uCA (time in which boosting of ICA cannot happen) then increse ICA
    if (m_ICA_last_boost_time < parameters.g_current_time - parameters.g_uCA) {
      
      // Increase IB and update IB last boost time
      m_ICA++;
      m_ICA_last_boost_time = parameters.g_current_time + modf(m_ICA_last_boost_time, &m_ICA_last_boost_time);
      
    }
    
    // Calculate a symptom success rate (phi) dependent on age and acquired immunity
    m_symptom_success_rate = parameters.g_phi0 * (parameters.g_phi1 + ((1 - parameters.g_phi1) / (1 + (pow(((m_ICA + m_ICM) / parameters.g_IC0), parameters.g_kC)))));
    
    // Calculate individuals InfectionStatus change 
    
    // Set outcome member probabilities
    // N -> D
    m_transition_probabilities[0] = m_symptom_success_rate * (1 - parameters.g_ft);
    // N -> T
    m_transition_probabilities[1] = m_symptom_success_rate *parameters.g_ft;
    // N -> A
    m_transition_probabilities[2] = 1 - m_symptom_success_rate;
    
    // Need to find those who are to be infected again who are D as they should not then go to being
    // asymptomatic trough the additonal infection. Going to T however is still fair.
    if (m_infection_state == DISEASED) {
      m_transition_probabilities[2] = 0;
    }
    
    // Set outcome probability sum 
    m_sum_transition_probabilities = m_transition_probabilities[0] + m_transition_probabilities[1] + m_transition_probabilities[2];
    
  }
  
  // Increase number of succesful bites
  m_number_of_succesful_bites++;
  
  //TODO: Add here a rpois that details the possibility that more than one stain is pushed - if there is you simply loop through the following three lines
  
  // Push the resultant state of infection
  m_infection_state_realisation_vector.emplace_back(m_transition_vector[sample1(m_transition_probabilities, m_sum_transition_probabilities)]);
  
  // Push the resultant state change time
  m_infection_time_realisation_vector.emplace_back(static_cast<int>(parameters.g_dur_E + parameters.g_current_time));
  
  // Allocate strains from mosquito
  if (mosquito.get_m_ruptured_oocyst_count() == 1)
  {
    m_infection_barcode_realisation_vector.emplace_back(
      Strain::generate_recombinant_barcode(
        mosquito.get_m_oocyst_barcode_male_vector(0),
        mosquito.get_m_oocyst_barcode_female_vector(0)
      )
    );
  }
  else
  {
    m_infection_barcode_realisation_vector.emplace_back(
      Strain::generate_recombinant_barcode(
        mosquito.get_m_oocyst_barcode_male_vector(runiform_int_1(1, mosquito.get_m_ruptured_oocyst_count()) - 1),
        mosquito.get_m_oocyst_barcode_female_vector(runiform_int_1(1, mosquito.get_m_ruptured_oocyst_count()) - 1)
      )
    );
  }
  
  // Set next event date as may have changed as a result of the bite
  set_m_day_of_next_event();
  
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// SCHEDULERS 
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Schedule person's infection state change
void Person::schedule_m_day_of_InfectionStatus_change(const Parameters &parameters)
{
  // Match infection state and schedule associated next state change
  switch (m_infection_state)
  {
  case DISEASED:
    m_day_of_InfectionStatus_change = rexpint1(1.0 / parameters.g_dur_D) + parameters.g_current_time + 1;
    break;
  case ASYMPTOMATIC:
    m_day_of_InfectionStatus_change = rexpint1(1.0 / parameters.g_dur_A) + parameters.g_current_time + 1;
    break;
  case SUBPATENT:
    m_day_of_InfectionStatus_change = rexpint1(1.0 / parameters.g_dur_U) + parameters.g_current_time + 1;
    break;
  case TREATED:
    m_day_of_InfectionStatus_change = rexpint1(1.0 / parameters.g_dur_T) + parameters.g_current_time + 1;
    break;
  case PROPHYLAXIS:
    m_day_of_InfectionStatus_change = rexpint1(1.0 / parameters.g_dur_P) + parameters.g_current_time + 1;
    break;
  default:
    assert(NULL && "Schedule Infection Status Change Error - person's infection status not D, A, U, T or P");
  break;
  }
  
}

// Schedule person's death day
void Person::schedule_m_day_of_death(const Parameters &parameters)
{
  // Throw if called on someone who is not SUSCEPTIBLE
  // assert(m_infection_state == SUSCEPTIBLE && "Death day schedule called for someone who is not susceptible");
  
  // Exppnential waiting time plus current day and 1 so not the same day
  m_day_of_death = rexpint1(1.0 / parameters.g_average_age) + parameters.g_current_time + 1;
}

// Schedule person's next strain clearance
void Person::schedule_m_day_of_strain_clearance(const Parameters &parameters)
{
  // Throw if this is called on someone with no strains
  assert(m_number_of_strains != 0 && "Tried to schedule a strain clearance day from an individual with no strains");
  
  // Work out potential new clearance day
  int possible_clearance_day = rexpint1(m_number_of_strains / parameters.g_dur_AU) + parameters.g_current_time + 1;
  
  // If the new clearance day is earlier then schedule it
  if (possible_clearance_day < m_day_of_strain_clearance || !m_day_of_strain_clearance || m_day_of_strain_clearance == parameters.g_current_time) {
    m_day_of_strain_clearance = possible_clearance_day;
  }
  
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// UPDATERS - functions that happen at the end of a day, i.e. ageing, dying, clearing strains etc. 
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Clear one strain from an individual
void Person::individual_strain_clearance()
{
  
  if (m_number_of_strains > 1)
  {
    
    // Draw a random number from 0 to number of strains
    m_temp_strain_to_be_deleted = runiform_int_1(0, m_number_of_strains - 1);
    
    // Swap the strain pointer and strain acquisition date at that position to the back
    std::swap(m_active_strains[m_temp_strain_to_be_deleted], m_active_strains.back());
    
    // Pop thus deleting the random strain pointer
    m_active_strains.pop_back();
    
    // Lower strain counter
    m_number_of_strains--;
    
    // Catch errors this may now cause in the next strain state change etc
    set_m_day_of_next_strain_state_change();
    
  }
  else
  {
    // Do nothing - if they are infected they should never clear all strains otherwise mosquitos will not pick up strains
  }
  
  // Throw if this is called on someone with no strains
  assert(m_number_of_strains != 0 && "Tried to clear a strain from an individual with no strains");
  
  // Set new day of next event
  set_m_day_of_next_event();
  
}

// Clear all strains from an individual
void Person::all_strain_clearance() {
  
  // Clear all strains
  
  // Throw if this is called on someone with no strains
  //assert(m_number_of_strains != 0 && "Tried to clear all strains from an individual with no strains");
  m_active_strains.clear();
  m_active_strain_contribution.clear();
  
  // Update numbers of strains and clearance dates
  m_number_of_strains = 0;
  m_number_of_realised_infections = 0;
  m_day_of_strain_clearance = 0;
  m_day_of_next_strain_state_change = std::numeric_limits<int>::max();
  m_more_than_one_strain_to_change_today_bool = false;
  
  // Clear all associated vectors
  m_infection_time_realisation_vector.clear();
  m_infection_state_realisation_vector.clear();
  m_infection_barcode_realisation_vector.clear();

  // Set new day of next event
  set_m_day_of_next_event();
  
}

// Kill person, i.e. reset age to 0, infections to 0, state to susceptible, immunities reset etc
void Person::die(const Parameters &parameters)
{
  // recover first
  recover(parameters);
  
  // Reset immunities, age
  m_IB = m_ID = m_ICA = 0;
  m_person_age = 0;
  
  // set their boost times back enough so that they are guaranteed to boost if bitten tomorrow
  m_IB_last_boost_time -= parameters.g_uB;
  m_ICA_last_boost_time -= parameters.g_uCA;
  m_ID_last_boost_time -= parameters.g_uD;
  
  // Set maternal immunity
  m_ICM = m_ICM_init = parameters.g_PM * parameters.g_mean_maternal_immunity;
  
  // Schedule death
  schedule_m_day_of_death(parameters);
  
  // Set new day of next event
  set_m_day_of_next_event();
}

// Recover to being susceptible, i.e. clearing all infections and strains and associated timings
void Person::recover(const Parameters &parameters)
{
  
  // Return to susceptible
  m_infection_state = SUSCEPTIBLE;
  
  // Clear strains
  if (m_number_of_strains > 0) {
    all_strain_clearance();
  }
  else 
  {
  // Clear waiting infection vectors
  m_infection_time_realisation_vector.clear();
  m_infection_state_realisation_vector.clear();
  m_infection_barcode_realisation_vector.clear();
  }
  
  // Reset other events to 0
  m_day_of_InfectionStatus_change = 0;
  
  // Update next event counter - due to nature of recovering the next event has to be their death
  m_day_of_next_event = m_day_of_death;
  
}

// Event handle, i.e. if any of person's death, state change, strain clearance or state change realisation days are today, respond accordingly
void Person::event_handle(const Parameters &parameters) {
  
  // If death occurs then no need to explore the other events
  if (m_day_of_death == m_day_of_next_event || m_person_age == parameters.g_max_age) {
    die(parameters);
  }
  else
  {
    // All other events could happen theoretically on the same day though so within same else block
    // Clear strain if time to do so
    if (m_day_of_strain_clearance == m_day_of_next_event) {
      
      individual_strain_clearance(); // clear strain
      if (m_number_of_strains > 1) {
        schedule_m_day_of_strain_clearance(parameters); // Schedule the next strain clearance
      }
      else {
        m_day_of_strain_clearance = 0; // If the strain cleared was the second last strain then there will be no clearnance date as we don't want to remove the last strain
      }
    }
    
    // Change Infection Status due to recoveries etc if time to do so
    if (m_day_of_InfectionStatus_change == m_day_of_next_event) {
      
      switch (m_infection_state)
      {
      case DISEASED:
        m_infection_state = ASYMPTOMATIC;
        schedule_m_day_of_InfectionStatus_change(parameters); // schedule next state change
        break;
      case ASYMPTOMATIC:
        m_infection_state = SUBPATENT;
        schedule_m_day_of_InfectionStatus_change(parameters); // schedule next state change
        break;
      case SUBPATENT:
        recover(parameters); // fully recover to susceptible
        break;
      case TREATED:
        m_infection_state = PROPHYLAXIS;
        schedule_m_day_of_InfectionStatus_change(parameters); // schedule next state change
        all_strain_clearance(); // When they are in prophylaxis we remove all the strains, as treated individuals still have strains
        break;
      case PROPHYLAXIS:
        recover(parameters); // fully recover to susceptible
        break;
      default:
        assert(NULL && "Event Infection Status Change Error - person's infection status not D, A, U, T or P");
      break;
      }
      
    }
    
    // Change strain state if time to do so
    if (m_day_of_next_strain_state_change == m_day_of_next_event) {
      
      // if we have flagged that more than one strain is changing strain state today then loop through all the strains
      if (m_more_than_one_strain_to_change_today_bool)
      {
        for (auto &strain : m_active_strains)
        {
          // if the strain is changing today then switch it accordingly
          if (strain.get_m_day_of_strain_infection_status_change() == parameters.g_current_time)
          {
            // We might end up looping through all the strains when a new treated strain is given to a human which has the same strain state change
            // day as a previous strain - in this case just pass over the treated strain as ultimatel
            if (strain.get_m_strain_infection_status() != Strain::TREATED) 
            {
              switch (strain.get_m_strain_infection_status())
              {
              case Strain::DISEASED:
                // If strain to change state is diseased it will become asymptomatic and its following state change must be the same as the humans
                strain.set_m_strain_infection_status(Strain::ASYMPTOMATIC);
                strain.set_m_day_of_strain_infection_status_change(m_day_of_InfectionStatus_change);
                break;
              case Strain::ASYMPTOMATIC:
                // If strain to change state is asymptomatic it will become subpatent and its following state change will not occur
                strain.set_m_strain_infection_status(Strain::SUBPATENT);
                strain.set_m_day_of_strain_infection_status_change(0);
                break;
              default:
                assert(NULL && "Strain state change equested on strain that is not diseasod of asymptomatic");
              break;
              }
            } 
          }
        }
        set_m_day_of_next_strain_state_change();
      }
      // if we have not flagged that there are multiple strains changing states today then just change the necessary strain
      else
      {
        switch (m_active_strains[m_temp_strain_to_next_change].get_m_strain_infection_status())
        {
        case Strain::DISEASED:
          // If strain to change state is diseased it will become asymptomatic and its following state change must be the same as the humans
          m_active_strains[m_temp_strain_to_next_change].set_m_strain_infection_status(Strain::ASYMPTOMATIC);
          m_active_strains[m_temp_strain_to_next_change].set_m_day_of_strain_infection_status_change(m_day_of_InfectionStatus_change);
          break;
        case Strain::ASYMPTOMATIC:
          // If strain to change state is asymptomatic it will become subpatent and its following state change will not occur
          m_active_strains[m_temp_strain_to_next_change].set_m_strain_infection_status(Strain::SUBPATENT);
          m_active_strains[m_temp_strain_to_next_change].set_m_day_of_strain_infection_status_change(0);
          break;
        default:
          assert(NULL && "Strain state change equested on strain that is not diseasod of asymptomatic");
        break;
        }
        set_m_day_of_next_strain_state_change();
      }
      
    }
    
    // Realise infection if time to do so
    // First check if there are any pending infections
    if (m_infection_time_realisation_vector.size() > m_number_of_realised_infections)
    {
      // Is the next pending infection for today
      if (m_infection_time_realisation_vector[m_number_of_realised_infections] == m_day_of_next_event)
      {
        // Reset this catch variable which ensiures multiple infections on the same day can be realised
        m_infection_realisation_empty_catch = 1;
        
        // Set infection state to first infection state and then clear that state realisation and time
        while (m_infection_realisation_empty_catch == 1)
        {
          // Assign infection state
          
          // If you are diseased already and the next infection would make you asymptomatic then remain diseased
          // If it were to make you treated then this is plausible
          
          if (m_infection_state == DISEASED)
          {
            if (m_infection_state_realisation_vector[m_number_of_realised_infections] == TREATED) {
              m_infection_state = TREATED;
            }
          }
          // If you are not diseased you can go to any of treated, diseased, and remain asymptomatic
          else  
          {
            m_infection_state = m_infection_state_realisation_vector[m_number_of_realised_infections];
          }
          
          // If infection state is treated then we clear all pending infections
          if (m_infection_state == TREATED)
          {
            
            // schedule next state change
            schedule_m_day_of_InfectionStatus_change(parameters); 
            
            // Push the pending strain
            m_active_strains.emplace_back(
              m_infection_barcode_realisation_vector[m_number_of_realised_infections],
              Strain::m_transition_vector[static_cast<int>(m_infection_state)],
              m_day_of_InfectionStatus_change
            );
            
            // We don't care if the new strain will state change earlier than the next state change as this would cause recovery
            
            // Increase the number of strains and realised infections
            m_number_of_strains++;
            m_number_of_realised_infections++;
            
            // Schedule new strain clearance day if there is more than one strain
            if (m_number_of_strains > 1)
            {
              schedule_m_day_of_strain_clearance(parameters);
            }

            // since they are being treated
            m_infection_time_realisation_vector.clear();
            m_infection_state_realisation_vector.clear();
            m_infection_barcode_realisation_vector.clear();
            
          }
          // Otherwise pop the time and state and schedule state change
          else
          {
            // Schedule new infection state change and pop the time and state queue
            schedule_m_day_of_InfectionStatus_change(parameters);

            // Push next strain
            m_active_strains.emplace_back(
              m_infection_barcode_realisation_vector[m_number_of_realised_infections],
              Strain::m_transition_vector[static_cast<int>(m_infection_state)],
              m_day_of_InfectionStatus_change
            );
            
            // If the new strain will change state earlier than anyother strain then update the state change event day
            if (m_day_of_InfectionStatus_change <= m_day_of_next_strain_state_change) {
              
              // if the next strain state change day is the same as the temp flag the bool
              if (m_day_of_next_strain_state_change == m_day_of_InfectionStatus_change)
              {
                m_more_than_one_strain_to_change_today_bool = true;
              }
              
              // Update the day of next strain state change and the position to the most recent strain i.e the end
              m_day_of_next_strain_state_change = m_day_of_InfectionStatus_change;
              m_temp_strain_to_next_change = m_number_of_strains;
              
            }
            
            // Increase the number of strains
            m_number_of_strains++;
            m_number_of_realised_infections++;
            
            // Schedule new strain clearance day if there is more than one strain
            if (m_number_of_strains > 1)
            {
              schedule_m_day_of_strain_clearance(parameters);
            }
            
          }
          
          // Little conditional loop to assess if there are pending infections still for today and if not to change the catch
          if (m_infection_time_realisation_vector.size() > m_number_of_realised_infections)
          {
            
            if (m_infection_time_realisation_vector[m_number_of_realised_infections] != parameters.g_current_time)
            {
              m_infection_realisation_empty_catch = 0;
            }
          }
          else { m_infection_realisation_empty_catch = 0; }
          
        }
        
      }
    }
    
    // Update the day of next event
    set_m_day_of_next_event();
  }
  
}

// Daily update function, i.e. increase ages, set maternal booleans, age dependent biting rate etc
double Person::update(Parameters &parameters)
{
  // Update age
  m_person_age++;
  m_age_dependent_biting_rate = 1 - (parameters.g_rho*exp(-m_person_age / parameters.g_a0));
  
  // Find out if pregnacy age
  if (m_person_age > (15 * 365) && m_person_age < (25 * 365)) {
    parameters.g_sum_maternal_immunity += m_ICA;
    parameters.g_total_mums++;
  }
  
  // Handle any events that are happening today
  if (m_day_of_next_event == parameters.g_current_time || m_person_age == parameters.g_max_age) {
    event_handle(parameters);
  }
  
  // Reset total bites and other counters designed to prevent unnecessary recalculations in cases of successive same day bites on an inidividual
  m_number_of_bites = 0;
  m_number_of_succesful_bites = 0;
  m_contribution_counter = 0;
  m_cA_counter = true;
  m_active_strain_contribution.clear();
  
  // Throw if the next event date is still today
  assert(m_day_of_next_event != parameters.g_current_time &&
    "Update function failed to handle event update as next event is still today");
  
  // Throw if the next event date is still today
  assert(m_day_of_strain_clearance != parameters.g_current_time &&
    "Update function failed to handle strain");
  
  // Throw if the next event date is still today
  assert(m_day_of_death != parameters.g_current_time &&
    "Update function failed to handle death");
  
  // Throw if the next event date is still today
  assert(m_day_of_next_strain_state_change != parameters.g_current_time &&
    "Update function failed to handle strain state change");
  
  // Throw if the next event date is still today
  assert(m_day_of_InfectionStatus_change != parameters.g_current_time &&
    "Update function failed to handle infection status change");
  
  // Throw if the next event date is less than today
  assert((m_day_of_next_event > parameters.g_current_time || m_day_of_next_event == 0) &&
    "Update function failed to handle infection status change");
  
  // Throw if the next event date is less than today
  assert((m_day_of_InfectionStatus_change > parameters.g_current_time || m_day_of_InfectionStatus_change == 0) &&
    "Update function failed to handle infection status change");
  
  // Throw if the next event date is less than today
  assert((m_day_of_strain_clearance > parameters.g_current_time || m_day_of_strain_clearance == 0) &&
    "Update function failed to handle infection status change");
  
  // Throw if the next event date is less than today
  assert((m_day_of_next_strain_state_change > parameters.g_current_time || m_day_of_next_strain_state_change == 0) &&
    "Update function failed to handle infection status change");
  
  return(m_age_dependent_biting_rate);
  
}


// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// LOGGERS - functions that report summaries on an individual
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int Person::log_daily_incidence(const Parameters &parameters) {
  
  // Are they considered to be a new case, i.e. S, A, U
  if (m_infection_state == SUSCEPTIBLE || m_infection_state == ASYMPTOMATIC || m_infection_state == SUBPATENT)
  {
    // Do they have any infections pending
    if (!m_infection_time_realisation_vector.empty())
    {
      // Is that infection pending for tomorrow (thuis we assume that this function is always called at the end of a day before the next day starts)
      if (m_infection_time_realisation_vector[m_number_of_realised_infections] == parameters.g_current_time + 1)
      {
        // Is that infection yielding a clinical case
        if (m_infection_state_realisation_vector[m_number_of_realised_infections] == TREATED || m_infection_state_realisation_vector[m_number_of_realised_infections] == DISEASED)
        {
          
          // Are they aged 0 - 5 years
          if (m_person_age <= (5 * 365))
          {
            return(2);
          }
          else
          {
            return(1);
          }
        }
      }
    }
  }
  
  return(0);
}

void Person::update_immunities_to_today(const Parameters &parameters) {
  
  m_IB *= exp((m_IB_last_calculated_time - parameters.g_current_time) / parameters.g_dB);
  // Update last calculated time
  m_IB_last_calculated_time = parameters.g_current_time;
  // First calculate what the current immunities should be given when it was last calculated
  m_ID *= exp((m_I_C_D_CM_last_calculated_time - parameters.g_current_time) / parameters.g_dB);
  m_ICA *= exp((m_I_C_D_CM_last_calculated_time - parameters.g_current_time) / parameters.g_dCA);
  m_ICM = m_ICM_init * exp(-m_person_age / parameters.g_dCM);
  
  // Update last calculated time
  m_I_C_D_CM_last_calculated_time = parameters.g_current_time;
}



