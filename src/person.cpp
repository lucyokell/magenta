#include "stdafx.h"
#include "person.h"
#include "probability.h"
#include <algorithm>
#include <numeric>
#include <cassert> // for error checking
#include "parameters.h"

// Declare the already initialised global paramters object
Parameters parameters;

// Default Constructor
Person::Person() :
  m_person_ID{ Person::s_person_ID_generator++ },
  m_individual_biting_rate{ set_m_individual_biting_rate(parameters.g_zeta_meanlog, parameters.g_zeta_sdlog) },
  m_person_age{ set_initial_m_person_age(parameters.g_average_age) },
  m_age_dependent_biting_rate{ set_m_age_dependent_biting_rate(parameters.g_rho, parameters.g_a0) }
  {
    // Allocate death day
    set_initial_m_day_of_death(parameters);
    
    // Increase static sum of psi required for normalisiing age dependent biting rates
    s_psi_sum += m_age_dependent_biting_rate;
    
    // Reserve storage space for vectors
    // TODO: Turn these into an input number
    m_active_strain_pointers.reserve(100);
    m_active_strain_acquisition_day.reserve(100);
  }

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// SETTERS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Set person's individual biting rate
double Person::set_m_individual_biting_rate(double zeta_meanlog, double zeta_sdlog) {
  
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

// Set day of next event
void Person::set_m_day_of_next_event() {
  
  // First if catches if person is completely susceptible, i.e. not with a pending infection as this is the initialisation criteria
  if (m_infection_time_realisation_queue.empty() && !m_day_of_strain_clearance && !m_day_of_InfectionStatus_change) {
    m_day_of_next_event = m_day_of_death;
  }
  else // find minimum that is not zero
  {
    // Start next event as death
    m_day_of_next_event = m_day_of_death;
    
    // Catch if queue is empty first
    if (!m_infection_time_realisation_queue.empty())
    {
      
      if (m_day_of_next_event > m_infection_time_realisation_queue.front() && m_infection_time_realisation_queue.front() > 0) {
        m_day_of_next_event = m_infection_time_realisation_queue.front();
      }
      
    }
    // Compare against strain clearance day
    if (m_day_of_next_event > m_day_of_strain_clearance && m_day_of_strain_clearance > 0) {
      m_day_of_next_event = m_day_of_strain_clearance;
    }
    
    // Compare against infection status change day
    if (m_day_of_next_event > m_day_of_InfectionStatus_change && m_day_of_InfectionStatus_change > 0) {
      m_day_of_next_event = m_day_of_InfectionStatus_change;
    }
  }
  
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// ALLOCATIONS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Allocate bite to person
// TODO: Give a mosquito pointer variable to this function
void Person::allocate_bite(const Parameters &parameters)
{
  // If we have already allocated a bite to this indivdual in this time step then we know we won't need to assess the immunity boosting again
  if (!m_number_of_bites) 
  {
    
    // If the last bite was less than g_uB (time in which boosting of IB cannot happen) then increse IB
    if (m_IB_last_boost_time < parameters.g_current_time - parameters.g_uB - m_immunity_boost_float) {
      
      // Increase IB and update IB last boost time
      m_IB++;
      m_IB_last_boost_time = parameters.g_current_time;
      
    }
    
    // First calculate what the current IB should be given when it was last calculated
    m_IB *= exp((m_IB_last_calculated_time - parameters.g_current_time) / parameters.g_dB);
    
    // Then update the biting success rate 
    m_biting_success_rate = parameters.g_b0 * (parameters.g_b1 + ((1 - parameters.g_b1) / (1 + (pow((m_IB / parameters.g_IB0), parameters.g_kB)))));
    
    // Update last calculated time
    m_IB_last_calculated_time = parameters.g_current_time;
  }
  
  // Work out if the bite has led to an infection
  if (rbinomial1(1, m_biting_success_rate)) {
    // Allocate mosquito pointer 
    // m_mosquito_pointers[m_number_of_bites] = mosquito_pointer;
    
    // Allocate infection
    allocate_infection(parameters);
  }
  
  // Increase number of bites 
  m_number_of_bites++;
  
  
}

// Allocate an infection to person, i.e. individuals who return >0 from allocate_force_of_infection()
void Person::allocate_infection(const Parameters &parameters)
{
  // Only need to calculate this once per day step
  if (!m_number_of_succesful_bites) 
  {
    
    // First calculate what the current immunities should be given when it was last calculated
    m_ID *= exp((m_I_C_D_CM_last_calculated_time - parameters.g_current_time) / parameters.g_dB);
    m_ICA *= exp((m_I_C_D_CM_last_calculated_time - parameters.g_current_time) / parameters.g_dCA);
    m_ICM *= exp((m_I_C_D_CM_last_calculated_time - parameters.g_current_time) / parameters.g_dCM);
    
    // Update last calculated time
    m_I_C_D_CM_last_calculated_time = parameters.g_current_time;
    
    // If the last bite was less than g_uD (time in which boosting of ID cannot happen) then increse ID
    if (m_ID_last_boost_time < parameters.g_current_time - parameters.g_uD - m_immunity_boost_float) {
      
      // Increase ID and update ID last boost time
      m_ID++;
      m_ID_last_boost_time = parameters.g_current_time;
      
    }
    
    // If the last bite was less than g_uCA (time in which boosting of ICA cannot happen) then increse ICA
    if (m_ICA_last_boost_time < parameters.g_current_time - parameters.g_uCA - m_immunity_boost_float) {
      
      // Increase IB and update IB last boost time
      m_ICA++;
      m_ICA_last_boost_time = parameters.g_current_time;
      
    }
    
    // Calculate a symptom success rate (phi) dependent on age and acquired immunity
    m_symptom_success_rate = parameters.g_phi0 * (parameters.g_phi1 + ((1 - parameters.g_phi1) / (1 + (pow(((m_ICA + m_ICM) / parameters.g_IC0), parameters.g_kC)))));
    
    // Calculate individuals InfectionStatus change 
    
    // TODO: So although i've caught cases where more than 1 bite may occur each day, I've not yet factored this into the state change
    // i.e. if 1st bite leads to treated, but a second bite could then cause them to go back to being diseased. Is this fair?
    
    // Set outcome member probabilities
    // N -> D
    m_transition_probabilities[0] = m_symptom_success_rate * (1 - parameters.g_ft);
    // N -> T
    m_transition_probabilities[1] = m_symptom_success_rate *parameters.g_ft;
    // N -> A
    m_transition_probabilities[2] = 1 - m_symptom_success_rate;
    
    // Need to find those who are to be infected again who are D as they should not then go to being
    // asymptomatic trough the additonal infection. Going to T however seems plausible.
    if (m_infection_state == DISEASED) {
      m_transition_probabilities[2] = 0;
    }
    
    // Set outcome probability sum 
    m_sum_transition_probabilities = m_transition_probabilities[0] + m_transition_probabilities[1] + m_transition_probabilities[2];
    
  }
  
  // Increase number of succesful bites
  m_number_of_succesful_bites++;
  
  // Push the resultant state of infection
  m_infection_state_realisation_queue.push(m_transition_vector[sample1(m_transition_probabilities, m_sum_transition_probabilities)]);
  
  // Push the resultant state change time
  // TODO: Chat if this is fair, otherwise this will have to become a different container that is sorted. 
  m_infection_time_realisation_queue.push(parameters.g_dur_E + parameters.g_current_time);
  
  // Allocate strains according to which mosquito bites were succesful
  // TODO: Pass in bite counter element from main so that we know which mosquito it is.
  // TODO: Allocate strain through pointer to mosquito barcode
  // m_infection_strain_realisation_queue.push(m_mosquito_pointers[bite counter element]->getStrain());
  
  // Set next event date as may have changed as a result of the bite
  set_m_day_of_next_event();
  
}

// Allocate a strain to individual
void Person::allocate_strain(const Parameters &parameters, Strain* strain_ptr)
{
  
  // TODO: Add strain pointer here
  
  // Increase the number of strains
  m_number_of_strains++;
  
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
  assert(m_infection_state == SUSCEPTIBLE && "Death day schedule called for someone who is not susceptible");
  
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
  /*
   TODO: Add this in when the strains fully come in
   
   if(m_number_of_strains > 1){
   
   // Draw a random number from 0 to number of strains
   int strain_to_be_deleted = runiform_int_1(0, m_number_of_strains);
   
   // Swap the strain pointer and strain acquisition date at that position to the back
   std::swap(m_active_strain_pointers[strain_to_be_deleted], m_active_strain_pointers.back());
   std::swap(m_active_strain_acquisition_day[strain_to_be_deleted], m_active_strain_acquisition_day.back());
   
   // Pop thus deleting the random strain pointer
   m_active_strain_pointers.pop_back();
   m_active_strain_acquisition_day.pop_back();
   
   } else {
   m_active_strain_pointers.clear();
   m_active_strain_acquisition_day.clear();
   }
   */
  
  // Throw if this is called on someone with no strains
  assert(m_number_of_strains != 0 && "Tried to clear a strain from an individual with no strains");
  m_number_of_strains--;
  
  // Set new day of next event
  set_m_day_of_next_event();
  
}

// Clear all strains from an individual
void Person::all_strain_clearance() {
  
  // Clear all strains
  // catch incase you naturally cleared strains before recovering
  // TODO: Check if this is something that should be overridden, i.e. if you have no strains and are still diseased you recover
  
  // Throw if this is called on someone with no strains
  assert(m_number_of_strains != 0 && "Tried to clear all strains from an individual with no strains");
  m_active_strain_pointers.clear();
  m_active_strain_acquisition_day.clear();
  
  // Update numbers of strains and clearance dates
  m_number_of_strains = 0;
  m_day_of_strain_clearance = 0;
  
  // Set new day of next event
  set_m_day_of_next_event();
}

// Kill person, i.e. reset age to 0, infections to 0, state to susceptible, immunities reset etc
void Person::die(const Parameters &parameters)
{
  // recover first
  recover(parameters);
  
  // Reset immunities, strains, age, boost times
  m_IB = m_ID = m_ICA = m_number_of_strains = m_person_age = m_ID_last_boost_time = m_IB_last_boost_time = m_ICA_last_boost_time = 0;
  
  // Set maternal immunity
  m_ICM = parameters.g_PM * s_mean_maternal_immunity * m_individual_biting_rate;
  
}

// Recover to being susceptible, i.e. clearing all infections and strains and associated timings
void Person::recover(const Parameters &parameter)
{
  
  // Return to susceptible
  m_infection_state = SUSCEPTIBLE;
  
  // Clear strains
  if (m_number_of_strains > 0) {
    all_strain_clearance();
  }
  
  // Clear waiting infection queue
  // TODO: Check if this queue is the correct clearing process
  // Apparently easiest way to clear a queue::
  m_infection_time_realisation_queue = std::queue<int>();
  m_infection_state_realisation_queue = std::queue<InfectionStatus>();
  m_infection_strain_realisation_queue = std::queue<Strain*>();
  
  // Clear strain vectors
  m_active_strain_pointers.clear();
  m_active_strain_acquisition_day.clear();
  
  // Reset other events to 0
  m_day_of_strain_clearance = 0;
  m_day_of_InfectionStatus_change = 0;
  
  // Schedule death
  schedule_m_day_of_death(parameters);
  
  // Update next event counter - due to nature of recovering the next event has to be their death
  m_day_of_next_event = m_day_of_death;
  
}

// Event handle, i.e. if any of person's death, state change, strain clearance or state change realisation days are today, respond accordingly
void Person::event_handle(const Parameters &parameters) {
  
  // If death occurs then no need to explore the other events
  if (m_day_of_death == m_day_of_next_event) {
    die(parameters);
  }
  else 
  {
    // All other events could happen theoretically on the same day though so within same else block
    // Clear strain if time to do so
    if (m_day_of_strain_clearance == m_day_of_next_event) {
      
      individual_strain_clearance(); // clear strain
      if (m_number_of_strains > 0) {
        schedule_m_day_of_strain_clearance(parameters); // Schedule the next strain clearance
      }
      else {
        m_day_of_strain_clearance = 0; // If the strain cleared was the last strain then there will be no clearnance date
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
        recover(parameters);
        break;
      case TREATED:
        m_infection_state = PROPHYLAXIS;
        schedule_m_day_of_InfectionStatus_change(parameters); // schedule next state change
        m_infection_time_realisation_queue = std::queue<int>(); // remove pending infections as assumend to be cleared with treatment
        m_infection_state_realisation_queue = std::queue<InfectionStatus>(); // TODO: Is this right biologically
        m_infection_strain_realisation_queue = std::queue<Strain*>();
        break;
      case PROPHYLAXIS:
        recover(parameters);
        break;
      default:
        assert(NULL && "Event Infection Status Change Error - person's infection status not D, A, U, T or P");
      break;
      }
      
    }
    
    // Realise infection if time to do so
    // Catch if there is no queue
    if (!m_infection_time_realisation_queue.empty()) 
    {
      if (m_infection_time_realisation_queue.front() == m_day_of_next_event) 
      {
        
        m_infection_realisation_empty_catch = 1;
        // Set infection state to first infection state and then clear that state realisation and time
        while (m_infection_realisation_empty_catch == 1)
        {
          m_infection_state = m_infection_state_realisation_queue.front();
          schedule_m_day_of_InfectionStatus_change(parameters);
          m_infection_state_realisation_queue.pop();
          m_infection_time_realisation_queue.pop();
          
          // Allocate active strain and time of acquisition to map
          // TODO: Add strain allocation in and removing of strain pointer in strain realisation vector
          // allocate_strain(parameters, m_infection_strain_realisation_queue.front());
          // m_infection_strain_realisation_queue.pop();
          
          m_number_of_strains++; // This will be removed as should be in allocate strain etc
          schedule_m_day_of_strain_clearance(parameters);
          
          if (!m_infection_time_realisation_queue.empty()) {
            if (m_infection_time_realisation_queue.front() == parameters.g_current_time) {}
            else { m_infection_realisation_empty_catch = 0; }
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
double Person::update(const Parameters &parameters)
{
  // Update age
  m_person_age++;
  
  // Find out if pregnacy age
  if (m_person_age > 20 * 365 && m_person_age < 21 * 365) {
    s_sum_maternal_immunity += m_ICA;
    s_total_mums++;
  }
  
  // Update age dependent biting rate and the sum
  s_psi_sum += m_age_dependent_biting_rate = (1 - (parameters.g_rho * exp(-m_person_age / parameters.g_a0)));
  
  // Handle any events that are happening today
  if (m_day_of_next_event == parameters.g_current_time) {
    event_handle(parameters);
  }
  
  // Reset total bites
  m_number_of_bites = 0;
  m_number_of_succesful_bites = 0;
  
  // Throw if the next event date is still today
  assert(m_day_of_next_event != parameters.g_current_time && 
    "Update function failed to handle event update as next event is still today");
  
  return(m_age_dependent_biting_rate);
  
}


// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// LOGGERS - functions that report summaries on an individual
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Caluclate daily incidence for population and 0-5 years, i.e. would they cause an incident case today
// Returns 2 if 0-5 and incident, 1 if just incident, and 0 if neither
int Person::log_daily_incidence(){
  
  // Are they considered to be a new case, i.e. S, A, U
  if(m_infection_state == SUSCEPTIBLE || m_infection_state == ASYMPTOMATIC || m_infection_state == SUBPATENT) 
  {
    // Do they have any infections pending
    if (!m_infection_time_realisation_queue.empty()) 
    {
      // Is that infection pending for today
      if (m_infection_time_realisation_queue.front() == parameters.g_current_time) 
      {
        // Is that infection yielding a clinical case
        if (m_infection_state_realisation_queue.front() == TREATED || m_infection_state_realisation_queue.front() == DISEASED) 
        {
          
          // Are they aged 0 - 5 years
          if(m_person_age <= (5*365))
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





