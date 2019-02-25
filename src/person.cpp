#include "person.h"

// Only Constructor
Person::Person(const Parameters &parameters) :
  
  m_person_age{ set_initial_m_person_age(parameters.g_average_age, parameters.g_max_age) },
  m_age_dependent_biting_rate{ set_m_age_dependent_biting_rate(parameters.g_rho, parameters.g_a0) },
  m_individual_biting_rate{ set_initial_m_individual_biting_rate(parameters.g_zeta_meanlog, parameters.g_zeta_sdlog) }
  {
    // Allocate death day
    set_initial_m_day_of_death(parameters);
    
    // Reserve storage space for vectors - trialing suggests the below should fit about right in terms of minimising vector resizing
    m_active_strains.reserve(static_cast<int>(100 * m_individual_biting_rate));
    m_resistant_strains.reserve(static_cast<int>(100 * m_individual_biting_rate));
    m_post_treatment_strains.reserve(static_cast<int>(100 * m_individual_biting_rate));
    
    m_active_strain_contribution.reserve(static_cast<int>(100 * m_individual_biting_rate));
    m_infection_time_realisation_vector.reserve(static_cast<int>(100 * m_individual_biting_rate));			// First pending infection time in position 0 to handle multiple infections times that have not been realised yet
    m_infection_state_realisation_vector.reserve(static_cast<int>(100 * m_individual_biting_rate));			// First pending infection state in position 0 to handle multiple infections states that have not been realised yet
    m_infection_barcode_realisation_vector.reserve(static_cast<int>(100 * m_individual_biting_rate));		// First pending infection barcode in position 0 to handle multiple infections states that have not been realised yet
    
    // if doing nmf then set this here
    if(parameters.g_nmf_flag){
      set_m_day_of_nmf(parameters);
    } 
    
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
    // work out the number of strains that are gametocytogenic, i.e. they were realised more than delay_gam time earlier
    for(int n = 0 ; n < m_number_of_strains ; n++){
      if(m_active_strains[n].get_m_day_of_strain_acquisition() < parameters.g_current_time - parameters.g_delay_gam){
        m_gametocytogenic_strains.emplace_back(n);
        m_gametocytogenic_infections++;
      }
    }
  }
  
  // If there are gametocytogenic infections then work out whether they led to onward infecion of the mosquito, otherwise return false
  if (m_gametocytogenic_infections)
  {
    
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
  else 
  {
    return(false);
  }
  
}

// Work out which strain barcodes are passed on as gametocytes
std::vector<boost::dynamic_bitset<>> Person::sample_two_barcodes(const Parameters &parameters)
{
  // Work out what the contribution for each strain is 
  if (m_contribution_counter == 0) {
    
    // loop over the smaller of the number of strains or gametocytogenic infections - need to do this min to handle when strains are cleared,
    // as it might occur that there are fewer strains than realised infections that could have led to gametocytogenic infections (if they had not been cleared)
    for (int n = 0; n < m_gametocytogenic_infections; n++)
    {
      
      // Match infection state and schedule associated next state change
      switch (m_active_strains[m_gametocytogenic_strains[n]].get_m_strain_infection_status())
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
      
      // if there are resistance elements then we nneed to create our modifier
      if(parameters.g_resistance_flag){
        m_active_strain_contribution.back() *= m_active_strains[m_gametocytogenic_strains[n]].relative_contribution(parameters);
      }
      
      if(parameters.g_vector_adaptation_flag){
        if(m_infection_state == TREATED){
        if(!m_active_strains[m_gametocytogenic_strains[n]].barcode_position(0)){
        m_active_strain_contribution.back() *= parameters.g_gametocyte_non_sterilisation;
        }
        }
      }
      
    }
    m_contribution_counter = 1;
    m_contribution_sum = std::accumulate(m_active_strain_contribution.begin(), m_active_strain_contribution.end(), 1.0);
  }
  
  // If individual only has one strain then catch for this as we won't need to draw a random sample for what strain is drawn
  if (m_number_of_strains == 1)
  {
    return(std::vector<boost::dynamic_bitset<>> { m_active_strains[0].get_m_barcode(), m_active_strains[0].get_m_barcode() });
  }
  else
  {
    // TODO: If we need selfing, then introduce effective selfing here, by making a m_temp_active_strain_contribution, for which the position that 
    // was drawn for the first barcode becomes x, such that p(selfing) = x/std::accumulate(m_temp_active_strain_contribution)
    return(std::vector<boost::dynamic_bitset<> > { m_active_strains[m_gametocytogenic_strains[sample1(m_active_strain_contribution, m_contribution_sum)]].get_m_barcode(), 
           m_active_strains[m_gametocytogenic_strains[sample1(m_active_strain_contribution, m_contribution_sum)]].get_m_barcode() });
  }
}

// Work out if late parasitological failure happened
bool Person::late_paristological_failure_boolean(const Parameters &parameters){
  
  // the default drug
  int drug_choice = parameters.g_drug_choice;
  int time_ago = 0;
  double prob_of_lpf = 0.0;
  
  // set up our post treatment vectors
  m_post_treatment_strains.clear();
  m_resistant_strains.clear();
  
  m_post_treatment_strains.reserve(m_number_of_strains);
  m_resistant_strains.reserve(m_number_of_strains);
  
  // are we doing mft, and if so what drug did they get this time
  if(parameters.g_mft_flag) {
    drug_choice = sample1(parameters.g_partner_drug_ratios, 1.0);  
  }
  
  //loop through the final m_numbers_of_last_passed as these have to be the relevant treatment strains
  for(int ts = 0; ts < m_number_of_strains ; ts++){
    
    prob_of_lpf = m_active_strains[ts].late_paristological_failure_prob(parameters, drug_choice);
    if(prob_of_lpf != 0.0) { 
      
      m_resistant_strains.emplace_back(m_active_strains[ts]);
      
      // if infection was more than dur_A ago then it is cleared for sure
      time_ago = ((parameters.g_current_time - m_active_strains[ts].get_m_day_of_strain_acquisition()) / parameters.g_dur_A);
      
      if( time_ago < 1 ) {
        if(rbernoulli1(prob_of_lpf * (1 - time_ago))) {
          m_post_treatment_strains.emplace_back(m_active_strains[ts]);
        }
      }

    }
  }
  
  return(m_post_treatment_strains.size());
  
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// SETTERS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Set person's individual biting rate
double Person::set_initial_m_individual_biting_rate(double zeta_meanlog, double zeta_sdlog) {
  
  double zeta = rlognorm1(zeta_meanlog, zeta_sdlog);
  while(zeta > 100) zeta = rlognorm1(zeta_meanlog, zeta_sdlog);
  return(zeta);
  
}

// Set person's age
int Person::set_initial_m_person_age(double average_age, int max_age) {
  
  int age = rexpint1(average_age);
  while(age > max_age) age = rexpint1(average_age);
  return(age);
  
}

// Set person's initial death day
void Person::set_initial_m_day_of_death(const Parameters &parameters)
{
  // Exppnential waiting time plus current day and 1 so not the same day
  m_day_of_death = rexpint1(parameters.g_average_age) + parameters.g_current_time + 1;
  
}

// Set person's relative chance of biting (psi)
double Person::set_m_age_dependent_biting_rate(double rho, double a0) {
  
  return(1 - (rho*exp(-m_person_age / a0)));
  
}

// Set day of next nmf event
void Person::set_m_day_of_nmf(const Parameters &parameters){
  
  while(parameters.g_nmf_age_brackets[m_nmf_age_band+1] < m_person_age){
    m_nmf_age_band++;
  }
  
  // Exppnential waiting time plus current day and 1 so not the same day
  m_day_of_nmf = rexpint1(parameters.g_mean_nmf_frequency[m_nmf_age_band]) + parameters.g_current_time + 1;
  
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
      m_temp_int = m_active_strains[s].get_m_day_of_strain_infection_status_change();
      if (m_temp_int <= m_day_of_next_strain_state_change && m_temp_int != 0)
      {
        
        // if the next strain state change day is the same as the temp flag the bool
        if (m_day_of_next_strain_state_change == m_temp_int)
        {
          m_more_than_one_strain_to_change_today_bool = true;
        }
        
        // keep earliest time and strain position 
        m_day_of_next_strain_state_change = m_temp_int;
        m_temp_strain_to_next_change = s;
        
      }
    }
  }
}


// Set day of next event
void Person::set_m_day_of_next_event() {
  
  // First if catches if person is completely susceptible, i.e. not with a pending infection as this is the initialisation criteria
  if (m_infection_time_realisation_vector.empty() && m_day_of_strain_clearance == 0 && m_day_of_InfectionStatus_change == 0 && m_day_of_nmf == 0) {
    m_day_of_next_event = m_day_of_death;
  }
  
  else // find minimum that is not zero
  {
    // Start next event as death
    m_day_of_next_event = m_day_of_death;
    
    // Compare against nmf
    if (m_day_of_next_event > m_day_of_nmf && m_day_of_nmf != 0) {
      m_day_of_next_event = m_day_of_nmf;
    }
    
    // Compare against strain state change
    if (m_day_of_next_event > m_day_of_next_strain_state_change && m_day_of_next_strain_state_change != 0) {
      m_day_of_next_event = m_day_of_next_strain_state_change;
    }
    
    // Catch if there are no pending infections first, i.e. if the number of realized infections is the same as the length of the infection realisation vector
    if (m_infection_time_realisation_vector.size() > static_cast<unsigned int>(m_number_of_realised_infections))
    {
      if (m_day_of_next_event > m_infection_time_realisation_vector[m_number_of_realised_infections] && m_infection_time_realisation_vector[m_number_of_realised_infections] != 0) {
        m_day_of_next_event = m_infection_time_realisation_vector[m_number_of_realised_infections];
      }
    }
    
    // Compare against strain clearance day
    /*
     if (m_day_of_next_event > m_day_of_strain_clearance && m_day_of_strain_clearance != 0) {
     m_day_of_next_event = m_day_of_strain_clearance;
     }
     */
    
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
void Person::allocate_bite(Parameters &parameters, Mosquito &mosquito)
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
  // Firstly if the human is treated or in prophylaxis or is prophylactic but mosquito has resistant strains and they are soon to recover.
  if (m_infection_state != TREATED || m_infection_state != PROPHYLAXIS || 
      (m_infection_state == PROPHYLAXIS && mosquito.check_resistance() && m_day_of_InfectionStatus_change > (parameters.g_current_time-5))
  )
  {
    
    // Random draw to see if the bite led to an infection
    
    if (rbernoulli1(m_biting_success_rate)) {
      
      // Allocate infection
      parameters.g_total_human_infections++;
      allocate_infection(parameters, mosquito);
    }
  }
  
  // Increase number of bites 
  m_number_of_bites++;
  
}

// Allocate an infection to person, i.e. individuals who return >0 from allocate_force_of_infection()
void Person::allocate_infection(Parameters &parameters, Mosquito &mosquito)
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
    m_transition_probabilities[1] = m_symptom_success_rate * parameters.g_ft;
    // N -> A
    m_transition_probabilities[2] = 1 - m_symptom_success_rate;
    
    // Set outcome probability sum 
    m_sum_transition_probabilities = m_transition_probabilities[0] + m_transition_probabilities[1] + m_transition_probabilities[2];
    
    // Draw what infection state they move to 
    m_temp_infection_state = m_transition_vector[sample1(m_transition_probabilities, m_sum_transition_probabilities)];
    
  }
  
  // Increase number of succesful bites
  m_number_of_succesful_bites++;
  
  // Allocate strains being passed on
  for(int cotransmission = 0; cotransmission < parameters.g_cotransmission_frequencies[parameters.g_cotransmission_frequencies_counter]; cotransmission++)
  {
    
    // if it is the first sporozoite then it is always taken. For successive sporozoites, we probabilistically decide based on their immunity
    if(cotransmission == 0 || rbernoulli1(m_biting_success_rate)) {
      
      if(!cotransmission){
        m_cotransmission_realisation_vector.emplace_back(false);
      } else {
        m_cotransmission_realisation_vector.back() = true;
        m_cotransmission_realisation_vector.emplace_back(true);
      }
      
      // Push the resultant state of infection
      m_infection_state_realisation_vector.emplace_back(m_temp_infection_state);
      
      // Push the resultant state change time
      m_infection_time_realisation_vector.emplace_back(static_cast<int>(parameters.g_dur_E + parameters.g_current_time));
      
      // Allocate strains from mosquito
      
      // if we are doing spatial then use the exported barcodes first - the human biting quueue is shuffled so distributed across humans fine.
      if(parameters.g_spatial_imported_human_infection_counter < parameters.g_spatial_total_imported_human_infections)
      {
        
        // asign the exported barcode and increase the count
        if(parameters.g_spatial_type == Parameters::METAPOPULATION)
        {
          m_infection_barcode_realisation_vector.emplace_back(parameters.g_spatial_imported_barcodes[parameters.g_spatial_imported_human_infection_counter]);
        } 
        else 
        {
          m_infection_barcode_realisation_vector.emplace_back(Strain::generate_next_barcode());  
        }
        
        parameters.g_spatial_imported_human_infection_counter++;
        
      }
      else 
      {
        
        m_infection_barcode_realisation_vector.emplace_back(mosquito.sample_sporozoite());
        
        // export a barcode if doing metapopulation spatial
        if(parameters.g_spatial_exported_barcode_counter < parameters.g_spatial_total_exported_barcodes)
        {
          parameters.g_spatial_exported_barcodes[parameters.g_spatial_exported_barcode_counter] = m_infection_barcode_realisation_vector.back();
          parameters.g_spatial_exported_barcode_counter++;  
        }
        
      }
      
      // are we simulating mutations
      if (parameters.g_mutation_flag){
        // are we still allocating mutations
        if(parameters.g_mutation_pos_allocator < parameters.g_num_loci){
          // were there any mutations at this position for the population
          if(parameters.g_mutations_today[parameters.g_mutation_pos_allocator]>0){
            m_infection_barcode_realisation_vector[m_infection_barcode_realisation_vector.size()-1][parameters.g_mutation_pos_allocator].flip();
            parameters.g_mutations_today[parameters.g_mutation_pos_allocator]--;
            if(parameters.g_mutations_today[parameters.g_mutation_pos_allocator] == 0){
              parameters.g_mutation_pos_allocator++;
            }
          } else {
            parameters.g_mutation_pos_allocator++;
          }
        }
      }
      
    }
    
  }
  
  // Increase cotransmission counter and catch for overflow
  if(++parameters.g_cotransmission_frequencies_counter == parameters.g_cotransmission_frequencies_size) {
    parameters.g_cotransmission_frequencies_counter = 0;
  }
  
  // Set next event date as may have changed as a result of the bite
  set_m_day_of_next_event();
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// SCHEDULERS & DRAWS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Draw day of infection state change 
int Person::draw_m_day_of_InfectionStatus_change(const Parameters &parameters)
{
  // Match infection state and schedule associated next state change
  switch (m_infection_state)
  {
  case DISEASED:
    return(rexpint1(parameters.g_dur_D) + parameters.g_current_time + 1);
  case ASYMPTOMATIC:
    return(rexpint1(parameters.g_dur_A) + parameters.g_current_time + 1);
  case SUBPATENT:
    return(rexpint1(parameters.g_dur_U) + parameters.g_current_time + 1);
  case TREATED:
    return(rexpint1(parameters.g_dur_T) + parameters.g_current_time + 1);
  case PROPHYLAXIS:
    return(rexpint1(parameters.g_dur_P) + parameters.g_current_time + 1);
  default:
    assert(NULL && "Schedule Infection Status Change Error - person's infection status not D, A, U, T or P");
  return(-1);
  break;
  }
}

// Schedule person's infection state change
void Person::schedule_m_day_of_InfectionStatus_change(const Parameters &parameters)
{
  m_day_of_InfectionStatus_change = draw_m_day_of_InfectionStatus_change(parameters);
}

// Schedule person's death day
void Person::schedule_m_day_of_death(const Parameters &parameters)
{
  // Throw if called on someone who is not SUSCEPTIBLE
  // assert(m_infection_state == SUSCEPTIBLE && "Death day schedule called for someone who is not susceptible");
  
  // Exppnential waiting time plus current day and 1 so not the same day
  m_day_of_death = rexpint1(parameters.g_average_age) + parameters.g_current_time + 1;
}

// Schedule person's next strain clearance#
// DEPRECATED - NOT USED ANY MORE
void Person::schedule_m_day_of_strain_clearance(const Parameters &parameters)
{
  // Throw if this is called on someone with no strains
  assert(m_number_of_strains != 0 && "Tried to schedule a strain clearance day from an individual with no strains");
  
  // Work out potential new clearance day
  int possible_clearance_day = rexpint1(m_number_of_strains * parameters.g_dur_AU) + parameters.g_current_time + 1;
  
  // If the new clearance day is earlier then schedule it
  if (possible_clearance_day < m_day_of_strain_clearance || !m_day_of_strain_clearance || m_day_of_strain_clearance == parameters.g_current_time) {
    m_day_of_strain_clearance = possible_clearance_day;
  }
  
  // FOR NOW
  m_day_of_strain_clearance = 0;
  
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
  m_day_last_treated = 0;
  
  // set their boost times back enough so that they are guaranteed to boost if bitten tomorrow
  m_IB_last_boost_time -= parameters.g_uB;
  m_ICA_last_boost_time -= parameters.g_uCA;
  m_ID_last_boost_time -= parameters.g_uD;
  
  // Set maternal immunity
  m_ICM = m_ICM_init = parameters.g_PM * parameters.g_mean_maternal_immunity;
  
  // Schedule death
  schedule_m_day_of_death(parameters);
  
  // Schedule nmf
  if(parameters.g_nmf_flag){
    set_m_day_of_nmf(parameters);
  }
  
  // Set new day of next event
  set_m_day_of_next_event();
}

// Recover to being susceptible, i.e. clearing all infections and strains and associated timings
void Person::recover(const Parameters &parameters)
{
  
  // Return to susceptible
  m_infection_state = SUSCEPTIBLE;
  
  // Clear strains
  all_strain_clearance();
  
  // Reset other events to 0
  m_day_of_InfectionStatus_change = 0;
  
  // Update next event counter - due to nature of recovering the next event has to be their death
  m_day_of_next_event = m_day_of_death;
  
}

// Treatment outcomes
void Person::treatment_outcome(const Parameters &parameters) {
  
  // are we doing resistance firstly
  if(parameters.g_resistance_flag) { 
    
    // are they currently slow parasite clearance
    if(m_slow_parasite_clearance_bool){
      
      m_slow_parasite_clearance_bool = false;
      m_treatment_outcome = SUCCESFULLY_TREATED;
      m_infection_state = PROPHYLAXIS;
      m_day_last_treated = parameters.g_current_time;
      schedule_m_day_of_InfectionStatus_change(parameters); // schedule next state change
      all_strain_clearance(); // When they are in prophylaxis we remove all the strains, as treated individuals still have strains
      
    } else {
      
      // are they actually infected, i.e. not just here because of nmf
      if(m_number_of_strains > 0) {
        
        // do they fail due to LPF
        if(late_paristological_failure_boolean(parameters)){
          late_paristological_failure(parameters);
          m_treatment_outcome = LPF;
        } else {
          // did no LPF happen because there were no resistant strains - straight to prophylaxis
          if(!m_resistant_strains.size()) {
            m_treatment_outcome = SUCCESFULLY_TREATED;
            m_infection_state = PROPHYLAXIS;
            m_day_last_treated = parameters.g_current_time;
            schedule_m_day_of_InfectionStatus_change(parameters); // schedule next state change
            all_strain_clearance(); // When they are in prophylaxis we remove all the strains, as treated individuals still have strains
          } else {
            // otherwise we consider a SPC
            slow_treatment_clearance(parameters);
          }
        } 
        
      } else {
        m_infection_state = PROPHYLAXIS;
        m_day_last_treated = parameters.g_current_time;
        schedule_m_day_of_InfectionStatus_change(parameters); // schedule next state change
      }
    }
  } else {
    m_treatment_outcome = SUCCESFULLY_TREATED;
    m_infection_state = PROPHYLAXIS;
    m_day_last_treated = parameters.g_current_time;
    schedule_m_day_of_InfectionStatus_change(parameters); // schedule next state change
    all_strain_clearance(); // When they are in prophylaxis we remove all the strains, as treated individuals still have strains
  }
  
}

// Late parasitological failure
void Person::late_paristological_failure(const Parameters &parameters) {
  
  // change their state to asymptomatic and then draw their state change day
  m_infection_state = ASYMPTOMATIC;
  schedule_m_day_of_InfectionStatus_change(parameters); // schedule next state change  
  
  // clear the active strains
  m_active_strains = m_post_treatment_strains;
  m_number_of_strains = m_post_treatment_strains.size();
  
  //loop through the new number of strains and draw their movement dates
  for(int ts = 0; ts < m_number_of_strains ; ts++){
    
    // Set strains to asymptomatic
    m_active_strains[ts].set_m_strain_infection_status(Strain::ASYMPTOMATIC);
    
    // Draw a potential time for this strain to move out of this state
    m_temp_int = draw_m_day_of_InfectionStatus_change(parameters);
    
    // If it's larger than the drawn human's state change day then set it equal to
    if (m_temp_int > m_day_of_InfectionStatus_change) {
      m_temp_int = m_day_of_InfectionStatus_change;
    }
    m_active_strains[ts].set_m_day_of_strain_infection_status_change(m_temp_int);
  }
}

// Slow parasite clearance
void Person::slow_treatment_clearance(const Parameters &parameters) {
  
  // clear the active strains and just leave the resistant ones
  m_active_strains = m_resistant_strains;
  m_number_of_strains = m_resistant_strains.size();
  
  // give a new end of treatment time
  m_day_of_InfectionStatus_change = rexpint1(parameters.g_dur_SPC) + parameters.g_current_time + 1;
  
  // set their slow parasite clearance flag
  m_slow_parasite_clearance_bool = true;
  
  //loop through the new number of strains and draw their new time to be after the end of treatment
  for(int ts = 0; ts < m_number_of_strains ; ts++){
    m_active_strains[ts].set_m_day_of_strain_infection_status_change(m_day_of_InfectionStatus_change + 1);
  }
  
}

// Seek treatment for nmf
void Person::seek_nmf_treatment(const Parameters &parameters){
  
  // first would they seek treatment
  if(rbernoulli1(parameters.g_ft)){
    
    // are they tested - if not we assume they are just given drugs in this situation
    if(rbernoulli1(parameters.g_prob_of_testing_nmf)) {
      // are they detected
      if(rbernoulli1(q_fun(parameters))){
        
        // then trigger treatment
        m_infection_state = TREATED;
        schedule_m_day_of_InfectionStatus_change(parameters);
        
        // Clear all pending infection vectors 
        m_infection_time_realisation_vector.clear();
        m_infection_state_realisation_vector.clear();
        m_infection_barcode_realisation_vector.clear();
      }
    } else {
      
      // then trigger treatment
      m_infection_state = TREATED;
      schedule_m_day_of_InfectionStatus_change(parameters);
      
      // Clear all pending infection vectors 
      m_infection_time_realisation_vector.clear();
      m_infection_state_realisation_vector.clear();
      m_infection_barcode_realisation_vector.clear();
      
    }
    
  }
  
  // set the new nmf day
  set_m_day_of_nmf(parameters);
  set_m_day_of_next_event();
  
}

// Event handle, i.e. if any of person's death, state change, strain clearance or state change realisation days are today, respond accordingly
void Person::event_handle(const Parameters &parameters) {
  
  // If death occurs then no need to explore the other events
  if (m_day_of_death == m_day_of_next_event || m_person_age >= parameters.g_max_age) {
    die(parameters);
  }
  else
  {
    
    // Handle nmf event if time to do so. If does this will update the day of infection status change 
    if(parameters.g_nmf_flag){
      if (m_day_of_nmf == m_day_of_next_event) {
        seek_nmf_treatment(parameters);
      }
    }
    // All other events could happen theoretically on the same day though so within same else block
    
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
        treatment_outcome(parameters); // handle their state of treatment accordingly
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
      if (true)
      {
        for (int n = 0 ; n < m_number_of_strains ; n++)
        {
          // if the strain is changing today then switch it accordingly
          if (m_active_strains[n].get_m_day_of_strain_infection_status_change() == parameters.g_current_time)
          {
            // We might end up looping through all the strains when a new treated strain is given to a human which has the same strain state change
            // day as a previous strain - in this case just pass over the treated strain
            if (m_active_strains[n].get_m_strain_infection_status() != Strain::TREATED)
            {
              switch (m_active_strains[n].get_m_strain_infection_status())
              {
              case Strain::DISEASED:
                // If strain to change state is diseased it will become asymptomatic and its following state change must be the same as the humans
                m_active_strains[n].set_m_strain_infection_status(Strain::ASYMPTOMATIC);
                m_active_strains[n].set_m_day_of_strain_infection_status_change(m_day_of_InfectionStatus_change);
                break;
              case Strain::ASYMPTOMATIC:
                // If strain to change state is asymptomatic it will become subpatent and its following state change will not occur
                m_active_strains[n].set_m_strain_infection_status(Strain::SUBPATENT);
                // If the human has only one strain then the state change for this strain should be the same, if not then random given U duration mean
                if (m_infection_state == SUBPATENT)
                {
                  m_active_strains[n].set_m_day_of_strain_infection_status_change(m_day_of_InfectionStatus_change);
                }
                else
                {
                  m_active_strains[n].set_m_day_of_strain_infection_status_change(rexpint1(parameters.g_dur_U) + parameters.g_current_time + 1);
                }
                break;
              case Strain::SUBPATENT:
                
                if(m_number_of_strains > 1) {
                  // If strain is subpatent to change then we clear it.
                  // Swap the strain pointer and strain acquisition date at that position to the back
                  std::swap(m_active_strains[n], m_active_strains.back());
                  
                  // Pop thus deleting the random strain pointer
                  m_active_strains.pop_back();
                  
                  // Lower strain counter and decrease n so that we check the strain we just put here
                  m_number_of_strains--;
                  n--;
                  if(m_number_of_strains==0) rcpp_out(parameters.g_h_quiet_print, "removed last strain\n!");
                }
                break;
              default:
                assert(NULL && "Strain state change equested on strain that is not diseasod or asymptomatic or subpatent");
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
          // If strain to change state is asymptomatic it will become subpatent and its following state change is assumed to be that of a subpatent infection
          m_active_strains[m_temp_strain_to_next_change].set_m_strain_infection_status(Strain::SUBPATENT);
          // If the human has only one strain then the state change for this strain should be the same, if not then random given U duration mean
          if(m_number_of_strains == 1)
          {
            m_active_strains[m_temp_strain_to_next_change].set_m_day_of_strain_infection_status_change(m_day_of_InfectionStatus_change);
          }
          else
          {
            m_active_strains[m_temp_strain_to_next_change].set_m_day_of_strain_infection_status_change(rexpint1(parameters.g_dur_U) + parameters.g_current_time + 1);
          }
          break;
        case Strain::SUBPATENT:
          
          if(m_number_of_strains > 1) {
            // If strain is subpatent to change then we clear it.
            // Swap the strain pointer and strain acquisition date at that position to the back
            std::swap(m_active_strains[m_temp_strain_to_next_change], m_active_strains.back());
            
            // Pop thus deleting the random strain pointer
            m_active_strains.pop_back();
            
            // Lower strain counter
            m_number_of_strains--;
          }
          break;
        default:
          assert(NULL && "Strain state change equested on strain that is not diseasod or asymptomatic or subpatent");
        break;
        }
        set_m_day_of_next_strain_state_change();
      }
      
    }
    
    // Realise infection if time to do so
    // First check if there are any pending infections
    if (m_infection_time_realisation_vector.size() > static_cast<unsigned int>(m_number_of_realised_infections))
    {
      // Is the next pending infection for today
      if (m_infection_time_realisation_vector[m_number_of_realised_infections] == m_day_of_next_event)
      {
        
        // Assign infection state
        
        // If you are diseased already and the next infection would make you asymptomatic then remain diseased
        // If it were to make you treated then this is plausible as maybe the extended duration of the fever that would result 
        // from an additional clinical disease may cause you to seek treatment
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
          
          // if you are now diseased this is a treatment failure
          if (m_infection_state == DISEASED){
            m_treatment_outcome = NOT_TREATED;
          }
        }
        
        // schedule next state change
        // Draw a potential time for this new resultant infection to change state
        // Reason it is a potential time is that if you are  in state A for example, 
        // and currently you would move to U in 100 days, it is not right that this additonal
        // infection could produce a change to U sooner, however, greater than does. 
        m_temp_int = draw_m_day_of_InfectionStatus_change(parameters);
        if (m_temp_int > m_day_of_InfectionStatus_change) {
          m_day_of_InfectionStatus_change = m_temp_int;
        } else {
          m_temp_int = m_day_of_InfectionStatus_change;
        }
        
        // Reset this catch variable which ensiures multiple infections on the same day can be realised
        m_infection_realisation_empty_catch = 1;
        
        // loop through if there are multiple sporozoites to pass on
        while (m_infection_realisation_empty_catch)
        {
          
          // if the person is asymptomatic we allow differeet sporozoites to have different parasitaemias
          // by allowing them to reach < 200p/uL at different times, but not larger than the drawn time for 
          // the humnan to move to state U.
          if (m_infection_state == Person::ASYMPTOMATIC) {
            if (m_infection_realisation_empty_catch != 1)
            {
              // Draw a potential time for this strain to move out of this state
              m_temp_int = draw_m_day_of_InfectionStatus_change(parameters);
              
              // If it's large than the drawn human's state change day then set it equal to
              if (m_temp_int > m_day_of_InfectionStatus_change) {
                m_temp_int = m_day_of_InfectionStatus_change;
              }
            } 
          } 
          
          // Push the pending strain
          m_active_strains.emplace_back(m_infection_barcode_realisation_vector[m_number_of_realised_infections],
                                        Strain::m_transition_vector[static_cast<int>(m_infection_state)],
                                                                   m_temp_int,
                                                                   parameters.g_current_time);
          
          // Increase the number of strains and realised infections
          m_number_of_strains++;
          m_number_of_realised_infections++;
          m_infection_realisation_empty_catch++;
          
          // If the new strain will change state earlier than anyother strain then update the state change event day
          if (m_temp_int <= m_day_of_next_strain_state_change) {
            
            // if the next strain state change day is the same as the temp flag the bool
            if (m_day_of_next_strain_state_change == m_temp_int)
            {
              m_more_than_one_strain_to_change_today_bool = true;
            }
            
            // Update the day of next strain state change and the position to the most recent strain i.e the end
            m_day_of_next_strain_state_change = m_temp_int;
            m_temp_strain_to_next_change = m_number_of_strains;
            
          }
          
          // Little conditional loop to assess if there are pending infections still for today and if not to change the catch
          if (m_infection_time_realisation_vector.size() > static_cast<unsigned int>(m_number_of_realised_infections))
          {
            if (m_infection_time_realisation_vector[m_number_of_realised_infections] != parameters.g_current_time)
            {
              m_infection_realisation_empty_catch = 0;
            }
          }
          else 
          { 
            m_infection_realisation_empty_catch = 0; 
          }
          
        }
        
        // If infection state is treated then we clear all pending infection times meaning that we 
        // don't allow a new infection occurring bettween today and when the individual moves from T->P 
        // to cause them to move out of T back to D/A.
        if (m_infection_state == TREATED)
        {
          // Clear all pending infection vectors 
          m_infection_time_realisation_vector.clear();
          m_infection_state_realisation_vector.clear();
          m_infection_barcode_realisation_vector.clear();
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
  m_age_dependent_biting_rate = 1 - (parameters.g_rho * exp(-m_person_age / parameters.g_a0));
  
  // Find out if pregnacy age
  if (m_person_age > (15 * 365) && m_person_age < (25 * 365)) {
    parameters.g_sum_maternal_immunity += m_ICA;
    parameters.g_total_mums++;
  }
  
  // reset this now so we know who has gone to be treated
  m_treatment_outcome = NOT_CLINICAL;
  
  // Handle any events that are happening today
  if (m_day_of_next_event == parameters.g_current_time || m_person_age >= parameters.g_max_age) {
    event_handle(parameters);
  }
  
  // Reset total bites and other counters designed to prevent unnecessary recalculations in cases of successive same day bites on an inidividual
  m_number_of_bites = 0;
  m_number_of_succesful_bites = 0;
  m_contribution_counter = 0;
  m_cA_counter = true;
  m_active_strain_contribution.clear();
  m_gametocytogenic_infections = 0;
  m_gametocytogenic_strains.clear();
  
  
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
    if (m_infection_time_realisation_vector.size() > static_cast<unsigned int>(m_number_of_realised_infections))
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



