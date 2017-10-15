#include "stdafx.h"
#include "mosquito.h"

// Non class constructor which will inialise a random strain
Mosquito::Mosquito(const Parameters &parameters) :
  m_day_of_death(rexpint1(parameters.g_mean_mosquito_age) + parameters.g_current_time + 1),
  m_day_of_next_blood_meal(runiform_int_1(1, 3) + parameters.g_current_time)
{
  // Reserve vector memory
  m_oocyst_rupture_time_vector.reserve(10);
  m_oocyst_barcode_male_vector.reserve(10);
  m_oocyst_barcode_female_vector.reserve(10);
  
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// ALLOCATIONS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Allocate male and female gameteocyte resulting froom biting infected human
void Mosquito::allocate_gametocytes(const Parameters &parameters, std::vector<barcode_t> gametocytes)
{
  // Push gametocyte barcodes and time of bursting
  
  m_oocyst_barcode_male_vector.emplace_back(gametocytes[0]);
  m_oocyst_barcode_female_vector.emplace_back(gametocytes[1]);
  m_oocyst_rupture_time_vector.emplace_back(parameters.g_current_time + static_cast<int>(parameters.g_delay_mos));
  
  // Update infection status if mosquito is susceptible
  if (m_mosquito_infection_state == SUSCEPTIBLE)
  {
    m_mosquito_infection_state = EXPOSED;
  }
  
  schedule_m_day_of_next_event();
  
};

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// SCHEDULERS 
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// Schedule mosquito's death day
// Also return true if the scheduled day is one such that we care about bites this mosquito takes, i.e.
// it will live long enought to become infectious
bool Mosquito::schedule_m_day_of_death(const Parameters &parameters)
{
  
  // Exponential waiting time plus current day and 1 so not the same day
  // If less than 11 also move their blood meal day  to after their death
  // and set next event to their death day as 
  // and return true so we don't bother resetting anything else about the mosquito in the die function
  m_day_of_death = rexpint1(parameters.g_mean_mosquito_age);
  
  if (m_day_of_death < 11) {
    m_day_of_next_blood_meal = parameters.g_current_time + 12;
    
    m_day_of_death += parameters.g_current_time + 1;
    m_day_of_next_event = m_day_of_death;
    return(false);
  }
  
  m_day_of_death += parameters.g_current_time + 1;
  return(true);
  
}

// Schedule mosquito's day of next blood meal
void Mosquito::schedule_m_day_of_blood_meal(const Parameters & parameters)
{
  
  // Current day plus 3 for the moment, but might change in future for multiple biting etc so 
  m_day_of_next_blood_meal = parameters.g_current_time + 3;
  
}

// Set day of next event
void Mosquito::schedule_m_day_of_next_event()
{
  
  // Start next event as death
  m_day_of_next_event = m_day_of_death;
  
  // Catch if queue is empty first
  // Check if the next pending oocyst rupture is the next event
  if (m_oocyst_rupture_time_vector.size() > m_ruptured_oocyst_count)
  {
    if (m_day_of_next_event > m_oocyst_rupture_time_vector[m_ruptured_oocyst_count] && m_oocyst_rupture_time_vector[m_ruptured_oocyst_count] > 0) {
      m_day_of_next_event = m_oocyst_rupture_time_vector[m_ruptured_oocyst_count];
    }
  }
  
  //}
  
  // Compare against strain clearance day
  if (m_day_of_next_event > m_day_of_next_blood_meal && m_day_of_next_blood_meal > 0) {
    m_day_of_next_event = m_day_of_next_blood_meal;
  }
  
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// UPDATERS 
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Kill mosquito, i.e. reset age to 0, infections to 0, state to susceptible, ino barcodes etc
// Returns true if mosquito lives long enough to surive incubation, otherwise false
bool Mosquito::die(const Parameters & parameters)
{
  
  // Reset whether mosquito is off season to false, as if we are killing th mosquito we must be considering it
  m_mosquito_off_season = false;
  
  // Schedule death day. If returns true, i.e. mosquito won't survive incubation, then return early as
  // no need to reset everything else
  if (!schedule_m_day_of_death(parameters)) {
    return (false);
  }
  
  // Reset mosquito's age and infection status changes
  m_ruptured_oocyst_count = 0;
  m_day_of_next_event = 0;
  m_mosquito_infected = false;
  
  // Clear vectors
  m_oocyst_rupture_time_vector.clear();
  m_oocyst_barcode_male_vector.clear();
  m_oocyst_barcode_female_vector.clear();
  
  // Make mosquito susceptible again
  m_mosquito_infection_state = SUSCEPTIBLE;
  
  // Schedule next blood meal
  // schedule_m_day_of_blood_meal(parameters);
  schedule_m_day_of_blood_meal(parameters);
  
  
  // Schedule next event
  schedule_m_day_of_next_event();
  
  return(true);
}


// Daily update function, i.e. check for sesonal mosquitoes and then handle any events
bool Mosquito::update(Parameters &parameters)
{
  // if the mosquito is off season, and there is a positive deficit, i.e. a surplus of mosquitoes then return early
  if (m_mosquito_off_season && parameters.g_mosquito_deficit > 0) {
    return(false);
  } 
  
  // if the mosquito is off season, and there is a deficit, i.e. a need of mosquitoes then return true early as will be biting today
  // DEPRECATED IF STATEMENT
  // This will never be reached as when a new mosquito first comes back on season we don't update it we simply copy an already on season one
  // in order to keep the correct mosquito popualtion infection. 
  if (m_mosquito_off_season && parameters.g_mosquito_deficit < 0) {
    parameters.g_mosquito_deficit++;
    return(die(parameters));
  }
  
  // Handle any events that are happening today
  if (m_day_of_next_event == parameters.g_current_time) {
    // if event handle returns true the mosquito is biting today 
    return(event_handle(parameters)); 
  }
  
  // if there was no event or seasonality then return false
  return(false);
  
}

// Event handle, i.e. mosquito death, state change, biting handle.
// If it returns true then the mosquito is biting today
bool Mosquito::event_handle(const Parameters &parameters)
{
  
  // If death occurs then no need to explore the other events
  if (m_day_of_death == m_day_of_next_event)
  {
    return(die(parameters));
  }
  else
  {
    // All other events could happen theoretically on the same day though so within same else block
    // Update blood meal day
    if (m_day_of_next_blood_meal == m_day_of_next_event)
    {
      
      // Schedule next blood meal
      schedule_m_day_of_blood_meal(parameters);
      
    }
    
    // Realise oocyst rupturing if time to do so
    // If there are more pending oocyst times than ruptures then see if the pending is happening today
    if (m_oocyst_rupture_time_vector.size() > m_ruptured_oocyst_count)
    {
      
      // if the vector is not for today, i.e. oocyst does not burst today
      if (m_oocyst_rupture_time_vector[m_ruptured_oocyst_count] == m_day_of_next_event)
      {
        
        // Update the catch - this is needed for when two oocysts might burst on the same day - this is not possible at the moment but might be incorporated later if we assume mosquitos
        // might produce more than 1 oocyst upon one feeding event
        // m_oocyst_realisation_empty_catch = 1;
        
        // while (m_oocyst_realisation_empty_catch == 1)
        // {
        
        // Make infection status infected
        m_mosquito_infection_state = INFECTED;
        m_mosquito_infected = true;
        
        // Update the burst oocyst count
        m_ruptured_oocyst_count++;
        
        
        // Update the catching - i.e. if there are no oocysts to rupture today then do not update the catch
        /*if (!m_pending_oocyst_time_queue.empty())
         {
         if (m_pending_oocyst_time_queue.front() != parameters.g_current_time)
         {
         m_oocyst_realisation_empty_catch = 0;
         }
         }
         else
         {
         m_oocyst_realisation_empty_catch = 0;
         }
         
         }*/
        
        
      }
    }
    
    // Update the day of next event
    schedule_m_day_of_next_event();
    
    // Throw if the next event is in the past unless mosquito is off season
    assert((m_day_of_next_event > parameters.g_current_time || m_mosquito_off_season) &&
      "Mosquitoes stuck in the past");
    
    // If mosquito is biting today return true
    if (m_day_of_next_blood_meal == (parameters.g_current_time + 3)) {
      return(true);
    }
    else
    {
      return(false);
    }
  }
  
  
  
}
