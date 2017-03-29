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
void Mosquito::schedule_m_day_of_death(const Parameters &parameters)
{
  // Throw if called on mosquito which is not SUSCEPTIBLE
  assert(m_mosquito_infection_state == SUSCEPTIBLE && "Death day schedule called for mosquito which is not susceptible");
  
  // Exponential waiting time plus current day and 1 so not the same day
  m_day_of_death = rexpint1(parameters.g_mean_mosquito_age) + parameters.g_current_time + 1;
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
void Mosquito::die(const Parameters & parameters)
{
  
  // TODO: Add seasonality hook here to only replace upon death if seasonality is allowing
  
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
  
  // Schedule next blood meal and death
  schedule_m_day_of_blood_meal(parameters);
  schedule_m_day_of_death(parameters);
  
  // Schedule next event
  schedule_m_day_of_next_event();
  
}

// Daily update function, i.e. check events through event_handle, and returns either -1 or moqsuito ID if mosquito is biting today
bool Mosquito::update(Parameters &parameters)
{
  
  // Handle any events that are happening today
  if (m_day_of_next_event == parameters.g_current_time) {
    event_handle(parameters);
  }
  
  // Throw if the next event date is still today
  assert(m_day_of_next_event != parameters.g_current_time &&
    "Update function failed to handle event update as next event is still today");
  
  // If mosquito is biting today return its ID otherwise return -1
  if (m_day_of_next_blood_meal == (parameters.g_current_time + 3)) {
    return(true);
  }
  else
  {
    return(false);
  }
  
}

// Event handle, i.e. mosquito death, state change, biting handle
void Mosquito::event_handle(const Parameters & parameters)
{
  
  // If death occurs then no need to explore the other events
  if (m_day_of_death == m_day_of_next_event) 
  {
    die(parameters);
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
  }
  
  
}
