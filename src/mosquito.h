//
//  MAGENTA
//  mosquito.h
//
//  Created: OJ on 27/02/2017
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Class describing a parasite strain.
//
// ---------------------------------------------------------------------------


#ifndef MOSQUITO_H
#define MOSQUITO_H

#include <vector>
#include <queue>
#include <bitset>
#include <random>
#include <cassert> // for error checking
#include "probability.h"
#include "strain.h"
#include "parameters.h"

class Person;

class Mosquito {
  
public:
  enum InfectionStatus
  {
    SUSCEPTIBLE,	// 0
    EXPOSED,		// 1
    INFECTED,		// 2
    NUMBER_OF_STATES = 3
  };
  
  bool m_mosquito_infected = false;					// bool for infection speed
  bool m_mosquito_biting_today = false;     // bool for bting day
  
private:
  
  int m_mosquito_ID;							// mosquito ID
  bool m_mosquito_off_season = false;			// is mosquito off season, i.e. currently "dead" due to lack of seasonality relatd carrying capacity
  InfectionStatus m_mosquito_infection_state;	// infection status associated with a strain
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // EVENT VARIABLES 
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  std::vector<int> m_oocyst_rupture_time_vector;	// First oocyst in mosquito at position 0 etc vector to store pending sporozoite times, i.e. day of biting + 10
  std::vector<boost::dynamic_bitset<> > m_oocyst_barcode_male_vector;	// First oocyst barcode male in mosquito at position 0 etc vector to handle pending sporozoite barcodes from male
  std::vector<boost::dynamic_bitset<> > m_oocyst_barcode_female_vector;	// First oocyst barcode male in mosquito at position 0 etc vector to handle pending sporozoite barcodes from female
  std::vector<boost::dynamic_bitset<> > m_generated_sporozoites_vector; // Sporozites that have already been simulated for this mosquito
  int m_generated_sporozoite_count = 0;						// Count of spz generated. 
  
  int m_day_of_death = 0;						// Mosquito's time of death
  int m_day_of_next_blood_meal = 0;			// Mosquito's day of next blood meal, i.e. day of current blood meal + 3
  int m_day_of_next_event = 0;					// Mosquito's closest event day
  
  uint8_t m_oocyst_realisation_empty_catch = 0;	// Variable that allows a check for empty queues when dealing with more than one infection realisation on a day
  uint8_t m_ruptured_oocyst_count = 0;			// counter for number of burst oocysts so that the related paretnal barcodes are allocated to the correct vector position
  
public:
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // CONSTRUCTORS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Default class constructor
  Mosquito() {};
  
  // Non class constructor which will inialise a random strain
  Mosquito(const Parameters &parameters);
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // GETTERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Get mosquito ID
  int get_m_mosquito_ID() { return(m_mosquito_ID); }
  
  // Get mosquito off season
  bool get_m_mosquito_off_season() { return(m_mosquito_off_season); }
  
  // Get mosquito infection status
  InfectionStatus get_m_mosquito_infection_state() { return(m_mosquito_infection_state); }
  
  // Get Mosquito's time of death
  int get_m_day_of_death() { return(m_day_of_death); }
  
  // Get Mosquito's time of next blood meal
  int get_m_day_of_next_blood_meal() { return(m_day_of_next_blood_meal); }
  
  // Get Mosquito's time of next event
  int get_m_day_of_next_event() { return(m_day_of_next_event); }
  
  // Get Mosquito's ruptured oocyst count
  int get_m_ruptured_oocyst_count() { return(m_ruptured_oocyst_count); }
  
  // Get Mosquito's ruptured oocyst count
  int get_m_generated_sporozoite_count() { return(m_generated_sporozoite_count); }
  
  // Get Mosquito's oocyst male barcode given a chosen oocyst count
  boost::dynamic_bitset<> get_m_oocyst_barcode_male_vector(int x) { return(m_oocyst_barcode_male_vector[x]); }
  
  // Get Mosquito's oocyst female barcode given a chosen oocyst count
  boost::dynamic_bitset<> get_m_oocyst_barcode_female_vector(int x) { return(m_oocyst_barcode_female_vector[x]); }
  
  // Get Mosquito's sporozoite x
  boost::dynamic_bitset<> get_m_generated_sporozoites_vector(int x) { return(m_generated_sporozoites_vector[x]); }
  
  // Get Mosquito's oocyst rupture time vector
  std::vector<int> get_m_oocyst_rupture_time_vector() { return(m_oocyst_rupture_time_vector); }
  
  // Get Mosquito's oocyst male barcode vector
  std::vector<boost::dynamic_bitset<>> get_m_oocyst_barcode_male_vector() { return(m_oocyst_barcode_male_vector); }
  
  // Get Mosquito's oocyst female barcode vector
  std::vector<boost::dynamic_bitset<>> get_m_oocyst_barcode_female_vector() { return(m_oocyst_barcode_female_vector); }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // SETTERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Set mosquito ID
  void set_m_mosquito_ID(int x) { m_mosquito_ID = x; }
  
  // Set mosquito's seasonality, i.e. is the mosquito off consideration due to seasonality
  void set_m_mosquito_off_season(bool x) { m_mosquito_off_season = x; }
  
  // Set strain infection status
  void set_m_mosquito_infection_state(InfectionStatus x) { m_mosquito_infection_state = x; }
  
  // Set Mosquito's time of death
  void set_m_day_of_death(int x) { m_day_of_death = x; }
  
  // Set Mosquito's time of next blood meal
  void set_m_day_of_next_blood_meal(int x) { m_day_of_next_blood_meal = x; }
  
  // Set Mosquito's time of next event
  void set_m_day_of_next_event(int x) { m_day_of_next_event = x; }
  
  // Set Mosquito's ruptured oocyst count
  void set_m_generated_sporozoite_count(int x) { m_generated_sporozoite_count = x; }
  
  // Set Mosquito's ruptured oocyst count
  void set_m_ruptured_oocyst_count(int x) { m_ruptured_oocyst_count = x; }
  
  // Set Mosquito's ruptured oocyst male barcode by emplacing to back
  void set_m_oocyst_barcode_male_vector(boost::dynamic_bitset<> x) { m_oocyst_barcode_male_vector.emplace_back(x); }
  
  // Set Mosquito's ruptured oocyst female barcode by emplacing to back
  void set_m_oocyst_barcode_female_vector(boost::dynamic_bitset<> x) { m_oocyst_barcode_female_vector.emplace_back(x); }
  
  // Set Mosquito's new sporozoite by emplacing to back
  void set_m_generated_sporozoites_vector(boost::dynamic_bitset<> x) { m_generated_sporozoites_vector.emplace_back(x); }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // CHECKERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Function to quickly check if mosquito is infected
  bool check_infection() {  return((m_mosquito_infection_state == INFECTED) ? true : false); }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // SEMI - SETTERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Set mosquito's oocyst barcode male vector
  void set_m_oocyst_barcode_male_vector(std::vector<boost::dynamic_bitset<>> x) { m_oocyst_barcode_male_vector = x;}
  
  // Set mosquito's oocyst barcode female vector
  void set_m_oocyst_barcode_female_vector(std::vector<boost::dynamic_bitset<>> x) { m_oocyst_barcode_female_vector = x;}
  
  // Set mosquito's oocyst rupture time vector
  void set_m_oocyst_rupture_time_vector(std::vector<int> x)  
  { 
    for (unsigned int i = 0; i < x.size(); i++) 
    {
      m_oocyst_rupture_time_vector.emplace_back(x[i]);
    }
    
  }
  
  // Set mosquito's oocyst barcode male vector
  void set_m_oocyst_barcode_male_vector_from_vector_of_vector_bool(std::vector<std::vector<bool> > x)
  {
    m_oocyst_barcode_male_vector.reserve(x.size());
    unsigned int temp_barcode_iterator = 0;
    
    for(unsigned int i = 0; i < x.size() ; i++){
      
      // fetch vector<bool> and turn into barcode
      for(temp_barcode_iterator = 0; temp_barcode_iterator < Parameters::g_barcode_length ; temp_barcode_iterator++ )
      {
        Strain::temp_barcode[temp_barcode_iterator] = x[i][temp_barcode_iterator];
      }
      m_oocyst_barcode_male_vector.emplace_back(Strain::temp_barcode);
    }
  }
  
  // Set mosquito's oocyst barcode male vector
  void set_m_oocyst_barcode_female_vector_from_vector_of_vector_bool(std::vector<std::vector<bool> > x)
  {
    m_oocyst_barcode_female_vector.reserve(x.size());
    unsigned int temp_barcode_iterator = 0;
    
    for(unsigned int i = 0; i < x.size() ; i++){
      
      // fetch vector<bool> and turn into barcode
      for(temp_barcode_iterator = 0; temp_barcode_iterator < Parameters::g_barcode_length ; temp_barcode_iterator++ )
      {
        Strain::temp_barcode[temp_barcode_iterator] = x[i][temp_barcode_iterator];
      }
      m_oocyst_barcode_female_vector.emplace_back(Strain::temp_barcode);
    }
  }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // ALLOCATIONS ETC
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Allocate male and female gameteocyte resulting froom biting infected human
  void allocate_gametocytes(const Parameters &parameters, 
                            boost::dynamic_bitset<> &gam_1,
                            boost::dynamic_bitset<> &gam_2);
  
  // Handle bite, i.e. allocating gametocytes and other housekeeping
  void handle_bite(Parameters &parameters, Person &person);
  
  // Allocate sporozoite arising from importation
  void allocate_imported_sporozoite(boost::dynamic_bitset<> x); 
  
  // Sample one sporozoite to be passed on to the human
  boost::dynamic_bitset<> sample_sporozoite();
    
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // SCHEDULERS 
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Schedule mosquitos's death day
  bool schedule_m_day_of_death(const Parameters &parameters);
  
  // Schedule mosquito's next strain clearance
  void schedule_m_day_of_blood_meal(Parameters &parameters);
  
  // Set day of next event
  void schedule_m_day_of_next_event();
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // UPDATERS - functions that happen at the end of a day, i.e. ageing, dying, event handling etc. 
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Kill mosquito, i.e. reset age to 0, infections to 0, state to susceptible, ino barcodes etc
  bool die(Parameters &parameters);
  
  // Daily update function, i.e. check events through event_handle, and returns either -1 or moqsuito ID if mosquito is biting today
  bool update(Parameters &parameters);
  
  // Event handle, i.e. mosquito death, state change, biting handle
  bool event_handle(Parameters &parameters);
  
  
};

#endif
