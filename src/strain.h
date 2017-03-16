//
//  MAGENTA
//  strain.h
//
//  Created: OJ on 13/01/2017
//	Most recent edits: OJ on 10/03/2017
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Class describing a parasite strain.
//
// ---------------------------------------------------------------------------


#ifndef STRAIN_H
#define STRAIN_H

#include <vector>
#include <queue>
#include <bitset>
#include <random>
#include "probability.h"

const int barcode_length = 24; // TODO: Check whether this is the best way to declare this - perhaps should be a parameter, thus not const
const int barcode_length_max_bits = static_cast<int>(pow(2, barcode_length));
using barcode_t = std::bitset<barcode_length>;


class Strain {
  
public:
  enum InfectionStatus
  {
    SUSCEPTIBLE,	// 0
    DISEASED,		// 1
    ASYMPTOMATIC,	// 2
    SUBPATENT,		// 3
    TREATED,		// 4
    PROPHYLAXIS,	// 5
    NUMBER_OF_STATES = 6
  };
  
  
  // Infection transition options. 
  const static std::vector<InfectionStatus> m_transition_vector;
  
private:
  
  barcode_t m_barcode;										// barcode sequence
  InfectionStatus m_strain_infection_status;					// infection status associated with a strain
  int m_day_of_strain_infection_status_change;				// day that strain would move infection status
  
  // TODO: Phenotype parameters to do here for resistance work
  
  
public:
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // CONSTRUCTORS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Default class constructor which will inialise a random strain
  Strain();
  
  // Known barcode constructor
  Strain(const barcode_t &barcode);
  
  // Known barcode and infection status constructor
  Strain(const barcode_t &barcode, const Strain::InfectionStatus &infectionStatus);
  
  // Known barcode, infection status and day of infection status change outcome constructor
  Strain(const barcode_t &barcode, const Strain::InfectionStatus &infectionStatus, const int &dayOfInfectionStatusChange);
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // GETTERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Get barcode
  barcode_t get_m_barcode() { return(m_barcode); }
  
  // Get strain infection status
  InfectionStatus get_m_strain_infection_status() { return(m_strain_infection_status); }
  
  // Get day of strain infection status
  int get_m_day_of_strain_infection_status_change() { return(m_day_of_strain_infection_status_change); }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // SETTERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Set random barcode based on barcode length only
  void set_random_barcode();
  
  // Set barcode pointer
  void set_m_barcode(barcode_t x) { m_barcode = x; }
  
  // Set strain infection status
  void set_m_strain_infection_status(InfectionStatus x) { m_strain_infection_status = x; }
  
  // Set strain infection status day change
  void set_m_day_of_strain_infection_status_change(int x) { m_day_of_strain_infection_status_change = x; }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // CLASS STATIC FUNCTIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Generate a random barcode
  static barcode_t generate_random_barcode();
  
  // Generate a random recombinant barcode given two barcodes
  static barcode_t generate_recombinant_barcode(barcode_t x, barcode_t y);
  
};

#endif
