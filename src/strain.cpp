#include "stdafx.h"
#include "strain.h"

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Default class constructor which will inialise a random strain
Strain::Strain() : 
  m_barcode(std::bitset<barcode_length>(runiform_int_1(1, barcode_length_max_bits)))
{
  
}

// Known barcode constructor
Strain::Strain(const barcode_t &barcode) :
  m_barcode(barcode)
{
}

// Known barcode and infection status constructor
Strain::Strain(const barcode_t &barcode, const Strain::InfectionStatus &infectionStatus) :
  m_barcode(barcode),
  m_strain_infection_status(infectionStatus)
{
}

// Known barcode, infection status and day of infection status change outcome constructor
Strain::Strain(const barcode_t &barcode, const Strain::InfectionStatus &infectionStatus, const int &dayOfInfectionStatusChange) :
  m_barcode(barcode),
  m_strain_infection_status(infectionStatus),
  m_day_of_strain_infection_status_change(dayOfInfectionStatusChange)
{
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// STATIC INITIALIZERS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Initialise static class const vectors
const std::vector<Strain::InfectionStatus> Strain::m_transition_vector{ SUSCEPTIBLE, DISEASED, ASYMPTOMATIC, SUBPATENT, TREATED, PROPHYLAXIS };

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// SEMI-SETTERS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void Strain::set_random_barcode()
{
  m_barcode = std::bitset<barcode_length>(runiform_int_1(1, barcode_length_max_bits));
}

// Generate a random recombinant barcode given two barcodes
barcode_t Strain::generate_random_barcode()
{
  return(std::bitset<barcode_length>(runiform_int_1(1, barcode_length_max_bits)));
}


// Generate a random recombinant barcode given two barcodes
barcode_t Strain::generate_recombinant_barcode(barcode_t x, barcode_t y)
{
  // find the different positions, then where these cross with a random other barcode which represents seperate segragation, i.e. 50% recombination, 
  // and then where these differ with x
  return( ( (x ^ y) & generate_random_barcode()) ^ x );
}

