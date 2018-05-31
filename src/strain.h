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
#include "parameters.h"
#include "util.h"



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

  // temporary barcodes for all purposes
  static boost::dynamic_bitset<> temp_barcode;
  static boost::dynamic_bitset<> temp_identity_barcode;
  static boost::dynamic_bitset<> temp_crossovers;
  
private:
  
  boost::dynamic_bitset<> m_barcode;									// barcode sequence
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
  Strain(boost::dynamic_bitset<> barcode);
  
  // Known barcode and infection status constructor
  Strain(boost::dynamic_bitset<> barcode, const Strain::InfectionStatus &infectionStatus);
  
  // Known barcode, infection status and day of infection status change outcome constructor
  Strain(boost::dynamic_bitset<> barcode, const Strain::InfectionStatus &infectionStatus, const int &dayOfInfectionStatusChange);
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // GETTERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Get barcode
  boost::dynamic_bitset<> get_m_barcode() { return(m_barcode); }
  
  // Get strain infection status
  InfectionStatus get_m_strain_infection_status() { return(m_strain_infection_status); }
  
  // Get day of strain infection status
  int get_m_day_of_strain_infection_status_change() { return(m_day_of_strain_infection_status_change); }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // SETTERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Set barcode pointer
  void set_m_barcode(boost::dynamic_bitset<> x) { m_barcode = x; }
  
  // Set strain infection status
  void set_m_strain_infection_status(InfectionStatus x) { m_strain_infection_status = x; }
  
  // Set strain infection status day change
  void set_m_day_of_strain_infection_status_change(int x) { m_day_of_strain_infection_status_change = x; }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // CLASS STATIC FUNCTIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Generate next barcode. For non IBD this is just a random. IBD its the next identity
  static boost::dynamic_bitset<> generate_next_barcode();
  
  // Generate a random barcode
  static boost::dynamic_bitset<> generate_random_barcode();
  
  // Generate a random ordinary barcode, i.e. 1 bit per loci 
  static boost::dynamic_bitset<> generate_random_ordinary_barcode();
  
  // Generate a random ordinary barcode, i.e. 1 bit per loci 
  static boost::dynamic_bitset<> generate_random_ibd_barcode();
  
  // Generate a random identity barcode, i.e. where the barcode represents num_loci * ibd_length
  static boost::dynamic_bitset<> generate_random_identity_barcode();

  // Generate the next identity barcode, i.e. if it was 0101 it will now be 111, if num_loci = 2
  static boost::dynamic_bitset<> generate_next_ibd_barcode();
  
  // Generate a random recombinant barcode given two barcodes
  static boost::dynamic_bitset<> generate_recombinant_barcode(boost::dynamic_bitset<> x, boost::dynamic_bitset<> y);
  
  // Generate a random barcode given probability of each SNP, i.e. PLAF
  static boost::dynamic_bitset<> generate_random_barcode_given_SNP_frequencies(std::vector<double> x);
  
  // Create a bitset by replicating each bit n times
  static boost::dynamic_bitset<> replicate_by_bit(boost::dynamic_bitset<> x, unsigned int n);
  
  
};

#endif
