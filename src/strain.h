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

template <class T>
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
  
  boost::dynamic_bitset<T> m_barcode;									// barcode sequence
  InfectionStatus m_strain_infection_status;					// infection status associated with a strain
  int m_day_of_strain_infection_status_change;				// day that strain would move infection status
  int m_day_of_acquisition;                           // day that the strain was acquired by an individual/mosquito
  
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
  
  // Known barcode, infection status and day of infection status change outcome constructor and day of acquisition
  Strain(boost::dynamic_bitset<> barcode, const Strain::InfectionStatus &infectionStatus, 
         const int &dayOfInfectionStatusChange, const int & day_of_acquisition);
  
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // GETTERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Get barcode
  boost::dynamic_bitset<T> get_m_barcode() { return(m_barcode); }
  
  // Get strain infection status
  InfectionStatus get_m_strain_infection_status() { return(m_strain_infection_status); }
  
  // Get day of strain infection status
  int get_m_day_of_strain_infection_status_change() { return(m_day_of_strain_infection_status_change); }
  
  // Get day of strain acquisition
  int get_m_day_of_strain_acquisition() { return(m_day_of_acquisition); }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // SETTERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Set barcode pointer
  void set_m_barcode(boost::dynamic_bitset<> x) { m_barcode = x; }
  
  // Set strain infection status
  void set_m_strain_infection_status(InfectionStatus x) { m_strain_infection_status = x; }
  
  // Set strain infection status day change
  void set_m_day_of_strain_infection_status_change(int x) { m_day_of_strain_infection_status_change = x; }
  
  // Set strain acquisition day
  void set_m_day_of_strain_acquisition(int x) { m_day_of_acquisition = x; }
  
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

  // Generate the next identity barcode, i.e. if it was 0101 it will now be 1111, if num_loci = 2
  static boost::dynamic_bitset<> generate_next_ibd_barcode();
  
  // Generate a random recombinant barcode given two barcodes
  static boost::dynamic_bitset<> generate_recombinant_barcode(boost::dynamic_bitset<> x, boost::dynamic_bitset<> y);
  
  // Generate a random barcode given probability of each SNP, i.e. PLAF
  static boost::dynamic_bitset<> generate_random_barcode_given_SNP_frequencies(std::vector<double> &x);
  
  // Create a bitset by replicating each bit n times
  static boost::dynamic_bitset<> replicate_by_bit(boost::dynamic_bitset<> x, unsigned int n);
  
  // Turns our ibd barcode into a vector of the ints making it up
  static std::vector<unsigned long> ibd_barcode_to_integer_vector(boost::dynamic_bitset<> x);
  
  // distances between one bitset and a vector range 
  static unsigned int distance_of_bitset_a_and_x(boost::dynamic_bitset<> a, 
                                                  std::vector<boost::dynamic_bitset<> >::const_iterator start, 
                                                  std::vector<boost::dynamic_bitset<> >::const_iterator end);
  
  // distances between one bitset and a vector range of bitsets
  static unsigned int ibd_distance_of_bitset_a_and_x(boost::dynamic_bitset<> a, 
                                                      std::vector<boost::dynamic_bitset<> >::const_iterator start, 
                                                      std::vector<boost::dynamic_bitset<> >::const_iterator end);
  
  // distances between one bitset and a vector of vectors of bitsets
  static double distance_of_bitset_a_and_vec_x(boost::dynamic_bitset<> a, 
                                                std::vector<std::vector< boost::dynamic_bitset<> > >::const_iterator start, 
                                                std::vector<std::vector< boost::dynamic_bitset<> > >::const_iterator end);
  
  // mean distance between all bitsets within a vector of bitsets
  static double distance_mean_within_bitsets(std::vector<boost::dynamic_bitset<> > x, unsigned int bl);

  // mean ibd distance between all bitsets within a vector of bitsets
  static double ibd_distance_mean_within_bitsets(std::vector<boost::dynamic_bitset<> > x, unsigned int bl);
  
};

#include "Strain.inl"

#endif
