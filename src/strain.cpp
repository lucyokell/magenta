#include "strain.h"

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Default class constructor which will inialise a random strain
Strain::Strain() :
  m_barcode(Strain::generate_next_barcode())
{
}

// Known barcode constructor
Strain::Strain(boost::dynamic_bitset<> barcode) :
  m_barcode(barcode)
{
}

// Known barcode and infection status constructor
Strain::Strain(boost::dynamic_bitset<> barcode, const Strain::InfectionStatus &infectionStatus) :
  m_barcode(barcode),
  m_strain_infection_status(infectionStatus)
{
}

// Known barcode, infection status and day of infection status change outcome constructor
Strain::Strain(boost::dynamic_bitset<> barcode, const Strain::InfectionStatus &infectionStatus, const int &dayOfInfectionStatusChange) :
  m_barcode(barcode),
  m_strain_infection_status(infectionStatus),
  m_day_of_strain_infection_status_change(dayOfInfectionStatusChange)
{
}

// Known barcode, infection status and day of infection status change outcome constructor and acquisition dat
Strain::Strain(boost::dynamic_bitset<> barcode, const Strain::InfectionStatus &infectionStatus, 
               const int &dayOfInfectionStatusChange, const int & day_of_acquisition) :
  m_barcode(barcode),
  m_strain_infection_status(infectionStatus),
  m_day_of_strain_infection_status_change(dayOfInfectionStatusChange),
  m_day_of_acquisition(day_of_acquisition)
{
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// STATIC INITIALIZERS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Initialise static class const vectors
const std::vector<Strain::InfectionStatus> Strain::m_transition_vector{ SUSCEPTIBLE, DISEASED, ASYMPTOMATIC, SUBPATENT, TREATED, PROPHYLAXIS };

// temp barcodes for speed
boost::dynamic_bitset<> Strain::temp_barcode(Parameters::g_barcode_length);
boost::dynamic_bitset<> Strain::temp_identity_barcode(Parameters::g_ibd_length);
boost::dynamic_bitset<> Strain::temp_crossovers(Parameters::g_num_loci);

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// STRAIN - UTILS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Generate next barcode. For non IBD this is just a random. IBD its the next identity
boost::dynamic_bitset<> Strain::generate_next_barcode()
{
  // Match barcode type
  switch (Parameters::g_barcode_type)
  {
  case Parameters::ORDINARY:
    return(generate_random_barcode_given_SNP_frequencies(Parameters::g_plaf));
  case Parameters::IBD:
    return(generate_next_ibd_barcode());
  default:
    Rcpp::stop("Unrecognised barcode_type");
  break;
  }
  
}

// Generate a random barcode fiven plaf
boost::dynamic_bitset<> Strain::generate_random_barcode()
{
  // Match barcode type
  switch (Parameters::g_barcode_type)
  {
  case Parameters::ORDINARY:
    return(generate_random_barcode_given_SNP_frequencies(Parameters::g_plaf));
  case Parameters::IBD:
    return(generate_random_ibd_barcode());
  default:
    Rcpp::stop("Unrecognised barcode_type");
  break;
  }
}

// Generate a random ordinary barcode, i.e. 1 bit per loci 
boost::dynamic_bitset<> Strain::generate_random_ordinary_barcode()
{
  for(unsigned int i = 0; i < Parameters::g_barcode_length ; i++){
    Strain::temp_barcode[i] = rbernoulli1(0.5);
  }
  return(Strain::temp_barcode);
}

// Generate a random identity barcode, i.e. bitset of length equal to ibd_length
boost::dynamic_bitset<> Strain::generate_random_identity_barcode()
{
  for(unsigned int i = 0; i < Parameters::g_ibd_length ; i++){
    Strain::temp_identity_barcode[i] = rbernoulli1(0.5);
  }
  
  return(Strain::temp_identity_barcode);
}

// Generate a random identity barcode, i.e. where the barcode represents num_loci * ibd_length
boost::dynamic_bitset<> Strain::generate_random_ibd_barcode()
{
  unsigned int count = 0;
 
  // create and fill the temp barcode
  for(unsigned int i = 0; i < Parameters::g_num_loci ; i++) {
    
    Strain::temp_identity_barcode = generate_random_identity_barcode();
    
    for(unsigned int j = 0; j < Parameters::g_ibd_length ; j++) {
      
      Strain::temp_barcode[count++] =  Strain::temp_identity_barcode[j];
    }
  }
  
  return(Strain::temp_barcode);
}

// Generate a random identity barcode, i.e. where the barcode represents num_loci * ibd_length
boost::dynamic_bitset<> Strain::generate_next_ibd_barcode()
{
  // create the next identity bitset
  Strain::temp_identity_barcode =  boost::dynamic_bitset<>(Parameters::g_ibd_length, Parameters::g_identity_id++);
  
  unsigned int count = 0;
  
  // then fill in th temp by repeating this
  for(unsigned int i = 0; i < Parameters::g_num_loci ; i++) {
    
    for(unsigned int j = 0; j < Parameters::g_ibd_length ; j++) {
      
      Strain::temp_barcode[count++] =  Strain::temp_identity_barcode[j];
    }
  }
  
  return(Strain::temp_barcode);
}

// Generate a crossovers bitset using the prob_crossovers, and is thus length of num_loci
boost::dynamic_bitset<> generate_crossover()
{
  for(unsigned int i = 0; i < Parameters::g_num_loci ; i++){
    Strain::temp_crossovers[i] = rbernoulli1(Parameters::g_prob_crossover[i]);
  }
  return(Strain::temp_crossovers);
}

// Generate a random recombinant barcode given two barcodes
boost::dynamic_bitset<> Strain::generate_recombinant_barcode(boost::dynamic_bitset<> x, boost::dynamic_bitset<> y)
{
  // find the different positions, then where these cross with another 
  // bitset created using the prob_crossover which represents chance of segragation. 
  // i.e. 0.5 equals indepedent segregation, 0 equals loci next to each other and
  // 0% recombination. Then where these differ with x. 
  // For the IBD style, our random barcode
  
  // Match barcode type
  switch (Parameters::g_barcode_type)
  {
  case Parameters::ORDINARY:
    return( ( (x ^ y) & generate_random_barcode_given_SNP_frequencies(Parameters::g_prob_crossover) ) ^ x );
  case Parameters::IBD:
    return( ( (x ^ y) & replicate_by_bit(generate_crossover(),Parameters::g_ibd_length) ) ^ x );
  default:
    Rcpp::stop("Unrecognised barcode_type");
  break;
  }

}

// Generate a random barcode given probability of each SNP, i.e. PLAF
boost::dynamic_bitset<> Strain::generate_random_barcode_given_SNP_frequencies(std::vector<double> x)
{
  // create and fill the temp barcode
  for(unsigned int i = 0; i < Parameters::g_barcode_length ; i++){
    Strain::temp_barcode[i] = rbernoulli1(x[i]);
  }
  
  // return barcode
  return(Strain::temp_barcode);
}

// Stretches a bitset by replicating each bit n times
boost::dynamic_bitset<> Strain::replicate_by_bit(boost::dynamic_bitset<> x, unsigned int n)
{
  
  unsigned int count = 0;
  
  // create and fill the temp barcode
  for(unsigned int i = 0; i < Parameters::g_num_loci ; i++) {
    
    for(unsigned int j = 0; j < n ; j++) {
      
      Strain::temp_barcode[count++] = x[i];
      
    }
  }
  
  // return barcode
  return(Strain::temp_barcode);
}

// PLAF test
// [[Rcpp::export]]
SEXP test_barcode_from_PLAF(Rcpp::NumericVector plaf, unsigned int n)
{
  Parameters::g_barcode_length = n;
  Strain::temp_barcode = boost::dynamic_bitset<>(Parameters::g_barcode_length);
  std::vector<double> plaf_c = Rcpp::as<std::vector<double> >(plaf);
  boost::dynamic_bitset<> b = Strain::generate_random_barcode_given_SNP_frequencies(plaf_c);
  return(bitset_to_sexp(b));
}

// IBD recombinant test
// [[Rcpp::export]]
SEXP test_recombinant_with_ibd(SEXP barcode_1,
                               SEXP barcode_2,
                               unsigned int bl, unsigned int nl,
                               unsigned int ib, Rcpp::NumericVector pc
                               )
{
  Parameters::g_barcode_length = bl;
  Parameters::g_num_loci = nl;
  Parameters::g_ibd_length = ib;
  Parameters::g_prob_crossover = Rcpp::as<std::vector<double> >(pc);
  Parameters::g_barcode_type = Parameters::IBD;
  
  Strain::temp_barcode = boost::dynamic_bitset<>(Parameters::g_barcode_length);
  Strain::temp_crossovers = boost::dynamic_bitset<>(Parameters::g_num_loci);
  
  boost::dynamic_bitset<> barcode_a = sexp_to_bitset(barcode_1, bl);
  boost::dynamic_bitset<> barcode_b = sexp_to_bitset(barcode_2, bl);
  Strain::temp_barcode = Strain::generate_recombinant_barcode(barcode_a, barcode_b);
  return(bitset_to_sexp(Strain::temp_barcode));
}

