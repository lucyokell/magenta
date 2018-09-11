#include "universe_utils.h"

// Create universe structure for all important variables
struct Universe {
  // Human storage
  std::vector<Person> population;
  std::vector<double> psi_vector;
  std::vector<double> zeta_vector;
  std::vector<double> pi_vector;
  // Mosquito storage
  std::vector<Mosquito> scourge;
  // Parameter storage
  Parameters parameters;
};

// q function for detectability
double q_fun(const Parameters &parms, unsigned int age, double ID) {
  
  double fd =  1 - ((1-parms.g_fD0) / (1 + pow((age/parms.g_aD),parms.g_gD)));
  return(parms.g_d1 + ((1-parms.g_d1) / (1 + pow((ID/parms.g_ID0),parms.g_kD)*fd)));
}


//' Returns the population's parasite genetics for ibd style summarised by pibd for given sample size and state
//'
//' @param paramList paramList containing statePtr, sample_size, and sample_states
//' @return list of population information
//' @export
// [[Rcpp::export]]
Rcpp::List population_get_genetics_ibd_df_n(Rcpp::List paramList) {
  
  // create universe and do parameter conversions
  Rcpp::XPtr<Universe> u_ptr = Rcpp::as<Rcpp::XPtr<Universe> > (paramList["statePtr"]);
  unsigned int sample_size = Rcpp::as<unsigned int>(paramList["sample_size"]);
  std::vector<unsigned int> sample_states = Rcpp::as<std::vector<unsigned int> >(paramList["sample_states"]);
  
  unsigned int pop_size = u_ptr->population.size();
  // handle if we are subsampling
  std::vector<int> pop_sample(pop_size);  
  std::iota(pop_sample.begin(), pop_sample.end(), 0);
  std::vector<unsigned int> all_states = {0, 1, 2, 3, 4, 5};
  
  // if we are subsampling states
  if(all_states != sample_states){
    pop_sample.clear();
    for(auto p : u_ptr->population) {
      if(std::find(sample_states.begin(), sample_states.end(), static_cast<unsigned int>(p.get_m_infection_state())) != sample_states.end()) {
        pop_sample.emplace_back(p.get_m_person_ID());
      }
    }
  }
  
  // if we are choosing a sample
  if(sample_size) {
    pop_sample = sample_without_replacement(pop_sample, sample_size);
  }
  pop_size = pop_sample.size();
  
  
  std::vector<double> pibd_d(pop_size, 0);
  std::vector<double> pibd_within_d(pop_size, 0);
  std::vector<double> pibd(pop_size, 0);
  std::vector<double> pibd_within(pop_size, 0);
  std::vector<unsigned int> ages(pop_size, 0);
  std::vector<unsigned int> last_treatment(pop_size, 0);
  std::vector<unsigned int> states(pop_size, 0);
  std::vector<std::vector<SEXP> > barcodes(pop_size);
  std::vector<std::vector<unsigned int> > barcode_states(pop_size);
  std::vector<std::vector<boost::dynamic_bitset<> > > bitsets(pop_size);
  std::vector<std::vector<boost::dynamic_bitset<> > > bitsets_d(pop_size);
  
  
  // temps
  double q_i;
  unsigned int el_i;
  std::vector<Strain> strains_i;
  std::vector<boost::dynamic_bitset<> > bitsets_i;
  std::vector<boost::dynamic_bitset<> > bitsets_d_i;
  std::vector<unsigned int> barcode_states_i;
  
  for (unsigned int el = 0; el < pop_size; el++) {
    
    // who are we sampling
    el_i = static_cast<unsigned int>(pop_sample[el]);
    
    // simple assignments first
    ages[el] = u_ptr->population[el_i].get_m_person_age();
    states[el] = static_cast<unsigned int>(u_ptr->population[el_i].get_m_infection_state());
    last_treatment[el] = u_ptr->parameters.g_current_time - u_ptr->population[el_i].get_m_day_last_treated();
    
    // what's their disease detectability
    q_i = q_fun(u_ptr->parameters, ages[el], u_ptr->population[el_i].get_m_ID());
    
    // clear temps
    bitsets_i.clear();
    bitsets_d_i.clear();
    strains_i.clear();
    barcode_states_i.clear();
    
    // get their strains
    strains_i = u_ptr->population[el_i].get_m_active_strains();
    
    // reserve space
    bitsets_i.reserve(strains_i.size());
    barcodes[el].reserve(strains_i.size());
    barcode_states_i.reserve(strains_i.size());
    

    
    // and then their barcode bitsets for all D and those detected in A
    for(auto s : strains_i) {
      if(s.get_m_strain_infection_status() == Strain::DISEASED || s.get_m_strain_infection_status() == Strain::TREATED) {
        bitsets_i.emplace_back(s.get_m_barcode());
        barcodes[el].emplace_back(bitset_to_sexp(bitsets_i.back()));
        bitsets_d_i.emplace_back(s.get_m_barcode());
        barcode_states_i.emplace_back(static_cast<unsigned int>(s.get_m_strain_infection_status()));
      } 
      if(s.get_m_strain_infection_status() == Strain::ASYMPTOMATIC){
        if(rbernoulli1(q_i)) {
          bitsets_i.emplace_back(s.get_m_barcode());
          barcodes[el].emplace_back(bitset_to_sexp(bitsets_i.back()));
          barcode_states_i.emplace_back(static_cast<unsigned int>(s.get_m_strain_infection_status()));
        }  
      }
    }
    
    // store these for later
    bitsets[el] = bitsets_i;
    bitsets_d[el] = bitsets_d_i;
    barcode_states[el] = barcode_states_i;
    
    // whats the pibd within the individuals bitsets that could be detected and within just those at high parasitaemia
    pibd_within[el] = Strain::ibd_distance_mean_within_bitsets(bitsets_i, u_ptr->parameters.g_barcode_length);
    pibd_within_d[el] = Strain::ibd_distance_mean_within_bitsets(bitsets_d_i, u_ptr->parameters.g_barcode_length);
    
  }
  
  // now work out the pibd within the population
  double pibd_i = 0.0;
  double pibd_d_i = 0.0;
  std::vector<int> contacts(pop_size,0);
  std::vector<int> contacts_d(pop_size,0);
  
  
  for(unsigned int i = 0; i < (pop_size-1); i++) {
    if(bitsets[i].size()){
      for(unsigned int j = i+1; j < pop_size; j++) {
        if(bitsets[j].size()) {
          pibd_i = pibd_d_i = 0;
          for(unsigned int k1 = 0; k1 < bitsets[i].size(); k1++){
            pibd_i = pibd_i + Strain::ibd_distance_of_bitset_a_and_x(bitsets[i][k1], bitsets[j].begin(), bitsets[j].end());
          }
          for(unsigned int k2 = 0; k2 < bitsets_d[i].size(); k2++){
            pibd_d_i = pibd_d_i + Strain::ibd_distance_of_bitset_a_and_x(bitsets_d[i][k2], bitsets_d[j].begin(), bitsets_d[j].end());
          }
          
          // running total of the number of bitset comparisons we have made for each individual
          contacts[i] += bitsets[i].size() * bitsets[j].size();
          contacts[j] += bitsets[i].size() * bitsets[j].size();
          contacts_d[i] += bitsets_d[i].size() * bitsets_d[j].size();
          contacts_d[j] += bitsets_d[i].size() * bitsets_d[j].size();
          
          // and the running total for the pibd
          pibd[i] += pibd_i;
          pibd[j] += pibd_i;
          pibd_d[i] += pibd_d_i;
          pibd_d[j] += pibd_d_i;
        }
      }
    }
  }
  
  // and now divide accordingly
  for(unsigned int c = 0; c < pop_size; c++) {
    pibd[c] /= (contacts[c] * u_ptr->parameters.g_num_loci);
    pibd_d[c] /= (contacts_d[c] * u_ptr->parameters.g_num_loci);
  }
  
  return(Rcpp::List::create(
      Rcpp::Named("pibd_within") = pibd_within,
      Rcpp::Named("pibd_within_d") = pibd_within_d, 
      Rcpp::Named("pibd") = pibd,
      Rcpp::Named("pibd_d") = pibd_d, 
      Rcpp::Named("age")= ages, 
      Rcpp::Named("state")=states,
      Rcpp::Named("barcodes")=barcodes,
      Rcpp::Named("barcode_states")=barcode_states));
  
  
}


//' Returns the population's parasite genetics summarised by coi for given sample size and state
//'
//' @param paramList paramList containing statePtr, sample_size, and sample_states
//' @return list of population information
//' @export
// [[Rcpp::export]]
Rcpp::List population_get_genetics_df_n(Rcpp::List paramList) {
  
  // create universe and do parameter conversions
  Rcpp::XPtr<Universe> u_ptr = Rcpp::as<Rcpp::XPtr<Universe> > (paramList["statePtr"]);
  unsigned int sample_size = Rcpp::as<unsigned int>(paramList["sample_size"]);
  std::vector<unsigned int> sample_states = Rcpp::as<std::vector<unsigned int> >(paramList["sample_states"]);
  
  unsigned int pop_size = u_ptr->population.size();
  
  // handle if we are subsampling
  std::vector<int> pop_sample(pop_size);  
  std::iota(pop_sample.begin(), pop_sample.end(), 0);
  std::vector<unsigned int> all_states = {0, 1, 2, 3, 4, 5};
  
  // if we are subsampling states
  if(all_states != sample_states){
    pop_sample.clear();
    for(unsigned int p = 0 ; p < u_ptr->population.size(); p++) {
      if(std::find(sample_states.begin(), sample_states.end(), static_cast<unsigned int>(u_ptr->population[p].get_m_infection_state())) != sample_states.end()) {
        pop_sample.emplace_back(p);
      }
    }
  }
  
  // if we are choosing a sample
  if(sample_size) {
    pop_sample = sample_without_replacement(pop_sample, sample_size);
    std::cout << "sampling" << std::endl;
  }
  pop_size = pop_sample.size();
  
  std::vector<unsigned int> coi(pop_size, 0);
  std::vector<unsigned int> coi_detected_micro(pop_size, 0);
  std::vector<unsigned int> moi(pop_size, 0);
  std::vector<unsigned int> ages(pop_size, 0);
  std::vector<unsigned int> last_treatment(pop_size, 0);
  std::vector<unsigned int> states(pop_size, 0);
  std::vector<std::vector<SEXP> > barcodes(pop_size);
  std::vector<std::vector<unsigned int> > barcode_states(pop_size);
  
  // temps
  double q_i;
  unsigned int i;
  std::vector<Strain> strains_i;
  std::vector<boost::dynamic_bitset<> > bitsets_u_i;
  std::vector<boost::dynamic_bitset<> > bitsets_not_u_i;
  std::vector<SEXP> barcodes_i;
  std::vector<unsigned int> barcode_states_i;
  
  rcpp_out(u_ptr->parameters.g_h_quiet_print, "preloop. pop_size = " + std::to_string(pop_size));
  
  for (unsigned int el = 0; el < pop_size; el++) {
    
    // who are we sampling
    i = pop_sample[el];
    rcpp_out(u_ptr->parameters.g_h_quiet_print, "uu_1: " + std::to_string(el) + " " + std::to_string(i));
    
    // simple assignments first
    ages[el] = u_ptr->population[i].get_m_person_age();
    states[el] = static_cast<unsigned int>(u_ptr->population[i].get_m_infection_state());
    last_treatment[el] = u_ptr->parameters.g_current_time - u_ptr->population[i].get_m_day_last_treated();
    
    //rcpp_out(u_ptr->parameters.g_h_quiet_print, "uu_2: " + std::to_string(el) + " " + std::to_string(i));
    
    // clear temps
    barcodes_i.clear();
    barcode_states_i.clear();
    bitsets_u_i.clear();
    bitsets_not_u_i.clear();
    strains_i.clear();
    
    //rcpp_out(u_ptr->parameters.g_h_quiet_print, "uu_3: " + std::to_string(el) + " " + std::to_string(i));
    
    // reserve space
    strains_i = u_ptr->population[i].get_m_active_strains();
    bitsets_u_i.reserve(strains_i.size());
    barcodes_i.reserve(strains_i.size());
    bitsets_not_u_i.reserve(strains_i.size());
    barcode_states_i.reserve(strains_i.size());
    
    rcpp_out(u_ptr->parameters.g_h_quiet_print, "uu_4: " + std::to_string(el) + " " + std::to_string(i));
    
    //rcpp_out(u_ptr->parameters.g_h_quiet_print, "uu_5: " + std::to_string(el) + " " + std::to_string(i));
    
    
    //rcpp_out(u_ptr->parameters.g_h_quiet_print, "uu_6: " + std::to_string(el) + " " + std::to_string(i));
    
    // for coi let's first get the barcodes for the non U state 
    for(unsigned int bu = 0; bu < strains_i.size(); bu++) {
      if(strains_i[bu].get_m_strain_infection_status() != Strain::SUBPATENT) {
        bitsets_not_u_i.emplace_back(strains_i[bu].get_m_barcode()); 
        barcode_states_i.emplace_back(static_cast<unsigned int> (strains_i[bu].get_m_strain_infection_status()));
        barcodes_i.emplace_back(bitset_to_sexp(bitsets_not_u_i.back())); 
      } 
      bitsets_u_i.emplace_back(strains_i[bu].get_m_barcode()); 
    }
    barcodes[el] = barcodes_i;
    barcode_states[el] = barcode_states_i;
    
    //rcpp_out(u_ptr->parameters.g_h_quiet_print, "uu_7: " + std::to_string(el) + " " + std::to_string(i));
    
    // and then the unique size
    q_i = q_fun(u_ptr->parameters, ages[el], u_ptr->population[i].get_m_ID());
    std::sort(bitsets_not_u_i.begin(), bitsets_not_u_i.end());
    coi_detected_micro[el] = ceil(q_i * (std::unique(bitsets_not_u_i.begin(), bitsets_not_u_i.end()) - bitsets_not_u_i.begin()));
    
    //rcpp_out(u_ptr->parameters.g_h_quiet_print, "uu_8: " + std::to_string(el) + " " + std::to_string(i));
    
    // similar for all
    std::sort(bitsets_u_i.begin(), bitsets_u_i.end());
    coi[el] = std::unique(bitsets_u_i.begin(), bitsets_u_i.end()) - bitsets_u_i.begin();
    
    rcpp_out(u_ptr->parameters.g_h_quiet_print, "uu_9: " + std::to_string(el) + " " + std::to_string(i));
    
  }
  
  rcpp_out(u_ptr->parameters.g_h_quiet_print, "postloop");
  
  return(Rcpp::List::create(
      Rcpp::Named("coi") = coi,
      Rcpp::Named("coi_detected_micro") = coi_detected_micro, 
      Rcpp::Named("age")= ages, 
      Rcpp::Named("state")=states,
      Rcpp::Named("last_treatment") = last_treatment, 
      Rcpp::Named("barcodes")=barcodes,
      Rcpp::Named("barcode_states")=barcode_states));
  
}