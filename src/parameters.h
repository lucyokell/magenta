//
//  MAGENTA
//  parameters.h
//
//  Created: OJ on 17/01/2017
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Class describing all global parameters.
//
// ---------------------------------------------------------------------------

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>
#include <vector>
#include <queue>
#include <bitset>
#include <boost/dynamic_bitset.hpp>

class Parameters {
  
public:
  
  int g_current_time;
  int g_calendar_day;
  std::vector<double> g_theta;
  double g_years;
  
  // spatial
  enum g_spatial_type_enum{
    NON, //0
    ISLAND, //1
    METAPOPULATION, //2
    NUMBER_OF_SPATIAL_TYPE_OPTIONS = 3
  };
  
  Parameters::g_spatial_type_enum g_spatial_type;
  std::vector<int> g_spatial_imported_cotransmission_frequencies;
  std::vector<int> g_spatial_imported_oocyst_frequencies;
  std::vector<int> g_spatial_exported_cotransmission_frequencies;
  std::vector<int> g_spatial_exported_oocyst_frequencies;
  unsigned int g_spatial_total_exported_barcodes;
  unsigned int g_spatial_total_exported_oocysts;
  unsigned int g_spatial_total_imported_barcodes;
  unsigned int g_spatial_total_imported_oocysts;
  unsigned int g_spatial_imported_barcode_counter;
  unsigned int g_spatial_imported_oocyst_counter;
  unsigned int g_spatial_exported_barcode_counter;
  unsigned int g_spatial_exported_oocyst_counter;
  std::vector<boost::dynamic_bitset<> > g_spatial_exported_barcodes;
  std::vector<boost::dynamic_bitset<> > g_spatial_imported_barcodes;
  std::vector<boost::dynamic_bitset<> > g_spatial_exported_oocysts;
  std::vector<boost::dynamic_bitset<> > g_spatial_imported_oocysts;
  
  // demographic parameters;
  unsigned int g_N;
  int g_max_age;
  int g_average_age;
  // epidemiological parameters;
  double g_EIR;
  double g_a0;
  double g_rho;
  double g_zeta_meanlog;
  double g_zeta_sdlog;
  double g_importation_human;
  double g_importation_mosquito;
  // diagnostic parameters;
  double g_ft;
  // entomological parameters;
  double g_mu0;
  double g_mean_mosquito_age;
  double g_beta_gradient;
  double g_beta_intercept;
  double g_ak;
  double g_Q0;
  // entomological time keeping variables
  int g_mosquito_deficit;
  unsigned int g_scourge_today;
  int g_mean_mv;
  std::vector<int> g_mosquito_next_biting_day_vector;
  int g_mosquito_biting_counter;
  int g_max_mosquito_biting_counter;
  // delays and durations;
  double g_delay_mos;
  double g_delay_gam;
  double g_dur_E;
  double g_dur_T;
  double g_dur_D;
  double g_dur_U;
  double g_dur_P;
  double g_dur_A;
  double g_dur_AU;
  // Pre-erythrocytic immunity parameters;
  double g_d1;
  double g_dID;
  double g_ID0;
  double g_kD;
  double g_uD;
  double g_aD;
  double g_fD0;
  double g_gD;
  // Blood stage immunity parameters;
  double g_b0;
  double g_b1;
  double g_dB;
  double g_IB0;
  double g_kB;
  double g_uB;
  // Acquired and maternal immunity parameters;
  double g_phi0;
  double g_phi1;
  double g_dCA;
  double g_IC0;
  double g_kC;
  double g_uCA;
  double g_PM;
  double g_dCM;
  double g_mean_maternal_immunity;
  double g_sum_maternal_immunity;
  int g_total_mums;
  // contributions to infectious reservoir by state and age;
  double g_gamma1;
  double g_cD;
  double g_cT;
  double g_cU;
  
  // barcode parameters
  static unsigned long long g_identity_id;
  static unsigned int g_num_loci;
  static unsigned int g_ibd_length;
  static unsigned int g_barcode_length;
  static std::vector<double> g_plaf;
  static std::vector<double> g_prob_crossover;
  enum g_barcode_type_enum{
    ORDINARY, //0
    IBD, //1
    NUMBER_OF_BARCODE_TYPE_OPTIONS = 2
  };
  static g_barcode_type_enum g_barcode_type;
  
  // Default Constructor
  Parameters();
  
};


#endif