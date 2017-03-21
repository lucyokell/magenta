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

#include "stdafx.h"
#include <math.h>
#include <vector>
#include <queue>
#include <bitset>

class Parameters {
  
public:
  
  int g_current_time;
  double g_years;
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
  // diagnostic parameters;
  double g_theta;
  double g_ft;
  // entomological parameters;
  double g_mu0;
  double g_beta_gradient;
  double g_beta_intercept;
  double g_ak;
  double g_Q0;
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
  
  // Default Constructor
  Parameters();
  
};


#endif