#ifndef __BBN_INJECTED_STRUCTURES_H__
#define __BBN_INJECTED_STRUCTURES_H__

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "bbn/common.h"

typedef double (*Spectrum)(double E, double  z, double E_0);
typedef std::map<std::string, std::string> map_parameters;
typedef std::map<std::string, Spectrum> map_spectrum;


class Structure_Spectrum{

public:

  int iterations;
  double redshift;
  double temperature;
  std::vector<double> Energy;
  std::vector<double> Spectrum;
  std::string spectrum_name;
  std::string species;


};



class Structure_Particle_Physics_Model{

public:

  double  M_x;
  double  E_0;
  double  zeta_x;
  double  tau_x;
  double  T_x;
  double  z_x;

};

class Structure_Spectrum_and_Precision_Parameters{

public:

  int number_iterations_photon;
  int number_iterations_electron;
  int z_step;
  int n_step;
  int eval_max;
  std::vector<int> weight;
  double divisor;
  double Gamma_Table_Size;
  double Energy_Table_Size;
  double E_min_table;
  std::string calculation_mode;
  std::string photon_spectrum_choice;
  std::string photon_spectrum_file_name;
  std::string electron_spectrum_choice;
  std::string electron_spectrum_file_name;
  std::string spectrum_mode;
  std::string inverse_compton_scattering;
  std::string double_photon_pair_creation;
  std::string pair_creation_in_nuclei;
  std::string compton_scattering;
  std::string photon_photon_diffusion;
  std::string check_energy_conservation;
  std::string integration_method;

  Spectrum Injected_Gamma_Spectrum;
  Spectrum Injected_Electron_Spectrum;


};

class Structure_Scan_Parameters_and_Results{

public:

  std::string nuclei;
  int i_min, i_max;
  int k_min, k_max;
  double Y_0, Y_min, Y_max;
  double tau_min;
  double tau_max;
  double zeta_min;
  double zeta_max;
  int tau_step;
  int zeta_step;
  std::vector<double> Results_scan_tau_x;
  std::vector<double> Results_scan_zeta_x;
  std::vector<double> Results_scan_Abundance;

};

class Structure_Output_Options{

  public:

  std::string results_files;
  std::string spectrum_files;
  std::string interaction_rate_files;
  std::string test_files;
  std::string task;
  std::string task_test;
  int EM_cascade_verbose;
  int BBN_constraints_verbose;
  int Input_verbose;
  int Output_verbose;
};

#endif // __BBN_INJECTED_STRUCTURES_H__
