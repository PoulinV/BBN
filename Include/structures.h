#ifndef __structures__
#define __structures__

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "common.h"

typedef double (*Spectrum)(double E, double  z, double E_0);
typedef std::map<std::string, std::string> map_parameters;
typedef std::map<std::string, Spectrum> map_spectrum;


struct Structure_Spectrum{

  int iterations;
  double redshift;
  std::vector<double> Energy;
  std::vector<double> Spectrum;
  std::string spectrum_name;
  std::string species;


};



struct Structure_Particle_Physics_Model{

  double  M_x;
  double  E_0;
  double  zeta_x;
  double  tau_x;
  double  T_x;
  double  z_x;

};

struct Structure_Spectrum_and_Precision_Parameters{

  int number_iterations_photon;
  int number_iterations_electron;
  int z_step;
  int n_step;
  std::string calculation_mode;
  std::string photon_spectrum_choice;
  std::string electron_spectrum_choice;
  std::string spectrum_mode;
  std::string inverse_compton_scattering;
  Spectrum Injected_Gamma_Spectrum;
  Spectrum Injected_Electron_Spectrum;

};

struct Structure_Scan_Parameters{

  std::string nuclei;
  double tau_min;
  double tau_max;
  double zeta_min;
  double zeta_max;
  int tau_step;
  int zeta_step;
};

struct Structure_Output_Options{
  std::string print_result;
  std::string results_files;
  std::string spectrum_files;
};

#endif
