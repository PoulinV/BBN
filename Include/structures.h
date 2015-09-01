#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>

#include "common.h"

typedef double (*Spectrum)(double E, double  z, double E_0);
typedef void (*Electron_Spectrum)(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                                  struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                                  struct Structure_Spectrum * Tmp_Electron_Spectrum);
struct Structure_Spectrum{

  int iterations;
  double redshift;
  vector<double> Energy;
  vector<double> Spectrum;
  string spectrum_name;
  string species;


};



struct Structure_Particle_Physics_Model{

  double  M_x;
  double  E_0;
  double  Zeta_x;
  double  tau_x;
  double  T_x;
  double  z_x;

};

struct Structure_Spectrum_and_Precision_Parameters{

  int number_iterations_photon;
  int number_iterations_electron;
  int z_step;
  int n_step;
  string calculation_mode;
  string photon_spectrum_choice;
  string electron_spectrum_choice;
  string spectrum_mode;
  string inverse_compton_scattering;
  Spectrum Injected_Gamma_Spectrum;
  Spectrum Injected_Electron_Spectrum;

};

struct Structure_Scan_Parameters{

  string nuclei;
  double tau_min;
  double tau_max;
  double zeta_min;
  double zeta_max;
  int tau_step;
  int zeta_step;
};

struct Structure_Output_Options{
  string print_result;
  string results_files;
  string spectrum_files;
};
