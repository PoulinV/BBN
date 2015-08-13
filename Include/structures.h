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

struct Structure_Gamma_Spectrum{

  int iterations;
  double redshift;
  vector<double> Gamma_Energy;
  vector<double> Gamma_Spectrum;



};

struct Structure_Electron_Spectrum{

  int iterations;
  double redshift;
  vector<double> Electron_Energy;
  vector<double> Electron_Spectrum;



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

  int iterations;
  int z_step;
  int n_step;
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
