#ifndef __tools__
#define __tools__


#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include "common.h"
/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

/*
*
*This module contains all the numerical tools necesserary to the code.
*Most of the algorithm comes from Numerical Recipes.
*
*/

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/


void Attribution_avec_correction(int &x);
void fill_two_tables_from_file(std::ifstream &file, std::vector<double> &vector_z, std::vector<double> &vector_y);
void linearint(std::vector<double> &xa, std::vector<double> &ya, int n, double x, double &y);
void polint(std::vector<double> &xa, std::vector<double> &ya, int n, double x, double &y, double &dy);
void print_spectrum(std::ostream &file, Structure_Spectrum * pt_Cascade_Spectrum, Structure_Particle_Physics_Model * pt_Particle_Model);
void print_spectrum_from_function(std::ostream &file, double (*func)(double,double,double), Structure_Particle_Physics_Model * pt_Particle_Model);
void print_spectrum_automatic_names(int number_files, Structure_Spectrum * pt_Cascade_Spectrum, Structure_Particle_Physics_Model * pt_Particle_Model);
void fill_structure_particle_physics_model(double M_x, double Zeta_x, double tau_x, Structure_Particle_Physics_Model * pt_Particle_Model);
void fill_structure_spectrum_and_precision_parameters(int number_iterations_photon,
                                                      int number_iterations_electron,
                                                      int z_step,
                                                      int n_step,
                                                      const std::string &calcutation_mode,
                                                      const std::string &photon_spectrum_choice,
                                                      const std::string &electron_spectrum_choice,
                                                      const std::string &spectrum_mode,
                                                      const std::string &inverse_compton_scattering,
                                                      double (*Gamma_Spectrum)(double, double, double),
                                                      double (*Electron_Spectrum)(double,double,double),
                                                      Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
void fill_structure_scan_parameters(const std::string &nuclei, double tau_min, double tau_max, double tau_step, double zeta_min, double zeta_max, double zeta_step, Structure_Scan_Parameters * pt_Scan_Parameters);
void fill_output_options(const std::string &print_result, const std::string &results_files, const std::string &spectrum_files, Structure_Output_Options * pt_Structure_Output_Options);
void check_energy_conservation(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                              Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                              Structure_Spectrum * pt_Gamma_Spectrum,
                              Structure_Spectrum * pt_Electron_Spectrum,
                               double &integrale);

#endif
