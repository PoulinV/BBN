#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
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
void fill_two_tables_from_file(ifstream &file, vector<double> &vector_z, vector<double> &vector_y);
void linearint(vector<double> &xa, vector<double> &ya, int n, double x, double &y);
void polint(vector<double> &xa, vector<double> &ya, int n, double x, double &y, double &dy);
void print_spectrum(ostream &file, struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum, struct Structure_Particle_Physics_Model * pt_Particle_Model);
void print_spectrum_from_function(ostream &file, double (*func)(double,double,double), struct Structure_Particle_Physics_Model * pt_Particle_Model);
void print_spectrum_automatic_names(int number_files, struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum, struct Structure_Particle_Physics_Model * pt_Particle_Model);
void Fill_Structure_Particle_Physics_Model(double M_x, double Zeta_x, double tau_x, struct Structure_Particle_Physics_Model * pt_Particle_Model);
void Fill_Structure_Spectrum_and_Precision_Parameters(int iterations,
                                                      int z_step,
                                                      int n_step,
                                                      string spectrum_choice,
                                                      string spectrum_mode,
                                                      string inverse_compton_scattering,
                                                      double (*Gamma_Spectrum)(double, double, double),
                                                      double (*Electron_Spectrum)(double, double, double),
                                                      struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
void Fill_Structure_Scan_Parameters(string nuclei, double tau_min, double tau_max, double tau_step, double zeta_min, double zeta_max, double zeta_step, struct Structure_Scan_Parameters * pt_Scan_Parameters);
void Fill_Output_Options(string print_result, string results_files, string spectrum_files, struct Structure_Output_Options * pt_Structure_Output_Options);

// double  trapzdb(double  (*func)(double,double,int,double,struct), double  a, double  b, double  d, int n, int i, struct Structure_Injected_Spectrum * pt_Injected_Spectrum);
// double  qsimp(double  (*func)(double,double,int,double,struct), double  a, double  b, double  d, int i, struct Structure_Injected_Spectrum * pt_Injected_Spectrum);
// double  trapzd_2(double  (*func)(double,double,int,double,struct), double  a, double  b, double  d, int n, int i, struct Structure_Injected_Spectrum * pt_Injected_Spectrum);
// double  qsimp_2(double  (*func)(double,double,int,double,struct), double  a, double  b, double  d, int i, struct Structure_Injected_Spectrum * pt_Injected_Spectrum);
// double  trapzd_3(double  (*func)(double,double,int,double,struct), double  a, double  b, double  d, int n, int i, struct Structure_Injected_Spectrum * pt_Injected_Spectrum);
// double  qsimp_4(double  (*func)(double,double,int,double,struct), double  a, double  b, double  d, int i, struct Structure_Injected_Spectrum * pt_Injected_Spectrum);
// double  trapzd_Z(double  (*func)(double,double,int,double,struct), double  a, double  b, double  d, int n, int i, struct Structure_Injected_Spectrum * pt_Injected_Spectrum);
// double  qsimp_Z(double  (*func)(double,double,int,double,struct), double  a, double  b, double  d, int i, struct Structure_Injected_Spectrum * pt_Injected_Spectrum);
// double  trapzd_E(double  (*func)(double,double,int,double,struct), double  a, double  b, double  d, int n, int i, struct Structure_Injected_Spectrum * pt_Injected_Spectrum);
// double  qsimp_E(double  (*func)(double,double,int,double,struct), double  a, double  b, double  d, int i, struct Structure_Injected_Spectrum * pt_Injected_Spectrum);
