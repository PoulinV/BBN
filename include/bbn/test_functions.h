#ifndef __BBN_TEST_FUNCTIONS_H__
#define __BBN_TEST_FUNCTIONS_H__

#include <iostream>
#include <fstream>
#include <sstream>

#include "bbn/common.h"
#include "bbn/structures.h"

double integrate_dsigma_phph(double E_MIN, double E_MAX, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,Structure_Output_Options * pt_Output_Options);
double integrate_dsigma_compton(double E_MIN, double E_MAX, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,Structure_Output_Options * pt_Output_Options);
double integrate_dsigma_NPC(double E_MIN,
                                double E_MAX,
                                double z,
                                Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                                Structure_Output_Options * pt_Output_Options);
double integrate_dsigma_pair_creation(double z,
        double E_ini,
        double E_max,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Output_Options * pt_Output_Options);
double integrate_dsigma_inverse_compton_electron_spectrum(double E_MIN,
                                double E_MAX,
                                double z,
                                Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                                Structure_Output_Options * pt_Output_Options);
double Function_Integrand_Spectre_Compton_times_bb_spectrum(double E_e, double E_gamma,  double E_gamma_bar);
double print_func_kawmor(double  z, double  E_0, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters, Structure_Output_Options * pt_Output_Options);
void check_energy_conservation(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                               Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                               Structure_Output_Options * pt_Output_Options,
                               Structure_Spectrum * pt_Gamma_Spectrum,
                               Structure_Spectrum * pt_Electron_Spectrum,
                               double &integrale);
double print_interaction_rate(double z,
                              double E_MIN,
                              double E_MAX,
                              Structure_Output_Options * pt_Output_Options,
                              Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
#endif // __BBN_TEST_FUNCTIONS_H__
