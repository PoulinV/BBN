#ifndef __BBN_EM_CASCADE_H__
#define __BBN_EM_CASCADE_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include "bbn/common.h"
#include "bbn/structures.h"


void integration_distribution_over_kernel(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Output_Options * pt_Output_Options,
        double E_i,
        double z,
        Structure_Spectrum * pt_Cascade_Spectrum,
        Structure_Spectrum * pt_Electron_Spectrum,
        const int step,
        const double Rate_electrons_E_e,
        const double Rate_photons_E_g,
        double &resultat_electrons,
        double &resultat_photons);
double compute_photons_rate(double E,
                           double z,
                           Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                           Structure_Output_Options * pt_Output_Options);
double compute_photons_kernel(double E_g,
                             double E_g_prime,
                             double z,
                             double Electron_Spectrum,
                             double Gamma_Spectrum,
                             Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                             Structure_Output_Options * pt_Output_Options);
double compute_electrons_kernel(double E_e,
                               double E_e_prime,
                               double z,
                               double Electron_Spectrum,
                               double Gamma_Spectrum,
                               Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                               Structure_Output_Options * pt_Output_Options);
double compute_electrons_rate(double E,
                             double z,
                             Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                             Structure_Output_Options * pt_Output_Options);
void Cascade_Spectrum_Reading_From_File(double z,
                                        Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                                        Structure_Spectrum * pt_Spectrum,
                                        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
void Cascade_Spectrum_Calculation(double z,
                                  Structure_Output_Options * pt_Output_Options,
                                  Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                                  Structure_Spectrum * pt_Cascade_Spectrum,
                                  Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
void Triangular_Spectrum(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                         Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                         Structure_Spectrum * pt_Cascade_Spectrum,
                         Structure_Spectrum * pt_Electron_Spectrum,
                         Structure_Output_Options * pt_Output_Options);








/*FUNCTIONS TO MODIFY
double Function_Integrand_Spectre_Compton(double E_e, double E_gamma,  double E_gamma_bar);

void Spectre_gamma_compton(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                           Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                           Structure_Spectrum * pt_Electron_Spectrum,
                           Structure_Spectrum * pt_Gamma_Spectrum);
void Spectrum_gamma_scattered(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                              Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                              Structure_Spectrum * pt_Input_Gamma_Spectrum,
                              Structure_Spectrum * pt_Output_Gamma_Spectrum);
void Spectre_electron_compton(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                              Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                              Structure_Spectrum * pt_Electron_Spectrum,
                              Structure_Spectrum * pt_Gamma_Spectrum);
void Spectrum_electron_scattered(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                                 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                                 Structure_Spectrum * pt_Input_Electron_Spectrum,
                                 Structure_Spectrum * pt_Output_Electron_Spectrum);
*/
#endif // __BBN_EM_CASCADE_H__
