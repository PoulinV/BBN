#ifndef __em_cascade__
#define __em_cascade__

#include <iostream>
#include <fstream>
#include <sstream>
#include "common.h"
#include "../Include/structures.h"

double dsigma_compton(double  x, double  z, double g);
double dsigma_phph(double  x, double  z,  double g);
double dsigma_NPC(double E_gamma, double z, double E_e);
double dsigma_pair_creation(double z,
													  double E_e,
													  double E_gamma,
													  Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
double rate_compton(double  x, double  z);
double rate_NPC(double  x, double  z);
double rate_gg_scattering(double  x, double  z);
double rate_pair_creation(double E, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
double rate_pair_creation_v2(double E_gamma, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
double Analytical_form_inverse_compton(double E_e, double E_gamma_bar, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
double integrator_simpson_rate_inverse_compton(double z,
																									 double E_ini,
																									 double E_max,
																									 double E_e,
																									 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
double print_interaction_rate(double z,
															double E_MIN,
															double E_MAX,
															Structure_Output_Options * pt_Output_Options,
													  	Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
void Cascade_Spectrum_Reading_From_File(double z,
																				Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                                        Structure_Spectrum * pt_Spectrum,
                                        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
void Cascade_Spectrum_Calculation(double z,
																	 Structure_Output_Options * pt_Output_Options,
																	 Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
																	 Structure_Spectrum * pt_Cascade_Spectrum,
																	 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
void Spectre_gamma_compton(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
													 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
													 Structure_Spectrum * pt_Electron_Spectrum,
													 Structure_Spectrum * pt_Gamma_Spectrum);
void Spectrum_gamma_scattered(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
													  	Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
														  Structure_Spectrum * pt_Input_Gamma_Spectrum,
													 	  Structure_Spectrum * pt_Output_Gamma_Spectrum);
void  Spectre_electron_compton(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
															 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
															 Structure_Spectrum * pt_Electron_Spectrum,
															 Structure_Spectrum * pt_Gamma_Spectrum);
void Spectrum_electron_scattered(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
																 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
																 Structure_Spectrum * pt_Input_Electron_Spectrum,
															 	 Structure_Spectrum * pt_Output_Electron_Spectrum);
void Triangular_Spectrum(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
												 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
												 Structure_Spectrum * pt_Cascade_Spectrum,
											 	 Structure_Spectrum * pt_Electron_Spectrum,
											 	 Structure_Output_Options * pt_Output_Options);
double Rate_Inverse_Compton(double E_e, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
double Function_Integrand_Spectre_Compton(double E_e, double E_gamma,  double E_gamma_bar);
double Function_Integrand_Spectre_Compton_version_q(double q, double E_e, double E_gamma_bar);
double function_integrand_pair_creation(double E_e, double E_gamma, double E_gamma_bar);
double integrand_rate_pair_creation(double s);
#endif
