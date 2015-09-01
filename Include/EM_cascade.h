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
double rate_compton(double  x, double  z);
double rate_NPC(double  x, double  z);
double rate_gg_scattering(double  x, double  z);
void Cascade_Spectrum_Reading_From_File(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
																				Structure_Spectrum * pt_Cascade_Spectrum,
																				double z,
																				int iterations);
void  Cascade_Spectrum_Calculation(double z,
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
											 	 Structure_Spectrum * pt_Electron_Spectrum);
double Rate_Inverse_Compton(double E_e, double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
double Function_Integrand_Spectre_Compton(long double E_gamma, double E_e, long double E_gamma_bar);
double Function_Integrand_Spectre_Compton_version_q(double q, double E_e, double E_gamma_bar);

#endif
