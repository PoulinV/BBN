#include "common.h"

double dsigma_compton(double  x, double  z, double g);
double dsigma_phph(double  x, double  z,  double g);
double gamma_compton(double  x, double  z);
double gamma_NPC(double  x, double  z);
double gamma_phph(double  x, double  z);
void Cascade_Spectrum_Reading_From_File(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
																				struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum,
																				double z,
																				int iterations);
void  Cascade_Spectrum_Calculation(double z,
																	 struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
																	 struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum,
																	 struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters);
void Spectre_gamma_compton(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
													 struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
													 struct Structure_Electron_Spectrum * pt_Electron_Spectrum,
													 struct Structure_Gamma_Spectrum * pt_Gamma_Spectrum);
void  Spectre_electron_compton(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
															 struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
															 struct Structure_Electron_Spectrum * pt_Electron_Spectrum,
															 struct Structure_Gamma_Spectrum * pt_Gamma_Spectrum);
