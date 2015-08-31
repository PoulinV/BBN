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
*This module has for purposes to take the injected gamma and/or e^+e^- spectra from your the model of interest.
*
*Some typical (most common) spectrum are already implemented : the "Universal Spectra" for high-energy electrons or photons injection,
*and the "Dirac spectrum" for below pair-production threshold injection. See e.g.  arXiv:1503.04852 for more details.
*
*/

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/





double Universal_Spectrum(double  E, double  z, double E_0);
double Dirac_Spectrum_After_One_Iteration(double  x, double  z, double E_0);
void Electron_dirac_spectrum_after_one_iteration(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
							                                   struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
							                                   struct Structure_Spectrum * pt_Electron_Spectrum);
void No_Electrons_Injected(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                           struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                           struct Structure_Spectrum * pt_Electron_Spectrum);
double No_Photons_Injected(double x, double z, double E_0);
