#ifndef __injected_spectrum__
#define __injected_spectrum__

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





double universal_spectrum(double  E, double  z, double E_0);
double photon_spectrum_after_one_iteration(double  x, double  z, double E_0);
double electron_dirac_spectrum_after_one_iteration(double  E, double  z, double E_0);
double no_electrons_injected(double x, double z, double E_0);
double no_photons_injected(double x, double z, double E_0);

#endif
