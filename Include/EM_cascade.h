#include "common.h"

double  dsigma_compton(double  x, double  z, double g);
double  dsigma_phph(double  x, double  z,  double g);
double  gamma_compton(double  x, double  z);
double  gamma_NPC(double  x, double  z);

double  gamma_phph(double  x, double  z);

double Dirac_Spectrum_After_One_Iteration(double  x, double  z, struct Structure_Particle_Physics_Model * pt_Particle_Model);
/*
void  Cascade_Spectrum_Calculation_From_File(ifstream &file, struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum, long n_step);
*/
void  Cascade_Spectrum_Calculation_From_Function(double (*func)(double,double,double),double z, struct Structure_Particle_Physics_Model * pt_Particle_Model, struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum, long n_step, int iterations);
