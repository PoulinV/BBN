#include "../Include/Injected_spectrum.h"
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



double Universal_Spectrum(double  E, double  z, double E_0){

	double  E_x = E_x_0/(1+z), E_c = E_c_0/(1+z);
	double  K_0 = E_0/(pow(E_x,2)*(2+log(E_c/E_x)));
	double  f;

	if(E < E_x) f = K_0*pow(E_x/E,1.5);
	else if(E > E_x && E < E_c) f =  K_0*pow(E_x/E,2);
	else {f = 0;}
	return f;
}
void Fill_Structure_Particle_Physics_Model(double M_x, double Zeta_x, double tau_x, struct Structure_Particle_Physics_Model * pt_Particle_Model){

	pt_Particle_Model->M_x = M_x;
	pt_Particle_Model->E_0 = M_x/2.;
	pt_Particle_Model->Zeta_x = Zeta_x;
	pt_Particle_Model->tau_x = tau_x;
	pt_Particle_Model->T_x	= pow(tau_x*(2*H_r),-0.5)*T_0;
	pt_Particle_Model->z_x = 	pow(tau_x*(2*H_r),-0.5) -1;

}
// void Initialize_Spectrum(double (*func)(double,double,struct),struct Structure_Particle_Model * pt_Particle_Model,struct Structure_Gamma_Spectrum * pt_Injected_Spectrum){
//
// 	pt_Cascade_Spectrum->Gamma_Energy.resize(Gamma_Table_Size);
//
//
// }
