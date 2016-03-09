#include "bbn/EM_cascade.h"
#include "bbn/injected_spectrum.h"
#include "bbn/structures.h"
#include "bbn/BBN_constraints.h"
#include "bbn/tools.h"
#include "bbn/test_functions.h"
#include "bbn/photons.h"
#include "bbn/electrons.h"
using namespace std;

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

/*
*This module has for purposes to take the injected gamma and/or e^+e^- spectra from your model of interest.
*
*Some typical (most common) spectrum are already implemented : the "Universal Spectra" for high-energy electrons or photons injection,
*and the "Dirac spectrum" for below pair-production threshold injection. See e.g.  arXiv:1503.04852 for more details.
*
*/

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/



double universal_spectrum(double  E, double  z, double E_0, Structure_Output_Options * pt_Output_Options)
{

    double  E_x = E_x_0/(1+z), E_c = E_c_0/(1+z);
    double E_phph = m_e*m_e/(T_0*(1+z));
    double  K_0 = E_0/(pow(E_x,2)*(2+log(E_c/E_x)));
    double  f;
    double Rate_photons_E_g = rate_NPC(E,z)+rate_compton(E,z);
    Rate_photons_E_g += rate_gg_scattering(E,z);
    //  Rate_photons_E_g = 1;
    if(E < E_x) {
        f = K_0*pow(E_x/E,1.5)/(Rate_photons_E_g);
    } else if(E > E_x && E < E_c) {
        f =  K_0*pow(E_x/E,2)/(Rate_photons_E_g);
    } else {
        f = 0;
    }
    if(pt_Output_Options->Input_verbose > 2)cout << " f = " << f << " E = "<< E <<" E_c = " << E_c << " z = " << z << " E_0 " << E_0 << endl;
    return f;
}

double no_electrons_injected(double x, double z, double E_0, Structure_Output_Options * pt_Output_Options)
{
    return 0;
}

double no_photons_injected(double x, double z, double E_0, Structure_Output_Options * pt_Output_Options)
{
    return 0;
}
double Dirac_Spectrum_After_One_Iteration(double  x, double  z, double E_0, Structure_Output_Options * pt_Output_Options)
{

    double E_c = E_c_0/(1+z), E_phph = m_e*m_e/(T_0*(1+z));

    double T = T_0*(1+z);
    double int_BB = 8./63.*pow(pi,4)*pow(T_0*(1+z),6);
    double f;
    double spectre_gamma_gamma,spectre_compton;
    if(x>E_0) {
        f=0;
    }
    else{
      if(x < E_phph) {
          spectre_gamma_gamma = dsigma_phph(E_0,z,x,pt_Output_Options);

      } else {
          spectre_gamma_gamma = 0;
      }
      spectre_compton = dsigma_compton(E_0,z,x,pt_Output_Options);
      f = spectre_gamma_gamma+spectre_compton;
      if(x==E_0)f+=1;
    }


    if(pt_Output_Options->Input_verbose > 2)cout << " f = " << f << " E = "<< x <<" E_phph = " << E_phph << " z = " << z << " E_0 " << E_0 << endl;
    return f;

}
// void Electron_dirac_spectrum_after_one_iteration(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
//         Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
//         Structure_Spectrum * pt_Electron_Spectrum)
// {
//     double z = pt_Electron_Spectrum->redshift;
//     double E_0 = pt_Particle_Physics_Model->E_0;
//     int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
//     double E_gamma_bb = 2.701*T_0*(1+z);
//     double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
//     double resultat = 0, Gamma_electron_Dirac = 0, E_e;
//     double dE = (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table)/ (double) (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1);
//     // cout << " E_gamma_bb = " << E_gamma_bb << " intbb = " << int_bb << endl;
//     Gamma_electron_Dirac = Rate_Inverse_Compton(E_0,z,pt_Spectrum_and_Precision_Parameters);
//     for(int i = 0 ; i < pt_Spectrum_and_Precision_Parameters->Energy_Table_Size ; i++) {
//         E_e = pt_Spectrum_and_Precision_Parameters->E_min_table + i*dE;
//         pt_Electron_Spectrum->Energy[i] = E_e;
//         resultat = 2*pi*r_e*r_e*m_e*m_e*int_bb/(E_gamma_bb*E_0*E_0*Gamma_electron_Dirac)*Function_Integrand_Spectre_Compton(E_0+E_gamma_bb-E_e,E_0,E_gamma_bb);
//         pt_Electron_Spectrum->Spectrum[i] = resultat/Rate_Inverse_Compton(E_e,z,pt_Spectrum_and_Precision_Parameters);
//         // cout << " E = " << E_e << " resultat = " << pt_Electron_Spectrum->Spectrum[i] << endl;
//     }
// }




/******************************************************************************************************************************************/
/*************************** A few functions that you need to slightly modify if you add an injection spectrum ****************************/
/***************************************** Those are the only functions you will have to modify ! *****************************************/
/******************************************************************************************************************************************/

void attribute_map_spectrum(map_spectrum &map_spectrum)
{
    map_spectrum["no_photons_injected"]=no_photons_injected;
    map_spectrum["no_electrons_injected"]=no_electrons_injected;
    map_spectrum["Dirac_Spectrum_After_One_Iteration"]=Dirac_Spectrum_After_One_Iteration;
    /*You need to add here your own spectrum of the type "Spectrum" as it is defined at the beginning of "structures.h", it means your_function(double E, double z, double E_0).*/
}
void check_name_spectrum(const string &value, string &error_value)
{
    if(value == "Dirac" || value == "none" || value == "universal" || value == "from_file" || value == "Dirac_Spectrum_After_One_Iteration") {
        error_value = "no";
    }
    /* If you want to add a spectrum, you need to add also a "error checking" line. Here is a template you can use.
    else if{value == "name_of_your_function"}{
    	error_value = "no";
    }
    */
    else {
        error_value = "yes";
    }
}
