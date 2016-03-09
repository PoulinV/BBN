#include "bbn/EM_cascade.h"
#include "bbn/injected_spectrum.h"
#include "bbn/structures.h"
#include "bbn/BBN_constraints.h"
#include "bbn/tools.h"
#include "bbn/test_functions.h"
#include "bbn/photons.h"
#include "bbn/electrons.h"


using namespace std;



 /**
  * @ detailed The EM_cascade module :
  * It contains all necessary functions for computing the cascade spectrum.
  */







void Cascade_Spectrum_Reading_From_File(double z,
                                        Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                                        Structure_Spectrum * pt_Spectrum,
                                        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters)
{
    ostringstream os;
    string name;
    double tmp_Energy, tmp_Spectrum;
    pt_Spectrum->redshift = z;

    if(pt_Spectrum->species == "photon") {
        if(pt_Spectrum_and_Precision_Parameters->calculation_mode=="iterative") {
            os << "output/Cascade_Spectrum_Folder/Spectrum_"<<pt_Spectrum->spectrum_name<<"_m" << pt_Particle_Physics_Model->M_x<<"_z"<< z <<"_" << pt_Spectrum_and_Precision_Parameters->number_iterations_photon <<"iterations.dat";
        } else if(pt_Spectrum_and_Precision_Parameters->calculation_mode=="triangular") {
            os << "output/Cascade_Spectrum_Folder/Spectrum_"<<pt_Spectrum->spectrum_name<<"_m" << pt_Particle_Physics_Model->M_x<<"_z"<< z <<"_" << "triangular.dat";
        }
    } else if(pt_Spectrum->species == "electron") {
        // if(pt_Spectrum_and_Precision_Parameters->calculation_mode=="iterative")  os << "output/Cascade_Spectrum_Folder/Spectrum_"<<pt_Spectrum->spectrum_name<<"m" << pt_Particle_Physics_Model->M_x<<"_z"<< z <<"_" << pt_Spectrum_and_Precision_Parameters->number_iterations_electron <<"iterations.dat";
        // else if(pt_Spectrum_and_Precision_Parameters->calculation_mode=="triangular")  os << "output/Cascade_Spectrum_Folder/Spectrum_"<<pt_Spectrum->spectrum_name<<"m" << pt_Particle_Physics_Model->M_x<<"_z"<< z <<"_" << "triangular.dat";
        if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_file_name=="automatic") {
            os << "output/Cascade_Spectrum_Folder/Spectrum_"<<pt_Spectrum->spectrum_name<<"_m" << pt_Particle_Physics_Model->M_x<<"_z"<< z <<"_" << "triangular.dat";
        } else {
            os << pt_Spectrum_and_Precision_Parameters->electron_spectrum_file_name;
        }
        // cout  << pt_Spectrum_and_Precision_Parameters->electron_spectrum_file_name;
    }
    name = os.str();
    ifstream file(name);

    if(file) {
        cout << "Importing file " << name << " in structure " << pt_Spectrum->spectrum_name << "." <<endl;
    } else {
        cout << "I couldn't recognize cascade spectrum file. Please check that it is present in the folder Cascade_Specrum_Folder with proper name : 'Spectrum_mXXX_zXXX_XXXiterations.dat' corresponding to the value of m, z and iterations you are using."<<endl;
        cout << "File I could read is " << name << endl;
        return;
    }
    while(file) {
        string line;
        getline(file, line);
        // stringstream is(ligne);
        // if(line[0] == '#' or line[0] == '\0') continue;
        file >> tmp_Energy >> tmp_Spectrum ;
        // cout << z << "  " << tmp_Energy << "  " << tmp_Spectrum<<endl;
        pt_Spectrum->Energy.push_back(tmp_Energy);
        pt_Spectrum->Spectrum.push_back(tmp_Spectrum);
        // pt_Spectrum->Spectrum.push_back(tmp_Spectrum*(rate_NPC(tmp_Energy,z)+rate_compton(tmp_Energy,z)+rate_gg_scattering(tmp_Energy,z)));

    }

    file.close();


}


double compute_electrons_kernel(double E_e,
                                double E_e_prime,
                                double z,
                                double Electron_Spectrum,
                                double Gamma_Spectrum,
                                Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                                Structure_Output_Options * pt_Output_Options){

  double result_electrons = 0.;
  double E_c = E_c_0/(1+z);
  double ICS_e = 0, NPC = 0, COM = 0, DP = 0;
  double gamma_e, gamma_prime;
  if(pt_Spectrum_and_Precision_Parameters->compton_scattering == "yes" && Gamma_Spectrum != 0) {
      COM = Gamma_Spectrum * dsigma_compton(E_e,z,(E_e+m_e-E_e_prime),pt_Output_Options);
  } else {
      COM = 0;
  }

  if(pt_Spectrum_and_Precision_Parameters->pair_creation_in_nuclei == "yes" && Gamma_Spectrum != 0) {
      NPC = Gamma_Spectrum * dsigma_NPC(E_e,z,E_e_prime,pt_Output_Options);
  } else {
      NPC = 0;
  }
  if(pt_Spectrum_and_Precision_Parameters->inverse_compton_scattering == "yes" && Electron_Spectrum != 0) {
      gamma_e = E_e/m_e;
      gamma_prime = E_e_prime/m_e;
      ICS_e = Electron_Spectrum * dsigma_inverse_compton_electron_spectrum(z, gamma_e,  gamma_prime,  pt_Spectrum_and_Precision_Parameters, pt_Output_Options);
  } else {
      ICS_e = 0;
  }

  if(pt_Spectrum_and_Precision_Parameters->double_photon_pair_creation=="yes" && E_e >= E_c && Gamma_Spectrum != 0) {
      DP = Gamma_Spectrum * dsigma_pair_creation_v2(z,E_e_prime,E_e,pt_Spectrum_and_Precision_Parameters,pt_Output_Options);
      if(DP<0) {
          DP = 0.;
      }
  } else {
      DP = 0.;
  }
  result_electrons += COM;
  result_electrons += NPC;
  result_electrons += ICS_e;
  result_electrons += DP;
  // cout << " (compute_electrons_kernel: ) " << result_electrons << endl;
  return result_electrons;
}

double compute_photons_kernel(double E_g,
                              double E_g_prime,
                              double z,
                              double Electron_Spectrum,
                              double Gamma_Spectrum,
                              Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                              Structure_Output_Options * pt_Output_Options){

  double result_photons = 0.;
  double PP= 0, CS= 0, ICS_g = 0;
  double E_c = E_c_0/(1+z);
  double E_phph = m_e*m_e/(T_0*(1+z));
  double gamma_e;

  if(	pt_Spectrum_and_Precision_Parameters->photon_photon_diffusion == "yes" &&  E_g < E_phph && Gamma_Spectrum != 0)	{
      PP = Gamma_Spectrum * (dsigma_phph(E_g,z,E_g_prime,pt_Output_Options));
  }
  else {
      PP = 0;
  }
  if(pt_Spectrum_and_Precision_Parameters->compton_scattering == "yes" && Gamma_Spectrum != 0){
      CS = Gamma_Spectrum * (dsigma_compton(E_g,z,E_g_prime,pt_Output_Options));
  } else {
      CS = 0;
  }
  if(pt_Spectrum_and_Precision_Parameters->inverse_compton_scattering == "yes" && Electron_Spectrum != 0) {
      gamma_e = E_g/m_e;

      // ICS_g = Electron_Spectrum*gamma_inverse_compton_analytical_v2(gamma_e,E_g_prime,z,pt_Output_Options)*m_e/(E_g_prime);
      ICS_g = Electron_Spectrum*gamma_inverse_compton_analytical(gamma_e,E_g_prime,z,3,pt_Spectrum_and_Precision_Parameters,pt_Output_Options);
      // ICS_g = Electron_Spectrum * dsigma_inverse_compton_electron_spectrum(z, gamma_e, gamma_e - E_g_prime/m_e,  pt_Spectrum_and_Precision_Parameters, pt_Output_Options);


  } else {
      ICS_g = 0;
  }
  result_photons += CS;
  result_photons += PP;
  result_photons += ICS_g;

  // cout << " (compute_photons_kernel: )  " << result_photons << endl;
  return result_photons;
}
double compute_photons_rate(double E,
                           double z,
                           Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                           Structure_Output_Options * pt_Output_Options){

double rate_PP= 0, rate_COM= 0, rate_DP= 0, rate_NP = 0, Rate_photons = 0;
double E_c = E_c_0/(1+z);
double E_phph = m_e*m_e/(T_0*(1+z));


                   if(pt_Spectrum_and_Precision_Parameters->pair_creation_in_nuclei == "yes") {
                       rate_NP= rate_NPC(E,z);
                   } else {
                       rate_NP = 0;
                   }
                   if(pt_Spectrum_and_Precision_Parameters->compton_scattering == "yes") {
                       rate_COM = rate_compton(E,z);
                   } else {
                       rate_COM = 0;
                   }
                   if(pt_Spectrum_and_Precision_Parameters->photon_photon_diffusion == "yes" && E < E_phph ) {
                       rate_PP=rate_gg_scattering(E,z);
                   }
                   else {
                       rate_PP = 0.;
                   }
                   if(pt_Spectrum_and_Precision_Parameters->double_photon_pair_creation=="yes" && E >= E_c) {
                       rate_DP=rate_pair_creation_v2(E,z,pt_Spectrum_and_Precision_Parameters);
                   } else {
                       rate_DP = 0.;
                   }

                   Rate_photons += rate_PP;
                   Rate_photons += rate_NP;
                   Rate_photons += rate_COM;
                   Rate_photons += rate_DP;

                   if(pt_Output_Options->EM_cascade_verbose > 2) {
                       #pragma omp critical(print)
                       {
                           cout << "(compute_photons_rate: ) at E = " << E << " rate_NP = " << rate_NP << " rate_COM = " << rate_COM << " rate_PP = " << rate_PP << " rate_DP = " << rate_DP << " tot = " << Rate_photons << endl;
                       }
                   }

                return Rate_photons;

}
double compute_electrons_rate(double E,
                              double z,
                              Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                              Structure_Output_Options * pt_Output_Options){


 // double Rate_electrons = integrator_simpson_rate_inverse_compton_v2(z,E/m_e,pt_Spectrum_and_Precision_Parameters,pt_Output_Options);
 double Rate_electrons = integrator_simpson_rate_inverse_compton(z,E,pt_Spectrum_and_Precision_Parameters,pt_Output_Options);

 if(pt_Output_Options->EM_cascade_verbose > 1) {
     #pragma omp critical(print)
     {
         cout << "(Rate electrons : ) at E = " << E << " tot = " << Rate_electrons << endl;
     }
 }

 return Rate_electrons;

}
void integration_distribution_over_kernel(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Output_Options * pt_Output_Options,
        double E_i,
        double z,
        Structure_Spectrum * pt_Cascade_Spectrum,
        Structure_Spectrum * pt_Electron_Spectrum,
        const int step,
        const int initial_step,
        const double Rate_electrons_E_e,
        const double Rate_photons_E_g,
        double &resultat_electrons,
        double &resultat_photons)
{

    double E_j, E_j_plus_1, E_j_minus_1, dE_j, dlogE_j;
    double E_switch = 5*E_c_0/(1+z);
    int switch_step = pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-pt_Spectrum_and_Precision_Parameters->Energy_Table_Size/10;
    double Rate_photons_E_j, Rate_electrons_E_j;
    int max_step = pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1;
    int num_steps = max_step-initial_step;
    double weight;
    // weight[step]=0.5;
    if(num_steps==1){
      // weight=0.5;// Trapezoidal rules with f(x0) = 0;
      weight=1.;//Open formula
    }
    else if(num_steps==2){
    /*Simpson's 1/3 rules with f(x0) = 0*/
    if(step == initial_step+1)weight=4./3.;
    else if(step == initial_step+2)weight=1./3.;
    /*Open formula*/
    // if(step == initial_step+1)weight=3./2.;
    // else if(step == initial_step+2)weight=-1./2.;
    }
    else if(num_steps==3){
      /*Simpson's 3/8 rules with f(x0) = 0*/
      if(step == initial_step+1 || step == initial_step+2)weight=9./8.;
      else if(step == initial_step+3)weight=3./8.;
      //Open formula*/
      // if(step == initial_step+1)weight = 23./12.;
      // else if(step == initial_step+2)weight=-16./12.;
      // else if(step == initial_step+3)weight=5./12.;
    }
    else if(num_steps==4){
      //Bode's rules with f(x0) = 0;
      if(step == initial_step+1 || step == initial_step+3)weight=64./45.;
      else if(step == initial_step+2 || step == initial_step+4)weight=24./45.;
      /*Open Formula*/
      // if(step == initial_step+1)weight = 55./24.;
      // else if(step == initial_step+2)weight=-59./24.;
      // else if(step == initial_step+3)weight=37./24.;
      // else if(step == initial_step+4)weight = -9./24.;
    }
    else{
      /*Taken from numerical recipe*/
      // if(step == initial_step+1 || step == max_step) weight=3./8.;
      // else if(step == initial_step+2 || step == max_step-1) weight=7./6.;
      // else if(step == initial_step+3 || step == max_step-2) weight=23./24.;
      // else weight=1.;
      /*Open Formula*/
      if(step==initial_step+1)weight=23./12.;
      else if(step==initial_step+2)weight=7./12.;
      else if(step==max_step-1)weight=13./12.;
      else if(step==max_step-2)weight=5./12.;
      else weight = 1.;
    }
    weight = 1. ;
    // if(step == max_step)weight=1;
    // else weight = 1.;
    // cout << "initial_step " << initial_step << " step " << step << " weight " << weight_integral[step] << " E_i " << E_i << endl;
    double E_c = E_c_0/(1+z), E_x = E_x_0/(1+z);
    double E_phph = m_e*m_e/(T_0*(1+z));
    double E_g = E_i, E_e = E_i;
    double PP= 0, CS= 0, ICS_g = 0, ICS_e= 0, NPC= 0, COM= 0, DP= 0, f_e = 0;
    double E_max = pt_Particle_Physics_Model->E_0, E_cmb_min, E_cmb_max, gamma_e;
    double E_s, x_j, gamma_prime;

    E_j = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table,(double) step/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1));
    // E_j_plus_1 = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table,((double) step+1)/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1));
    // E_j_minus_1 = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table,((double) step-1)/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1));
    dlogE_j = (log(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table))/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1);
    // dlogE_j = 2*(log(E_max)-log(E_i))/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1-initial_step);
    // E_j = exp(log(E_i)+(step-initial_step)*dlogE_j/2.);
    // if(step == (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1))dlogE_j;

    if(step<switch_step){
    	E_max = E_switch;
    	E_j = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table,(double) step/(switch_step));
      dlogE_j = (log(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table))/(switch_step);
    }
    else{
    	E_j = E_switch*pow(E_max/E_switch,(double) (step-switch_step)/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1-switch_step));
      dlogE_j = (log(E_max/E_switch))/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1-switch_step);
    }
    // cout << "dlogE_j " << dlogE_j << " E_j " << E_j << endl;

    // if(step==(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1))dE_j = E_j - E_j_minus_1;
    // else dE_j = (E_j_plus_1 - E_j_minus_1)/2.;
    if(Rate_electrons_E_e!=0) {
      resultat_electrons += dlogE_j*weight/Rate_electrons_E_e*E_j*compute_electrons_kernel(E_j,
                                                                            E_i,
                                                                            z,
                                                                            pt_Electron_Spectrum->Spectrum[step],
                                                                            pt_Cascade_Spectrum->Spectrum[step],
                                                                            pt_Spectrum_and_Precision_Parameters,
                                                                            pt_Output_Options);

    } else {
        resultat_electrons += 0;
    }

    if(pt_Output_Options->EM_cascade_verbose > 1) {
        cout <<"(Scattering electrons : ) at E = " << E_e << " E_j = " <<  E_j <<  " COM = " << COM << " NPC = " << NPC << " ICS_e = " << ICS_e << " DP = " << DP << endl;
    }

    if(Rate_photons_E_g!=0) {
        resultat_photons += dlogE_j*weight/Rate_photons_E_g*E_j*compute_photons_kernel(E_j,
                                                                        E_i,
                                                                        z,
                                                                        pt_Electron_Spectrum->Spectrum[step],
                                                                        pt_Cascade_Spectrum->Spectrum[step],
                                                                        pt_Spectrum_and_Precision_Parameters,
                                                                        pt_Output_Options);

    } else {
        resultat_photons += 0;
    }
    if(resultat_electrons<0) resultat_electrons = 0;
    if(resultat_photons<0) resultat_photons = 0;
    if(pt_Output_Options->EM_cascade_verbose > 1) {
        cout <<"(Scattering photons : ) at E = " << E_g << " E_j = " <<  E_j <<  " PP = " << PP << " CS = " << CS << " ICS_g = " << ICS_g << "f_e " << f_e << endl;
    }
}


void Triangular_Spectrum(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                         Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                         Structure_Spectrum * pt_Cascade_Spectrum,
                         Structure_Spectrum * pt_Electron_Spectrum,
                         Structure_Output_Options * pt_Output_Options)
{


    double E_e_minus_1, E_e_plus_1, E_j, E_j_minus_1, E_j_plus_1, dE_j, dE;
    double E_e, E_g, f_e, check;
    double z = pt_Cascade_Spectrum->redshift;

    Structure_Spectrum Tmp_Electron_Spectrum;
    Structure_Spectrum Tmp_Photon_Spectrum;
    Tmp_Electron_Spectrum.spectrum_name = "cascade_electrons_";
    Tmp_Electron_Spectrum.redshift = z;
    Tmp_Electron_Spectrum.Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
    Tmp_Electron_Spectrum.Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);


    Tmp_Photon_Spectrum.spectrum_name = "cascade_photons_";
    Tmp_Photon_Spectrum.redshift = z;
    Tmp_Photon_Spectrum.Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
    Tmp_Photon_Spectrum.Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);



    double E_c = E_c_0/(1+z), E_x = E_x_0/(1+z);
    double E_phph = m_e*m_e/(T_0*(1+z));
    double E_gamma_bb = 2.701*T_0*(1+z), E_e_ICS, E_g_lim;
    double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
    double E_0 = pt_Particle_Physics_Model->E_0;
    int reset = 0;
    int switch_step = pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-pt_Spectrum_and_Precision_Parameters->Energy_Table_Size/10;
    double integrale=E_0;
    double tau_x = pt_Particle_Physics_Model->tau_x;
    double Rate_photons_E_j = 0, Rate_photons_E_g = 0, Rate_electrons_E_j = 0, Rate_electrons_E_e = 0, Rate_electrons_E_e_2;
    double rate_PP= 0, rate_COM= 0, rate_DP= 0, rate_NP = 0;
    double n_x = n_y_0*pow(1+z,3)*exp(-1/(2*H_r*pow(1+z,2)*tau_x))/E_0;
    double E_max = E_0, E_cmb_min, E_cmb_max, gamma_e;
    double E_s, x_j, gamma_prime;
    double E_switch = 5*E_c;
    double PP= 0, CS= 0, ICS_g = 0, ICS_e= 0, NPC= 0, COM= 0, DP= 0;
    double error;
    if(pt_Output_Options->EM_cascade_verbose > 1) {
        #pragma omp critical(print)
        {
            cout << "E_c = " << E_c << "E_max = " << E_max << " E_phph = " << E_phph <<  endl;
        }
    }

    for(int i = (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1); i>=0 ; i--) {
        Rate_photons_E_g = 0;
        E_e = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table,(double) i/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1));
        E_e_plus_1 = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table,((double) i+1)/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-switch_step));
        E_e_minus_1 = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table,((double) i-1)/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-switch_step));
        if(i<switch_step){
        	E_max = E_switch;
        	E_e = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table,(double) i/(switch_step));
        	E_e_minus_1 = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table,((double) i-1)/(switch_step));
        }
        else{
        	E_e = E_switch*pow(E_max/E_switch,(double) (i-switch_step)/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1-switch_step));
        	E_e_minus_1 = E_switch*pow(E_max/E_switch,((double) (i-switch_step-1))/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1-switch_step));
        }
        E_g = E_e;
        // if(i==(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1))dE = E_e - E_e_minus_1;
        // else dE = (E_e_plus_1 - E_e_minus_1)/2.;
        dE = E_e - E_e_minus_1;
        // error = pow(2*dE,3)/12.;
        pt_Electron_Spectrum->Energy[i] = E_e;
        pt_Cascade_Spectrum->Energy[i] = E_g;
        Tmp_Electron_Spectrum.Energy[i] = E_e;
        Tmp_Electron_Spectrum.Spectrum[i] = 0;
        Tmp_Photon_Spectrum.Energy[i] = E_g;
        Tmp_Photon_Spectrum.Spectrum[i] = 0;

        #pragma omp parallel sections // starts a new team
        {

            {
                Rate_electrons_E_e = compute_electrons_rate(E_g,z,pt_Spectrum_and_Precision_Parameters,pt_Output_Options);
            }


            #pragma omp section
            {
                Rate_photons_E_g = compute_photons_rate(E_g,z,pt_Spectrum_and_Precision_Parameters,pt_Output_Options);
            }

        }
        if(i==(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1)) {
            if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice == "Dirac") {
                pt_Electron_Spectrum->Spectrum[i]=1./(dE);
            }
            if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "Dirac") {
                pt_Cascade_Spectrum->Spectrum[i]=1./(dE);
                // pt_Cascade_Spectrum->Spectrum[i]=1./(dE);
                // pt_Electron_Spectrum->Spectrum[i] =  dE * pt_Cascade_Spectrum->Spectrum[i] * 1/4.*pi*r_e*r_e*pow(m_e,4)*dsigma_pair_creation(z,0.99*E_e,E_g,pt_Spectrum_and_Precision_Parameters,pt_Output_Options)/(E_g*E_g*E_g)/Rate_electrons_E_e;
            }
            pt_Electron_Spectrum->Spectrum[i]/=Rate_electrons_E_e;
            pt_Cascade_Spectrum->Spectrum[i]/=Rate_photons_E_g;
        }


            {
                int end = pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;
                #pragma omp parallel for ordered schedule(dynamic)
                for(int j = i+1 ; j < end ; j ++) {
                    double resultat_photons = 0 ,resultat_electrons = 0;
                    integration_distribution_over_kernel(pt_Particle_Physics_Model,
                                                         pt_Spectrum_and_Precision_Parameters,
                                                         pt_Output_Options,
                                                         E_e,
                                                         z,
                                                         pt_Cascade_Spectrum,
                                                         pt_Electron_Spectrum,
                                                         j,
                                                         i,
                                                         Rate_electrons_E_e,
                                                         Rate_photons_E_g,
                                                         resultat_electrons,
                                                         resultat_photons);

                    #pragma omp critical(dataupdate)
                    {
                        pt_Cascade_Spectrum->Spectrum[i] += resultat_photons;
                        Tmp_Photon_Spectrum.Spectrum[i] += resultat_photons*Rate_photons_E_g;
                        pt_Electron_Spectrum->Spectrum[i] += resultat_electrons;
                        Tmp_Electron_Spectrum.Spectrum[i] += resultat_electrons*Rate_electrons_E_e;
                    }
                }
            }


        if(pt_Output_Options->EM_cascade_verbose > 0) {
            #pragma omp critical(print)
            {
                cout << " ******************* at z = " << z << " E = "<< pt_Electron_Spectrum->Energy[i] << " pt_Electron_Spectrum->Spectrum[i] = " << pt_Electron_Spectrum->Spectrum[i] << " pt_Cascade_Spectrum->Spectrum[i] = " << pt_Cascade_Spectrum->Spectrum[i]  << " i = " << i << " ******************* " << endl;
                // cout << " ******************* at z = " << z << " E = "<< pt_Electron_Spectrum->Energy[i] << " " << Dirac_Spectrum_After_One_Iteration(E_g,z,E_0,pt_Spectrum_and_Precision_Parameters) << " ******************* " << endl;
            }
        }



    }
    if(pt_Spectrum_and_Precision_Parameters->check_energy_conservation == "yes") {
        check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Output_Options,&Tmp_Photon_Spectrum,&Tmp_Electron_Spectrum,integrale);

        for(int i = (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1); i>=0 ; i--) {

            pt_Cascade_Spectrum->Spectrum[i]=pt_Cascade_Spectrum->Spectrum[i]*E_0/integrale;
            Tmp_Photon_Spectrum.Spectrum[i]=Tmp_Photon_Spectrum.Spectrum[i]*E_0/integrale;
        }
        check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Output_Options,&Tmp_Photon_Spectrum,&Tmp_Electron_Spectrum,integrale);
    }

    else {
        cout << "*** No energy check energy conservation check requested ***" << endl;
    }
}








void  Cascade_Spectrum_Calculation(double z,
                                   Structure_Output_Options * pt_Output_Options,
                                   Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                                   Structure_Spectrum * pt_Cascade_Spectrum,
                                   Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters)
{
    double dE, E_i, E_e, E_gamma;
    double resultat, integrale;
    double E_c = E_c_0/(1+z);
    Structure_Spectrum Electron_Spectrum;
    Structure_Spectrum Inverse_Compton_Spectrum;
    Structure_Spectrum Diffused_Gamma_Spectrum;
    Structure_Spectrum Diffused_Electron_Spectrum;
    Structure_Spectrum Compton_Electron_Spectrum;
    Structure_Spectrum Tmp_Electron_Spectrum;
    double E_0 = pt_Particle_Physics_Model->E_0;
    double E_max = E_0;

    /*****Initialization****/
    pt_Cascade_Spectrum->Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
    pt_Cascade_Spectrum->Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
    pt_Cascade_Spectrum->species = "photon";
    pt_Cascade_Spectrum->spectrum_name = "total_photon_";
    pt_Cascade_Spectrum->redshift = z;

    Inverse_Compton_Spectrum.Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
    Inverse_Compton_Spectrum.Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
    Inverse_Compton_Spectrum.species = "photon";
    Inverse_Compton_Spectrum.spectrum_name = "ICS_";
    Inverse_Compton_Spectrum.redshift=z;

    Diffused_Gamma_Spectrum.Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
    Diffused_Gamma_Spectrum.Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
    Diffused_Gamma_Spectrum.species = "photon";
    Diffused_Gamma_Spectrum.spectrum_name = "diffused_photon_";
    Diffused_Gamma_Spectrum.redshift=z;

    Electron_Spectrum.Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
    Electron_Spectrum.Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
    Electron_Spectrum.species="electron";
    Electron_Spectrum.spectrum_name = "total_electrons_";
    Electron_Spectrum.redshift=z;

    Diffused_Electron_Spectrum.Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
    Diffused_Electron_Spectrum.Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
    Diffused_Electron_Spectrum.species = "electron";
    Diffused_Electron_Spectrum.spectrum_name = "diffused_electron_";
    Diffused_Electron_Spectrum.redshift=z;

    Compton_Electron_Spectrum.Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
    Compton_Electron_Spectrum.Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
    Compton_Electron_Spectrum.species = "electron";
    Compton_Electron_Spectrum.spectrum_name = "compton_electron_";
    Compton_Electron_Spectrum.redshift=z;

    Tmp_Electron_Spectrum.Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
    Tmp_Electron_Spectrum.Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
    Tmp_Electron_Spectrum.species = "electron";
    Tmp_Electron_Spectrum.spectrum_name = "tmp_electron_";
    Tmp_Electron_Spectrum.redshift=z;




    /*****************Start of the computation****************/

    if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "universal") {
      if(E_0>1.5*E_c) {
          E_max = 1.5*E_c;
      }
        for(int i=0; i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size; i++) {
            E_i = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table,(double) i/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1));
            pt_Cascade_Spectrum->Energy[i]=E_i;
            pt_Cascade_Spectrum->Spectrum[i]=universal_spectrum(E_i,z,E_0,pt_Output_Options);
            Electron_Spectrum.Energy[i]=E_i;
            Electron_Spectrum.Spectrum[i]=0;
        }
        pt_Cascade_Spectrum->spectrum_name="universal_photon_";
        if(pt_Spectrum_and_Precision_Parameters->check_energy_conservation == "yes") {
            check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Output_Options,pt_Cascade_Spectrum,&Electron_Spectrum,integrale);
        } else {
            if(pt_Output_Options->EM_cascade_verbose > 1) {
                cout << "*** No energy check energy conservation check requested ***" << endl;
            }
        }
        print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, pt_Cascade_Spectrum, pt_Particle_Physics_Model);

    } else {
        if(E_c <= E_0 && pt_Spectrum_and_Precision_Parameters->calculation_mode == "simplified") {


            for(int i=0; i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size; i++) {
                E_i = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table,(double) i/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1));
                pt_Cascade_Spectrum->Energy[i]=E_i;
                pt_Cascade_Spectrum->Spectrum[i]=universal_spectrum(E_i,z,E_0,pt_Output_Options);
                Electron_Spectrum.Energy[i]=E_i;
                Electron_Spectrum.Spectrum[i]=0;
            }
            if(pt_Spectrum_and_Precision_Parameters->check_energy_conservation == "yes") {
                check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Output_Options,pt_Cascade_Spectrum,&Electron_Spectrum,integrale);
            } else {
                if(pt_Output_Options->EM_cascade_verbose > 1) {
                    cout << "*** No energy check energy conservation check requested ***" << endl;
                }
            }
            if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing") {
                if(pt_Output_Options->EM_cascade_verbose>1) {
                    cout <<" I will now print the spectrum in files." << endl;
                }
                print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
            }
        }


        else {



            if(pt_Spectrum_and_Precision_Parameters->calculation_mode == "triangular") {

                if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice == "from_file") {
                    Electron_Spectrum.Energy.resize(0);
                    Electron_Spectrum.Spectrum.resize(0);
                    Cascade_Spectrum_Reading_From_File(z,pt_Particle_Physics_Model,&Electron_Spectrum,pt_Spectrum_and_Precision_Parameters);
                } else{
                  for(int i=0; i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size; i++) {
                      E_i = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_max/pt_Spectrum_and_Precision_Parameters->E_min_table,(double) i/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1));
                      Electron_Spectrum.Energy[i]=E_i;
                      Electron_Spectrum.Spectrum[i]=pt_Spectrum_and_Precision_Parameters->Injected_Electron_Spectrum(E_i,z,E_0,pt_Output_Options);
                  }
                }


                if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "from_file") {
                  pt_Cascade_Spectrum->Energy.resize(0);
                  pt_Cascade_Spectrum->Spectrum.resize(0);
                  Cascade_Spectrum_Reading_From_File(z,pt_Particle_Physics_Model,pt_Cascade_Spectrum,pt_Spectrum_and_Precision_Parameters);
                } else{
                    for(int i=0; i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size; i++) {
                        E_i = pt_Spectrum_and_Precision_Parameters->E_min_table*pow(E_0/pt_Spectrum_and_Precision_Parameters->E_min_table,(double) i/(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1));
                        pt_Cascade_Spectrum->Energy[i]=E_i;
                        pt_Cascade_Spectrum->Spectrum[i]=pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum(E_i,z,E_0,pt_Output_Options);
                    }
                }

                Electron_Spectrum.spectrum_name = "total_electrons_";
                pt_Cascade_Spectrum->spectrum_name = "total_photons_";
                Triangular_Spectrum(pt_Particle_Physics_Model,
                                    pt_Spectrum_and_Precision_Parameters,
                                    pt_Cascade_Spectrum,
                                    &Electron_Spectrum,
                                    pt_Output_Options);
                print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, &Electron_Spectrum, pt_Particle_Physics_Model);
                print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
            }
        }
    }



}








/*************************************FONCTIONS A MODIFIER***********************************/

/*
double Function_Integrand_Spectre_Compton(double E_e, double E_gamma,  double E_gamma_bar)
{
    double f;
    double Gamma_e = 4*E_gamma_bar*E_e/(m_e*m_e);
    double q = E_gamma/(Gamma_e*(E_e-E_gamma));
    double n = E_gamma_bar*E_gamma/(m_e*m_e);

    if(q>=0. && q<=1.) {
        f = 2*q*log(q)
            +(1+2*q)*(1-q)
            +2*n*q*(1-q);
    }
// 	if(q>=0. && q<=1.) {
// 			f = 2*q*log(q)
// 					+(1+2*q)*(1-q)
// 			  	+pow(Gamma_e*q,2)/(2*(1-Gamma_e*q))*(1-q);
// }
    else {
        f=0;
    }
    if(f<0) {
        // cout << "(Function_Integrand_Spectre_Compton : ) Eg= "<<E_gamma<<" Ee =" << E_e << " f = " << f << " q=" << q <<" Gamma_e = " << Gamma_e<< " Ecmb="<< E_gamma_bar<< endl;
        f=0;
    }
    // cout << " q = " << q << " Gamma_e = " << Gamma_e << " f = " << f << " E_gamma_bar =" << E_gamma_bar<<endl;
    return f/E_gamma_bar;
}


void  Spectre_electron_pair_creation(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
													 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
													 Structure_Spectrum * pt_Gamma_Spectrum,
													 Structure_Spectrum * pt_Electron_Spectrum){
			double E_e;
			double E1, E2, E3, E4, E5, E6, E7, f1, f2, f3, f4, f5, f6, f7;
			double dE = (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table)/ (double) (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1);
		 	for(int j =0; j<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;j++){
		 		resultat=0;
		 		Gamma_electron=0;
		 		E_e = pt_Spectrum_and_Precision_Parameters->E_min_table+j*dE;
		 		dE_2 = (pt_Particle_Physics_Model->E_0 - (E_e))/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
		 		h = dE_2/6.;

		 		for(int i=0;i<n_step-1;i++){


		 			if(i==0){
		 				E1=E_e;
		 				}
		 			else{
		 				E1=E7;
		 			}

		 			E2=E1 + h;
		 			E3=E1 + 2*h;
		 			E4=E1 + 3*h;
		 			E5=E1 + 4*h;
		 			E6=E1 + 5*h;
		 			E7=E1 + 6*h;
		 			if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E1, f1);
		 			else f1=0;
		 			if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E2, f2);
		 			else f2=0;
		 			if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E3, f3);
		 			else f3=0;
		 			if(E4<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E4, f4);
		 			else f4=0;
		 			if(E5<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E5, f5);
		 			else f5=0;
		 			if(E6<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E6, f6);
		 			else f6=0;
		 			if(E7<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E7, f7);
		 			else f7=0;

		 			f1*=(function_integrand_electrom_pair_creation(E_e,E1,E_gamma_bb))/pow(E1,3);
		 			f2*=(function_integrand_electrom_pair_creation(E_e,E2,E_gamma_bb))/pow(E2,3);
		 			f3*=(function_integrand_electrom_pair_creation(E_e,E3,E_gamma_bb))/pow(E3,3);
		 			f4*=(function_integrand_electrom_pair_creation(E_e,E4,E_gamma_bb))/pow(E4,3);
		 			f5*=(function_integrand_electrom_pair_creation(E_e,E5,E_gamma_bb))/pow(E5,3);
		 			f6*=(function_integrand_electrom_pair_creation(E_e,E6,E_gamma_bb))/pow(E6,3);
		 			f7*=(function_integrand_electrom_pair_creation(E_e,E7,E_gamma_bb))/pow(E7,3);
		 			resultat += dE_2/840. * (41*f1+216*f2+27*f3+272*f4+27*f5+216*f6+41*f7);

		 		}
				pt_Electron_Spectrum->s
			}

}*/


// void  Spectre_electron_compton(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
// 													 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
// 													 Structure_Spectrum * pt_Gamma_Spectrum,
// 													 Structure_Spectrum * pt_Electron_Spectrum){
//
//  	double z = pt_Electron_Spectrum->redshift;
// 	int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
// 	double E_0 = pt_Particle_Physics_Model->E_0;
// 	double E_gamma_bb = 2.701*T_0*(1+z);
// 	double E_x = E_x_0/(1+z), E_c = E_c_0/(1+z);
// 	double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
// 	double resultat = 0, Gamma_electron = 0;
// 	double dE, dE_2, h;
// 	double E1, E2, E3, E4, E5, E6, E7, f1, f2, f3, f4, f5, f6, f7;
// 	double E_gamma, E_e;
// 	double E_gamma_1, E_gamma_2, E_gamma_3;
// 	double resultat_monochromatique_1, resultat_monochromatique_2;
// 	dE = (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table)/ (double) (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1);
// 	for(int j =0; j<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;j++){
// 		resultat=0;
// 		Gamma_electron=0;
// 		E_e = pt_Spectrum_and_Precision_Parameters->E_min_table+j*dE;
// 		dE_2 = (pt_Particle_Physics_Model->E_0 - (E_e))/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
// 		h = dE_2/6.;
//
// 		for(int i=0;i<n_step-1;i++){
//
//
// 			if(i==0){
// 				E1=E_e;
// 				}
// 			else{
// 				E1=E7;
// 			}
//
// 			E2=E1 + h;
// 			E3=E1 + 2*h;
// 			E4=E1 + 3*h;
// 			E5=E1 + 4*h;
// 			E6=E1 + 5*h;
// 			E7=E1 + 6*h;
// 			if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E1, f1);
// 			else f1=0;
// 			if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E2, f2);
// 			else f2=0;
// 			if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E3, f3);
// 			else f3=0;
// 			if(E4<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E4, f4);
// 			else f4=0;
// 			if(E5<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E5, f5);
// 			else f5=0;
// 			if(E6<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E6, f6);
// 			else f6=0;
// 			if(E7<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E7, f7);
// 			else f7=0;
//
// 			f1*=(dsigma_compton(E1,z,(E1+m_e-E_e),pt_Output_Options)+dsigma_NPC(E1+m_e,z,E_e,pt_Output_Options));
// 			f2*=(dsigma_compton(E2,z,(E2+m_e-E_e),pt_Output_Options)+dsigma_NPC(E2+m_e,z,E_e,pt_Output_Options));
// 			f3*=(dsigma_compton(E3,z,(E3+m_e-E_e),pt_Output_Options)+dsigma_NPC(E3+m_e,z,E_e,pt_Output_Options));
// 			f4*=(dsigma_compton(E4,z,(E4+m_e-E_e),pt_Output_Options)+dsigma_NPC(E4+m_e,z,E_e,pt_Output_Options));
// 			f5*=(dsigma_compton(E5,z,(E5+m_e-E_e),pt_Output_Options)+dsigma_NPC(E5+m_e,z,E_e,pt_Output_Options));
// 			f6*=(dsigma_compton(E6,z,(E6+m_e-E_e),pt_Output_Options)+dsigma_NPC(E6+m_e,z,E_e,pt_Output_Options));
// 			f7*=(dsigma_compton(E7,z,(E7+m_e-E_e),pt_Output_Options)+dsigma_NPC(E7+m_e,z,E_e,pt_Output_Options));
// 			resultat += dE_2/840. * (41*f1+216*f2+27*f3+272*f4+27*f5+216*f6+41*f7);
//
// 		}
// 		Gamma_electron=Rate_Inverse_Compton(E_e,z,pt_Spectrum_and_Precision_Parameters);
// 		if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice=="Dirac")resultat += (dsigma_compton(E_0,z,(E_0+m_e-E_e),pt_Output_Options)+dsigma_NPC(E_0+m_e,z,E_0,pt_Output_Options))/(rate_NPC(E_0,z)+rate_compton(E_0,z)+rate_gg_scattering(E_0,z)) ;
// 		pt_Electron_Spectrum->Energy[j]=E_e;
// 		pt_Electron_Spectrum->Spectrum[j]=resultat;
// 		pt_Electron_Spectrum->Spectrum[j]/=(Gamma_electron);
// 		// cout << " E = " << pt_Electron_Spectrum->Energy[j] << " resultat = " << pt_Electron_Spectrum->Spectrum[j] << endl;
//
// 	}
//
// }
// void Spectrum_electron_scattered(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
// 																 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
// 																 Structure_Spectrum * pt_Input_Electron_Spectrum,
// 															 	 Structure_Spectrum * pt_Output_Electron_Spectrum){
//
// 	  	double z = pt_Input_Electron_Spectrum->redshift;
// 		 	int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
// 		 	double E_0 = pt_Particle_Physics_Model->E_0;
// 		 	double E_gamma_bb = 2.701*T_0*(1+z);
// 		 	double E_x = E_x_0/(1+z), E_c = E_c_0/(1+z);
// 		 	double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
// 		 	double resultat = 0, Gamma_electron = 0, Gamma_electron_Dirac=0, Gamma_electron_test=0;
// 		 	double dE, dE_2, h;
// 		 	double E1, E2, E3, E4, E5, E6, E7, F1, F2, F3, F4, F5, F6, F7, g1, g2, g3, g4, g5, g6, g7;
// 		 	double E_gamma, E_e;
// 		 	double E_gamma_1, E_gamma_2, E_gamma_3;
// 		 	double resultat_monochromatique_1, resultat_monochromatique_2;
// 		 	dE = (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table)/ (double) (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1);
//
//
//
// 		 	for(int j =0; j<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;j++){
// 		 		resultat=0;
// 		 		Gamma_electron=0;
// 				E_e = pt_Spectrum_and_Precision_Parameters->E_min_table+j*dE;
// 				n_step = pt_Spectrum_and_Precision_Parameters->n_step;
// 				if(E_e<15){
// 					n_step *=50;
// 				}
//
// 				dE_2 = (pt_Particle_Physics_Model->E_0 - (E_e))/ (double) (n_step-1);
//
// 				h = dE_2/6.;
//
// 				for(int i=0;i<n_step-1;i++){
//
// 					if(i==0){
// 						E1=E_e;
// 						}
// 					else{
// 						E1=E7;
// 					}
//
// 					E2=E1 + h;
// 					E3=E1 + 2*h;
// 					E4=E1 + 3*h;
// 					E5=E1 + 4*h;
// 					E6=E1 + 5*h;
// 					E7=E1 + 6*h;
//
// 					if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Electron_Spectrum->Energy, pt_Input_Electron_Spectrum->Spectrum, pt_Input_Electron_Spectrum->Energy.size(), E1, g1);
//           else g1=0;
//           if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Electron_Spectrum->Energy, pt_Input_Electron_Spectrum->Spectrum, pt_Input_Electron_Spectrum->Energy.size(), E2, g2);
//           else g2=0;
//           if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Electron_Spectrum->Energy, pt_Input_Electron_Spectrum->Spectrum, pt_Input_Electron_Spectrum->Energy.size(), E3, g3);
//           else g3=0;
//           if(E4<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Electron_Spectrum->Energy, pt_Input_Electron_Spectrum->Spectrum, pt_Input_Electron_Spectrum->Energy.size(), E4, g4);
//           else g4=0;
//           if(E5<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Electron_Spectrum->Energy, pt_Input_Electron_Spectrum->Spectrum, pt_Input_Electron_Spectrum->Energy.size(), E5, g5);
//           else g5=0;
//           if(E6<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Electron_Spectrum->Energy, pt_Input_Electron_Spectrum->Spectrum, pt_Input_Electron_Spectrum->Energy.size(), E6, g6);
//           else g6=0;
//           if(E7<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Electron_Spectrum->Energy, pt_Input_Electron_Spectrum->Spectrum, pt_Input_Electron_Spectrum->Energy.size(), E7, g7);
//           else g7=0;
//
//           F1=g1*Function_Integrand_Spectre_Compton(E1+E_gamma_bb-E_e,E1, E_gamma_bb)/((E1)*(E1));
//           F2=g2*Function_Integrand_Spectre_Compton(E2+E_gamma_bb-E_e,E2, E_gamma_bb)/((E2)*(E2));
//           F3=g3*Function_Integrand_Spectre_Compton(E3+E_gamma_bb-E_e,E3, E_gamma_bb)/((E3)*(E3));
//           F4=g4*Function_Integrand_Spectre_Compton(E4+E_gamma_bb-E_e,E4, E_gamma_bb)/((E4)*(E4));
//           F5=g5*Function_Integrand_Spectre_Compton(E5+E_gamma_bb-E_e,E5, E_gamma_bb)/((E5)*(E5));
//           F6=g6*Function_Integrand_Spectre_Compton(E6+E_gamma_bb-E_e,E6, E_gamma_bb)/((E6)*(E6));
//           F7=g7*Function_Integrand_Spectre_Compton(E7+E_gamma_bb-E_e,E7, E_gamma_bb)/((E7)*(E7));
//
//           resultat+= dE_2/840. * (41*F1+216*F2+27*F3+272*F4+27*F5+216*F6+41*F7);
//
// 					}
//
// 		 		pt_Output_Electron_Spectrum->Energy[j]=E_e;
// 				// cout << " E = " << pt_Output_Electron_Spectrum->Energy[j] << " resultat before constant = " << resultat << endl;
//
// 		 		Gamma_electron=Rate_Inverse_Compton(E_e,z,pt_Spectrum_and_Precision_Parameters);
// 				pt_Output_Electron_Spectrum->Spectrum[j]=resultat*2*pi*r_e*r_e*m_e*m_e*int_bb/(E_gamma_bb);
// 				// cout << " E = " << pt_Output_Electron_Spectrum->Energy[j] << " resultat spectre diffuse avant = " << pt_Output_Electron_Spectrum->Spectrum[j] << " Gamma_electron = " << Gamma_electron << " Gamma_test = " << Gamma_electron_test<<  endl;
//
// 				pt_Output_Electron_Spectrum->Spectrum[j]/=(Gamma_electron);
// 				// cout << " E = " << pt_Output_Electron_Spectrum->Energy[j] << " resultat spectre diffuse = " << pt_Output_Electron_Spectrum->Spectrum[j] << endl;
// 		}
//
// }
//
//
//
// void Spectre_gamma_compton(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
// 													 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
// 													 Structure_Spectrum * pt_Electron_Spectrum,
// 													 Structure_Spectrum * pt_Gamma_Spectrum){
//
// 	double z = pt_Gamma_Spectrum->redshift;
//
// 	int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
// 	double E_gamma_bb = 2.701*T_0*(1+z);
// 	double E_x = E_x_0/(1+z), E_c = E_c_0/(1+z);
// 	double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
// 	double resultat = 0;
// 	double dE, dE_2, h;
// 	double E1, E2, E3, E4, E5, E6, E7, f1, f2, f3, f4, f5, f6, f7, F1, F2, F3, F4, F5, F6, F7;
// 	double E_gamma, E_e;
// 	double q_min=0.0001, q1, q2, q3, dq;
// 	double Gamma_electron=0;
// 	double E_0=pt_Particle_Physics_Model->E_0;
// 	dq = (1-q_min)/ (double) (n_step-1);
// 	dE = (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table)/ (double) (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1);
//
// 		for(int j =0; j<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;j++){
// 			resultat=0;
// 			E_gamma = pt_Spectrum_and_Precision_Parameters->E_min_table + j*dE;
// 			dE_2 = (pt_Particle_Physics_Model->E_0 - (E_gamma+m_e))/ (double) (n_step-1);
// 			h = dE_2/6.;
// 			for(int i=0;i<n_step-1;i++){
// 				if(i==0){
// 					E1=E_gamma+m_e;
// 				}
// 				else{
// 					E1=E7;
// 				}
//
//
// 				E2=E1 + h;
// 				E3=E1 + 2*h;
// 				E4=E1 + 3*h;
// 				E5=E1 + 4*h;
// 				E6=E1 + 5*h;
// 				E7=E1 + 6*h;
// 				if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E1, f1);
// 				else f1=0;
// 				if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E2, f2);
// 				else f2=0;
// 				if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E3, f3);
// 				else f3=0;
// 				if(E4<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E4, f4);
// 				else f4=0;
// 				if(E5<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E5, f5);
// 				else f5=0;
// 				if(E6<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E6, f6);
// 				else f6=0;
// 				if(E7<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E7, f7);
// 				else f7=0;
// 				F1 = Function_Integrand_Spectre_Compton(E_gamma,E1, E_gamma_bb);
// 				F2 = Function_Integrand_Spectre_Compton(E_gamma,E2, E_gamma_bb);
// 				F3 = Function_Integrand_Spectre_Compton(E_gamma,E3, E_gamma_bb);
// 				F4 = Function_Integrand_Spectre_Compton(E_gamma,E4, E_gamma_bb);
// 				F5 = Function_Integrand_Spectre_Compton(E_gamma,E5, E_gamma_bb);
// 				F6 = Function_Integrand_Spectre_Compton(E_gamma,E6, E_gamma_bb);
// 				F7 = Function_Integrand_Spectre_Compton(E_gamma,E7, E_gamma_bb);
//
//
// 				f1*=2*F1/(E1*E1);
// 				f2*=2*F2/(E2*E2);
// 				f3*=2*F3/(E3*E3);
// 				f4*=2*F4/(E4*E4);
// 				f5*=2*F5/(E5*E5);
// 				f6*=2*F6/(E6*E6);
// 				f7*=2*F7/(E7*E7);
// 				// if(f1!=0){cout << "E1 = " << E1<< " f1 = " << f1 << endl;
// 				// }
// 				// cout << "E = " << pt_Spectrum_and_Precision_Parameters->E_min_table+j*dE<<"F1 = "<< F1 << " f1 apres = " << f1 << endl;
//
// 				resultat+= dE_2/840. * (41*f1+216*f2+27*f3+272*f4+27*f5+216*f6+41*f7);
//
// 			}
// 			pt_Gamma_Spectrum->Energy[j]=E_gamma ;
// 			pt_Gamma_Spectrum->Spectrum[j]=resultat*2*pi*r_e*r_e*m_e*m_e*int_bb/E_gamma_bb/(rate_NPC(E_gamma,z)+rate_compton(E_gamma,z)+rate_gg_scattering(E_gamma,z));
// 			if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice=="Dirac"){
// 				Gamma_electron = Rate_Inverse_Compton(E_0,z,pt_Spectrum_and_Precision_Parameters);
// 				pt_Gamma_Spectrum->Spectrum[j]+=2*pi*r_e*r_e*m_e*m_e*2*int_bb*Function_Integrand_Spectre_Compton(E_gamma,E_0,E_gamma_bb)/(E_gamma_bb*Gamma_electron*E_0*E_0)/(rate_NPC(E_gamma,z)+rate_compton(E_gamma,z)+rate_gg_scattering(E_gamma,z));
// 			}
// 			// cout << "E = "  << pt_Spectrum_and_Precision_Parameters->E_min_table+j*dE<< "resultat = " << 	pt_Gamma_Spectrum->Spectrum[j] << endl;
//
// 		}
//
//
// }
// void Spectrum_gamma_scattered(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
// 																 Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
// 																 Structure_Spectrum * pt_Input_Gamma_Spectrum,
// 															 	 Structure_Spectrum * pt_Output_Gamma_Spectrum){
//
//
//  double dE, dE_2, h, h2;
//  double E1, E2, E3, E4, E5, E6, E7, f1, f2, f3, f4, f5, f6, f7,E_gamma;
//  double resultat ;
//  double z = pt_Input_Gamma_Spectrum->redshift;
//  dE = (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table)/ (double) (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1);
//
// 	for(int i=0;i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;i++){
// 		E_gamma = pt_Spectrum_and_Precision_Parameters->E_min_table+i*dE;
// 		resultat = 0;
// 					for(int j=0; j<pt_Spectrum_and_Precision_Parameters->n_step-1;j++){
// 						dE_2 = (pt_Particle_Physics_Model->E_0 - (E_gamma))/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
// 						h2 = dE_2/6.;
//
// 						if(j==0){
// 							E1=E_gamma;
// 							}
// 						else{
// 							E1=E7;
// 						}
//
// 						E2=E1 + h2;
// 						E3=E1 + 2*h2;
// 						E4=E1 + 3*h2;
// 						E5=E1 + 4*h2;
// 						E6=E1 + 5*h2;
// 						E7=E1 + 6*h2;
// 						if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Gamma_Spectrum->Energy, pt_Input_Gamma_Spectrum->Spectrum, pt_Input_Gamma_Spectrum->Energy.size(), E1, f1);
// 						else f1=0;
// 						if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Gamma_Spectrum->Energy, pt_Input_Gamma_Spectrum->Spectrum, pt_Input_Gamma_Spectrum->Energy.size(), E2, f2);
// 						else f2=0;
// 						if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Gamma_Spectrum->Energy, pt_Input_Gamma_Spectrum->Spectrum, pt_Input_Gamma_Spectrum->Energy.size(), E3, f3);
// 						else f3=0;
// 						if(E4<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Gamma_Spectrum->Energy, pt_Input_Gamma_Spectrum->Spectrum, pt_Input_Gamma_Spectrum->Energy.size(), E4, f4);
// 						else f4=0;
// 						if(E5<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Gamma_Spectrum->Energy, pt_Input_Gamma_Spectrum->Spectrum, pt_Input_Gamma_Spectrum->Energy.size(), E5, f5);
// 						else f5=0;
// 						if(E6<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Gamma_Spectrum->Energy, pt_Input_Gamma_Spectrum->Spectrum, pt_Input_Gamma_Spectrum->Energy.size(), E6, f6);
// 						else f6=0;
// 						if(E7<pt_Particle_Physics_Model->E_0)linearint(pt_Input_Gamma_Spectrum->Energy, pt_Input_Gamma_Spectrum->Spectrum, pt_Input_Gamma_Spectrum->Energy.size(), E7, f7);
// 						else f7=0;
// 						// f1*=(dsigma_compton(E1,z,E_gamma));
// 						// f2*=(dsigma_compton(E2,z,E_gamma));
// 						// f3*=(dsigma_compton(E3,z,E_gamma));
// 						// f4*=(dsigma_compton(E4,z,E_gamma));
// 						// f5*=(dsigma_compton(E5,z,E_gamma));
// 						// f6*=(dsigma_compton(E6,z,E_gamma));
// 						// f7*=(dsigma_compton(E7,z,E_gamma));
// 						f1*=(dsigma_phph(E1,z,E_gamma)+dsigma_compton(E1,z,E_gamma));
// 						f2*=(dsigma_phph(E2,z,E_gamma)+dsigma_compton(E2,z,E_gamma));
// 						f3*=(dsigma_phph(E3,z,E_gamma)+dsigma_compton(E3,z,E_gamma));
// 						f4*=(dsigma_phph(E4,z,E_gamma)+dsigma_compton(E4,z,E_gamma));
// 						f5*=(dsigma_phph(E5,z,E_gamma)+dsigma_compton(E5,z,E_gamma));
// 						f6*=(dsigma_phph(E6,z,E_gamma)+dsigma_compton(E6,z,E_gamma));
// 						f7*=(dsigma_phph(E7,z,E_gamma)+dsigma_compton(E7,z,E_gamma));
//
// 						resultat += dE_2/840. * (41*f1+216*f2+27*f3+272*f4+27*f5+216*f6+41*f7);
// 						// cout << "Egamma = " << E_gamma << " E7 = " << E7 << " resultat = " << resultat<< " j = " << j << " i = " << i << endl;
// 					}
// 					pt_Output_Gamma_Spectrum->Energy[i]=E_gamma;
// 					pt_Output_Gamma_Spectrum->Spectrum[i]=resultat/(rate_NPC(E_gamma,z)+rate_compton(E_gamma,z)+rate_gg_scattering(E_gamma,z));
// 				}
// }









/**************************************MODE DE CALCUL ITERATIF EN DEVELOPPEMENT********************************/
// 	else if(pt_Spectrum_and_Precision_Parameters->calculation_mode == "iterative"){
// 		/********First step : compute initial ICS spectrum from the electon spectrum injected**********/
//
// 		for(int i=0;i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;i++){
// 				E1=pt_Spectrum_and_Precision_Parameters->E_min_table+i*dE;
// 				Tmp_Electron_Spectrum.Energy[i]=E1;
// 				Tmp_Electron_Spectrum.Spectrum[i]=pt_Spectrum_and_Precision_Parameters->Injected_Electron_Spectrum(E1,z,pt_Particle_Physics_Model->E_0);
// 				Electron_Spectrum.Energy[i]=E1;
// 				Electron_Spectrum.Spectrum[i]=Tmp_Electron_Spectrum.Spectrum[i];
// 		}
//
// 	for(int i=0;i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;i++){
// 			E1=pt_Spectrum_and_Precision_Parameters->E_min_table+i*dE;
// 			pt_Cascade_Spectrum->Energy[i]=E1;
// 			pt_Cascade_Spectrum->Spectrum[i]=0;
// 			// Electron_Spectrum.Spectrum[i]=Tmp_Electron_Spectrum.Spectrum[i];
// 	}
//
// 	// print_spectrum_automatic_names(0, &Electron_Spectrum, pt_Particle_Physics_Model);
//
// 	if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice != "none"){
//
// 				for(int i = 0; i<pt_Spectrum_and_Precision_Parameters->number_iterations_electron;i++){
// 						cout << " iteration electrons : " << i+1 << endl;
// 						Spectrum_electron_scattered(pt_Particle_Physics_Model,
// 																			pt_Spectrum_and_Precision_Parameters,
// 																			&Electron_Spectrum,
// 																			&Diffused_Electron_Spectrum);
// 						for(int j=0;j<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;j++){
// 						E1=pt_Spectrum_and_Precision_Parameters->E_min_table+j*dE;
// 						Electron_Spectrum.Energy[j]=E1;
// 						Electron_Spectrum.Spectrum[j]=Tmp_Electron_Spectrum.Spectrum[j]+Diffused_Electron_Spectrum.Spectrum[j];
// 						cout <<"E = " << Electron_Spectrum.Energy[j] << " Tmp_Electron_Spectrum = " << Tmp_Electron_Spectrum.Spectrum[j] << " Diffused_Electron_Spectrum = " << Diffused_Electron_Spectrum.Spectrum[j]<<endl;
//
// 						}
// 						print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, &Electron_Spectrum, pt_Particle_Physics_Model);
//
// 				}
//
// 				if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
// 						Diffused_Electron_Spectrum.spectrum_name = "electron_diffused_from_dirac_";
// 						print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, &Diffused_Electron_Spectrum, pt_Particle_Physics_Model);
// 				}
// 				Spectre_gamma_compton(pt_Particle_Physics_Model,
// 															pt_Spectrum_and_Precision_Parameters,
// 															&Electron_Spectrum,
// 															&Inverse_Compton_Spectrum);
// 				if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
// 						Inverse_Compton_Spectrum.spectrum_name = "ICS_from_e_injection_";
// 						if(pt_Output_Options->EM_cascade_verbose>1)cout <<" I will now print the spectrum in files." << endl;
// 						print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, &Inverse_Compton_Spectrum, pt_Particle_Physics_Model);
// 				}
//
//
//
// 	}
// 	/*********Second step : compute the gamma spectrum from gg->gg, ge->ge, gN->Nee for a certain number of iterations.*********/
//
// 		for(int i=0;i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;i++){
// 					pt_Cascade_Spectrum->spectrum_name="Initial_photon_spectrum_";
// 					E1=pt_Spectrum_and_Precision_Parameters->E_min_table+i*dE;
// 					pt_Cascade_Spectrum->Energy[i]=E1;
// 					if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice=="none")pt_Cascade_Spectrum->Spectrum[i]=0;
// 					else pt_Cascade_Spectrum->Spectrum[i]=pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum(E1,z,pt_Particle_Physics_Model->E_0)/(rate_NPC(E1,z)+rate_compton(E1,z)+rate_gg_scattering(E1,z));
// 					pt_Cascade_Spectrum->Spectrum[i]+=Inverse_Compton_Spectrum.Spectrum[i];
// 			}
// 			// check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,&Electron_Spectrum,integrale);
//
// 			if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
// 					pt_Cascade_Spectrum->spectrum_name = "Cascade_";
// 					if(pt_Output_Options->EM_cascade_verbose>1)cout <<" I will now print the spectrum in files." << endl;
// 					print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
// 			}
//
//
//
// 		for(int k = 0; k<pt_Spectrum_and_Precision_Parameters->number_iterations_photon;k++){
// 					if(pt_Output_Options->EM_cascade_verbose>1)cout<<"iteration : " << k+1 << endl;
//
//
// 			Spectrum_gamma_scattered(pt_Particle_Physics_Model,
// 															 pt_Spectrum_and_Precision_Parameters,
// 															 pt_Cascade_Spectrum,
// 															 &Diffused_Gamma_Spectrum);
// 			 for(int i=0;i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;i++){
// 					E_gamma = pt_Spectrum_and_Precision_Parameters->E_min_table+i*dE;
// 					pt_Cascade_Spectrum->Spectrum[i] = Diffused_Gamma_Spectrum.Spectrum[i] + Inverse_Compton_Spectrum.Spectrum[i];
// 					if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice!="none")pt_Cascade_Spectrum->Spectrum[i]+=pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum(E_gamma,z,pt_Particle_Physics_Model->E_0)/(rate_NPC(E_gamma,z)+rate_compton(E_gamma,z)+rate_gg_scattering(E_gamma,z));
// 				}
//
// 				// check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,&Electron_Spectrum,integrale);
//
// 				if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
// 					pt_Cascade_Spectrum->spectrum_name = "Cascade_";
// 					print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
// 				}
//
// 	}
// 	// check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,&Electron_Spectrum,integrale);
//
// 	/**********Third Step : Compute the associated electron spectrum and the gamma spectrum from ICS.**********/
// 	if(pt_Spectrum_and_Precision_Parameters->inverse_compton_scattering == "yes"){
// 		// pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice = "from_cascade";
//
// 		Spectre_electron_compton(pt_Particle_Physics_Model,
// 														 pt_Spectrum_and_Precision_Parameters,
// 														 pt_Cascade_Spectrum,
// 														 &Compton_Electron_Spectrum);
//
// 	for(int j=0; j<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size; j++){
// 		Tmp_Electron_Spectrum.Energy[j]=Compton_Electron_Spectrum.Energy[j];
// 		Tmp_Electron_Spectrum.Spectrum[j]=Compton_Electron_Spectrum.Spectrum[j];
// 	}
// 	for(int i = 0 ; i < pt_Spectrum_and_Precision_Parameters->number_iterations_electron; i++){
// 		cout << "iteration : " << i+1 << endl;
// 	 Spectrum_electron_scattered(pt_Particle_Physics_Model,
// 															pt_Spectrum_and_Precision_Parameters,
// 															&Tmp_Electron_Spectrum,
// 															&Diffused_Electron_Spectrum);
// 	for(int j=0; j<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size; j++){
// 		cout <<"E = " << Electron_Spectrum.Energy[j] << " Compton_Electron_Spectrum = " << Compton_Electron_Spectrum.Spectrum[j] << " Diffused_Electron_Spectrum = " << Diffused_Electron_Spectrum.Spectrum[j]<<endl;
// 		Tmp_Electron_Spectrum.Spectrum[j]= Compton_Electron_Spectrum.Spectrum[j]+Diffused_Electron_Spectrum.Spectrum[j];
// 		cout << " Electron_Spectrum = " << Electron_Spectrum.Spectrum[j] << endl;
// 	}
// 	print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, &Tmp_Electron_Spectrum, pt_Particle_Physics_Model);
//
// }
// for(int j=0; j<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size; j++){
// 	Electron_Spectrum.Energy[j]=Tmp_Electron_Spectrum.Energy[j];
// 	Electron_Spectrum.Spectrum[j]+=Tmp_Electron_Spectrum.Spectrum[j];
// }
// print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, &Electron_Spectrum, pt_Particle_Physics_Model);
//
// 	if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice == "Dirac"){
// 		pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice = "from_cascade";
// 	}
// 		Spectre_gamma_compton(pt_Particle_Physics_Model,
// 													pt_Spectrum_and_Precision_Parameters,
// 													&Tmp_Electron_Spectrum,
// 													&Inverse_Compton_Spectrum);
// 	if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice == "from_cascade"){
// 		pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice = "Dirac";
// 	}
// 		for(int l = 0; l<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size ; l++){
// 			cout <<"pt_Cascade_Spectrum->Spectrum[l] = " << pt_Cascade_Spectrum->Spectrum[l] << " Inverse_Compton_Spectrum.Spectrum[l] = " << Inverse_Compton_Spectrum.Spectrum[l] << endl;
// 			pt_Cascade_Spectrum->Spectrum[l]+=Inverse_Compton_Spectrum.Spectrum[l];
// 		}
//
// 		if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
// 			Compton_Electron_Spectrum.spectrum_name = "compton_electron_";
// 			print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, &Compton_Electron_Spectrum, pt_Particle_Physics_Model);
// 		}
// 		if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
// 			pt_Cascade_Spectrum->spectrum_name = "total_";
// 			print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
// 		}
//
// 		// check_energy_conservation(pt_Particle_Physics_Model,pt_Spectrum_and_Precision_Parameters,pt_Cascade_Spectrum,&Electron_Spectrum,integrale);
//
// 		if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
// 			Inverse_Compton_Spectrum.spectrum_name = "ICS_from_cascade_";
// 			print_spectrum(pt_Output_Options,pt_Spectrum_and_Precision_Parameters, &Inverse_Compton_Spectrum, pt_Particle_Physics_Model);
// 		}
// 	}
// 	}
