#include "../Include/tools.h"
#include "../Include/EM_cascade.h"


void print_spectrum(ostream &file, struct Structure_Spectrum * pt_Cascade_Spectrum, struct Structure_Particle_Physics_Model * pt_Particle_Model){


  double dE = (pt_Particle_Model->E_0 - E_min)/Gamma_Table_Size;
  double z = pt_Cascade_Spectrum->redshift;
  double E = E_min;
  int i = 0;
  while(i<Gamma_Table_Size){
    file << E << "   " << pt_Cascade_Spectrum->Spectrum[i]/(rate_NPC(E,z)+rate_compton(E,z)+rate_gg_scattering(E,z)) << endl;
    i++;
    E+=dE;
  }

}
void print_spectrum_from_function(ostream &file, double (*func)(double,double,double),double z, struct Structure_Particle_Physics_Model * pt_Particle_Model){


  double dE = (pt_Particle_Model->E_0 - E_min)/Gamma_Table_Size;
  double E = E_min;
  int i = 0;
  while(i<Gamma_Table_Size){
    file << E << "   " << (*func)(E,z,pt_Particle_Model->E_0)/(rate_NPC(E,z)+rate_compton(E,z)+rate_gg_scattering(E,z)) << endl;
    i++;
    E+=dE;
  }

}
void print_spectrum_automatic_names(int iterations, struct Structure_Spectrum * pt_Cascade_Spectrum, struct Structure_Particle_Physics_Model * pt_Particle_Model){
  ostringstream os;
  string name;


  double dE = (pt_Particle_Model->E_0 - E_min)/pt_Cascade_Spectrum->Energy.size();
  double z = pt_Cascade_Spectrum->redshift;
  double E = E_min;
  int i = 0;

	  mkdir("Cascade_Spectrum_Folder", 01777);
    os << "Cascade_Spectrum_Folder/Spectrum_"<<pt_Cascade_Spectrum->spectrum_name<<"m" << pt_Particle_Model->M_x<<"_z"<< z <<"_" << iterations <<"iterations.dat";
    name = os.str();
    ofstream file(name);
    if(verbose>1)cout << "Printing in file " << name << endl;

    if(pt_Cascade_Spectrum->species=="photon"){
    while(i<Gamma_Table_Size){
      file << pt_Cascade_Spectrum->Energy[i] << "   " << pt_Cascade_Spectrum->Spectrum[i] << endl;
      // file << pt_Cascade_Spectrum->Energy[i] << "   " << pt_Cascade_Spectrum->Spectrum[i]/(rate_NPC(pt_Cascade_Spectrum->Energy[i],z)+rate_compton(pt_Cascade_Spectrum->Energy[i],z)+rate_gg_scattering(pt_Cascade_Spectrum->Energy[i],z)) << endl;
      // if(i<Gamma_Table_Size-1)file << pt_Cascade_Spectrum->Energy[i] << "   " << pt_Cascade_Spectrum->Spectrum[i]/(rate_NPC(pt_Cascade_Spectrum->Energy[i],z)+rate_compton(pt_Cascade_Spectrum->Energy[i],z)+rate_gg_scattering(pt_Cascade_Spectrum->Energy[i],z)) << endl;
      // if(i==Gamma_Table_Size-1)file << pt_Cascade_Spectrum->Energy[i] << "   " << pt_Cascade_Spectrum->Spectrum[i]/(rate_NPC(pt_Cascade_Spectrum->Energy[i],z)+rate_compton(pt_Cascade_Spectrum->Energy[i],z)+rate_gg_scattering(pt_Cascade_Spectrum->Energy[i],z))+1/(rate_NPC(pt_Particle_Model->E_0,z)+rate_compton(pt_Particle_Model->E_0,z)+rate_gg_scattering(pt_Particle_Model->E_0,z)) << endl;
      i++;
    }
  }
    if(pt_Cascade_Spectrum->species=="electron"){
      while(i<Electron_Table_Size){
        file << pt_Cascade_Spectrum->Energy[i] << "   " << pt_Cascade_Spectrum->Spectrum[i] << endl;
        i++;
      }
    }

}

void linearint(vector<double> &xa, vector<double> &ya, int n, double x, double &y){
				int i,m,ns=1;
				float x0,x1,y0,y1;
				double dift,dif;
				if(isnan(xa[0])==0){
          dif=sqrt(pow(x-xa[0],2));
          y = ya[0];
        }
        else{
          cout << "careful, you're initial point isn't defined !" << endl;
          return;
        }
				// cout << "dif = " << dif << endl;
				for (i=1;i<n;i++) {
						if ( (dift=sqrt(pow(x-xa[i],2))) < dif) {
						ns=i;
						dif=dift;
						y0=ya[i];
						x0=xa[i];
						y=ya[i];
						}
						// if(ns!=i)break;
				}
        if(i<n-1){
          x1=xa[ns+1];
          y1=ya[ns+1];
          if(x1!=x0)y=y0+(x-x0)/(x1-x0)*fabs(y1-y0);
          // cout << " y = " << y << " x - x0 = " << x-x0 << endl;
          else y = y0;
          if(isnan(y)==1){cout << "y is a nan ! I automatically put it to 0." << endl;
          y=0;
          }
        }
}

void Fill_Structure_Particle_Physics_Model(double M_x, double Zeta_x, double tau_x, struct Structure_Particle_Physics_Model * pt_Particle_Model){

	pt_Particle_Model->M_x = M_x;
	pt_Particle_Model->E_0 = M_x/2.;
	pt_Particle_Model->Zeta_x = Zeta_x;
	pt_Particle_Model->tau_x = tau_x;
	pt_Particle_Model->T_x	= pow(tau_x*(2*H_r),-0.5)*T_0;
	pt_Particle_Model->z_x = 	pow(tau_x*(2*H_r),-0.5) -1;

}

void Fill_Structure_Scan_Parameters(string nuclei, double tau_min, double tau_max, double tau_step, double zeta_min, double zeta_max, double zeta_step, struct Structure_Scan_Parameters * pt_Scan_Parameters){

  pt_Scan_Parameters->nuclei = nuclei;
	pt_Scan_Parameters->tau_min = tau_min;
	pt_Scan_Parameters->tau_max = tau_max;
	pt_Scan_Parameters->tau_step = tau_step;
  pt_Scan_Parameters->zeta_min = zeta_min;
	pt_Scan_Parameters->zeta_max = zeta_max;
	pt_Scan_Parameters->zeta_step = zeta_step;

}

void Fill_Structure_Spectrum_and_Precision_Parameters(int number_iterations_photon,
                                                      int number_iterations_electron,
                                                      int z_step,
                                                      int n_step,
                                                      string photon_spectrum_choice,
                                                      string electron_spectrum_choice,
                                                      string spectrum_mode,
                                                      string inverse_compton_scattering,
                                                      double (*Gamma_Spectrum)(double, double, double),
                                                      void (*Electron_Spectrum)(struct Structure_Particle_Physics_Model *, struct Structure_Spectrum_and_Precision_Parameters *, struct Structure_Spectrum *),
                                                      struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){

	pt_Spectrum_and_Precision_Parameters->number_iterations_photon = number_iterations_photon;
  pt_Spectrum_and_Precision_Parameters->number_iterations_electron = number_iterations_electron;
	pt_Spectrum_and_Precision_Parameters->z_step = z_step;
	pt_Spectrum_and_Precision_Parameters->n_step = n_step;
  pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice = photon_spectrum_choice;
  pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice = electron_spectrum_choice;
  pt_Spectrum_and_Precision_Parameters->spectrum_mode = spectrum_mode;
  pt_Spectrum_and_Precision_Parameters->inverse_compton_scattering = inverse_compton_scattering;
  pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum = (*Gamma_Spectrum);
  pt_Spectrum_and_Precision_Parameters->Injected_Electron_Spectrum = (*Electron_Spectrum);
}
void check_energy_conservation(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                               struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                               struct Structure_Spectrum * pt_Gamma_Spectrum,
                               struct Structure_Spectrum * pt_Electron_Spectrum,
                               double &integrale){

double dE = (pt_Particle_Physics_Model->E_0 - E_min) / (double) (pt_Spectrum_and_Precision_Parameters->n_step - 1);
double h = dE/6.;
double dE_2, h2;
double E1, E2, E3, E4, E5, E6, E7, f1, f2, f3, f4, f5, f6, f7, g1, g2, g3, g4, g5, g6, g7, E_gamma, E_e;
double resultat_1 = 0, resultat_2 = 0, resultat_3 = 0, resultat_4 = 0, resultat_5 = 0;
double F1, F2, F3, F4, F5, F6, F7;
double z = pt_Gamma_Spectrum->redshift;
vector<double> Gamma_Spectrum_Integrated_Over_Kernel;
vector<double> Gamma_Spectrum_Integrated_Over_Kernel_energy;
vector<double> Electron_Spectrum_Integrated_Over_Kernel;
vector<double> Electron_Spectrum_Integrated_Over_Kernel_energy;
double E_gamma_bb = 2.701*T_0*(1+z);
double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
double dE_3 = (pt_Particle_Physics_Model->E_0 - E_min) / (double) (pt_Spectrum_and_Precision_Parameters->z_step - 1);


if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice=="universal"){
  E_gamma = 50;
        for(int i=0;i<pt_Spectrum_and_Precision_Parameters->n_step-1;i++){
          if(i==0){
            E1=E_min;
            }
          else{
            E1=E7;
          }

          E2=E1 + h;
          E3=E2 + h;
          E4=E3 + h;
          E5=E4 + h;
          E6=E5 + h;
          E7=E6 + h;
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

          f1*=E1;
          f2*=E2;
          f3*=E3;
          f4*=E4;
          f5*=E5;
          f6*=E6;
          f7*=E7;
          resultat_1 += dE/840. * (41*f1+216*f2+27*f3+272*f4+27*f5+216*f6+41*f7);

        }

        integrale = resultat_1;
}
else{
for(int i=0;i<pt_Spectrum_and_Precision_Parameters->z_step;i++){
  E_gamma = E_min+i*dE_3;
  E_e = E_gamma;
  resultat_1 = 0;
  resultat_2 = 0;
  resultat_3 = 0;
  resultat_4 = 0;
  double resultat_5 = 0;
  dE_2 = (pt_Particle_Physics_Model->E_0 - (E_gamma))/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
  h2 = dE_2/6.;

        for(int j=0; j<pt_Spectrum_and_Precision_Parameters->n_step-1;j++){

          if(j==0){
            E1=E_gamma;
            }
          else{
            E1=E7;
          }

          E2=E1 + h2;
          E3=E1 + 2*h2;
          E4=E1 + 3*h2;
          E5=E1 + 4*h2;
          E6=E1 + 5*h2;
          E7=E1 + 6*h2;
          /********** First : Integrate photon spectrum over gamma->gamma kernel ******/
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

          F1=f1*(dsigma_phph(E1,z,E_gamma)+dsigma_compton(E1,z,E_gamma));
          F2=f2*(dsigma_phph(E2,z,E_gamma)+dsigma_compton(E2,z,E_gamma));
          F3=f3*(dsigma_phph(E3,z,E_gamma)+dsigma_compton(E3,z,E_gamma));
          F4=f4*(dsigma_phph(E4,z,E_gamma)+dsigma_compton(E4,z,E_gamma));
          F5=f5*(dsigma_phph(E5,z,E_gamma)+dsigma_compton(E5,z,E_gamma));
          F6=f6*(dsigma_phph(E6,z,E_gamma)+dsigma_compton(E6,z,E_gamma));
          F7=f7*(dsigma_phph(E7,z,E_gamma)+dsigma_compton(E7,z,E_gamma));
          // cout << " E1 = " << E1 << endl;

          resultat_1 += dE_2/840. * (41*F1+216*F2+27*F3+272*F4+27*F5+216*F6+41*F7);
          if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "Dirac" && j==pt_Spectrum_and_Precision_Parameters->n_step-2){
            resultat_1+=(dsigma_phph(pt_Particle_Physics_Model->E_0,z,E_gamma)+dsigma_compton(pt_Particle_Physics_Model->E_0,z,E_gamma))/(rate_NPC(pt_Particle_Physics_Model->E_0,z)+rate_compton(pt_Particle_Physics_Model->E_0,z)+rate_gg_scattering(pt_Particle_Physics_Model->E_0,z));
          }
          /******** Second : Integrate electron spectrum over e->gamma kernel ********/


          if(E1+m_e<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E1+m_e, g1);
          else g1=0;
          if(E2+m_e<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E2+m_e, g2);
          else g2=0;
          if(E3+m_e<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E3+m_e, g3);
          else g3=0;
          if(E4+m_e<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E4+m_e, g4);
          else g4=0;
          if(E5+m_e<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E5+m_e, g5);
          else g5=0;
          if(E6+m_e<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E6+m_e, g6);
          else g6=0;
          if(E7+m_e<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E7+m_e, g7);
          else g7=0;


          F1=2*g1*Function_Integrand_Spectre_Compton(E_gamma,E1+m_e, E_gamma_bb)/((E1+m_e)*(E1+m_e));
          F2=2*g2*Function_Integrand_Spectre_Compton(E_gamma,E2+m_e, E_gamma_bb)/((E2+m_e)*(E2+m_e));
          F3=2*g3*Function_Integrand_Spectre_Compton(E_gamma,E3+m_e, E_gamma_bb)/((E3+m_e)*(E3+m_e));
          F4=2*g4*Function_Integrand_Spectre_Compton(E_gamma,E4+m_e, E_gamma_bb)/((E4+m_e)*(E4+m_e));
          F5=2*g5*Function_Integrand_Spectre_Compton(E_gamma,E5+m_e, E_gamma_bb)/((E5+m_e)*(E5+m_e));
          F6=2*g6*Function_Integrand_Spectre_Compton(E_gamma,E6+m_e, E_gamma_bb)/((E6+m_e)*(E6+m_e));
          F7=2*g7*Function_Integrand_Spectre_Compton(E_gamma,E7+m_e, E_gamma_bb)/((E7+m_e)*(E7+m_e));
          // cout << "E7 "<<E7<<" F7  = " << F7 << endl;

          resultat_2 += dE_2/840. * (41*F1+216*F2+27*F3+272*F4+27*F5+216*F6+41*F7);
          if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice == "Dirac" && j==pt_Spectrum_and_Precision_Parameters->n_step-2){
            resultat_2+=2*Function_Integrand_Spectre_Compton(E_gamma,pt_Particle_Physics_Model->E_0, E_gamma_bb)/((pt_Particle_Physics_Model->E_0)*(pt_Particle_Physics_Model->E_0))/Rate_Inverse_Compton(pt_Particle_Physics_Model->E_0,z,pt_Spectrum_and_Precision_Parameters);
          }
          /*******Third : Integrate photon spectrum over gamma->e kernel*******/
          F1=f1*(dsigma_compton(E1,z,(E1+m_e-E_e))+dsigma_NPC(E1+m_e,z,E_e));
          F2=f2*(dsigma_compton(E2,z,(E2+m_e-E_e))+dsigma_NPC(E2+m_e,z,E_e));
          F3=f3*(dsigma_compton(E3,z,(E3+m_e-E_e))+dsigma_NPC(E3+m_e,z,E_e));
          F4=f4*(dsigma_compton(E4,z,(E4+m_e-E_e))+dsigma_NPC(E4+m_e,z,E_e));
          F5=f5*(dsigma_compton(E5,z,(E5+m_e-E_e))+dsigma_NPC(E5+m_e,z,E_e));
          F6=f6*(dsigma_compton(E6,z,(E6+m_e-E_e))+dsigma_NPC(E6+m_e,z,E_e));
          F7=f7*(dsigma_compton(E7,z,(E7+m_e-E_e))+dsigma_NPC(E7+m_e,z,E_e));
          resultat_3 += dE_2/840. * (41*F1+216*F2+27*F3+272*F4+27*F5+216*F6+41*F7);
          if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "Dirac" && j==pt_Spectrum_and_Precision_Parameters->n_step-2){
            resultat_3+=(dsigma_compton(pt_Particle_Physics_Model->E_0,z,(pt_Particle_Physics_Model->E_0+m_e-E_e))+dsigma_NPC(pt_Particle_Physics_Model->E_0+m_e,z,pt_Particle_Physics_Model->E_0))/(rate_NPC(pt_Particle_Physics_Model->E_0,z)+rate_compton(pt_Particle_Physics_Model->E_0,z)+rate_gg_scattering(pt_Particle_Physics_Model->E_0,z));
          }


          /******** Fourth : Integrate electron spectrum over e->e kernel ********/

          if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E1, g1);
          else g1=0;
          if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E2, g2);
          else g2=0;
          if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E3, g3);
          else g3=0;
          if(E4<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E4, g4);
          else g4=0;
          if(E5<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E5, g5);
          else g5=0;
          if(E6<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E6, g6);
          else g6=0;
          if(E7<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E7, g7);
          else g7=0;

          F1=g1*Function_Integrand_Spectre_Compton(E1+E_gamma_bb-E_e,E1, E_gamma_bb)/((E1)*(E1));
          F2=g2*Function_Integrand_Spectre_Compton(E2+E_gamma_bb-E_e,E2, E_gamma_bb)/((E2)*(E2));
          F3=g3*Function_Integrand_Spectre_Compton(E3+E_gamma_bb-E_e,E3, E_gamma_bb)/((E3)*(E3));
          F4=g4*Function_Integrand_Spectre_Compton(E4+E_gamma_bb-E_e,E4, E_gamma_bb)/((E4)*(E4));
          F5=g5*Function_Integrand_Spectre_Compton(E5+E_gamma_bb-E_e,E5, E_gamma_bb)/((E5)*(E5));
          F6=g6*Function_Integrand_Spectre_Compton(E6+E_gamma_bb-E_e,E6, E_gamma_bb)/((E6)*(E6));
          F7=g7*Function_Integrand_Spectre_Compton(E7+E_gamma_bb-E_e,E7, E_gamma_bb)/((E7)*(E7));

          resultat_4+= dE_2/840. * (41*F1+216*F2+27*F3+272*F4+27*F5+216*F6+41*F7);
          if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice == "Dirac" && j==pt_Spectrum_and_Precision_Parameters->n_step-2){
            resultat_4+=Function_Integrand_Spectre_Compton(pt_Particle_Physics_Model->E_0+E_gamma_bb-E_e,pt_Particle_Physics_Model->E_0, E_gamma_bb)/((pt_Particle_Physics_Model->E_0)*(pt_Particle_Physics_Model->E_0))/Rate_Inverse_Compton(pt_Particle_Physics_Model->E_0,z,pt_Spectrum_and_Precision_Parameters);
          }

        }

        resultat_2*=2*pi*r_e*r_e*m_e*m_e*int_bb/E_gamma_bb;
        resultat_4*=2*pi*r_e*r_e*m_e*m_e*int_bb/E_gamma_bb;


        // cout << "Egamma = " << E_gamma << " resultat = " << resultat_1   << "resultat 2 =  "<< resultat_2<< " resultat 3 = " << resultat_3 << "resultat 4 = " << resultat_4<<endl;
        // if(pt_Spectrum_and_Precision_Parameters->inverse_compton_scattering=="no" && pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice=="none"){
        //   resultat_2 = 0;
        //   resultat_3 = 0;
        //   resultat_4 = 0;
        // }
        // if(pt_Spectrum_and_Precision_Parameters->inverse_compton_scattering=="no" && pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice=="none"){
        //   resultat_1 = 0;
        //   resultat_2 = 0;
        //   resultat_3 = 0;
        // }
        Gamma_Spectrum_Integrated_Over_Kernel_energy.push_back(E_gamma);
        Gamma_Spectrum_Integrated_Over_Kernel.push_back(resultat_1+resultat_2);
        Electron_Spectrum_Integrated_Over_Kernel_energy.push_back(E_e);
        Electron_Spectrum_Integrated_Over_Kernel.push_back(resultat_3+resultat_4);
  }
  resultat_1 = 0;
  resultat_2 = 0;
  resultat_3 = 0;
  resultat_4 = 0;
        for(int i=0;i<pt_Spectrum_and_Precision_Parameters->n_step-1;i++){
          if(i==0){
            E1=E_min;
            }
          else{
            E1=E7;
          }

          E2=E1 + h;
          E3=E2 + h;
          E4=E3 + h;
          E5=E4 + h;
          E6=E5 + h;
          E7=E6 + h;
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

          f1*=E1*(rate_NPC(E1,z)+rate_compton(E1,z)+rate_gg_scattering(E1,z));
          f2*=E2*(rate_NPC(E2,z)+rate_compton(E2,z)+rate_gg_scattering(E2,z));
          f3*=E3*(rate_NPC(E3,z)+rate_compton(E3,z)+rate_gg_scattering(E3,z));
          f4*=E4*(rate_NPC(E4,z)+rate_compton(E4,z)+rate_gg_scattering(E4,z));
          f5*=E5*(rate_NPC(E5,z)+rate_compton(E5,z)+rate_gg_scattering(E5,z));
          f6*=E6*(rate_NPC(E6,z)+rate_compton(E6,z)+rate_gg_scattering(E6,z));
          f7*=E7*(rate_NPC(E7,z)+rate_compton(E7,z)+rate_gg_scattering(E7,z));
          resultat_1 += dE/840. * (41*f1+216*f2+27*f3+272*f4+27*f5+216*f6+41*f7);
          if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "Dirac" &&  i==pt_Spectrum_and_Precision_Parameters->n_step-2){
            resultat_1+=pt_Particle_Physics_Model->E_0;
          }

          if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E1, f1);
          else f1=0;
          if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E2, f2);
          else f2=0;
          if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E3, f3);
          else f3=0;
          if(E4<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E4, f4);
          else f4=0;
          if(E5<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E5, f5);
          else f5=0;
          if(E6<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E6, f6);
          else f6=0;
          if(E7<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E7, f7);
          else f7=0;

          f1*=E1*(Rate_Inverse_Compton(E1,z,pt_Spectrum_and_Precision_Parameters));
          f2*=E2*(Rate_Inverse_Compton(E2,z,pt_Spectrum_and_Precision_Parameters));
          f3*=E3*(Rate_Inverse_Compton(E3,z,pt_Spectrum_and_Precision_Parameters));
          f4*=E4*(Rate_Inverse_Compton(E4,z,pt_Spectrum_and_Precision_Parameters));
          f5*=E5*(Rate_Inverse_Compton(E5,z,pt_Spectrum_and_Precision_Parameters));
          f6*=E6*(Rate_Inverse_Compton(E6,z,pt_Spectrum_and_Precision_Parameters));
          f7*=E7*(Rate_Inverse_Compton(E7,z,pt_Spectrum_and_Precision_Parameters));
          resultat_2 += dE/840. * (41*f1+216*f2+27*f3+272*f4+27*f5+216*f6+41*f7);
          if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice == "Dirac" &&  i==pt_Spectrum_and_Precision_Parameters->n_step-2){
            resultat_2+=pt_Particle_Physics_Model->E_0;
          }

          if(E1<pt_Particle_Physics_Model->E_0)linearint(Gamma_Spectrum_Integrated_Over_Kernel_energy, Gamma_Spectrum_Integrated_Over_Kernel, Gamma_Spectrum_Integrated_Over_Kernel_energy.size(), E1, g1);
          else g1=0;
          if(E2<pt_Particle_Physics_Model->E_0)linearint(Gamma_Spectrum_Integrated_Over_Kernel_energy, Gamma_Spectrum_Integrated_Over_Kernel, Gamma_Spectrum_Integrated_Over_Kernel_energy.size(), E2, g2);
          else g2=0;
          if(E3<pt_Particle_Physics_Model->E_0)linearint(Gamma_Spectrum_Integrated_Over_Kernel_energy, Gamma_Spectrum_Integrated_Over_Kernel, Gamma_Spectrum_Integrated_Over_Kernel_energy.size(), E3, g3);
          else g3=0;
          if(E4<pt_Particle_Physics_Model->E_0)linearint(Gamma_Spectrum_Integrated_Over_Kernel_energy, Gamma_Spectrum_Integrated_Over_Kernel, Gamma_Spectrum_Integrated_Over_Kernel_energy.size(), E4, g4);
          else g4=0;
          if(E5<pt_Particle_Physics_Model->E_0)linearint(Gamma_Spectrum_Integrated_Over_Kernel_energy, Gamma_Spectrum_Integrated_Over_Kernel, Gamma_Spectrum_Integrated_Over_Kernel_energy.size(), E5, g5);
          else g5=0;
          if(E6<pt_Particle_Physics_Model->E_0)linearint(Gamma_Spectrum_Integrated_Over_Kernel_energy, Gamma_Spectrum_Integrated_Over_Kernel, Gamma_Spectrum_Integrated_Over_Kernel_energy.size(), E6, g6);
          else g6=0;
          if(E7<pt_Particle_Physics_Model->E_0)linearint(Gamma_Spectrum_Integrated_Over_Kernel_energy, Gamma_Spectrum_Integrated_Over_Kernel, Gamma_Spectrum_Integrated_Over_Kernel_energy.size(), E7, g7);
          else g7=0;

          g1*=E1;
          g2*=E2;
          g3*=E3;
          g4*=E4;
          g5*=E5;
          g6*=E6;
          g7*=E7;

          resultat_3+=dE/840. * (41*g1+216*g2+27*g3+272*g4+27*g5+216*g6+41*g7);


          if(E1<pt_Particle_Physics_Model->E_0)linearint(Electron_Spectrum_Integrated_Over_Kernel_energy, Electron_Spectrum_Integrated_Over_Kernel, Electron_Spectrum_Integrated_Over_Kernel_energy.size(), E1, g1);
          else g1=0;
          if(E2<pt_Particle_Physics_Model->E_0)linearint(Electron_Spectrum_Integrated_Over_Kernel_energy, Electron_Spectrum_Integrated_Over_Kernel, Electron_Spectrum_Integrated_Over_Kernel_energy.size(), E2, g2);
          else g2=0;
          if(E3<pt_Particle_Physics_Model->E_0)linearint(Electron_Spectrum_Integrated_Over_Kernel_energy, Electron_Spectrum_Integrated_Over_Kernel, Electron_Spectrum_Integrated_Over_Kernel_energy.size(), E3, g3);
          else g3=0;
          if(E4<pt_Particle_Physics_Model->E_0)linearint(Electron_Spectrum_Integrated_Over_Kernel_energy, Electron_Spectrum_Integrated_Over_Kernel, Electron_Spectrum_Integrated_Over_Kernel_energy.size(), E4, g4);
          else g4=0;
          if(E5<pt_Particle_Physics_Model->E_0)linearint(Electron_Spectrum_Integrated_Over_Kernel_energy, Electron_Spectrum_Integrated_Over_Kernel, Electron_Spectrum_Integrated_Over_Kernel_energy.size(), E5, g5);
          else g5=0;
          if(E6<pt_Particle_Physics_Model->E_0)linearint(Electron_Spectrum_Integrated_Over_Kernel_energy, Electron_Spectrum_Integrated_Over_Kernel, Electron_Spectrum_Integrated_Over_Kernel_energy.size(), E6, g6);
          else g6=0;
          if(E7<pt_Particle_Physics_Model->E_0)linearint(Electron_Spectrum_Integrated_Over_Kernel_energy, Electron_Spectrum_Integrated_Over_Kernel, Electron_Spectrum_Integrated_Over_Kernel_energy.size(), E7, g7);
          else g7=0;

          g1*=E1;
          g2*=E2;
          g3*=E3;
          g4*=E4;
          g5*=E5;
          g6*=E6;
          g7*=E7;

          resultat_4+=dE/840. * (41*g1+216*g2+27*g3+272*g4+27*g5+216*g6+41*g7);

          // cout << " E7 = " << E7 << "resultat_1 " << resultat_1 << " resultat_2 = " << resultat_2 << "resultat_3 = " << resultat_3 <<  " i = " << i <<  endl;

        }
        integrale = (resultat_1-resultat_3)+(resultat_2-resultat_4);
      }
      if(verbose>1){cout << "The total energy contained in " <<   pt_Gamma_Spectrum->spectrum_name  << "spectrum is " << resultat_1 << " - " << resultat_3 << " = " << resultat_1-resultat_3 << " MeV";
        cout << " and in " << pt_Electron_Spectrum->spectrum_name << "spectrum is " << resultat_2 << " - " << resultat_4 << " = " << resultat_2-resultat_4 << " MeV ";
        cout << "for a total of "<< integrale <<" MeV, you had injected " << pt_Particle_Physics_Model->E_0 << " MeV." << endl;
      }

}
void Fill_Output_Options(string print_result, string results_files, string spectrum_files, struct Structure_Output_Options * pt_Structure_Output_Options){

  pt_Structure_Output_Options->print_result = print_result;
  pt_Structure_Output_Options->results_files = results_files;
  pt_Structure_Output_Options->spectrum_files = spectrum_files;
}
