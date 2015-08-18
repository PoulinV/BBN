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

void Fill_Structure_Spectrum_and_Precision_Parameters(int iterations,
                                                      int z_step,
                                                      int n_step,
                                                      string photon_spectrum_choice,
                                                      string electron_spectrum_choice,
                                                      string spectrum_mode,
                                                      string inverse_compton_scattering,
                                                      double (*Gamma_Spectrum)(double, double, double),
                                                      double (*Electron_Spectrum)(double, double, double),
                                                      struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){

	pt_Spectrum_and_Precision_Parameters->iterations = iterations;
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
                               struct Structure_Spectrum * pt_Cascade_Spectrum,
                               double &integrale){

double dE = (pt_Particle_Physics_Model->E_0 - E_min) / (double) (pt_Spectrum_and_Precision_Parameters->n_step - 1);
double h = dE/6.;
double E1, E2, E3, E4, E5, E6, E7, f1, f2, f3, f4, f5, f6, f7, resultat = 0;
double z = pt_Cascade_Spectrum->redshift;
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
          if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Cascade_Spectrum->Energy, pt_Cascade_Spectrum->Spectrum, pt_Cascade_Spectrum->Energy.size(), E1, f1);
          else f1=0;
          if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Cascade_Spectrum->Energy, pt_Cascade_Spectrum->Spectrum, pt_Cascade_Spectrum->Energy.size(), E2, f2);
          else f2=0;
          if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Cascade_Spectrum->Energy, pt_Cascade_Spectrum->Spectrum, pt_Cascade_Spectrum->Energy.size(), E3, f3);
          else f3=0;
          if(E4<pt_Particle_Physics_Model->E_0)linearint(pt_Cascade_Spectrum->Energy, pt_Cascade_Spectrum->Spectrum, pt_Cascade_Spectrum->Energy.size(), E4, f4);
          else f4=0;
          if(E5<pt_Particle_Physics_Model->E_0)linearint(pt_Cascade_Spectrum->Energy, pt_Cascade_Spectrum->Spectrum, pt_Cascade_Spectrum->Energy.size(), E5, f5);
          else f5=0;
          if(E6<pt_Particle_Physics_Model->E_0)linearint(pt_Cascade_Spectrum->Energy, pt_Cascade_Spectrum->Spectrum, pt_Cascade_Spectrum->Energy.size(), E6, f6);
          else f6=0;
          if(E7<pt_Particle_Physics_Model->E_0)linearint(pt_Cascade_Spectrum->Energy, pt_Cascade_Spectrum->Spectrum, pt_Cascade_Spectrum->Energy.size(), E7, f7);
          else f7=0;
          f1*=E1*(rate_NPC(E1,z)+rate_compton(E1,z)+rate_gg_scattering(E1,z));
          f2*=E2*(rate_NPC(E2,z)+rate_compton(E2,z)+rate_gg_scattering(E2,z));
          f3*=E3*(rate_NPC(E3,z)+rate_compton(E3,z)+rate_gg_scattering(E3,z));
          f4*=E4*(rate_NPC(E4,z)+rate_compton(E4,z)+rate_gg_scattering(E4,z));
          f5*=E5*(rate_NPC(E5,z)+rate_compton(E5,z)+rate_gg_scattering(E5,z));
          f6*=E6*(rate_NPC(E6,z)+rate_compton(E6,z)+rate_gg_scattering(E6,z));
          f7*=E7*(rate_NPC(E7,z)+rate_compton(E7,z)+rate_gg_scattering(E7,z));
          // f1*=E1;
          // f2*=E2;
          // f3*=E3;
          // f4*=E4;
          // f5*=E5;
          // f6*=E6;
          // f7*=E7;

          resultat += dE/840. * (41*f1+216*f2+27*f3+272*f4+27*f5+216*f6+41*f7);

        }
        integrale = resultat;
      if(verbose>1)cout << "The total energy contained in " <<   pt_Cascade_Spectrum->spectrum_name  << "spectrum is " << resultat << " MeV, you had injected " << pt_Particle_Physics_Model->E_0 << " MeV." << endl;
}
void Fill_Output_Options(string print_result, string results_files, string spectrum_files, struct Structure_Output_Options * pt_Structure_Output_Options){

  pt_Structure_Output_Options->print_result = print_result;
  pt_Structure_Output_Options->results_files = results_files;
  pt_Structure_Output_Options->spectrum_files = spectrum_files;
}
