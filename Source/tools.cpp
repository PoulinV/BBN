#include "../Include/tools.h"
#include "../Include/EM_cascade.h"


void print_spectrum(ostream &file, struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum, struct Structure_Particle_Physics_Model * pt_Particle_Model){


  double dE = (pt_Particle_Model->E_0 - E_min)/Gamma_Table_Size;
  double z = pt_Cascade_Spectrum->redshift;
  double E = E_min;
  int i = 0;
  while(i<Gamma_Table_Size){
    file << E << "   " << pt_Cascade_Spectrum->Gamma_Spectrum[i]/(gamma_NPC(E,z)+gamma_compton(E,z)+gamma_phph(E,z)) << endl;
    i++;
    E+=dE;
  }

}
void print_spectrum_from_function(ostream &file, double (*func)(double,double,double),double z, struct Structure_Particle_Physics_Model * pt_Particle_Model){


  double dE = (pt_Particle_Model->E_0 - E_min)/Gamma_Table_Size;
  double E = E_min;
  int i = 0;
  while(i<Gamma_Table_Size){
    file << E << "   " << (*func)(E,z,pt_Particle_Model->E_0)/(gamma_NPC(E,z)+gamma_compton(E,z)+gamma_phph(E,z)) << endl;
    i++;
    E+=dE;
  }

}
void print_spectrum_automatic_names(int iterations, struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum, struct Structure_Particle_Physics_Model * pt_Particle_Model){
  ostringstream os;
  string name;


  double dE = (pt_Particle_Model->E_0 - E_min)/pt_Cascade_Spectrum->Gamma_Energy.size();
  double z = pt_Cascade_Spectrum->redshift;
  double E = E_min;
  int i = 0;

	  mkdir("Cascade_Spectrum_Folder", 01777);
    os << "Cascade_Spectrum_Folder/Spectrum_m" << pt_Particle_Model->M_x<<"_z"<< z <<"_" << iterations <<"iterations.dat";
    cout << " T = " << T_0*(1+z) << endl;
    name = os.str();
    ofstream file(name);
    if(verbose>1)cout << "Printing in file " << name << endl;

    while(i<Gamma_Table_Size){
      if(i<Gamma_Table_Size-1)file << pt_Cascade_Spectrum->Gamma_Energy[i] << "   " << pt_Cascade_Spectrum->Gamma_Spectrum[i]/(gamma_NPC(pt_Cascade_Spectrum->Gamma_Energy[i],z)+gamma_compton(pt_Cascade_Spectrum->Gamma_Energy[i],z)+gamma_phph(pt_Cascade_Spectrum->Gamma_Energy[i],z)) << endl;
      if(i==Gamma_Table_Size-1)file << pt_Cascade_Spectrum->Gamma_Energy[i] << "   " << pt_Cascade_Spectrum->Gamma_Spectrum[i]/(gamma_NPC(pt_Cascade_Spectrum->Gamma_Energy[i],z)+gamma_compton(pt_Cascade_Spectrum->Gamma_Energy[i],z)+gamma_phph(pt_Cascade_Spectrum->Gamma_Energy[i],z))+1/(gamma_NPC(pt_Particle_Model->E_0,z)+gamma_compton(pt_Particle_Model->E_0,z)+gamma_phph(pt_Particle_Model->E_0,z)) << endl;
      i++;

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

void Fill_Output_Options(string print_result, string results_files, string spectrum_files, struct Structure_Output_Options * pt_Structure_Output_Options){

  pt_Structure_Output_Options->print_result = print_result;
  pt_Structure_Output_Options->results_files = results_files;
  pt_Structure_Output_Options->spectrum_files = spectrum_files;
}
