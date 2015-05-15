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
void print_spectrum_automatic_names(int number_files, struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum, struct Structure_Particle_Physics_Model * pt_Particle_Model){
  ostringstream os;
  string name;


  double dE = (pt_Particle_Model->E_0 - E_min)/Gamma_Table_Size;
  double z = pt_Cascade_Spectrum->redshift;
  double E = E_min;
  int i = 0;


    os << "Output/file_" << number_files << ".dat";
    name = os.str();
    ofstream file(name);

    while(i<Gamma_Table_Size){
      file << E << "   " << pt_Cascade_Spectrum->Gamma_Spectrum[i]/(gamma_NPC(E,z)+gamma_compton(E,z)+gamma_phph(E,z)) << endl;
      i++;
      E+=dE;
    }




}

void linearint(vector<double> &xa, vector<double> &ya, int n, double x, double &y){
				int i,m,ns=1;
				float x0,x1,y0,y1;
				double dift,dif;
				// cout << "xa[1]" << xa[1] << endl;
				dif=sqrt(pow(x-xa[1],2));
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
        if(i!=n-1){
          x1=xa[ns+1];
          y1=ya[ns+1];
          if(x1-x0!=0.0)y=y0+(x-x0)/(x1-x0)*fabs(y1-y0);

        }
}
