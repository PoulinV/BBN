#define FUNC(x,a,i,g,E_O,Z_x,z_0) ((*func)(x,a,i,g,E_O,Z_x,z_0))
#include <iostream>
using namespace std;
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#define EPS 1.0e-1
#define JMAX 8
#include <sstream>
#include "../Include/Abondancev5_2.h"
            /******Programme*******/
/******Calcul graphe tau-zeta pour Deuterium et Helium, spectre standard******/
int main()
{

  double Emin[21];
  Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ; Emin[19] = 5.49 ; Emin[20] = 7.72 ;
  //~ for(int i=0;i<14;i++)Emin[i]+=0.;
  double T;
  double T_0 = 2.7255*0.862*pow(10,-10);
  double H_0 = 2.187*pow(10,-18);
  double H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5), Z_x = pow(10,-6) , n_y_0=3.154*pow(10,-30) , E_0 ;
  double B ;
  double z,s_Li,s_Be,s_He,s_De,A,dY_Li,dY_Be,dY_De,dY_He,A_theo;

  double tau, Ex , Ec, s,s2;
  double pi=3.14159;

  ofstream output_FonctionAna("../TableIntegraleSpectre/results_quatre_iterations_Produc2H_30MeV_corrected.dat");
  ifstream file("../OutputExtra/Correction4iterations_30MeV.dat");
  fill_table_from_file(file,vector_E_correction,vector_correction);
  for(int i; i<vector_E_correction.size();i++)  cout<< vector_E_correction[i] << "  " << vector_correction[i] << endl;
  //~ ofstream output_Be("results_Graphe_Be.dat");
  double M_x,E_max;
  // M_x=300;
  // E_0=M_x/2;
  E_0 = 30;
  T=1*pow(10,-1);

  tau=pow(10,4);double tau_moy=pow(10,6);double tau_moy2=pow(10,8);double tau_moy3=pow(10,10);double tau_max=pow(10,12);
  //~ double T_min = 0.511*0.511/(1.6*22);
  double T_max = pow(tau_max*(2*H_r),-0.5)*T_0;
  //~ double T_max = 0.511*0.511/(2.2*22);
  //~ tau=pow(T_min/T_0,-2)/(2*H_r);tau_max=pow(T_max/T_0,-2)/(2*H_r);
  cout << " tau min = " << tau << " tau max = "<< tau_max << endl;
  double z_max = 35000;
  double h_moy = 100;
  double h_moy2 = 2*h_moy;
    double h_moy3 = 3*h_moy;
  double h_max = 4*h_moy;

  double Delta = (tau_moy-tau)/h_moy;
  double Delta2 = (tau_moy2-tau_moy)/h_moy;
  double Delta3 = (tau_moy3-tau_moy2)/h_moy;
    double Delta4 = (tau_max-tau_moy2)/h_moy;
  cout << "Delta 1 = " << Delta << " Delta 2 = " << Delta2 << " Delta 3 = " << Delta3 << endl;
  for(int h=1;h<h_max;h++)
  {
    if(h==h_moy)Delta=Delta2;
    if(h==h_moy2)Delta=Delta3;
        if(h==h_moy3)Delta=Delta4;
    if(h>1)tau+=Delta;
    T=pow(0.2*tau*(2*H_r),-0.5)*T_0;
    z = T/T_0 - 1;
      Ex = 0.261121/(80*T_0*(1+z)), Ec = 0.261121/(22*T_0*(1+z));

    cout << "****************h = "<<h<<"  tau = "<<tau<<"***************"<<endl;
    cout << "Ex = " << Ex << " Ec = " << Ec << " T= " << T << " tau = " << tau <<" z = "<< z << endl;
    s=0;

      for(int i=17;i<18;i++)
      {
                if(Ec>Emin[i]){
                  //  s+=qsimp_Z(func_standard,Emin[i],Ec,z,i,E_0,0,z);}
                  //  s+=qsimp_Z(func_quatre_interactions_Produc_2H,Emin[i],E_0,z,i,E_0,0,z);}
									// if(i==17)s+=2*qsimp_Z(func_quatre_interactions_produc2H,Emin[i],E_0,z,i,E_0,0,z);
									if(i==17)s+=qsimp_Z(func_quatre_interactions_produc2H,Emin[i],E_0,z,i,E_0,0,z);}
          //~ cout << "s + = " << s << endl;
      }
      cout << " E_0 = " << E_0 << " s = " << s << endl;
//			M_x=60;
//			E_0=M_x/2;
//			s=0;
//			for(int i=17;i<19;i++)
//			{
//				if(Ec>Emin[i]){
//					// s+=qsimp(func_standard,Emin[i],Ec,z,i,E_0);}
//					if(i==17)s+=2*qsimp(func_standard,Emin[i],Ec,z,i,E_0);
//					if(i==18)s+=qsimp(func_standard,Emin[i],Ec,z,i,E_0);}
//						//~ cout << "s + = " << s << endl;
//					}
//					cout << " E_0 = " << E_0 << " s = " << s*5000/E_0 << endl;
      output_FonctionAna << log10(z) << "  "  << log10(s)<< endl;
  }

  output_FonctionAna.close();

  //~ output_Be.close();
  return 0;
}
