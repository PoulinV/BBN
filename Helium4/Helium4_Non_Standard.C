#define FUNC(x,a,i,g,E_O,Z_x,z_0) ((*func)(x,a,i,g,E_O,Z_x,z_0))
#include <iostream>
using namespace std;
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#define EPS 1.0e-2
#define JMAX 10
#include <sstream>
#include "../Include/Abondancev5_2.h"
            /******Programme*******/
/******Calcul graphe tau-zeta pour Deuterium, spectre dirac et Destruction uniquement******/
/*****ATTENTION - vérifier Abondancev4.h pour savoir quel spectre est utilisé*****/
int main()
{

  float Emin[19];
  Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
  //~ for(int i=0;i<14;i++)Emin[i]+=0.;
  float T;
  double T_0 = 2.7255*0.862*pow(10,-10);
  double H_0 = 2.187*pow(10,-18);
  double H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5), Z_x = pow(10,-6) , n_y_0=3.154*pow(10,-30) , E_0 ;
  double B ;
  double z,z_0,s_Li,s_De,s_He,A,dY_Li,dY_De,dY_4He,A_theo,Y_0 = 2.63*pow(10,-5);
  double omega_b = 0.02273;
  double tau_n = 880.3, tau_n_0 = 880.3;
  double Y_4He_0 = 0.248*pow(omega_b/0.02273,0.39)*pow(tau_n/tau_n_0,0.72), Y_4He_min = 0.2368; //0.2368
  cout << " Y_4He_0 = " << Y_4He_0 << endl;
  float tau, Ex , Ec, s,s2;
  double pi=3.14159;
  int writing_variable;
  double Z_x_old;
  // ofstream output_De("../Output4He/results_Graphe_He_standard_sans_reinj_corrected.dat");
  ofstream output_De("Output4He/results_Graphe_He_non_standard_avec_reinj_30MeV_30juil2015.dat");
  ifstream file0("../TableIntegraleSpectre/results_standard_Destruc4He.dat");

  ifstream file("../TableIntegraleSpectre/results_Reinj_Monochromatique_Destruc4He_30MeV.dat");
  ifstream file2("../TableIntegraleSpectre/results_deux_iterations_Destruc4He_30MeV.dat");
  ifstream file3("../TableIntegraleSpectre/results_trois_iterations_Destruc4He_30MeV.dat");
  ifstream file4("../TableIntegraleSpectre/results_quatre_iterations_Destruc4He_30MeV.dat");
  // ofstream output_De("Output4He/results_Graphe_He_non_standard_avec_reinj_70MeV_trois_iterations.dat");
  // ifstream file("TableIntegraleSpectre/results_Reinj_Monochromatique_Destruc4He_70MeV.dat");
  // ifstream file2("TableIntegraleSpectre/results_deux_iterations_Destruc4He_70MeV.dat");
  // ifstream file3("TableIntegraleSpectre/results_trois_iterations_Destruc4He_70MeV.dat");
  fill_table_from_file(file0,vector_z_s4He_destruc_standard,vector_s4He_destruc_standard);

  fill_table_from_file(file,vector_z_s4He_destruc,vector_s4He_destruc);
  // for(int i=0;i<vector_z.size();i++)output_Check << vector_z[i]<< " "<< vector_s4He_destruc[i] << endl;
  fill_table_from_file(file2,vector_z_s4He_destruc_deux_iterations,vector_s4He_destruc_deux_iterations);
  fill_table_from_file(file3,vector_z_s4He_destruc_trois_iterations,vector_s4He_destruc_trois_iterations);
  fill_table_from_file(file4,vector_z_s4He_destruc_quatre_iterations,vector_s4He_destruc_quatre_iterations);
  // for(int i=0;i<vector_z_s4He_destruc_quatre_iterations.size();i++)cout << vector_z_s4He_destruc_quatre_iterations[i]<< " "<< vector_s4He_destruc_quatre_iterations[i] << endl;
  int x =0;
  Attribution_avec_correction(x);
  double M_x,E_max;
  M_x=60;
  E_0=M_x/2;
  //~ T=1*pow(10,-1);

  double tau_max=pow(10,10);double tau_moy=pow(10,6);double tau_moy2=pow(10,8);
  double T_min = 0.511*0.511/(22*E_0);
  //~ double T_max = pow(tau_max*(2*H_r),-0.5)*T_0;
  //~ double T_max = 0.511*0.511/(2.2*22);
  // tau=5*pow(T_min/T_0,-2)/(2*H_r);

  tau=pow(10,4);
  //~ tau_max=pow(T_max/T_0,-2)/(2*H_r);
  // cout << " tau min = " << tau << " tau max = "<< tau_max << endl;
  double z_max = 35000;
  double h_moy= 100;
  double h_moy2 = 2*h_moy;
  double h_max = 3*h_moy;
  // double h_max = 3*h_moy;
  double Delta = (tau_moy-tau)/h_moy;
  double Delta2 = (tau_moy2-tau_moy)/h_moy;
  double Delta3 = (tau_max-tau_moy2)/h_moy;
  // cout << "Delta 1 = " << Delta << " Delta 2 = " << Delta2 << " Delta 3 = " << Delta3 << endl;
  for(int h=0;h<h_max;h++)
  {
    writing_variable=0;
    if(h==h_moy){Delta=Delta2;}
    if(h==h_moy2){Delta=Delta3;}
    if(h>1)tau+=Delta;
    T=pow(0.2*tau*(2*H_r),-0.5)*T_0;
    z = T/T_0 - 1;
    // if(h==0)z_0 = z;
    Ex = 0.261121/(80*T), Ec = 0.261121/(22*T);
    cout << "****************h = "<<h<<"  tau = "<<tau<<"***************"<<endl;
    cout << "Ex = " << Ex << " Ec = " << Ec << " T= " << T << " tau = " << tau <<" z = "<< z << endl;
    s_Li = 0; s_De = 0;
    s_De = 0; s_He = 0;

    s_He=qsimp(K_Perte_4He,z_max,z,z,15,E_0,Z_x,z);
      // s_He=qsimp_Z(K_perte_non_standard_4He,z_max,z,z,15,E_0,Z_x,z);
      // cout <<z_0 << " "<< z << "  "<< K_perte_non_standard_4He( z,  z, 15,  15.,  E_0, Z_x, z_0) << endl;
// K_perte_non_standard_4He( z,  z, 15,  15.,  E_0, Z_x, z_0);

    for(int l=10;l>1;l--)
    {

      for(double k=1;k<10;k+=0.01)
      {
        Z_x=k*pow(10,-l);
            //  cout << Z_x << endl;
        B=Z_x*n_y_0/(E_0*H_r*tau);
        dY_4He=s_He*B;
        if(exp(-1*dY_4He)*Y_4He_0<Y_4He_min && writing_variable == 0 )
        {
          // cout << " Y_De = " << exp(-1*dY_4He)*Y_0 << " tau =  " << tau << " Z_x " << Z_x << endl;
          output_De <<tau <<"  "<<  Z_x << "  " << endl;
          writing_variable = 1;
          break;
        }
      }
      if(writing_variable == 1)break;
    }
  }
  output_De.close();
  // output_Check.close();

  return 0;
}
