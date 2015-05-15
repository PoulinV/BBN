#define FUNC(x,a,i,g,E_O,Z_x,z) ((*func)(x,a,i,g,E_O,Z_x,z))
#include <iostream>
using namespace std;
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#define EPS 1.0e-1
#define JMAX 10
#include <sstream>
#include "Abondancev5_3He.h"
            /******Programme*******/
/******Calcul graphe tau-zeta pour Helium-3, spectre standard******/
/*****ATTENTION - vérifier Abondancev5_2.h pour savoir quel spectre est utilisé*****/
int main()
{

  double Emin[19];
  Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
  //~ for(int i=0;i<14;i++)Emin[i]+=0.;
  double T;
  double T_0 = 2.7255*0.862*pow(10,-10);
  double H_0 = 2.187*pow(10,-18);
  double H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5), Z_x = pow(10,-6) , n_y_0=3.154*pow(10,-30) , E_0 ;
  double B ;
  double z,s_Li,s_Be,s_He,s_De,A,dY_Li,dY_3_He,dY_De,dY_He,A_theo;
  double omega_b = 0.02207;
  double tau_n = 880.3, tau_n_0 = 880.3;
  double Y_0 = 1.02*pow(10,-5)*pow(omega_b/0.02273,-0.59)*pow(tau_n/tau_n_0,0.15), Y_0_De =  2.53*pow(10,-5)*pow(omega_b/0.02273,-1.62)*pow(tau_n/tau_n_0,0.41);
  double tau, Ex , Ec, s,s2;
  double pi=3.14159;
  double Y_De_max_cyburt =1*pow(10,-4), Y_De_max_standard = 3.48*pow(10,-5),Y_De_min_cyburt=1*pow(10,-5),Y_De_min_standard = 2.56*pow(10,-5) ;

  //~ ofstream output_De("results_Graphe_De.dat");
  //~ ofstream output_He("results_Graphe_He_test.dat");
  // ofstream output_3He("results_Graphe_3He_sans_reinj.dat");
  ofstream output_Best_Constraints("Best_Constraints_70MeV_avec_reinjection_two_iterations.dat");
  // ofstream output_Best_Constraints("Best_Constraints_spectre_universel_new.dat");
  double M_x,E_max;
  M_x=140;
  E_0=M_x/2;
  //~ T=1*pow(10,-1);

  tau=pow(10,4);double tau_moy=2*pow(10,6);double tau_moy2=pow(10,8);double tau_max=pow(10,10);
  //~ double T_min = 0.511*0.511/(1.6*22);
  //~ double T_max = pow(tau_max*(2*H_r),-0.5)*T_0;
  //~ double T_max = 0.511*0.511/(2.2*22);
  //~ tau=pow(T_min/T_0,-2)/(2*H_r);
  //~ tau_max=pow(T_max/T_0,-2)/(2*H_r);
  cout << " tau min = " << tau << " tau max = "<< tau_max << " E_0 = "<< E_0<<endl;
  //~ double z_max = 35000;
  double z_max = 10000;
  double h_moy = 100;
  double h_moy2 = 2*h_moy;
  double h_max = 3*h_moy;
  //~ double Delta = (tau_max-tau)/h_moy;
  double Delta = (tau_moy-tau)/h_moy;
  double Delta2 = (tau_moy2-tau_moy)/h_moy;
  double Delta3 = (tau_max-tau_moy2)/h_moy;
  cout << "Delta 1 = " << Delta << " Delta 2 = " << Delta2 << " Delta 3 = " << Delta3 << endl;
  //~ tau=tau_moy;
  double count;
  for(int h=1;h<h_max;h++)
  {
    if(h==h_moy)Delta=Delta2;
    if(h==h_moy2)Delta=Delta3;
    if(h>1)tau+=Delta;
    T=pow(0.2*tau*(2*H_r),-0.5)*T_0;
    z = T/T_0 - 1;
    count=0;

    cout << "****************h = "<<h<<"  tau = "<<tau<<"***************"<<endl;
    cout << " T= " << T << " tau = " << tau <<" z = "<< z << endl;






    for(int l=12;l>1;l--)
    {
      for(double k=1;k<10;k+=0.1)
      {
        Z_x=k*pow(10,-l);
        if(tau>pow(10,6)){A = qsimp(K_Perte_3He_non_standard,z_max,z,z,16,E_0,Z_x,z);
        B = qsimp(S_Gain_3He_non_standard,z_max,z,z,16,E_0,Z_x,z);}
        // A = qsimp(K_Perte_3He_sans_reinj,z_max,z,z,16,E_0,Z_x,z);
        // B = qsimp(S_Gain_3He_sans_reinj,z_max,z,z,16,E_0,Z_x,z);
        // cout << " A = "<< A << " B = " << B << endl;
      //~
        //~ for(int n=1;n<2;n++)B = trapzd(S_Gain,z_max,z,z,2,14,E_0,Z_x,z);
                  //~ test = trapzd(func_x,z_max,z,z,n,14,E_0,Z_x,z);
                  //~ }
        //~ cout << " A = " << A << " B = " << B << endl;
        //~ cout << "Gain = " << S_Gain(z,z,14,z,E_0,Z_x,z)<< endl;
        dY_3_He=exp(-Z_x*n_y_0/(E_0*H_r*tau)*A)*(Y_0 + Z_x*n_y_0/(E_0*H_r*tau)*B);

        if(tau<pow(10,7)){A = qsimp(K_Perte_non_standard,z_max,z,z,14,E_0,Z_x,z);
        B = qsimp(S_Gain_non_standard,z_max,z,z,14,E_0,Z_x,z);}

        // A = qsimp(K_Perte_sans_reinj,z_max,z,z,14,E_0,Z_x,z);
        // B = qsimp(S_Gain_sans_reinj,z_max,z,z,14,E_0,Z_x,z);
        dY_De=exp(-Z_x*n_y_0/(E_0*H_r*tau)*A)*(Y_0_De + Z_x*n_y_0/(E_0*H_r*tau)*B);
        // dY_3_He=(Y_0 +Z_x*n_y_0/(E_0*H_r*tau)*B);
      //  cout << "dY3He = " << dY_3_He << endl;
        if(dY_3_He>1.5*pow(10,-5) && count==0)
        {
          output_Best_Constraints << Z_x <<"  "<< tau  << endl;
          count++;
          cout << "Y3He = " << dY_3_He << endl;
          // cout << " s = " << s << "dY3He = " << dY_3_He << endl;
        }
        if(dY_De>Y_De_max_standard && count==0)
        {
          output_Best_Constraints << Z_x <<"  "<< tau  << endl;
          count++;
          cout << "YDe = " << dY_De << endl;
        }
        if(dY_De<Y_De_min_standard&& count == 0)
        {
          output_Best_Constraints << Z_x <<"  "<< tau  << endl;
          count++;
          cout << "YDe = " << dY_De << endl;
        }

        if(count==1)break;
      }
      if(count==1)break;
    }
  }

  //~ output_De.close();
  //~ output_He.close();
  output_Best_Constraints.close();
  return 0;
}
