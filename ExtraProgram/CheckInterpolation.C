#define FUNC(x,a,i,g,E_O,Z_x,z_0) ((*func)(x,a,i,g,E_O,Z_x,z_0))
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <iostream>
#include <string>
using namespace std;
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#define EPS 1.0e-1
#define JMAX 10
#include <sstream>
#include "../Include/Abondancev5_2.h"
/******Programme*******/
/******Calcul graphe tau-zeta pour Deutérium AVEC PRODUCTION via Helium, spectre non standard******/
/*****ATTENTION - vérifier Abondancev5_2.h pour savoir quel spectre est utilisé*****/
/****** Pour les intervalles des abondances, voir arxiv:1307.6955*****/
int main()
{
    double Y_De_max_cyburt =1*pow(10,-4), Y_De_max_standard = 3.48*pow(10,-5),Y_De_min_cyburt=1*pow(10,-5),Y_De_min_standard = 2.56*pow(10,-5) ;
    float Emin[19];
    Emin[0] = 2.47 ; Emin[1] = 7.25 ; Emin[2] = 10.95 ;  Emin[3] = 8.725 ; Emin[4] = 9.98 ; Emin[5] = 22.28 ; Emin[6] = 23.05 ; Emin[7] = 1.59; Emin[8] = 5.61 ; Emin[9] = 9.31 ; Emin[10] = 7.08 ; Emin[11] = 10.68 ; Emin[12] = 22.17 ; Emin[13] = 21.4 ; Emin[14] = 2.23 ; Emin[15] = 19.82 ; Emin[16] = 20.58 ; Emin[17] = 23.85 ; Emin[18] = 26.08 ;
    //~ for(int i=0;i<14;i++)Emin[i]+=0.;
    float T;
    double T_0 = 2.7255*0.862*pow(10,-10);
    double H_0 = 2.187*pow(10,-18);
    double H_r = 2.187*pow(10,-18)*pow((1+7/8*pow(4/11,4/3)*3.046)*5.46*pow(10,-5),0.5), Z_x = pow(10,-6) , n_y_0=3.154*pow(10,-30) , E_0 ;
    double omega_b = 0.02207;
    double tau_n = 880.3, tau_n_0 = 880.3;
    double z,count1,count2,count3,Y_0 = 2.53*pow(10,-5)*pow(omega_b/0.02273,-1.62)*pow(tau_n/tau_n_0,0.41),Y_De,A,B;
    cout << " Y_2H_0 = " << Y_0 << endl;

    float tau, Ex , Ec, s,s2;
    double pi=3.14159;

    ofstream output_Check("../OutputExtra/Interpolation.dat");
    ifstream file("../TableIntegraleSpectre/results_Reinj_Monochromatique_Destruc4He_4MeV.dat");
    ifstream file2("../TableIntegraleSpectre/results_deux_iterations_Destruc4He_4MeV.dat");
    ifstream file3("../TableIntegraleSpectre/results_trois_iterations_Destruc4He_4MeV.dat");
    ifstream file4("../TableIntegraleSpectre/results_Reinj_Monochromatique_Destruc2H_4MeV.dat");
    ifstream file5("../TableIntegraleSpectre/results_deux_iterations_Destruc2H_4MeV.dat");
    ifstream file6("../TableIntegraleSpectre/results_trois_iterations_Destruc2H_4MeV.dat");
    ifstream file7("../TableIntegraleSpectre/results_Reinj_Monochromatique_Produc2H_4MeV.dat");
    ifstream file8("../TableIntegraleSpectre/results_deux_iterations_Produc2H_4MeV.dat");
    ifstream file9("../TableIntegraleSpectre/results_trois_iterations_Produc2H_4MeV.dat");

    double M_x,E_max;
    M_x=8;
    E_0=M_x/2;
    //~ T=1*pow(10,-1);

      fill_table_from_file(file,vector_z_s4He_destruc,vector_s4He_destruc);
      // for(int i=0;i<vector_z.size();i++)output_Check << vector_z[i]<< " "<< vector_s4He_destruc[i] << endl;
      fill_table_from_file(file2,vector_z_s4He_destruc_deux_iterations,vector_s4He_destruc_deux_iterations);
      fill_table_from_file(file3,vector_z_s4He_destruc_trois_iterations,vector_s4He_destruc_trois_iterations);
      // for(int i=0;i<vector_z_s4He_destruc_trois_iterations.size();i++)output_Check << vector_z_s4He_destruc_trois_iterations[i]<< " "<< vector_s4He_destruc_trois_iterations[i] << endl;
        fill_table_from_file(file4,vector_z_s2H_destruc,vector_s2H_destruc);
        // for(int i=0;i<vector_z.size();i++)output_Check << vector_z[i]<< " "<< vector_s2H_destruc[i] << endl;
        fill_table_from_file(file5,vector_z_s2H_destruc_deux_iterations,vector_s2H_destruc_deux_iterations);
        fill_table_from_file(file6,vector_z_s2H_destruc_trois_iterations,vector_s2H_destruc_trois_iterations);
          fill_table_from_file(file7,vector_z_s2H_produc,vector_s2H_produc);
          // for(int i=0;i<vector_z.size();i++)output_Check << vector_z[i]<< " "<< vector_s2H_produc[i] << endl;
          fill_table_from_file(file8,vector_z_s2H_produc_deux_iterations,vector_s2H_produc_deux_iterations);
          fill_table_from_file(file9,vector_z_s2H_produc_trois_iterations,vector_s2H_produc_trois_iterations);
          for(int i=0;i<vector_z_s2H_destruc.size();i++)cout << vector_z_s2H_destruc[i]<< " "<< vector_s2H_destruc[i] << endl;


    tau=pow(10,4);double tau_moy=pow(10,6);double tau_moy2=pow(10,8);double tau_max=pow(10,10);
    //~ double T_min = 0.511*0.511/(1.6*22);
    double T_max = pow(tau_max*(2*H_r),-0.5)*T_0;
    //~ double T_max = 0.511*0.511/(2.2*22);
    //~ tau=pow(T_min/T_0,-2)/(2*H_r);tau_max=pow(T_max/T_0,-2)/(2*H_r);
    cout << " tau min = " << tau << " tau max = "<< tau_max << endl;
    double z_max = 35000;
    double h_moy = 100;
    double h_moy2 = 2*h_moy;
    double h_max = 3*h_moy;

    double Delta = (tau_moy-tau)/h_moy;
    double Delta2 = (tau_moy2-tau_moy)/h_moy;
    double Delta3 = (tau_max-tau_moy2)/h_moy;
    cout << "Delta 1 = " << Delta << " Delta 2 = " << Delta2 << " Delta 3 = " << Delta3 << endl;
  double s_2H_produc=0, s_2H_destruc=0, s_4He_destruc =0;
    double  y, dy;

    for(int h=0;h<h_max;h++)
    {
        if(h==h_moy){Delta=Delta2;tau=tau_moy;}
        if(h==h_moy2){Delta=Delta3;tau=tau_moy2;}
        if(h>1)tau+=Delta;
        T=pow(0.2*tau*(2*H_r),-0.5)*T_0;
        z = T/T_0 - 1;
        count1 = 0;count2=0;count3=0;
        Ex = 0.261121/(80*T_0*(1+z)), Ec = 0.261121/(22*T_0*(1+z));

        cout << "****************h = "<<h<<"  tau = "<<tau<<"***************"<<endl;
        // linearint(vector_z_s4He_destruc_trois_iterations,vector_s4He_destruc_trois_iterations,vector_z_s4He_destruc_trois_iterations.size(),log10(z),y,dy);
        // s_4He_destruc=y;
        linearint(vector_z_s2H_destruc_trois_iterations,vector_s2H_destruc_trois_iterations,vector_z_s2H_destruc_trois_iterations.size(),log10(z),y,dy);
        s_2H_destruc=y;
        // linearint(vector_z_s2H_produc_trois_iterations,vector_s2H_produc_trois_iterations,vector_z_s2H_produc_trois_iterations.size(),log10(z),y,dy);
        // s_2H_produc=y;
        // output_Check  << log10(z) << " " << s_4He_destruc << " " << s_2H_destruc << " " << s_2H_produc << endl;
        output_Check  << log10(z) << " " << s_4He_destruc << " " << s_2H_destruc  << endl;
          }
      output_Check.close();

    return 0;
}
