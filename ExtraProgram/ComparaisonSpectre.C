#define FUNC(x,a,i,g,E_O,Z_x,z_0) ((*func)(x,a,i,g,E_O,Z_x,z_0))

#include <iostream>
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

    float T;
    double T_0 = 2.7255*0.862*pow(10,-10);

    double z,z_0,E_0;

    float tau, Ex , Ec, s,s2;
    double pi=3.14159;

    ofstream output_De("../OutputExtra/SpectreTotaliterations_high_precision.dat");
    // ofstream output_Correction("../OutputExtra/Correction4iterations_70MeV.dat");
    // ofstream output_De("../OutputExtra/SpectreTotal_30MeV.dat");
    // ofstream output_De("../OutputExtra/Check_jmax.dat");
    //
    ifstream file("../OutputExtra/Spectre4iterations_high_precision.dat");
    ifstream file2("../OutputExtra/Spectre5iterations_high_precision.dat");
    ifstream file3("../OutputExtra/SpectreTotal7Iterations.dat");
    fill_table_from_file(file,vector_E,vector_Spectre);
    fill_table_from_file(file2,vector_E_cinq_iterations,vector_Spectre_cinq_iterations);
    fill_table_from_file(file3,vector_E_sept_iterations,vector_Spectre_sept_iterations);
  double y, dy;
  // for(double E=1;E<70;E++){linearint(vector_E,vector_Spectre,vector_E.size(),E,y,dy);
  // cout<< E << "   " << y << endl;}
  // for(int i; i<vector_E_cinq_iterations.size();i++)  cout<< vector_E_cinq_iterations[i] << "  " << vector_Spectre_cinq_iterations[i] << endl;
        // ofstream output_De("results_Graphe_De_avec_production_avec_reinjection_30MeV_deux_iterations.dat");
    //~ ofstream output_He("results_Graphe_He.dat");
    //~ ofstream output_Be("results_Graphe_Be.dat");
    double M_x,E_max;
    // double y, dy;
    M_x=140;
    E_0=M_x/2;
    T = pow(10,-4);
    z = T/T_0 -1.;
    z_0 = 1*z;
    cout << " z = " << z << endl;
    //~ tau=tau_moy;
    //~ double test;
    int h=0;
    double E = 1;
    double DeltaE = 1;
    // output_De << " E " << "         spectre une iteration  " << "       spectre deux iterations "<< "       spectre trois iterations  " << endl;
    double a,b,c=0;
    while(E<=E_0){
    h++;

    /**************Print Only one spectra*************/
    // a = integrale_une_interaction_monochromatique(E, z, 1, 1., E_0, 1., z_0);
    // // output_De << E << "     "<< a ;
    // a+= integrale_deux_interactions(E, z, 1, 1., E_0, 1., z_0);
    // // output_De << "     " << a;
    // a+= integrale_trois_interactions(E, z, 1, 1., E_0, 1., z_0);
    // // output_De << "     " << a;
    // a= integrale_quatre_interactions(E, z, 1, 1., E_0, 1., z_0);
    // a= integrale_cinq_interactions(E, z, 1, 1., E_0, 1., z_0);
    // output_De  << E << "     " << a << endl;
    /*************************************************/
    /***Integrale ponderee par les sections efficaces***/
    // a = func_non_standard(E, z, 14, 1., E_0, 1., z_0);
    // output_De<< E  << "     " << a ;
    // a += func_une_interaction_monochromatique(E, z, 14, 1., E_0, 1., z_0);
    // output_De << "     "<< a ;
    // a+= func_deux_interactions(E, z, 14, 1., E_0, 1., z_0);
    // output_De << "     " << a;
    // a+= func_trois_interactions(E, z, 14, 1., E_0, 1., z_0);
    // output_De << "     " << a;
    // a+= func_quatre_interactions(E, z, 14, 1., E_0, 1., z_0);
    // output_De << "     " << a;
    // a+= func_cinq_interactions(E, z, 14, 1., E_0, 1., z_0);
    // output_De << "     " << a;
    // a+= func_six_interactions(E, z, 14, 1., E_0, 1., z_0);
    // output_De << "     " << a;
    // a+= func_sept_interactions(E, z, 14, 1., E_0, 1., z_0);
    /***************************************************/
    /****************Print all spectra****************************/
    output_De << E << "     " << integrale_spectre_standard(E, z, 1, 1., E_0, 1., z_0);
    if(E-E_0>-0.5)a = 1./(gamma_NPC(E, z, 1, 1., E_0, 1., z_0)+gamma_compton(E, z, 1, 1., E_0, 1., z_0)+gamma_phph(E, z, 1, 1., E_0, 1., z_0));
    else a = 0;
    cout << " a = " << a << endl;
    // output_De  << E << "     " << a ;
    a += integrale_une_interaction_monochromatique(E, z, 1, 1., E_0, 1., z_0);
    output_De << "     "<< a ;
    a+= integrale_deux_interactions(E, z, 1, 1., E_0, 1., z_0);
    output_De << "     " << a;
    a+= integrale_trois_interactions(E, z, 1, 1., E_0, 1., z_0);
    output_De << "     " << a;
    // output_De << "     " << a;
    a+= integrale_quatre_interactions(E, z, 1, 1., E_0, 1., z_0);
    output_De<< "     " << a ;
    b=a;
    a+= integrale_cinq_interactions(E, z, 1, 1., E_0, 1., z_0);
    output_De<< "     " << a ;
    // output_De << "     " << a ;
    a+= integrale_six_interactions(E, z, 1, 1., E_0, 1., z_0);
    output_De << "     " << a ;
    a+= integrale_sept_interactions(E, z, 1, 1., E_0, 1., z_0);
    // if(E-E_0>-0.5)a += 1./(gamma_NPC(E, z, 1, 1., E_0, 1., z_0)+gamma_compton(E, z, 1, 1., E_0, 1., z_0)+gamma_phph(E, z, 1, 1., E_0, 1., z_0));
    output_De << "     " << a ;
    // // //  b=integrale_spectre_gamma_compton(E, z, 1, 1., E_0, 1., z_0);
    // // //  output_De << E <<"     " << b ;
    c=Spectre_Gamma_compton_avec_iterations(E, z, 1, 1., E_0, 1., z_0)/(gamma_NPC(E, z, 1, 1., E_0, 1., z_0)+gamma_compton(E, z, 1, 1., E_0, 1., z_0)+gamma_phph(E, z, 1, 1., E_0, 1., z_0)) ;
    a+=c;
    output_De << "     " << a ;
    if(integrale_quatre_interactions(E, z, 1, 1., E_0, 1., z_0)>0)b=a/b;
    else b=0;
    // output_Correction << E << "     " << b << "     " << b*integrale_quatre_interactions(E, z, 1, 1., E_0, 1., z_0)<< endl;
    output_De << "     " << c << endl;
    /*******************************************************************************/


    cout << " E " << E<< endl;
    // if(h==21){DeltaE=DeltaE*10;
    //   h=1;
    // }
    E+=DeltaE;
  }
  // s=qsimp(spectre_dirac,0.00001,E_0,z,z,E_0,1., z_0);
  // cout << " int spectre_dirac = " << s << endl;
  // s=qsimp(integrale_une_interaction_monochromatique,0.00001,E_0,z,z,E_0,1., z_0);
  // cout << " int Spectre_une_interaction_monochromatique = " << s << endl;
  // s=qsimp(integrale_deux_interactions,0.00001,E_0,z,z,E_0,1., z_0);
  // cout << " Spectre_deux_interactions = " << s << endl;
  // s=qsimp(integrale_trois_interactions,0.00001,E_0,z,z,E_0,1., z_0);
  // cout << " Spectre_trois_interactions = " << s << endl;
  // s=qsimp(integrale_spectre_standard,0.00001,E_0,z,z,E_0,1., z_0);
  // cout << " integrale_spectre_standard = " << s << endl;


    output_De.close();
    // output_Correction.close();
    //~ output_He.close();
    //~ output_Be.close();
    return 0;
}
