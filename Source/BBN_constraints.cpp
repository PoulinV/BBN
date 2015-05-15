#include "../Include/EM_cascade.h"
#include "../Include/Injected_spectrum.h"


void Constraints_from_destruction_only(string nuclei, double z_initial, double z_final, int z_step, struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum, struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,  long n_step){
double Y_0, Y_min;
int i_min,i_max;
double tau = pt_Particle_Physics_Model->tau_x;
double Z_x = pt_Particle_Physics_Model->Zeta_x;
double E_0 = pt_Particle_Physics_Model->E_0;
vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section;
vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_redshift;
double resultat;
double E1,E2,E3;
double f1,f2,f3;
double z1,z2,z3;
double cross_sections_1,cross_sections_2,cross_sections_3;
double z,h,log10_dz,dE;
double B,s_4He,Abundance;

  if(nuclei =="4He"){
    cout<<"I could recognize nuclei, it's 4He."<<endl;
    Y_0 = Y_4He_0;
    Y_min = Y_4He_Min;
    i_min = 15;
    i_max = 18;
  }
  else{
    cout<<"I couldn't recognize nuclei."<<endl;
    return;
    }

  log10_dz=(log10(z_initial)-log10(z_final))/(double) z_step;
  cout<<"log10dz="<<log10_dz<<endl;
  dE = (pt_Particle_Physics_Model->E_0 - E_min)/ (double) n_step;
  h = dE/2.;
  for(int j = 0; j<=z_step;j++){
    if(j==1)cout<<"I start generating cascade spectrum and convolute them with cross sections."<<endl;
    z=pow(10,log10(z_initial)-j*log10_dz);
    cout<<"redshift = " << z << endl;

    Cascade_Spectrum_Calculation_From_Function(Dirac_Spectrum_After_One_Iteration, z, pt_Particle_Physics_Model, pt_Cascade_Spectrum, 1000, 0);
    resultat = 0;


    for(int k = 0;k<=n_step;k++){
      cross_sections_1 = 0;
      cross_sections_2 = 0;
      cross_sections_3 = 0;
      if(k==0){
        E1=E_min;
        }
      else{
        E1=E3;
      }

      E2=E1 + h;
      E3=E2 + h;
      for(int i = i_min; i<=i_max;i++){
        cross_sections_1 += cross_section(E1,i);
        cross_sections_2 += cross_section(E2,i);
        cross_sections_3 += cross_section(E3,i);
      }
        linearint(pt_Cascade_Spectrum->Gamma_Energy, pt_Cascade_Spectrum->Gamma_Spectrum, pt_Cascade_Spectrum->Gamma_Energy.size(), E1, f1);
        linearint(pt_Cascade_Spectrum->Gamma_Energy, pt_Cascade_Spectrum->Gamma_Spectrum, pt_Cascade_Spectrum->Gamma_Energy.size(), E2, f2);
        linearint(pt_Cascade_Spectrum->Gamma_Energy, pt_Cascade_Spectrum->Gamma_Spectrum, pt_Cascade_Spectrum->Gamma_Energy.size(), E3, f3);
        f1 *=cross_sections_1;
        f2 *=cross_sections_2;
        f3 *=cross_sections_3;
        resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
        cout << " resultat integrale E = " << resultat << endl;

    }

    Cascade_Spectrum_Integrated_Over_Cross_Section.push_back(log10(resultat));
    Cascade_Spectrum_Integrated_Over_Cross_Section_redshift.push_back(log10(z));

  }
  cout <<"I'm done generating spectrum! I start to integrate over z."<<endl;

  log10_dz=(log10(z_initial)-log10(z_final))/(double) n_step;
  h=log10_dz/2.;
  resultat=0;
  for(int i = 0; i<n_step;i++){

    if(i==0){
      z1=z_initial;
      }
    else{
      z1=z3;
    }

    z2=pow(10,log10(z1) - h);
    z3=pow(10,log10(z2) - h);

      linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift, Cascade_Spectrum_Integrated_Over_Cross_Section, Cascade_Spectrum_Integrated_Over_Cross_Section.size(), log10(z1), f1);
      cout<<"redshift 1= " << z1 <<"interpolation = "<<f1<< endl;
      linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift, Cascade_Spectrum_Integrated_Over_Cross_Section, Cascade_Spectrum_Integrated_Over_Cross_Section.size(), log10(z2), f2);
      cout<<"redshift 2= " << z2 <<"interpolation = "<<f2<< endl;
      linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift, Cascade_Spectrum_Integrated_Over_Cross_Section, Cascade_Spectrum_Integrated_Over_Cross_Section.size(), log10(z3), f3);
      cout<<"redshift 3= " << z3 <<"interpolation = "<<f3<< endl;
      f1=pow(10,f1);
      f2=pow(10,f2);
      f3=pow(10,f3);
      f1*=exp(-1./(2*H_r*tau*(z1+1)*(z1+1)));
      f2*=exp(-1./(2*H_r*tau*(z2+1)*(z2+1)));
      f3*=exp(-1./(2*H_r*tau*(z3+1)*(z3+1)));

      resultat += pow(10,h) * (f1/3. + 4.*f2/3. + f3/3.);
      cout << " resultat integrale z = " << resultat << endl;

  }


  B=Z_x*n_y_0/(E_0*H_r*tau);
  Abundance=exp(-resultat*B)*Y_0;

  cout << "The final abundance = " << Abundance << endl;
  if(Abundance<Y_min)cout<<"your model is ruled out"<<endl;
  else cout << "your model is fine for BBN !"<<endl;

}
