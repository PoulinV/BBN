#include "../Include/EM_cascade.h"
#include "../Include/Injected_spectrum.h"
int l=0;
void Spectrum_and_cross_sections_convolution(struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum,
                                             struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                                             int i_min,
                                             int i_max,
                                             double &resultat,
                                             double z,
                                             long n_step,
                                             string spectrum_choice){
  double tau = pt_Particle_Physics_Model->tau_x;
  double Z_x = pt_Particle_Physics_Model->Zeta_x;
  double E_0 = pt_Particle_Physics_Model->E_0;
  double E_c ;
  double E1,E2,E3;
  double f1=0,f2=0,f3=0;
  double z1,z2,z3;
  double cross_sections_1,cross_sections_2,cross_sections_3;
  double h,log10_dz,dE;
      if(E_0>n_step)n_step *=100;
      resultat = 0;
      E_c = E_c_0/(1+z);
      dE = (E_0 - E_min)/ (double) n_step;
      h = dE/2.;

      for(int k = 0;k<=n_step;k++){
        f1=0,f2=0,f3=0;
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
          /**** If you are including 2H production, the channel 17 has a multiplicity of 2, it is now taken into account ****/
          if(i_min == 17 && i_max == 18 && i == i_min){
            cross_sections_1 += cross_section(E1,i);
            cross_sections_2 += cross_section(E2,i);
            cross_sections_3 += cross_section(E3,i);
          }
        }

          if(spectrum_choice == "Dirac"){
            if(E_0<E_c)linearint(pt_Cascade_Spectrum->Gamma_Energy, pt_Cascade_Spectrum->Gamma_Spectrum, pt_Cascade_Spectrum->Gamma_Energy.size(), E1, f1);
            else if (E_0>E_c)f1 = Universal_Spectrum(E1,z,pt_Particle_Physics_Model->E_0);
            if(E1>E_c)f1=0;
            if(E_0<E_c)linearint(pt_Cascade_Spectrum->Gamma_Energy, pt_Cascade_Spectrum->Gamma_Spectrum, pt_Cascade_Spectrum->Gamma_Energy.size(), E2, f2);
            else if (E_0>E_c)f2 = Universal_Spectrum(E2,z,pt_Particle_Physics_Model->E_0);
            if(E2>E_c)f2=0;
            if(E_0>E_c && E3<E_c)linearint(pt_Cascade_Spectrum->Gamma_Energy, pt_Cascade_Spectrum->Gamma_Spectrum, pt_Cascade_Spectrum->Gamma_Energy.size(), E3, f3);
            else if (E_0>E_c)f3 = Universal_Spectrum(E3,z,pt_Particle_Physics_Model->E_0);
            if(E3>E_c)f3=0;
          }
          else if(spectrum_choice == "universal"){
          if(E1 <= E_c)f1 = Universal_Spectrum(E1,z,pt_Particle_Physics_Model->E_0);
          if(E2 <= E_c)f2 = Universal_Spectrum(E2,z,pt_Particle_Physics_Model->E_0);
          if(E3 <= E_c)f3 = Universal_Spectrum(E3,z,pt_Particle_Physics_Model->E_0);
          }
          // if(z==462770 && l==0){
          //   cout << z << " " << E1 << " " << f1/(gamma_NPC(E1,z)+gamma_compton(E1,z)+gamma_phph(E1,z)) << " E_c = " << E_c << endl;
          //   l++;}
          f1 *=cross_sections_1;
          f1 /=(gamma_NPC(E1,z)+gamma_compton(E1,z)+gamma_phph(E1,z));
          f2 *=cross_sections_2;
          f2 /=(gamma_NPC(E2,z)+gamma_compton(E2,z)+gamma_phph(E2,z));
          f3 *=cross_sections_3;
          f3 /=(gamma_NPC(E3,z)+gamma_compton(E3,z)+gamma_phph(E3,z));
          resultat += h * (f1/3. + 4.*f2/3. + f3/3.);

          if(spectrum_choice == "Dirac" && k==n_step && E_c>pt_Particle_Physics_Model->E_0){
            // cout << "works"<<endl;

            for(int i = i_min; i<=i_max;i++){
              resultat+=cross_section(pt_Particle_Physics_Model->E_0,i)/(gamma_NPC(pt_Particle_Physics_Model->E_0,z)+gamma_compton(pt_Particle_Physics_Model->E_0,z)+gamma_phph(pt_Particle_Physics_Model->E_0,z));
              if(i_min == 17 && i_max == 18 && i == i_min)resultat+=cross_section(pt_Particle_Physics_Model->E_0,i)/(gamma_NPC(pt_Particle_Physics_Model->E_0,z)+gamma_compton(pt_Particle_Physics_Model->E_0,z)+gamma_phph(pt_Particle_Physics_Model->E_0,z));
            }
          }

      }


    // cout << "resultat = "<<resultat << endl;
}
//
// void Check_model_from_destruction_only(string nuclei,
//                                        struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum,
//                                        struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
//                                        double &Abundance,
//                                        double z_initial,
//                                        double z_final,
//                                        int z_step,
//                                        long n_step,
//                                        int iterations,
//                                        string spectrum_choice,
//                                        string spectrum_mode){
//
// double Y_0, Y_min, Y_max;
// int i_min,i_max, k_min, k_max;
// double resultat;
// double tau_x = pt_Particle_Physics_Model->tau_x;
// double zeta_x = pt_Particle_Physics_Model->Zeta_x;
// double E_0 = pt_Particle_Physics_Model->E_0;
// double E_c;
// double B;
// double f1, f2, f3, z1, z2, z3, h, dz, z;
// double log10_dz = (log10(z_initial)-log10(z_final))/(double) z_step;
// double E1, dE;
// vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei;
// vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei;
//
//
//
//
// Check_nuclei(nuclei, i_min, i_max, k_min, k_max, Y_min, Y_max, Y_0);
// for(int j = 0; j<=z_step;j++){
//
//   z=pow(10,log10(z_initial)-j*log10_dz);
//   E_c = E_c_0/(1+z);
//   if(verbose>1)cout<<"redshift = " << z << " still " << z_step-j << " to go " << endl;
//
//   if(spectrum_choice == "Dirac"){
//       if(spectrum_mode=="writing" || spectrum_mode == "nothing")Cascade_Spectrum_Calculation(spectrum_choice, z, pt_Particle_Physics_Model, pt_Cascade_Spectrum, n_step, iterations, spectrum_mode);
//       else if(spectrum_mode=="reading"){
//
//         if(E_c <= pt_Particle_Physics_Model->E_0 ){
//           Cascade_Spectrum_Calculation("universal", z, pt_Particle_Physics_Model, pt_Cascade_Spectrum, n_step, iterations, spectrum_mode);
//           }
//
//       else Cascade_Spectrum_Reading_From_File(pt_Particle_Physics_Model,pt_Cascade_Spectrum, z, iterations);
//       }
//    }
//    else if(spectrum_choice == "universal"){
//      Cascade_Spectrum_Calculation(spectrum_choice, z, pt_Particle_Physics_Model, pt_Cascade_Spectrum, n_step, iterations, spectrum_mode);
//    }
//
//   Spectrum_and_cross_sections_convolution(pt_Cascade_Spectrum, pt_Particle_Physics_Model, i_min, i_max, resultat, z, n_step, spectrum_choice);
//
//   Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.push_back(log10(resultat));
//   Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei.push_back(log10(z));
//  }
//
//  if(verbose>1)cout <<"I'm done generating spectrum! I start to integrate over z."<<endl;
//
//  dz=(z_initial-z_final)/(double) n_step;
//  h=dz/2;
//  resultat=0;
//
//  for(int i = 0; i<=n_step;i++){
//
//    if(i==0){
//      z1=z_initial;
//      }
//    else{
//      z1=z3;
//    }
//
//
//    z2=z1-h;
//    z3=z2-h;
//
//      linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.size(), log10(z1), f1);
//      // cout<<"redshift 1= " << z1 <<"interpolation = "<<f1<< endl;
//      if(f1<0)f1=0;
//      linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.size(), log10(z2), f2);
//      // cout<<"redshift 2= " << z2 <<"interpolation = "<<f2<< endl;
//      if(f2<0)f2=0;
//
//      linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.size(), log10(z3), f3);
//      // cout<<"redshift 3= " << z3 <<"interpolation = "<<f3<< endl;
//      if(f3<0)f3=0;
//
//      f1=pow(10,f1);
//      f2=pow(10,f2);
//      f3=pow(10,f3);
//      f1*=exp(-1./(2*H_r*tau_x*(z1+1)*(z1+1)));
//      f2*=exp(-1./(2*H_r*tau_x*(z2+1)*(z2+1)));
//      f3*=exp(-1./(2*H_r*tau_x*(z3+1)*(z3+1)));
//
//      resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
//      // cout << " resultat integrale z = " << resultat << endl;
//
//  }
//
//
//   B=zeta_x*n_y_0/(E_0*H_r*tau_x);
//   if(verbose>1)cout << "zeta_x = " << zeta_x << "resultat = " << resultat << " B = " << B << endl;
//   Abundance=exp(-resultat*B);
//   Abundance*=Y_0;
//
//     cout << "The final abundance = " << Abundance << endl;
//     if(Abundance<Y_min){
//       cout<<"your model is ruled out"<<endl;
//     }
//     else cout << "your model is fine for BBN !"<<endl;
//
//   Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei.clear();
//   Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.clear();
//
// }


void Compute_Constraints_from_destruction_only(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                                               struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                                               struct Structure_Scan_Parameters * pt_Scan_Parameters,
                                               struct Structure_Output_Options * pt_Output_Options){

  struct Structure_Gamma_Spectrum Cascade_Spectrum;
  vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei;
  vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei;

  /******** This step attributes locally the values of precision and scan parameters ********/
  double tau_min, tau_max, tau_step, zeta_min, zeta_max, zeta_step, z_step, n_step, iterations;
  tau_min = pt_Scan_Parameters->tau_min;
  tau_max = pt_Scan_Parameters->tau_max;
  tau_step = pt_Scan_Parameters->tau_step;
  zeta_min = pt_Scan_Parameters->zeta_min;
  zeta_max = pt_Scan_Parameters->zeta_max;
  zeta_step = pt_Scan_Parameters->zeta_step;
  n_step = pt_Spectrum_and_Precision_Parameters->n_step;
  z_step = pt_Spectrum_and_Precision_Parameters->z_step;
  iterations = pt_Spectrum_and_Precision_Parameters->iterations;
  /******************************************************************************************/

  double z_initial =  (5*(pow(2*H_r*tau_min,-0.5)-1)), z_final =  35000, z;
  int i_min, i_max, k_min, k_max;
  double zeta_x, tau_x;
  double log10_dtau = (log10(tau_max)-log10(tau_min))/(double) tau_step;
  double log10_dZ = (log10(zeta_max)-log10(zeta_min))/(double) zeta_step;
  double log10_dz = (log10(z_initial)-log10(z_final))/(double) z_step;
  double Y_0, Y_min, Y_max;
  double Abundance, resultat, B, E_c;
  double f1, f2, f3, z1, z2, z3, h, dz;
  double E1, dE;

    ostringstream os;
    string name;
    if(pt_Spectrum_and_Precision_Parameters->spectrum_choice=="universal")os << "Output/Results_destruc_only_"<< pt_Scan_Parameters->nuclei << "_m"<<pt_Particle_Physics_Model->M_x<<"MeV_Universal_Spectrum.dat";
    else if(pt_Spectrum_and_Precision_Parameters->spectrum_choice=="Dirac")os << "Output/Results_destruc_only_"<< pt_Scan_Parameters->nuclei << "_m"<<pt_Particle_Physics_Model->M_x<<"MeV_Dirac_Spectrum_"<<pt_Spectrum_and_Precision_Parameters->iterations<<"iterations.dat";
    name = os.str();
    ofstream file(name);
    if(file){
      if(verbose>1)cout<<" I could open file " << name << endl;
    }
    else{
      if(verbose>1)cout<<" I couldn't open file, i'm shutting down." << endl;
      return;
    }


 Check_nuclei(pt_Scan_Parameters->nuclei, i_min, i_max, k_min, k_max, Y_min, Y_max, Y_0);

 dE= (pt_Particle_Physics_Model->E_0 - E_min)/ (double) (Gamma_Table_Size-1);
 cout << " z initial = " << z_initial << " z final = " << z_final << endl;
 for(int j = 0; j<=z_step;j++){

   z=pow(10,log10(z_initial)-j*log10_dz);
   E_c = E_c_0/(1+z);
   if(verbose>1)cout<<"redshift = " << z << " still " << z_step-j << " to go " << endl;

   if(pt_Spectrum_and_Precision_Parameters->spectrum_choice == "Dirac"){
       if(pt_Spectrum_and_Precision_Parameters->spectrum_mode=="writing" || pt_Spectrum_and_Precision_Parameters->spectrum_mode == "nothing"){
         Cascade_Spectrum_Calculation(z,
                                    pt_Particle_Physics_Model,
                                    &Cascade_Spectrum,
                                    pt_Spectrum_and_Precision_Parameters);
       }
       else if(pt_Spectrum_and_Precision_Parameters->spectrum_mode=="reading"){
         if(E_c <= pt_Particle_Physics_Model->E_0 ){
           Cascade_Spectrum.Gamma_Energy.resize(Gamma_Table_Size);
           Cascade_Spectrum.Gamma_Spectrum.resize(Gamma_Table_Size);
           for(int i=0;i<Gamma_Table_Size;i++){
             Cascade_Spectrum.Gamma_Energy[i]=E_min+i*dE;
             Cascade_Spectrum.Gamma_Spectrum[i]=Universal_Spectrum(E_min+i*dE,z,pt_Particle_Physics_Model->E_0);

            }
         }


      Cascade_Spectrum_Reading_From_File(pt_Particle_Physics_Model,
                                               &Cascade_Spectrum,
                                               z,
                                               pt_Spectrum_and_Precision_Parameters->iterations);
     }

    }
    // else if(spectrum_choice == "universal"){
    //   Cascade_Spectrum_Calculation(spectrum_choice, z, &Particle_Physics_Model, &Cascade_Spectrum, n_step, iterations, spectrum_mode);
    // }
    // for(int i = 0 ; i < Gamma_Table_Size ; i ++){
    //   cout << z << "  " << Cascade_Spectrum.Gamma_Energy[i] << "  " << Cascade_Spectrum.Gamma_Spectrum[i]/(gamma_NPC(Cascade_Spectrum.Gamma_Energy[i],z)+gamma_compton(Cascade_Spectrum.Gamma_Energy[i],z)+gamma_phph(Cascade_Spectrum.Gamma_Energy[i],z))<<endl;
    // }
   Spectrum_and_cross_sections_convolution(&Cascade_Spectrum,
                                           pt_Particle_Physics_Model,
                                           i_min,
                                           i_max,
                                           resultat,
                                           z,
                                           n_step,
                                           pt_Spectrum_and_Precision_Parameters->spectrum_choice);

   Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.push_back(log10(resultat));
   Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei.push_back(log10(z));
   Cascade_Spectrum.Gamma_Energy.clear();
   Cascade_Spectrum.Gamma_Spectrum.clear();
  }
  for(int dtau = 0 ; dtau <= tau_step ; dtau++){
    tau_x = pow(10,log10(tau_min)+log10_dtau*dtau);
    if((dtau==tau_step) && (tau_x!=tau_max))cout<<"erreur : probleme de pas logarithmique en tau"<<endl;

    if(verbose>0)cout << "Current lifetime analysed : " << tau_x << endl;
    Fill_Structure_Particle_Physics_Model(pt_Particle_Physics_Model->M_x, zeta_min, tau_x, pt_Particle_Physics_Model); // MANDATORY STEP

    z_initial = 5*pt_Particle_Physics_Model->z_x;
    // z_final = z_min;
    dz=(z_initial-z_final)/(double) n_step;
    h=dz/2;
    resultat= 0;
    for(int i = 0; i<=n_step;i++){

      if(i==0){
        z1=z_initial;
        }
      else{
        z1=z3;
      }

      z2=z1-h;
      z3=z2-h;
      // for(int l = 0;l<Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.size();l++)cout << Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei[l] << "   " << Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei[l] <<  endl;

        linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.size(), log10(z1), f1);
        // cout<<"redshift 1= " << z1 <<"interpolation = "<<f1<< endl;
        if(f1<0)f1=0;
        linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.size(), log10(z2), f2);
        // cout<<"redshift 2= " << z2 <<"interpolation = "<<f2<< endl;
        if(f2<0)f2=0;

        linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.size(), log10(z3), f3);
        // cout<<"redshift 3= " << z3 <<"interpolation = "<<f3<< endl;
        if(f3<0)f3=0;

        f1=pow(10,f1);
        f2=pow(10,f2);
        f3=pow(10,f3);
        f1*=exp(-1./(2*H_r*tau_x*(z1+1)*(z1+1)));
        f2*=exp(-1./(2*H_r*tau_x*(z2+1)*(z2+1)));
        f3*=exp(-1./(2*H_r*tau_x*(z3+1)*(z3+1)));

        // resultat += pow(10,h) * (f1/3. + 4.*f2/3. + f3/3.);
        resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
      }
      for(int dZ = 0 ; dZ <= zeta_step ; dZ++){
      zeta_x = pow(10,log10(zeta_min)+log10_dZ*dZ);
      // cout << "tau_x = "<< tau_x<<" zeta_x = " << zeta_x << "resultat = " << resultat << "E0" << pt_Particle_Physics_Model->E_0<<endl;
      // B=zeta_x*n_y_0/(H_r*tau_x);
      B=zeta_x*n_y_0/(pt_Particle_Physics_Model->E_0*H_r*tau_x);
      Abundance=exp(-resultat*B);
      Abundance*=Y_0;
      cout << " B = " << B << "Abundance = "<<Abundance <<"Y0 = " << Y_0 << "Y_min = " << Y_min << endl;

      if(pt_Output_Options->print_result == "yes"){
        if(Abundance<Y_min){
          cout << "The final abundance = " << Abundance << endl;
          cout<<"your model is ruled out"<<endl;
        }
      }
      if(Abundance < Y_min){
        if(tau_x==tau_max)cout << "Exporting results in " << name << "."<<endl;

        file << tau_x << "  " << zeta_x << "  " << Abundance << endl;
        break;
      }
    }
  }

    Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei.clear();
    Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.clear();
}

//
//
// void Check_model_from_destruction_and_production(string nuclei,
//                                                  struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum,
//                                                  struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
//                                                  double &Abundance,
//                                                  struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters
//                                                  ){
//
//   double z,dz,h;
//   double z1, z2, z3;
//   vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei;
//   vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He;
//   vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei;
//   vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei;
//   vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He;
//   vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei;
//   vector<double> Integration_over_z_source_term_redshift;
//   vector<double> Integration_over_z_source_term;
//   double z_initial = pt_Spectrum_and_Precision_Parameters->z_initial,
//   z_final = pt_Spectrum_and_Precision_Parameters->z_final,
//   z_step = pt_Spectrum_and_Precision_Parameters->z_step,
//   n_step = pt_Spectrum_and_Precision_Parameters->n_step,
//   iterations = pt_Spectrum_and_Precision_Parameters->iterations;
//   int i_min, i_max;
//   int j_min, j_max;
//   int k_min, k_max;
//   int l_min, l_max;
//   double zeta_x=pt_Particle_Physics_Model->Zeta_x, tau_x=pt_Particle_Physics_Model->tau_x;
//   double E_c ;
//   double log10_dz=(log10(z_initial)-log10(z_final))/(double) z_step;
//   double Y_0, Y_min, Y_max;
//   double K_0, K_min, K_max;
//   double resultat;
//   double resultat_destruc_nuclei, resultat_source_term ;
//   double B;
//   double f1, f2, f3, g1, g2, g3;
//
//
//
//     // Fill_Structure_Particle_Physics_Model(M_x, zeta_x, tau_x, &Particle_Physics_Model); // MANDATORY STEP
//
//     /*****************************************************************************************************************************************************************************/
//     /**** First step : Compute destruction of the nuclei for several time and store it in a table Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei. ****/
//     /*********************************************************** It will be needed to integrate over in step four.****************************************************************/
//     /*****************************************************************************************************************************************************************************/
//     /******************************************************************************************************************************************************************************/
//     /**** Second step : Compute destruction of the 4He for several time and store it in a table resultat_destruc. It is the father nuclei that will enter the production term. ****/
//     /************************************************************* The evolution of its abundance needs to be followed. ***********************************************************/
//     /******************************************************************************************************************************************************************************/
//     /********************************************************************************************************************************************************************************/
//     /**** Third step : Compute production of the nuclei through 4He destruction. The difference between this step and the previous one is that only the channels destructing 4He ****/
//     /********************************************************************* and producing the nuclei are now open. *******************************************************************/
//     /********************************************************************************************************************************************************************************/
//
//     Check_nuclei(nuclei, i_min, i_max, k_min, k_max, Y_min, Y_max, Y_0);
//     Check_nuclei("4He", j_min, j_max, l_min, l_max, K_min, K_max, K_0);
//     if(verbose>1)cout << " z initial = " << z_initial << " z final = " << z_final << endl;
//     for(int j = 0; j<=z_step;j++){
//
//       z=pow(10,log10(z_initial)-j*log10_dz);
//       E_c = E_c_0/(1+z);
//       if(verbose>1)cout<<"redshift = " << z << " still " << z_step-j << " to go " << endl;
//
//       if(spectrum_choice == "Dirac"){
//           if(spectrum_mode=="writing" || spectrum_mode == "nothing")Cascade_Spectrum_Calculation(z,
//                                                                                                  pt_Particle_Physics_Model,
//                                                                                                  &Cascade_Spectrum,
//                                                                                                  pt_Spectrum_and_Precision_Parameters);
//           else if(spectrum_mode=="reading"){
//
//             if(E_c <= Particle_Physics_Model->E_0 ){
//               for(int i=0;i<Gamma_Table_Size;i++){
//                 pt_Cascade_Spectrum->Gamma_Energy[i]=E_min+i*dE;
//                 pt_Cascade_Spectrum->Gamma_Spectrum[i]=Universal_Spectrum(E_min+i*dE,z,pt_Particle_Physics_Model->E_0);
//               }
//             }
//
//           else Cascade_Spectrum_Reading_From_File(pt_Particle_Physics_Model,
//                                                   &Cascade_Spectrum,
//                                                   z,
//                                                   pt_Spectrum_and_Precision_Parameters->iterations);
//           }
//        }
//        else if(spectrum_choice == "universal"){
//          Cascade_Spectrum_Calculation(z,
//                                       pt_Particle_Physics_Model,
//                                       &Cascade_Spectrum,
//                                       pt_Spectrum_and_Precision_Parameters);
//        }
//
//       Spectrum_and_cross_sections_convolution(pt_Cascade_Spectrum, pt_Particle_Physics_Model, i_min, i_max, resultat, z, n_step, pt_Spectrum_and_Precision_Parameters->spectrum_mode);
//
//       Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.push_back(log10(resultat));
//       Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei.push_back(log10(z));
//       // cout << z << "  " << resultat  << endl;
//       // for(int l = 0;l<Gamma_Table_Size;l++)cout << pt_Cascade_Spectrum->Gamma_Energy[l] << "   " << pt_Cascade_Spectrum->Gamma_Spectrum[l] <<  endl;
//
//
//       Spectrum_and_cross_sections_convolution(pt_Cascade_Spectrum, pt_Particle_Physics_Model, j_min, j_max, resultat, z, n_step, pt_Spectrum_and_Precision_Parameters->spectrum_mode);
//
//       Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He.push_back(log10(resultat));
//       Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He.push_back(log10(z));
//       // cout << z << "  " << resultat  << endl;
//       // for(int l = 0;l<Gamma_Table_Size;l++)cout << pt_Cascade_Spectrum->Gamma_Energy[l] << "   " << pt_Cascade_Spectrum->Gamma_Spectrum[l] <<  endl;
//
//       Spectrum_and_cross_sections_convolution(pt_Cascade_Spectrum, pt_Particle_Physics_Model, k_min, k_max, resultat, z, n_step, pt_Spectrum_and_Precision_Parameters->spectrum_mode);
//
//       Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei.push_back(log10(resultat));
//       Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei.push_back(log10(z));
//       // cout << z << "  " << resultat << endl;
//       // for(int l = 0;l<Gamma_Table_Size;l++)cout << pt_Cascade_Spectrum->Gamma_Energy[l] << "   " << pt_Cascade_Spectrum->Gamma_Spectrum[l] <<  endl;
//
//     }
//
//   /********************************************************************************************************************************************************************************/
//   /****************** Fourth step : Every integral that is independant of zeta_x is done one and for all. The source term involved a double integration over z. *******************/
//   /************************************************* To do so, every step in z is stored in the table Integration_over_z_source_term. *********************************************/
//   /********************************************************************************************************************************************************************************/
//
// if(verbose>1)cout <<"I'm done convoluting the spectrum with cross sections! I start to integrate over z."<<endl;
//
// // log10_dz=(log10(z_initial)-log10(z_final))/(double) n_step;
// dz=(z_initial-z_final)/(double) n_step;
// h=dz/2;
// // h=log10_dz/log10(2.);
// resultat_destruc_nuclei= 0;
// resultat_source_term = 0;
// for(int i = 0; i<=n_step;i++){
//
//   if(i==0){
//     z1=z_initial;
//     }
//   else{
//     z1=z3;
//   }
//
//   // z2=pow(10,log10(z1) - h);
//   // z3=pow(10,log10(z2) - h);
//   z2=z1-h;
//   z3=z2-h;
//   // for(int l = 0;l<Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.size();l++)cout << Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei[l] << "   " << Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei[l] <<  endl;
//
//     //First integral : Destruction of the nuclei
//     linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.size(), log10(z1), f1);
//     // cout<<"redshift 1= " << z1 <<"interpolation = "<<f1<< endl;
//     if(f1<0)f1=0;
//     linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.size(), log10(z2), f2);
//     // cout<<"redshift 2= " << z2 <<"interpolation = "<<f2<< endl;
//     if(f2<0)f2=0;
//
//     linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.size(), log10(z3), f3);
//     // cout<<"redshift 3= " << z3 <<"interpolation = "<<f3<< endl;
//     if(f3<0)f3=0;
//
//     f1=pow(10,f1);
//     f2=pow(10,f2);
//     f3=pow(10,f3);
//     f1*=exp(-1./(2*H_r*tau_x*(z1+1)*(z1+1)));
//     f2*=exp(-1./(2*H_r*tau_x*(z2+1)*(z2+1)));
//     f3*=exp(-1./(2*H_r*tau_x*(z3+1)*(z3+1)));
//
//     // resultat += pow(10,h) * (f1/3. + 4.*f2/3. + f3/3.);
//     resultat_destruc_nuclei += h * (f1/3. + 4.*f2/3. + f3/3.);
//     // cout << " resultat_destruc_nuclei_1 integrale z = " << resultat_destruc_nuclei_1 << endl;
//
//
//     //Second integral : Destruction of the nuclei
//     linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He.size(), log10(z1), g1);
//     // cout<<"redshift 1= " << z1 <<"interpolation = "<<g1<< endl;
//     if(g1<0)g1=0;
//     linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He.size(), log10(z2), g2);
//     // cout<<"redshift 2= " << z2 <<"interpolation = "<<g2<< endl;
//     if(g2<0)g2=0;
//
//     linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He.size(), log10(z3), g3);
//     // cout<<"redshift 3= " << z3 <<"interpolation = "<<g3<< endl;
//     if(g3<0)g3=0;
//
//     g1=pow(10,g1);
//     g2=pow(10,g2);
//     g3=pow(10,g3);
//     g1*=exp(-1./(2*H_r*tau_x*(z1+1)*(z1+1)));
//     g2*=exp(-1./(2*H_r*tau_x*(z2+1)*(z2+1)));
//     g3*=exp(-1./(2*H_r*tau_x*(z3+1)*(z3+1)));
//
//
//     // resultat += pow(10,h) * (f1/3. + 4.*f2/3. + f3/3.);
//     resultat_source_term += h * ((g1-f1)/3. + 4.*(g2-f2)/3. + (g3-f3)/3.);
//     // cout << " resultat_source_term = " << resultat_source_term << " z = "<< z3 << " g3 = " << g3 << " f3 = " << f3 << " exp(-1./(2*H_r*tau_x*(z3+1)*(z3+1))) = " << exp(-1./(2*H_r*tau_x*(z3+1)*(z3+1))) << endl;
//     Integration_over_z_source_term_redshift.push_back(z3);
//     Integration_over_z_source_term.push_back(resultat_source_term);
//
//
//
//   }
//   // for(int l = 0;l<Integration_over_z_source_term_redshift.size();l++)cout << Integration_over_z_source_term_redshift[l] << "   " << Integration_over_z_source_term[l] <<  endl;
//
//     /********************************************************************************************************************************************************************************/
//     /********************* Last step : We can now perform integrals that are not independant of zeta. Since we stored previous results in tables or variable, ***********************/
//     /********************************************************* there is only one integral per zeta_x left. It is big time gain. *****************************************************/
//     /********************************************************************************************************************************************************************************/
//
//     // for(int dZ = 0 ; dZ <= zeta_step ; dZ++){
//     //   zeta_x = pow(10,log10(zeta_min)+log10_dZ*dZ);
//       B=zeta_x*n_y_0/(pt_Particle_Physics_Model->E_0*H_r*tau_x);
//       if(verbose>1)cout << "zeta_x = " << zeta_x << "resultat = " << resultat_destruc_nuclei << " B = " << B << endl;
//
//
//       if(verbose>1)cout <<"It's the last step : I start to integrate over z for the value of zeta."<<endl;
//
//       // log10_dz=(log10(z_initial)-log10(z_final))/(double) n_step;
//
//       // h=log10_dz/log10(2.);
//       resultat=0;
//       for(int i = 0; i<=n_step;i++){
//
//         if(i==0){
//           z1=z_initial;
//           }
//         else{
//           z1=z3;
//         }
//
//         // z2=pow(10,log10(z1) - h);
//         // z3=pow(10,log10(z2) - h);
//         z2=z1-h;
//         z3=z2-h;
//
//           linearint(Integration_over_z_source_term_redshift, Integration_over_z_source_term, Integration_over_z_source_term_redshift.size(), z1, g1);
//           // cout<<"redshift 1= " << z1 <<"interpolation = "<<g1<< endl;
//           // if(g1<0)g1=0;
//           linearint(Integration_over_z_source_term_redshift, Integration_over_z_source_term, Integration_over_z_source_term_redshift.size(), z2, g2);
//           // cout<<"redshift 2= " << z2 <<"interpolation = "<<g2<< endl;
//           // if(g2<0)g2=0;
//           linearint(Integration_over_z_source_term_redshift, Integration_over_z_source_term, Integration_over_z_source_term_redshift.size(), z3, g3);
//           // cout<<"redshift 3= " << z3 <<"interpolation = "<<g3<< endl;
//           // if(g3<0)g3=0;
//
//           g1=exp(-B*g1);
//           g2=exp(-B*g2);
//           g3=exp(-B*g3);
//           // cout<<"redshift 1= " << z1 <<"interpolation = "<<g1<< endl;
//
//           // cout<<"redshift 2= " << z2 <<"interpolation = "<<g2<< endl;
//
//           // cout<<"redshift 3= " << z3 <<"interpolation = "<<g3<< endl;
//
//
//           linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei.size(), log10(z1), f1);
//           // cout<<"redshift 1= " << z1 <<"interpolation = "<<f1<< endl;
//           if(f1<0)f1=0;
//           linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei.size(), log10(z2), f2);
//           // cout<<"redshift 2= " << z2 <<"interpolation = "<<f2<< endl;
//           if(f2<0)f2=0;
//
//           linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei.size(), log10(z3), f3);
//           // cout<<"redshift 3= " << z3 <<"interpolation = "<<f3<< endl;
//           if(f3<0)f3=0;
//           f1=pow(10,f1);
//           f2=pow(10,f2);
//           f3=pow(10,f3);
//           f1*=exp(-1./(2*H_r*tau_x*(z1+1)*(z1+1)));
//           f2*=exp(-1./(2*H_r*tau_x*(z2+1)*(z2+1)));
//           f3*=exp(-1./(2*H_r*tau_x*(z3+1)*(z3+1)));
//           f1*=g1*K_0;
//           f2*=g2*K_0;
//           f3*=g3*K_0;
//
//           // resultat += pow(10,h) * (f1/3. + 4.*f2/3. + f3/3.);
//           resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
//           // cout << " resultat integrale z = " << resultat << endl;
//
//       }
//       if(verbose>1)cout << " resultat source term = " << resultat << " Y_0 = " << Y_0 << endl;
//       Abundance=exp(-resultat_destruc_nuclei*B)*(Y_0+B*resultat);
//       // Abundance=exp(-resultat_destruc_nuclei*B)*(Y_0);
//
//
//         if(Abundance<Y_min){
//           cout<<"Your model is ruled out !"<<endl;
//           cout <<"You have overdestructed "<<nuclei<<" : Y_final = " << Abundance << " < Y_min = " << Y_min << "." << endl;
//         }
//         else if(Abundance>Y_max){
//           cout<<"Your model is ruled out !"<<endl;
//           cout <<"You have overproduced "<<nuclei<<" : Y_final = " << Abundance << " > Y_max = " << Y_max << "." << endl;
//         }
//         else {
//           cout << "The final abundance = " << Abundance << "." << endl;
//           cout << "Your model is fine for BBN !"<<endl;
//         }
//
//
//     // }
//     Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.clear();
//     Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He.clear();
//     Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei.clear();
//     Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei.clear();
//     Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He.clear();
//     Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei.clear();
//     Integration_over_z_source_term_redshift.clear();
//     Integration_over_z_source_term.clear();
// }
void Compute_constraints_from_destruction_and_production(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                                                         struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                                                         struct Structure_Scan_Parameters * pt_Scan_Parameters,
                                                         struct Structure_Output_Options * pt_Output_Options){

  double z,dz,h;
  double z1, z2, z3;
  vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei;
  vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He;
  vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei;
  vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei;
  vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He;
  vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei;
  vector<double> Integration_over_z_source_term_redshift;
  vector<double> Integration_over_z_source_term;

  struct Structure_Gamma_Spectrum Cascade_Spectrum;
  /******** This step attributes locally the values of precision and scan parameters ********/
  double tau_min, tau_max, tau_step, zeta_min, zeta_max, zeta_step, z_step, n_step, iterations;
  tau_min = pt_Scan_Parameters->tau_min;
  tau_max = pt_Scan_Parameters->tau_max;
  tau_step = pt_Scan_Parameters->tau_step;
  zeta_min = pt_Scan_Parameters->zeta_min;
  zeta_max = pt_Scan_Parameters->zeta_max;
  zeta_step = pt_Scan_Parameters->zeta_step;
  n_step = pt_Spectrum_and_Precision_Parameters->n_step;
  z_step = pt_Spectrum_and_Precision_Parameters->z_step;
  iterations = pt_Spectrum_and_Precision_Parameters->iterations;
  /******************************************************************************************/
  int i_min, i_max;
  int j_min, j_max;
  int k_min, k_max;
  int l_min, l_max;
  double zeta_x, tau_x;
  // double M_x = pt_Particle_Physics_Model->M_x;
  double E1, dE;
  double z_initial =  (5*(pow(2*H_r*tau_min,-0.5)-1)), z_final =  35000;
  double log10_dz=(log10(z_initial)-log10(z_final))/(double) z_step;
  double log10_dtau = (log10(tau_max)-log10(tau_min))/(double) tau_step;
  double log10_dZ = (log10(zeta_max)-log10(zeta_min))/(double) zeta_step;
  double Y_0, Y_min, Y_max;
  double K_0, K_min, K_max;
  double Abundance, resultat, B, E_c;
  double resultat_destruc_nuclei, resultat_source_term ;
  double f1, f2, f3, g1, g2, g3;







    /*****************************************************************************************************************************************************************************/
    /**** First step : Compute destruction of the nuclei for several time and store it in a table Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei. ****/
    /*********************************************************** It will be needed to integrate over in step four.****************************************************************/
    /*****************************************************************************************************************************************************************************/
    /******************************************************************************************************************************************************************************/
    /**** Second step : Compute destruction of the 4He for several time and store it in a table resultat_destruc. It is the father nuclei that will enter the production term. ****/
    /************************************************************* The evolution of its abundance needs to be followed. ***********************************************************/
    /******************************************************************************************************************************************************************************/
    /********************************************************************************************************************************************************************************/
    /**** Third step : Compute production of the nuclei through 4He destruction. The difference between this step and the previous one is that only the channels destructing 4He ****/
    /********************************************************************* and producing the nuclei are now open. *******************************************************************/
    /********************************************************************************************************************************************************************************/
    ostringstream os;
    string name;
    os << "Output/Results_destruc_and_production_"<< pt_Scan_Parameters->nuclei << "_m"<<pt_Particle_Physics_Model->M_x<<"MeV.dat";
    name = os.str();
    ofstream file(name);
    if(file){
      if(verbose>1)cout<<" I could open file " << name << endl;
    }
    else{
      if(verbose>1)cout<<" I couldn't open file, i'm shutting down." << endl;
      return;
    }
    Check_nuclei(pt_Scan_Parameters->nuclei, i_min, i_max, k_min, k_max, Y_min, Y_max, Y_0);
    Check_nuclei("4He", j_min, j_max, l_min, l_max, K_min, K_max, K_0);
    if(verbose>1)cout << " z initial = " << z_initial << " z final = " << z_final << endl;
    for(int j = 0; j<=z_step;j++){

      z=pow(10,log10(z_initial)-j*log10_dz);
      E_c = E_c_0/(1+z);
      if(verbose>1)cout<<"redshift = " << z << " still " << z_step-j << " to go " << endl;


      if(pt_Spectrum_and_Precision_Parameters->spectrum_choice == "Dirac"){
          if(pt_Spectrum_and_Precision_Parameters->spectrum_mode=="writing" || pt_Spectrum_and_Precision_Parameters->spectrum_mode == "nothing")
          Cascade_Spectrum_Calculation(z,
                                      pt_Particle_Physics_Model,
                                      &Cascade_Spectrum,
                                      pt_Spectrum_and_Precision_Parameters);
          else if(pt_Spectrum_and_Precision_Parameters->spectrum_mode=="reading"){

            if(E_c <= pt_Particle_Physics_Model->E_0 ){
              Cascade_Spectrum.Gamma_Energy.resize(Gamma_Table_Size);
              Cascade_Spectrum.Gamma_Spectrum.resize(Gamma_Table_Size);
              for(int i=0;i<Gamma_Table_Size;i++){
                Cascade_Spectrum.Gamma_Energy[i]=E_min+i*dE;
                Cascade_Spectrum.Gamma_Spectrum[i]=Universal_Spectrum(E_min+i*dE,z,pt_Particle_Physics_Model->E_0);
              }
            }

          else Cascade_Spectrum_Reading_From_File(pt_Particle_Physics_Model,
                                                  &Cascade_Spectrum,
                                                  z,
                                                  pt_Spectrum_and_Precision_Parameters->iterations);
          }
       }
        // else if(spectrum_choice == "universal"){
        //   Cascade_Spectrum_Calculation(spectrum_choice, z, &Particle_Physics_Model, &Cascade_Spectrum, n_step, iterations, spectrum_mode);
        // }

       Spectrum_and_cross_sections_convolution(&Cascade_Spectrum, pt_Particle_Physics_Model, i_min, i_max, resultat, z, n_step, pt_Spectrum_and_Precision_Parameters->spectrum_choice);


      Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.push_back(log10(resultat));
      Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei.push_back(log10(z));
      // cout << z << "  " << resultat  << endl;
      // for(int l = 0;l<Gamma_Table_Size;l++)cout << Cascade_Spectrum.Gamma_Energy[l] << "   " << Cascade_Spectrum.Gamma_Spectrum[l] <<  endl;


      Spectrum_and_cross_sections_convolution(&Cascade_Spectrum, pt_Particle_Physics_Model, j_min, j_max, resultat, z, n_step, pt_Spectrum_and_Precision_Parameters->spectrum_choice);

      Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He.push_back(log10(resultat));
      Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He.push_back(log10(z));
      // cout << z << "  " << resultat  << endl;
      // for(int l = 0;l<Gamma_Table_Size;l++)cout << Cascade_Spectrum.Gamma_Energy[l] << "   " << Cascade_Spectrum.Gamma_Spectrum[l] <<  endl;

      Spectrum_and_cross_sections_convolution(&Cascade_Spectrum, pt_Particle_Physics_Model, k_min, k_max, resultat, z, n_step, pt_Spectrum_and_Precision_Parameters->spectrum_choice);

      Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei.push_back(log10(resultat));
      Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei.push_back(log10(z));
      // cout << z << "  " << resultat << endl;
      // for(int l = 0;l<Gamma_Table_Size;l++)cout << Cascade_Spectrum.Gamma_Energy[l] << "   " << Cascade_Spectrum.Gamma_Spectrum[l] <<  endl;
      Cascade_Spectrum.Gamma_Energy.clear();
      Cascade_Spectrum.Gamma_Spectrum.clear();
    }


  // for(int l = 0;l<Integration_over_z_source_term_redshift.size();l++)cout << Integration_over_z_source_term_redshift[l] << "   " << Integration_over_z_source_term[l] <<  endl;
  /********************************************************************************************************************************************************************************/
  /****************** Fourth step : Every integral that is independant of zeta_x is done once and for all. The source term involved a double integration over z. *******************/
  /************************************************* To do so, every step in z is stored in the table Integration_over_z_source_term. *********************************************/
  /********************************************************************************************************************************************************************************/

if(verbose>1)cout <<"I'm done convoluting the spectrum with cross sections! I start to integrate over z."<<endl;

    for(int dtau = 0 ; dtau <= tau_step ; dtau++){
      tau_x = pow(10,log10(tau_min)+log10_dtau*dtau);
      if((dtau==tau_step) && (tau_x!=tau_max))cout<<"erreur : probleme de pas logarithmique en tau"<<endl;

      if(verbose>0)cout << "Current lifetime analysed : " << tau_x << endl;

      z_initial = 5*pt_Particle_Physics_Model->z_x;
      // z_final = z_min;

    // log10_dz=(log10(z_initial)-log10(z_final))/(double) n_step;
    dz=(z_initial-z_final)/(double) n_step;
    h=dz/2;
    // h=log10_dz/log10(2.);
    resultat_destruc_nuclei= 0;
    resultat_source_term = 0;
    for(int i = 0; i<n_step;i++){

      if(i==0){
        z1=z_initial;
        }
      else{
        z1=z3;
      }

      // z2=pow(10,log10(z1) - h);
      // z3=pow(10,log10(z2) - h);
      z2=z1-h;
      z3=z2-h;
      // for(int l = 0;l<Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.size();l++)cout << Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei[l] << "   " << Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei[l] <<  endl;

        //First integral : Destruction of the nuclei
        linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.size(), log10(z1), f1);
        // cout<<"redshift 1= " << z1 <<"interpolation = "<<f1<< endl;
        if(f1<0)f1=0;
        linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.size(), log10(z2), f2);
        // cout<<"redshift 2= " << z2 <<"interpolation = "<<f2<< endl;
        if(f2<0)f2=0;

        linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.size(), log10(z3), f3);
        // cout<<"redshift 3= " << z3 <<"interpolation = "<<f3<< endl;
        if(f3<0)f3=0;

        f1=pow(10,f1);
        f2=pow(10,f2);
        f3=pow(10,f3);
        f1*=exp(-1./(2*H_r*tau_x*(z1+1)*(z1+1)));
        f2*=exp(-1./(2*H_r*tau_x*(z2+1)*(z2+1)));
        f3*=exp(-1./(2*H_r*tau_x*(z3+1)*(z3+1)));

        // resultat += pow(10,h) * (f1/3. + 4.*f2/3. + f3/3.);
        resultat_destruc_nuclei += h * (f1/3. + 4.*f2/3. + f3/3.);
        // cout << " resultat_destruc_nuclei_1 integrale z = " << resultat_destruc_nuclei << " h = " << h <<  endl;


        //Second integral : Destruction of the nuclei
        linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He.size(), log10(z1), g1);
        // cout<<"redshift 1= " << z1 <<"interpolation = "<<g1<< endl;
        if(g1<0)g1=0;
        linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He.size(), log10(z2), g2);
        // cout<<"redshift 2= " << z2 <<"interpolation = "<<g2<< endl;
        if(g2<0)g2=0;

        linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He.size(), log10(z3), g3);
        // cout<<"redshift 3= " << z3 <<"interpolation = "<<g3<< endl;
        if(g3<0)g3=0;

        g1=pow(10,g1);
        g2=pow(10,g2);
        g3=pow(10,g3);
        g1*=exp(-1./(2*H_r*tau_x*(z1+1)*(z1+1)));
        g2*=exp(-1./(2*H_r*tau_x*(z2+1)*(z2+1)));
        g3*=exp(-1./(2*H_r*tau_x*(z3+1)*(z3+1)));


        // resultat += pow(10,h) * (f1/3. + 4.*f2/3. + f3/3.);
        resultat_source_term += h * ((g1-f1)/3. + 4.*(g2-f2)/3. + (g3-f3)/3.);
        // cout << " resultat_source_term = " << resultat_source_term << " z = "<< z3 << " g3 = " << g3 << " f3 = " << f3 << " exp(-1./(2*H_r*tau_x*(z3+1)*(z3+1))) = " << exp(-1./(2*H_r*tau_x*(z3+1)*(z3+1))) << endl;
        Integration_over_z_source_term_redshift.push_back(z3);
        Integration_over_z_source_term.push_back(resultat_source_term);



      }

          /********************************************************************************************************************************************************************************/
          /********************* Last step : We can now perform integrals that are not independant of zeta. Since we stored previous results in tables or variable, ***********************/
          /********************************************************* there is only one integral per zeta_x left. It is big time gain. *****************************************************/
          /********************************************************************************************************************************************************************************/


    for(int dZ = 0 ; dZ <= zeta_step ; dZ++){
      zeta_x = pow(10,log10(zeta_min)+log10_dZ*dZ);
      Fill_Structure_Particle_Physics_Model(pt_Particle_Physics_Model->M_x, zeta_x, tau_x, pt_Particle_Physics_Model); // MANDATORY STEP
      B=zeta_x*n_y_0/(pt_Particle_Physics_Model->E_0*H_r*tau_x);
      // if(verbose>1)cout << "zeta_x = " << zeta_x << "resultat = " << resultat_destruc_nuclei << " B = " << B << endl;


      resultat=0;
      for(int i = 0; i<n_step;i++){

        if(i==0){
          z1=z_initial;
          }
        else{
          z1=z3;
        }

        // z2=pow(10,log10(z1) - h);
        // z3=pow(10,log10(z2) - h);
        z2=z1-h;
        z3=z2-h;

          linearint(Integration_over_z_source_term_redshift, Integration_over_z_source_term, Integration_over_z_source_term_redshift.size(), z1, g1);
          // cout<<"redshift 1= " << z1 <<"interpolation = "<<g1<< endl;
          // if(g1<0)g1=0;
          linearint(Integration_over_z_source_term_redshift, Integration_over_z_source_term, Integration_over_z_source_term_redshift.size(), z2, g2);
          // cout<<"redshift 2= " << z2 <<"interpolation = "<<g2<< endl;
          // if(g2<0)g2=0;
          linearint(Integration_over_z_source_term_redshift, Integration_over_z_source_term, Integration_over_z_source_term_redshift.size(), z3, g3);
          // cout<<"redshift 3= " << z3 <<"interpolation = "<<g3<< endl;
          // if(g3<0)g3=0;

          g1=exp(-B*g1);
          g2=exp(-B*g2);
          g3=exp(-B*g3);
          // cout<<"redshift 1= " << z1 <<"interpolation = "<<g1<< endl;

          // cout<<"redshift 2= " << z2 <<"interpolation = "<<g2<< endl;

          // cout<<"redshift 3= " << z3 <<"interpolation = "<<g3<< endl;


          linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei.size(), log10(z1), f1);
          // cout<<"redshift 1= " << z1 <<"interpolation = "<<f1<< endl;
          if(f1<0)f1=0;
          linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei.size(), log10(z2), f2);
          // cout<<"redshift 2= " << z2 <<"interpolation = "<<f2<< endl;
          if(f2<0)f2=0;

          linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei.size(), log10(z3), f3);
          // cout<<"redshift 3= " << z3 <<"interpolation = "<<f3<< endl;
          if(f3<0)f3=0;
          f1=pow(10,f1);
          f2=pow(10,f2);
          f3=pow(10,f3);
          f1*=exp(-1./(2*H_r*tau_x*(z1+1)*(z1+1)));
          f2*=exp(-1./(2*H_r*tau_x*(z2+1)*(z2+1)));
          f3*=exp(-1./(2*H_r*tau_x*(z3+1)*(z3+1)));
          f1*=g1*K_0;
          f2*=g2*K_0;
          f3*=g3*K_0;

          // resultat += pow(10,h) * (f1/3. + 4.*f2/3. + f3/3.);
          resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
          // cout << " resultat integrale z = " << resultat << endl;

      }
      // if(verbose>1)cout << " resultat source term = " << resultat << " Y_0 = " << Y_0 << endl;
      Abundance=exp(-resultat_destruc_nuclei*B)*(Y_0+B*resultat);
      // Abundance=exp(-resultat_destruc_nuclei*B)*(Y_0);

      if(pt_Output_Options->print_result == "yes"){
        if(Abundance<Y_min){
          cout<<"Your model is ruled out !"<<endl;
          cout <<"You have overdestructed "<<pt_Scan_Parameters->nuclei<<" : Y_final = " << Abundance << " < Y_min = " << Y_min << "." << endl;
        }
        else if(Abundance>Y_max){
          cout<<"Your model is ruled out !"<<endl;
          cout <<"You have overproduced "<<pt_Scan_Parameters->nuclei<<" : Y_final = " << Abundance << " > Y_max = " << Y_max << "." << endl;
        }
        // else {
        //   cout << "The final abundance = " << Abundance << "." << endl;
        //   cout << "Your model is fine for BBN !"<<endl;
        // }
      }
      if(Abundance < Y_min || Abundance > Y_max){
        if(tau_x==tau_max)cout << "Exporting results in " << name << "."<<endl;
        file << tau_x << "  " << zeta_x << "  " << Abundance << endl;
        break;
      }
    }
    }
    Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.clear();
    Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He.clear();
    Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei.clear();
    Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei.clear();
    Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He.clear();
    Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei.clear();
    Integration_over_z_source_term_redshift.clear();
    Integration_over_z_source_term.clear();
}
double  cross_section(double  x, int i)
{



/*******fonctions a integrer*******/

/****Processus 0 = 7Li(y,t)4He****/
/****Processus 1 = 7Li(y,n)6Li****/
/****Processus 2 = 7Li(y,2np)4He****/
/****Processus 3 = 7Li(y,2H)4He+n****/
/****Processus 4 = 7Li(y,p)6He, 6He(p,n)6Li***/
/****Processus 5 = 7Li(y,3H)3H+p****/
/****Processus 6 = 7Li(y,3H)3He+n****/
/****Processus 7 = 7Be(y,3He)4He****/
/****Processus 8 = 7Be(y,p)6Li****/
/****Processus 9 = 7Be(y,2pn)4He****/
/****Processus 10 = 7Be(y,2H)4He+p****/
/****Processus 11 = 7Be(y,n)6Be, 6Be(n,4He)2p+n****/
/****Processus 12 = 7Be(y,3He)3He+n****/
/****Processus 13 = 7Be(y,3H)3He+p****/
/****Processus 14 = d(y,n)p****/
/****Processus 15 = 4He(y,p)t****/
/****Processus 16 = 4He(y,n)3He****/
/****Processus 17 = 4He(y,d)d****/
/****Processus 18 = 4He(y,np)d****/
/****Processus 19 = 3He(y,p)d****/
/****Processus 20 = 3He(y,np)p****/

	double  y=0;

	if(i == 0)
	{
	double  Q = 2.467032;
	//y = (0.105*2371*pow(x,-2)*exp(-2.5954*pow(x-Q,-0.5))*exp(-2.056*(x-Q))*(1+2.2875*pow(x-Q,2)-1.1798*pow(x-Q,3)+2.5279*pow(x-Q,4)))*pow(x,d);
	if(x>=Q)y = 0.057*931.434*pow(x,-2)*exp(-2.59*pow(x-Q,-0.5));
	}
	else if(i == 1)
	{
	double  Q = 7.249962;
	if(x>=Q)y = (0.176*pow(Q,1.51)*pow(x-Q,0.49)*pow(x,-2)+1205*pow(Q,5.5)*pow(x-Q,5)*pow(x,-10.5)+0.06/(1+pow((x-7.46)/0.188,2)));
	}
	else if(i == 2)
	{
	double  Q = 10.948850;
	if(x>=Q)y = 122*pow(Q,4)*pow(x-Q,3)*pow(x,-7);
	}
	else if(i == 3)
	{
	double  Q0 = 8.725;
	double  Q1 = 23;
	if(x>=Q0)y = 3.8 * pow(Q0,2.3)*(x-Q0)/pow(x,3.3);
	if(x>=Q1) y += 2.1*pow(Q1,1.5)*(x-Q1)/pow(x,2.5);

	}
	else if(i == 4)
	{
	double  Q = 9.98;
	if(x>=Q)y = 10.8 *pow(Q,2)*pow(x-Q,1.2)*(x-Q)/pow(x,3.2);
	}
	else if(i == 5)
	{
	double  Q = 22.28;
	if(x>=Q)y = 1.44*pow(10.,3)*pow(Q,20)*pow(x-Q,2.4)/pow(x,22.4);
	}
	else if(i == 6)
	{
	double  Q = 23.05;
	if(x>=Q)y = 1.44*pow(10.,3)*pow(Q,20)*pow(x-Q,2.4)/pow(x,22.4);
	}
	else if(i == 7)
	{
		double  Q = 1.586627;
		//y =  (0.504*2371*pow(x,-2)*exp(-5.1909*pow(x-Q,-0.5))*exp(-0.548*(x-Q))*(1-0.428*pow(x-Q,2)+0.543*pow(x-Q,3)-0.115*pow(x-Q,4)));
		if(x>=Q)y = 0.26*931.434*pow(x,-2)*exp(-5.19*pow(x-Q,-0.5));
		else y=0;
	}
	else if(i == 8)
	{
		double  Q = 5.605794;
		if(x>=Q)y = (32.6*pow(Q,10)*pow((x-Q),2)*pow(x,-12)+2.27*pow(10.,6)*pow(Q,8.8335)*pow((x-Q),13)*pow(x,-21.8335));
	}
	else if(i == 9)
	{
		double  Q = 9.30468;
		if(x>=Q)y = 133*pow(Q,4)*pow(x-Q,3)*pow(x,-7);
	}
	else if(i == 10)
	{
	double  Q0 = 7.08;
	double  Q1 = 23;
	if(x>=Q0)y = 3.8 * pow(Q0,2.3)*(x-Q0)/pow(x,3.3);
	if(x>=Q1) y += 2.1*pow(Q1,1.5)*(x-Q1)/pow(x,2.5);
	}
	else if(i == 11)
	{
	double  Q = 10.68;
	if(x>=Q)y = 10.8 *pow(Q,2)*pow(x-Q,1.2)*(x-Q)/pow(x,3.2);
	}
	else if(i == 12)
	{
	double  Q = 22.17;
	if(x>=Q)y = 1.44*pow(10.,3)*pow(Q,20)*pow(x-Q,2.4)/pow(x,22.4);
	}
	else if(i == 13)
	{
	double  Q = 21.4;
	if(x>=Q)y = 1.44*pow(10.,3)*pow(Q,20)*pow(x-Q,2.4)/pow(x,22.4);
	}
	else if(i == 14)
	{
	double  Q = 2.224573;
	if(x>=Q)y = 18.75*(pow(pow(Q*(x-Q),0.5)/x,3)+0.007947*pow(pow(Q*(x-Q),0.5)/x,2)*pow(pow(Q,0.5)-pow(0.037,0.5),2)/(x-Q+0.037));
	}
	else if(i == 15)
	{
	double  Q = 19.813852;
	//~ y = 128.9*pow(Q,4.524)*pow(x-Q,2.512)/pow(x,4.524+2.512);
	if(x>=Q)y = 19.5*pow(Q,3.5)*pow(x-Q,1)/pow(x,4.5);
	}
	else if(i == 16)
	{
	double  Q = 20.577615;
	//~ y = 31.68*pow(Q,3.663)*pow(x-Q,1.580)/pow(x,3.663+1.580);
	if(x>=Q)y = 17.1*pow(Q,3.5)*pow(x-Q,1.)/pow(x,4.5);
	}
	else if(i == 17)
	{
	double  Q = 23.846527;
	if(x>=Q)y = 10.7*pow(Q,10.2)*pow(x-Q,3.4)/pow(x,13.6);
	}
	else if(i == 18)
	{
	double  Q = 26.0711;
	if(x>=Q)y = 21.7*pow(Q,4.0)*pow(x-Q,3.0)/pow(x,7.0);
	}
	else if(i == 19)
	{
	double  Q = 5.483485;
	if(x>=Q)y = 8.88*pow(Q,1.75)*pow(x-Q,1.65)/pow(x,3.4);
	}
	else if(i == 20)
	{
	double  Q = 7.718058;
	if(x>=Q)y = 16.7*pow(Q,1.95)*pow(x-Q,2.3)/pow(x,4.25);
	}
	y=y*2.569*pow(10.,-6);//Conversion mb en Mev^-2
	return y;

}


void Check_nuclei(string nuclei,
                  int &i_min,
                  int &i_max,
                  int &k_min,
                  int &k_max,
                  double &Y_min,
                  double &Y_max,
                  double &Y_0){

    if(nuclei =="4He"){
      if(verbose>1)cout<<"I could recognize nuclei, it's 4He."<<endl;
      Y_0 = Y_4He_0;
      Y_min = Y_4He_Min;
      Y_max = Y_4He_Max;
      if(verbose>1)cout << " Y_4He_0 = " << Y_0 << " Y_4He_min = " << Y_min << " Y_4He_max = " << Y_max << endl;
      //Destruction cross-sections
      i_min = 15;
      i_max = 18;
      //Production cross-sections
      k_min = 0;
      k_max = 0; //4He cannot be produced in significant quantities by the non-thermal BBN !
    }
    else if(nuclei=="3He"){
      if(verbose>1)cout<<"I could recognize nuclei, it's 3He."<<endl;
      Y_0 = Y_3He_0;
      Y_min = Y_3He_Min;
      Y_max = Y_3He_Max;
      if(verbose>1)cout << " Y_3He_0 = " << Y_0 << " Y_3He_min = " << Y_min << " Y_3He_max = " << Y_max << endl;
      //Destruction cross-sections
      i_min = 19;
      i_max = 20;
      //Production cross-sections
      k_min = 15;
      k_max = 16;
    }
    else if(nuclei=="2H"){
      if(verbose>1)cout<<"I could recognize nuclei, it's 2H."<<endl;
      Y_0 = Y_2H_0;
      Y_min = Y_2H_Min;
      Y_max = Y_2H_Max;
      if(verbose>1)cout << " Y_2H_0 = " << Y_0 << " Y_2H_min = " << Y_min << " Y_2H_max = " << Y_max << endl;
      //Destruction cross-sections
      i_min = 14;
      i_max = 14;
      //Production cross-sections
      k_min = 17;
      k_max = 18;
    }
    else if(nuclei=="7Li"){
      if(verbose>1)cout<<"I could recognize nuclei, it's 7Li."<<endl;
      Y_0 = Y_7Li_0;
      Y_min = Y_7Li_Min;
      Y_max = Y_7Li_Max;

      if(verbose>1)cout << " Y_7li_0 = " << Y_0 << " Y_7Li_min = " << Y_min << " Y_7Li_max = " << Y_max << endl;
      //Destruction cross-sections
      i_min = 0;
      i_max = 6;
      //Production cross-sections
      k_min = 0;
      k_max = 0;
    }
    else if(nuclei=="7Be"){
      if(verbose>1)cout<<"I could recognize nuclei, it's 7Be."<<endl;
      Y_0 = Y_7Be_0;
      Y_min = Y_7Be_Min;
      Y_max = Y_7Be_Max;

      if(verbose>1)cout << " Y_7Be_0 = " << Y_0 << " Y_7Be_min = " << Y_min << " Y_7Be_max = " << Y_max << endl;
      //Destruction cross-sections
      i_min = 7;
      i_max = 13;
      //Production cross-sections
      k_min = 0;
      k_max = 0;
    }
    else{
      if(verbose>0){cout<<"I couldn't recognize nuclei."<<endl;
      cout << "The string 'nuclei' needs to be a word from this list : '4He', '3He', '2H', '7Li, '7Be'." << endl;
      cout << "Please check. " << endl;
      }
      return;
      }

}
