#include "bbn/EM_cascade.h"
#include "bbn/injected_spectrum.h"
#include "bbn/structures.h"
#include "bbn/BBN_constraints.h"
#include "bbn/tools.h"
#include "bbn/test_functions.h"

using namespace std;


static void Spectrum_and_cross_sections_convolution(Structure_Spectrum * pt_Cascade_Spectrum,
        Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        const int i_min,
        const int i_max,
        double &resultat,
        double z)
{
    double tau = pt_Particle_Physics_Model->tau_x;
    double Z_x = pt_Particle_Physics_Model->zeta_x;
    double E_0 = pt_Particle_Physics_Model->E_0;
    double E_max;
    double E_c = E_c_0/(1+z);;
    double E[pt_Spectrum_and_Precision_Parameters->eval_max];
    double f[pt_Spectrum_and_Precision_Parameters->eval_max];
    double cross_sections[pt_Spectrum_and_Precision_Parameters->eval_max];
    int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
    double h,dE;
    int y;
    if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "universal" || pt_Spectrum_and_Precision_Parameters->calculation_mode== "simplified"){
      E_max = E_0;
      if(E_0 > 1.5*E_c)E_max = 1.5*E_c;
    }
    else{
      E_max = E_0;
    }
    // E_max = E_0;
    // if(E_0 > 1.5*E_c) {
    //     E_max = 1.5*E_c;
    // }
    resultat = 0;
    dE = (E_max - pt_Spectrum_and_Precision_Parameters->E_min_table)/ (double) n_step;
    y = 0;
    // while(dE>pt_Spectrum_and_Precision_Parameters->E_min_table){
    //   dE/=10.;
    //   y++;
    // }
    h = dE/(pt_Spectrum_and_Precision_Parameters->eval_max-1);

    // cout << " dE = " << dE  << " y " << y << endl;
    // ds = (E*E_gamma_bb/(m_e*m_e) - 1)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);

    for(int i=0; i<pow(10,y)*n_step; i++) {
        for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {
            cross_sections[eval]=0.;
        }

        // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
        for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {

            E[eval]=pt_Spectrum_and_Precision_Parameters->E_min_table+eval*h+i*dE;
            linearint(pt_Cascade_Spectrum->Energy, pt_Cascade_Spectrum->Spectrum, pt_Cascade_Spectrum->Energy.size(), E[eval], f[eval]);
            for(int j = i_min; j<=i_max; j++) {
                cross_sections[eval] += cross_section(E[eval],j);
                /**** If you are including 2H production, the channel 17 has a multiplicity of 2, it is now taken into account ****/
                if(i_min == 17 && i_max == 18 && j == i_min) {
                    cross_sections[eval] += cross_section(E[eval],j);
                }
            }
            resultat += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval]*cross_sections[eval];
            // cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
        }
    }

    if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "Dirac" && E_c>pt_Particle_Physics_Model->E_0 && pt_Spectrum_and_Precision_Parameters->calculation_mode == "iterative") {
        // cout << "works"<<endl;

        for(int i = i_min; i<=i_max; i++) {
            resultat+=cross_section(pt_Particle_Physics_Model->E_0,i)/(rate_NPC(pt_Particle_Physics_Model->E_0,z)+rate_compton(pt_Particle_Physics_Model->E_0,z)+rate_gg_scattering(pt_Particle_Physics_Model->E_0,z));
            if(i_min == 17 && i_max == 18 && i == i_min) {
            resultat+=cross_section(pt_Particle_Physics_Model->E_0,i)/(rate_NPC(pt_Particle_Physics_Model->E_0,z)+rate_compton(pt_Particle_Physics_Model->E_0,z)+rate_gg_scattering(pt_Particle_Physics_Model->E_0,z));
            }
        }
    }
    // cout << "resultat " << resultat << endl;
}

static void  Compute_Constraints_from_destruction_only_loop(const int step,
        Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Scan_Parameters_and_Results * pt_Scan_Parameters_and_Results,
        Structure_Output_Options * pt_Output_Options,
        vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei,
        vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei)
{


    double z_step, n_step;
    double tau_min;
    double z_initial, z_final, z_x;
    double  E_c, z;
    int i_min, i_max;

    i_min = pt_Scan_Parameters_and_Results->i_min;
    i_max = pt_Scan_Parameters_and_Results->i_max;

    z_step = pt_Spectrum_and_Precision_Parameters->z_step;
    n_step = pt_Spectrum_and_Precision_Parameters->n_step;
    tau_min = pt_Scan_Parameters_and_Results->tau_min;
    z_initial =  (5*(pow(2*H_r*tau_min,-0.5)-1));
    z_final =  Z_MIN_SCAN;
    z = z_initial*pow(z_final/z_initial,(double) step/(z_step));

    double resultat;

    Structure_Spectrum Cascade_Spectrum;
    Cascade_Spectrum.species="photon";
    Cascade_Spectrum.spectrum_name="total_photons";

    //  z=pow(10,log10(z_initial)-j*log10_dz);aa
    E_c = E_c_0/(1+z);
    if(pt_Output_Options->BBN_constraints_verbose > 0)  {
        #pragma omp critical(print)
        {
            cout<<"redshift = " << z << " still " << z_step-step << " points to go." << endl;

        }
    }

    if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "universal") {
        Cascade_Spectrum_Calculation(z,
                                     pt_Output_Options,
                                     pt_Particle_Physics_Model,
                                     &Cascade_Spectrum,
                                     pt_Spectrum_and_Precision_Parameters);
    } else {
        if(pt_Spectrum_and_Precision_Parameters->spectrum_mode=="writing" || pt_Spectrum_and_Precision_Parameters->spectrum_mode == "nothing") {
            Cascade_Spectrum_Calculation(z,
                                         pt_Output_Options,
                                         pt_Particle_Physics_Model,
                                         &Cascade_Spectrum,
                                         pt_Spectrum_and_Precision_Parameters);
        } else if(pt_Spectrum_and_Precision_Parameters->spectrum_mode=="reading") {
            //    if(E_c <= pt_Particle_Physics_Model->E_0 ){
            //      Cascade_Spectrum.Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
            //      Cascade_Spectrum.Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
            //      for(int i=0;i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;i++){
            //        Cascade_Spectrum.Energy[i]=pt_Spectrum_and_Precision_Parameters->E_min_table+i*dE;
            //        Cascade_Spectrum.Spectrum[i]=universal_spectrum(pt_Spectrum_and_Precision_Parameters->E_min_table+i*dE,z,pt_Particle_Physics_Model->E_0);
            //
            //       }
            //    }
            //
            //
            // else
            Cascade_Spectrum_Reading_From_File(z,
                                               pt_Particle_Physics_Model,
                                               &Cascade_Spectrum,
                                               pt_Spectrum_and_Precision_Parameters);
            if(pt_Output_Options->BBN_constraints_verbose>2) {
                #pragma omp critical(print)
                {
                    cout << "redshift " << z << endl;
                    for(int i = 0; i < Cascade_Spectrum.Spectrum.size(); i++)
                    {
                        cout << "read spectrum from file = "<< Cascade_Spectrum.Spectrum[i]  << " energy = " << Cascade_Spectrum.Energy[i]<<endl;
                    }
                }
            }
        }

    }

    Spectrum_and_cross_sections_convolution(&Cascade_Spectrum,
                                            pt_Particle_Physics_Model,
                                            pt_Spectrum_and_Precision_Parameters,
                                            i_min,
                                            i_max,
                                            resultat,
                                            z);
    #pragma omp critical(dataupdate)
    {
        Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei[step]=log10(resultat);
        Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei[step]=log10(z);
    }

    if(pt_Output_Options->BBN_constraints_verbose>1) {
        cout << "results convolution spectrum and cross section = "<< resultat << " redshift " << z << endl;
    }
    Cascade_Spectrum.Energy.clear();
    Cascade_Spectrum.Spectrum.clear();



}

static void  compute_constraints_from_destruction_and_production_loop(const int step,
        Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Scan_Parameters_and_Results * pt_Scan_Parameters_and_Results,
        Structure_Scan_Parameters_and_Results * pt_Destruction_4He,
        Structure_Output_Options * pt_Output_Options,
        vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei,
        vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei,
        vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He,
        vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He,
        vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei,
        vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei)
{


    double z_step, n_step;
    double tau_min;
    double z_initial, z_final, z_x;
    double  E_c, z;
    int i_min, i_max, j_min, j_max, k_min, k_max;

    i_min = pt_Scan_Parameters_and_Results->i_min;
    i_max = pt_Scan_Parameters_and_Results->i_max;
    k_min = pt_Scan_Parameters_and_Results->k_min;
    k_max = pt_Scan_Parameters_and_Results->k_max;
    j_min = pt_Destruction_4He->i_min;
    j_max = pt_Destruction_4He->i_max;

    z_step = pt_Spectrum_and_Precision_Parameters->z_step;
    n_step = pt_Spectrum_and_Precision_Parameters->n_step;
    tau_min = pt_Scan_Parameters_and_Results->tau_min;
    z_initial =  (5*(pow(2*H_r*tau_min,-0.5)-1));
    z_final =  Z_MIN_SCAN;
    z = z_initial*pow(z_final/z_initial,(double) step/(z_step));

    double resultat;

    Structure_Spectrum Cascade_Spectrum;
    Cascade_Spectrum.species="photon";
    Cascade_Spectrum.spectrum_name="total_photons";

    //  z=pow(10,log10(z_initial)-j*log10_dz);aa
    E_c = E_c_0/(1+z);
    if(pt_Output_Options->BBN_constraints_verbose > 0)  {
        #pragma omp critical(print)
        {
            cout<<"redshift = " << z << " still " << z_step-step << " points to go." << endl;

        }
    }

    if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "universal") {
        Cascade_Spectrum_Calculation(z,
                                     pt_Output_Options,
                                     pt_Particle_Physics_Model,
                                     &Cascade_Spectrum,
                                     pt_Spectrum_and_Precision_Parameters);
    } else {
        if(pt_Spectrum_and_Precision_Parameters->spectrum_mode=="writing" || pt_Spectrum_and_Precision_Parameters->spectrum_mode == "nothing") {
            Cascade_Spectrum_Calculation(z,
                                         pt_Output_Options,
                                         pt_Particle_Physics_Model,
                                         &Cascade_Spectrum,
                                         pt_Spectrum_and_Precision_Parameters);
        } else if(pt_Spectrum_and_Precision_Parameters->spectrum_mode=="reading") {
            //    if(E_c <= pt_Particle_Physics_Model->E_0 ){
            //      Cascade_Spectrum.Energy.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
            //      Cascade_Spectrum.Spectrum.resize(pt_Spectrum_and_Precision_Parameters->Energy_Table_Size);
            //      for(int i=0;i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;i++){
            //        Cascade_Spectrum.Energy[i]=pt_Spectrum_and_Precision_Parameters->E_min_table+i*dE;
            //        Cascade_Spectrum.Spectrum[i]=universal_spectrum(pt_Spectrum_and_Precision_Parameters->E_min_table+i*dE,z,pt_Particle_Physics_Model->E_0);
            //
            //       }
            //    }
            //
            //
            // else
            Cascade_Spectrum_Reading_From_File(z,
                                               pt_Particle_Physics_Model,
                                               &Cascade_Spectrum,
                                               pt_Spectrum_and_Precision_Parameters);
            if(pt_Output_Options->BBN_constraints_verbose>2) {
                #pragma omp critical(print)
                {
                    cout << "redshift " << z << endl;
                    for(int i = 0; i < Cascade_Spectrum.Spectrum.size(); i++)
                    {
                        cout << "read spectrum from file = "<< Cascade_Spectrum.Spectrum[i]  << " energy = " << Cascade_Spectrum.Energy[i]<<endl;
                    }
                }
            }
        }

    }

    Spectrum_and_cross_sections_convolution(&Cascade_Spectrum,
                                            pt_Particle_Physics_Model,
                                            pt_Spectrum_and_Precision_Parameters,
                                            i_min,
                                            i_max,
                                            resultat,
                                            z);
    #pragma omp critical(dataupdate)
    {
        Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei[step]=log10(resultat);
        Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei[step]=log10(z);
        if(pt_Output_Options->BBN_constraints_verbose>1) {
            cout << "results convolution spectrum and cross section = "<< resultat << " redshift " << z << endl;
        }
    }




    Spectrum_and_cross_sections_convolution(&Cascade_Spectrum,
                                            pt_Particle_Physics_Model,
                                            pt_Spectrum_and_Precision_Parameters,
                                            j_min,
                                            j_max,
                                            resultat,
                                            z);

    #pragma omp critical(dataupdate)
    {
        Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He[step]=log10(resultat);
        Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He[step]=log10(z);
        if(pt_Output_Options->BBN_constraints_verbose>1) {
            cout << "results convolution spectrum and cross section destruction 4He = "<< resultat << "redshift" << z << endl;
        }
    }


    Spectrum_and_cross_sections_convolution(&Cascade_Spectrum,
                                            pt_Particle_Physics_Model,
                                            pt_Spectrum_and_Precision_Parameters,
                                            k_min,
                                            k_max,
                                            resultat,
                                            z);
    #pragma omp critical(dataupdate)
    {
        Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei[step]=log10(resultat);
        Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei[step]=log10(z);
        if(pt_Output_Options->BBN_constraints_verbose>1) {
            cout << "results convolution spectrum and cross section production nuclei = "<< resultat << "redshift" << z << endl;
        }
    }

    Cascade_Spectrum.Energy.clear();
    Cascade_Spectrum.Spectrum.clear();



}

void Compute_Constraints_from_destruction_only(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Scan_Parameters_and_Results * pt_Scan_Parameters_and_Results,
        Structure_Output_Options * pt_Output_Options)
{



    vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei;
    vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei;

    /******** This step attributes locally the values of precision and scan parameters ********/
    double tau_min, tau_max, tau_step, zeta_min, zeta_max, zeta_step, z_step, n_step, number_iterations_photon;
    tau_min = pt_Scan_Parameters_and_Results->tau_min;
    tau_max = pt_Scan_Parameters_and_Results->tau_max;
    tau_step = pt_Scan_Parameters_and_Results->tau_step;
    zeta_min = pt_Scan_Parameters_and_Results->zeta_min;
    zeta_max = pt_Scan_Parameters_and_Results->zeta_max;
    zeta_step = pt_Scan_Parameters_and_Results->zeta_step;
    n_step = pt_Spectrum_and_Precision_Parameters->n_step;
    z_step = pt_Spectrum_and_Precision_Parameters->z_step;
    number_iterations_photon = pt_Spectrum_and_Precision_Parameters->number_iterations_photon;
    /******************************************************************************************/
    Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.resize(z_step+1);
    Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei.resize(z_step+1);

    double z_initial =  (5*(pow(2*H_r*tau_min,-0.5)-1)), z_final =  Z_MIN_SCAN, z_x;
    double zeta_x, tau_x;
    double log10_dtau = (log10(tau_max)-log10(tau_min))/(double) tau_step;
    double log10_dZ = (log10(zeta_max)-log10(zeta_min))/(double) zeta_step;
    double log10_dz = (log10(z_initial)-log10(z_final))/(double) z_step;
    double Y_0, Y_min, Y_max;
    double Abundance, resultat, B, E_c;
    double f[pt_Spectrum_and_Precision_Parameters->eval_max], z_array[pt_Spectrum_and_Precision_Parameters->eval_max], y, h, dz;
    double f1,f2,f3,z1,z2,z3;
    double E1, dE;


    Check_nuclei(pt_Scan_Parameters_and_Results);
    Y_0 = pt_Scan_Parameters_and_Results->Y_0;
    Y_min = pt_Scan_Parameters_and_Results->Y_min;
    Y_max = pt_Scan_Parameters_and_Results->Y_max;

    dE= (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table)/ (double) (pt_Spectrum_and_Precision_Parameters->Energy_Table_Size-1);

    if(pt_Output_Options->BBN_constraints_verbose>1) {
        cout << " z initial = " << z_initial << " z final = " << z_final << endl;
    }
    if(pt_Output_Options->BBN_constraints_verbose > 0) {
        cout << "I start generating spectrum for each redshift. You asked for " << z_step << " points." << endl;
    }
    {
        int end = z_step;
        #pragma omp parallel for ordered schedule(dynamic)
        for(int j = 0; j<=end; j++) {

            Compute_Constraints_from_destruction_only_loop(j,
                    pt_Particle_Physics_Model,
                    pt_Spectrum_and_Precision_Parameters,
                    pt_Scan_Parameters_and_Results,
                    pt_Output_Options,
                    Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei,
                    Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei);


        }
    }

    if(pt_Output_Options->BBN_constraints_verbose > 0) {
        cout <<"I'm done convoluting the spectrum with cross sections! I start to integrate over z."<<endl;
    }

    for(int dtau = 0 ; dtau <= tau_step ; dtau++) {
        // tau_x = pow(10,log10(tau_min)+log10_dtau*dtau);
        tau_x = tau_min*pow(tau_max/tau_min,(double) dtau/(tau_step));

        if((dtau==tau_step) && (tau_x!=tau_max)) {
            cout<<"erreur : probleme de pas logarithmique en tau"<<endl;
        }

        if(pt_Output_Options->BBN_constraints_verbose > 0) {
            cout << "Current lifetime analysed : " << tau_x << endl;
        }
        z_x = pow(tau_x*(2*H_r),-0.5)-1;
        z_initial = 5*z_x;

        // z_final = z_min;
        dz=(z_initial-z_final)/(double) n_step;

        y = 0;
        while(dz>z_initial) {
            dz/=10.;
            y++;
        }
        h = dz/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
        resultat = 0.;
        // cout << " dE = " << dE  << " y " << y << endl;
        // ds = (E*E_gamma_bb/(m_e*m_e) - 1)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);

        for(int i=0; i<pow(10,y)*n_step; i++) {

            // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
            for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {
                if(eval == 0) {
                    if(i==0)	{
                        z_array[eval]=z_initial;
                    } else {
                        z_array[eval]=z_array[pt_Spectrum_and_Precision_Parameters->eval_max-1];
                    }
                } else {
                    z_array[eval]=z_array[0]-eval*h;
                }


                linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.size(), log10(z_array[eval]), f[eval]);
                // cout<<"redshift 1= " << z1 <<"interpolation = "<<f1<< endl;
                if(f[eval]<0) {
                    f[eval]=0;
                }
                f[eval]=pow(10,f[eval]);
                f[eval]*=exp(-1./(2*H_r*tau_x*(z_array[eval]+1)*(z_array[eval]+1)));

                resultat += dz/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
                // cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
            }
        }


        for(int dzeta = 0 ; dzeta <= zeta_step ; dzeta++) {
            // zeta_x = pow(10,log10(zeta_min)+log10_dZ*dzeta);
            zeta_x = zeta_min*pow(zeta_max/zeta_min,(double) dzeta/(zeta_step));

            // cout << "tau_x = "<< tau_x<<" zeta_x = " << zeta_x << "resultat = " << resultat << "E0" << pt_Particle_Physics_Model->E_0<<endl;
            // B=zeta_x*n_y_0/(H_r*tau_x);
            B=zeta_x*n_y_0/(pt_Particle_Physics_Model->E_0*H_r*tau_x);
            Abundance=exp(-resultat*B);
            Abundance*=Y_0;
            // cout << " B = " << B << "Abundance = "<<Abundance <<"Y0 = " << Y_0 << "Y_min = " << Y_min << endl;


            if(Abundance < Y_min || Abundance > Y_max) {
                if(pt_Output_Options->BBN_constraints_verbose>1) {
                    cout << "The final abundance = " << Abundance << endl;
                }
                pt_Scan_Parameters_and_Results->Results_scan_tau_x.push_back(tau_x);
                pt_Scan_Parameters_and_Results->Results_scan_zeta_x.push_back(zeta_x);
                pt_Scan_Parameters_and_Results->Results_scan_Abundance.push_back(Abundance);
                break;
            }
        }
    }


    print_results_scan(pt_Output_Options,
                       pt_Spectrum_and_Precision_Parameters,
                       pt_Scan_Parameters_and_Results,
                       pt_Particle_Physics_Model);

    Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei.clear();
    Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.clear();
}


static void Compute_constraints_from_destruction_and_production_loop_2(const int step,
        Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
        Structure_Scan_Parameters_and_Results * pt_Scan_Parameters_and_Results,
        Structure_Scan_Parameters_and_Results * pt_Destruction_4He,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Output_Options * pt_Output_Options,
        vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei,
        vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He,
        vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei,
        vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei,
        vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He,
        vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei)
{


  double zeta_x, tau_x, z_x;
  /******** This step attributes locally the values of precision and scan parameters ********/
  double tau_min, tau_max, tau_step, zeta_min, zeta_max, zeta_step, z_step, n_step, number_iterations_photon;
  tau_min = pt_Scan_Parameters_and_Results->tau_min;
  tau_max = pt_Scan_Parameters_and_Results->tau_max;
  tau_step = pt_Scan_Parameters_and_Results->tau_step;
  zeta_min = pt_Scan_Parameters_and_Results->zeta_min;
  zeta_max = pt_Scan_Parameters_and_Results->zeta_max;
  zeta_step = pt_Scan_Parameters_and_Results->zeta_step;
  n_step = pt_Spectrum_and_Precision_Parameters->n_step;
  z_step = pt_Spectrum_and_Precision_Parameters->z_step;
  number_iterations_photon = pt_Spectrum_and_Precision_Parameters->number_iterations_photon;
  /******************************************************************************************/
  double z_initial =  (5*(pow(2*H_r*tau_min,-0.5)-1)), z_final =  Z_MIN_SCAN;
  double Y_0, Y_min, Y_max;
  double K_0, K_min, K_max;
  double dz,h;
  double f[pt_Spectrum_and_Precision_Parameters->eval_max], g[pt_Spectrum_and_Precision_Parameters->eval_max], z_array[pt_Spectrum_and_Precision_Parameters->eval_max], y;
  double resultat_destruc_nuclei, resultat_source_term ;
  double Abundance, resultat, B, E_c;

  Y_0 = pt_Scan_Parameters_and_Results->Y_0;
  Y_min = pt_Scan_Parameters_and_Results->Y_min;
  Y_max = pt_Scan_Parameters_and_Results->Y_max;

  K_0 = pt_Destruction_4He->Y_0;
  K_min = pt_Destruction_4He->Y_min;
  K_max = pt_Destruction_4He->Y_max;


  vector<double> Integration_over_z_source_term_redshift;
  vector<double> Integration_over_z_source_term;
  // tau_x = pow(10,log10(tau_min)+log10_dtau*dtau);
  tau_x = tau_min*pow(tau_max/tau_min,(double) step/(tau_step));
  #pragma omp critical(print)
  {
      if((step==tau_step) && (tau_x!=tau_max)) {
          cout<<"erreur : probleme de pas logarithmique en tau"<<endl;
      }

      if(pt_Output_Options->BBN_constraints_verbose > 0) {
          cout << "Current lifetime analysed : " << tau_x << endl;
      }
  }
  z_x = pow(tau_x*(2*H_r),-0.5)-1;
  z_initial = 5*z_x;
  // z_final = z_min;

  // log10_dz=(log10(z_initial)-log10(z_final))/(double) n_step;
  dz=(z_initial-z_final)/(double) n_step;
  y = 0;
  while(dz>z_initial) {
      dz/=10.;
      y++;
  }
  h = dz/(pt_Spectrum_and_Precision_Parameters->eval_max-1);

  // cout << " dE = " << dE  << " y " << y << endl;
  // ds = (E*E_gamma_bb/(m_e*m_e) - 1)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
  resultat_destruc_nuclei= 0;
  resultat_source_term = 0;
  for(int i=0; i<pow(10,y)*n_step; i++) {

      // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
      for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {
          if(eval == 0) {
              if(i==0)	{
                  z_array[eval]=z_initial;
              } else {
                  z_array[eval]=z_array[pt_Spectrum_and_Precision_Parameters->eval_max-1];
              }
          } else {
              z_array[eval]=z_array[0]-eval*h;
          }


          linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.size(), log10(z_array[eval]), f[eval]);
          // cout<<"redshift 1= " << z1 <<"interpolation = "<<f1<< endl;
          if(f[eval]<0) {
              f[eval]=0;
          }
          f[eval]=pow(10,f[eval]);
          f[eval]*=exp(-1./(2*H_r*tau_x*(z_array[eval]+1)*(z_array[eval]+1)));

          resultat_destruc_nuclei += dz/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
          #pragma omp critical(print)
          {
            if(pt_Output_Options->BBN_constraints_verbose>2) {
                cout << " resultat integrale destruc nuclei over z = " << resultat_destruc_nuclei << endl;
            }
          }
          linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He, Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He.size(), log10(z_array[eval]), g[eval]);
          if(g[eval]<0) {
              g[eval]=0;
          }
          g[eval]=pow(10,g[eval]);
          g[eval]*=exp(-1./(2*H_r*tau_x*(z_array[eval]+1)*(z_array[eval]+1)));

          resultat_source_term += dz/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*(g[eval]-f[eval]);

          // cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
      }
      Integration_over_z_source_term_redshift.push_back(z_array[pt_Spectrum_and_Precision_Parameters->eval_max-1]);
      Integration_over_z_source_term.push_back(resultat_source_term);
      #pragma omp critical(print)
      {
        if(pt_Output_Options->BBN_constraints_verbose>2) {
            cout << " result integrale over z source term = " << resultat_source_term << endl;
        }
      }

  }


  /********************************************************************************************************************************************************************************/
  /********************* Last step : We can now perform integrals that are not independant of zeta. Since we stored previous results in tables or variable, ***********************/
  /********************************************************* there is only one integral per zeta_x left. It is big time gain. *****************************************************/
  /********************************************************************************************************************************************************************************/

  for(int dzeta = 0 ; dzeta <= zeta_step ; dzeta++) {
      // zeta_x = pow(10,log10(zeta_min)+log10_dZ*dZ);
      zeta_x = zeta_min*pow(zeta_max/zeta_min,(double) dzeta/(zeta_step));
      B=zeta_x*n_y_0/(pt_Particle_Physics_Model->E_0*H_r*tau_x);
      // if(pt_Output_Options->BBN_constraints_verbose>1)cout << "zeta_x = " << zeta_x << "resultat = " << resultat_destruc_nuclei << " B = " << B << endl;

      // cout << " dE = " << dE  << " y " << y << endl;
      // ds = (E*E_gamma_bb/(m_e*m_e) - 1)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
      resultat = 0;
      for(int i=0; i<pow(10,y)*n_step; i++) {

          // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
          for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++) {
              if(eval == 0) {
                  if(i==0)	{
                      z_array[eval]=z_initial;
                  } else {
                      z_array[eval]=z_array[pt_Spectrum_and_Precision_Parameters->eval_max-1];
                  }
              } else {
                  z_array[eval]=z_array[0]-eval*h;
              }


              linearint(Integration_over_z_source_term_redshift, Integration_over_z_source_term, Integration_over_z_source_term_redshift.size(), z_array[eval], g[eval]);
              g[eval]=exp(-B*g[eval]);
              linearint(Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei, Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei.size(), log10(z_array[eval]), f[eval]);

              if(f[eval]<0) {
                  f[eval]=0;
              }
              f[eval]=pow(10,f[eval]);
              f[eval]*=exp(-1./(2*H_r*tau_x*(z_array[eval]+1)*(z_array[eval]+1)));
              f[eval]*=g[eval]*K_0;

              resultat += dz/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
              // cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
          }
          #pragma omp critical(print)
          {
            if(pt_Output_Options->BBN_constraints_verbose>2) {
                cout << " resultat integrale z = " << resultat << endl;
            }
          }
      }
      // if(pt_Output_Options->BBN_constraints_verbose>1)cout << " resultat source term = " << resultat << " Y_0 = " << Y_0 << endl;
      Abundance=exp(-resultat_destruc_nuclei*B)*(Y_0+B*resultat);
      // Abundance=exp(-resultat_destruc_nuclei*B)*(Y_0);


      if(Abundance < Y_min || Abundance > Y_max) {
        #pragma omp critical(print)
        {
          if(pt_Output_Options->BBN_constraints_verbose>1) {
              cout << "The final abundance = " << Abundance << endl;
          }
          pt_Scan_Parameters_and_Results->Results_scan_tau_x[step]=tau_x;
          pt_Scan_Parameters_and_Results->Results_scan_zeta_x[step]=zeta_x;
          pt_Scan_Parameters_and_Results->Results_scan_Abundance[step]=Abundance;
        }
        break;

      }

  }
  Integration_over_z_source_term_redshift.clear();
  Integration_over_z_source_term.clear();
}



void Compute_constraints_from_destruction_and_production(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Scan_Parameters_and_Results * pt_Scan_Parameters_and_Results,
        Structure_Output_Options * pt_Output_Options)
{

    double z,dz,h;
    vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei;
    vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He;
    vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei;
    vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei;
    vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He;
    vector<double> Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei;
    vector<double> Integration_over_z_source_term_redshift;
    vector<double> Integration_over_z_source_term;
    Structure_Scan_Parameters_and_Results Destruction_4He;
    Destruction_4He.nuclei = "4He";
    Structure_Spectrum Cascade_Spectrum;
    Cascade_Spectrum.species="photon";
    /******** This step attributes locally the values of precision and scan parameters ********/
    double tau_min, tau_max, tau_step, zeta_min, zeta_max, zeta_step, z_step, n_step, number_iterations_photon;
    tau_min = pt_Scan_Parameters_and_Results->tau_min;
    tau_max = pt_Scan_Parameters_and_Results->tau_max;
    tau_step = pt_Scan_Parameters_and_Results->tau_step;
    zeta_min = pt_Scan_Parameters_and_Results->zeta_min;
    zeta_max = pt_Scan_Parameters_and_Results->zeta_max;
    zeta_step = pt_Scan_Parameters_and_Results->zeta_step;
    n_step = pt_Spectrum_and_Precision_Parameters->n_step;
    z_step = pt_Spectrum_and_Precision_Parameters->z_step;
    number_iterations_photon = pt_Spectrum_and_Precision_Parameters->number_iterations_photon;
    /******************************************************************************************/


    Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.resize(z_step+1);
    Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He.resize(z_step+1);
    Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei.resize(z_step+1);
    Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei.resize(z_step+1);
    Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He.resize(z_step+1);
    Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei.resize(z_step+1);

    pt_Scan_Parameters_and_Results->Results_scan_tau_x.resize(tau_step+1);
    pt_Scan_Parameters_and_Results->Results_scan_zeta_x.resize(tau_step+1);
    pt_Scan_Parameters_and_Results->Results_scan_Abundance.resize(tau_step+1);

    int i_min, i_max;
    int j_min, j_max;
    int k_min, k_max;
    int l_min, l_max;
    double zeta_x, tau_x, z_x;
    // double M_x = pt_Particle_Physics_Model->M_x;
    double E1, dE;
    double z_initial =  (5*(pow(2*H_r*tau_min,-0.5)-1)), z_final =  Z_MIN_SCAN;
    double log10_dz=(log10(z_initial)-log10(z_final))/(double) z_step;
    double log10_dtau = (log10(tau_max)-log10(tau_min))/(double) tau_step;
    double log10_dZ = (log10(zeta_max)-log10(zeta_min))/(double) zeta_step;
    double Y_0, Y_min, Y_max;
    double K_0, K_min, K_max;
    double Abundance, resultat, B, E_c;
    double f[pt_Spectrum_and_Precision_Parameters->eval_max], g[pt_Spectrum_and_Precision_Parameters->eval_max], z_array[pt_Spectrum_and_Precision_Parameters->eval_max], y;
    double resultat_destruc_nuclei, resultat_source_term ;








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


    Check_nuclei(pt_Scan_Parameters_and_Results);
    Check_nuclei(&Destruction_4He);




    if(pt_Output_Options->BBN_constraints_verbose > 0) {
        cout << "I start generating spectrum for each redshift and convolute them with cross sections. You asked for " << z_step << " points." << endl;
    }
    int end = z_step;
    {
      #pragma omp parallel for ordered schedule(dynamic)
      for(int j = 0; j<=end; j++) {
        compute_constraints_from_destruction_and_production_loop(j,
                pt_Particle_Physics_Model,
                pt_Spectrum_and_Precision_Parameters,
                pt_Scan_Parameters_and_Results,
                &Destruction_4He,
                pt_Output_Options,
                Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei,
                Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei,
                Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He,
                Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He,
                Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei,
                Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei);
      }
    }


    /********************************************************************************************************************************************************************************/
    /****************** Fourth step : Every integral that is independant of zeta_x is done once and for all. The source term involved a double integration over z. *******************/
    /************************************************* To do so, every step in z is stored in the table Integration_over_z_source_term. *********************************************/
    /********************************************************************************************************************************************************************************/

    if(pt_Output_Options->BBN_constraints_verbose > 0) {
        cout <<"I'm done convoluting the spectrum with cross sections! I start to integrate over z."<<endl;
    }
    end = tau_step;
    {
      #pragma omp parallel for ordered schedule(dynamic)

      for(int dtau = 0 ; dtau <= end ; dtau++) {
      Compute_constraints_from_destruction_and_production_loop_2(dtau,
              pt_Particle_Physics_Model,
              pt_Scan_Parameters_and_Results,
              &Destruction_4He,
              pt_Spectrum_and_Precision_Parameters,
              pt_Output_Options,
              Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei,
              Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He,
              Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei,
              Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei,
              Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He,
              Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei);
        }
      }
    print_results_scan(pt_Output_Options,
                       pt_Spectrum_and_Precision_Parameters,
                       pt_Scan_Parameters_and_Results,
                       pt_Particle_Physics_Model);

    Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei.clear();
    Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He.clear();
    Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei.clear();
    Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei.clear();
    Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He.clear();
    Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei.clear();

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

    if(i == 0) {
        double  Q = 2.467032;
        //y = (0.105*2371*pow(x,-2)*exp(-2.5954*pow(x-Q,-0.5))*exp(-2.056*(x-Q))*(1+2.2875*pow(x-Q,2)-1.1798*pow(x-Q,3)+2.5279*pow(x-Q,4)))*pow(x,d);
        if(x>=Q) {
            y = 0.057*931.434*pow(x,-2)*exp(-2.59*pow(x-Q,-0.5));
        }
    } else if(i == 1) {
        double  Q = 7.249962;
        if(x>=Q) {
            y = (0.176*pow(Q,1.51)*pow(x-Q,0.49)*pow(x,-2)+1205*pow(Q,5.5)*pow(x-Q,5)*pow(x,-10.5)+0.06/(1+pow((x-7.46)/0.188,2)));
        }
    } else if(i == 2) {
        double  Q = 10.948850;
        if(x>=Q) {
            y = 122*pow(Q,4)*pow(x-Q,3)*pow(x,-7);
        }
    } else if(i == 3) {
        double  Q0 = 8.725;
        double  Q1 = 23;
        if(x>=Q0) {
            y = 3.8 * pow(Q0,2.3)*(x-Q0)/pow(x,3.3);
        }
        if(x>=Q1) {
            y += 2.1*pow(Q1,1.5)*(x-Q1)/pow(x,2.5);
        }

    } else if(i == 4) {
        double  Q = 9.98;
        if(x>=Q) {
            y = 10.8 *pow(Q,2)*pow(x-Q,1.2)*(x-Q)/pow(x,3.2);
        }
    } else if(i == 5) {
        double  Q = 22.28;
        if(x>=Q) {
            y = 1.44*pow(10.,3)*pow(Q,20)*pow(x-Q,2.4)/pow(x,22.4);
        }
    } else if(i == 6) {
        double  Q = 23.05;
        if(x>=Q) {
            y = 1.44*pow(10.,3)*pow(Q,20)*pow(x-Q,2.4)/pow(x,22.4);
        }
    } else if(i == 7) {
        double  Q = 1.586627;
        //y =  (0.504*2371*pow(x,-2)*exp(-5.1909*pow(x-Q,-0.5))*exp(-0.548*(x-Q))*(1-0.428*pow(x-Q,2)+0.543*pow(x-Q,3)-0.115*pow(x-Q,4)));
        if(x>=Q) {
            y = 0.26*931.434*pow(x,-2)*exp(-5.19*pow(x-Q,-0.5));
        } else {
            y=0;
        }
    } else if(i == 8) {
        double  Q = 5.605794;
        if(x>=Q) {
            y = (32.6*pow(Q,10)*pow((x-Q),2)*pow(x,-12)+2.27*pow(10.,6)*pow(Q,8.8335)*pow((x-Q),13)*pow(x,-21.8335));
        }
    } else if(i == 9) {
        double  Q = 9.30468;
        if(x>=Q) {
            y = 133*pow(Q,4)*pow(x-Q,3)*pow(x,-7);
        }
    } else if(i == 10) {
        double  Q0 = 7.08;
        double  Q1 = 23;
        if(x>=Q0) {
            y = 3.8 * pow(Q0,2.3)*(x-Q0)/pow(x,3.3);
        }
        if(x>=Q1) {
            y += 2.1*pow(Q1,1.5)*(x-Q1)/pow(x,2.5);
        }
    } else if(i == 11) {
        double  Q = 10.68;
        if(x>=Q) {
            y = 10.8 *pow(Q,2)*pow(x-Q,1.2)*(x-Q)/pow(x,3.2);
        }
    } else if(i == 12) {
        double  Q = 22.17;
        if(x>=Q) {
            y = 1.44*pow(10.,3)*pow(Q,20)*pow(x-Q,2.4)/pow(x,22.4);
        }
    } else if(i == 13) {
        double  Q = 21.4;
        if(x>=Q) {
            y = 1.44*pow(10.,3)*pow(Q,20)*pow(x-Q,2.4)/pow(x,22.4);
        }
    } else if(i == 14) {
        double  Q = 2.224573;
        if(x>=Q) {
            y = 18.75*(pow(pow(Q*(x-Q),0.5)/x,3)+0.007947*pow(pow(Q*(x-Q),0.5)/x,2)*pow(pow(Q,0.5)-pow(0.037,0.5),2)/(x-Q+0.037));
        }
    } else if(i == 15) {
        double  Q = 19.813852;
        //~ y = 128.9*pow(Q,4.524)*pow(x-Q,2.512)/pow(x,4.524+2.512);
        if(x>=Q) {
            y = 19.5*pow(Q,3.5)*pow(x-Q,1)/pow(x,4.5);
        }
    } else if(i == 16) {
        double  Q = 20.577615;
        //~ y = 31.68*pow(Q,3.663)*pow(x-Q,1.580)/pow(x,3.663+1.580);
        if(x>=Q) {
            y = 17.1*pow(Q,3.5)*pow(x-Q,1.)/pow(x,4.5);
        }
    } else if(i == 17) {
        double  Q = 23.846527;
        if(x>=Q) {
            y = 10.7*pow(Q,10.2)*pow(x-Q,3.4)/pow(x,13.6);
        }
    } else if(i == 18) {
        double  Q = 26.0711;
        if(x>=Q) {
            y = 21.7*pow(Q,4.0)*pow(x-Q,3.0)/pow(x,7.0);
        }
    } else if(i == 19) {
        double  Q = 5.483485;
        if(x>=Q) {
            y = 8.88*pow(Q,1.75)*pow(x-Q,1.65)/pow(x,3.4);
        }
    } else if(i == 20) {
        double  Q = 7.718058;
        if(x>=Q) {
            y = 16.7*pow(Q,1.95)*pow(x-Q,2.3)/pow(x,4.25);
        }
    }
    y=y*2.569*pow(10.,-6);//Conversion mb en Mev^-2
    return y;

}


static void Check_nuclei(Structure_Scan_Parameters_and_Results * pt_Scan_Parameters_and_Results)
{

    if(pt_Scan_Parameters_and_Results->nuclei =="4He") {
        pt_Scan_Parameters_and_Results->Y_0 = Y_4He_0;
        pt_Scan_Parameters_and_Results->Y_min = Y_4He_Min;
        pt_Scan_Parameters_and_Results->Y_max = Y_4He_Max;
        //Destruction cross-sections
        pt_Scan_Parameters_and_Results->i_min = 15;
        pt_Scan_Parameters_and_Results->i_max = 18;
        //Production cross-sections
        pt_Scan_Parameters_and_Results->k_min = 0;
        pt_Scan_Parameters_and_Results->k_max = 0; //4He cannot be produced in significant quantities by the non-thermal BBN !
    } else if(pt_Scan_Parameters_and_Results->nuclei=="3He") {
        pt_Scan_Parameters_and_Results->Y_0 = Y_3He_0;
        pt_Scan_Parameters_and_Results->Y_min = Y_3He_Min;
        pt_Scan_Parameters_and_Results->Y_max = Y_3He_Max;
        //Destruction cross-sections
        pt_Scan_Parameters_and_Results->i_min = 19;
        pt_Scan_Parameters_and_Results->i_max = 20;
        //Production cross-sections
        pt_Scan_Parameters_and_Results->k_min = 15;
        pt_Scan_Parameters_and_Results->k_max = 16;
    } else if(pt_Scan_Parameters_and_Results->nuclei=="2H") {
        pt_Scan_Parameters_and_Results->Y_0 = Y_2H_0;
        pt_Scan_Parameters_and_Results->Y_min = Y_2H_Min;
        pt_Scan_Parameters_and_Results->Y_max = Y_2H_Max;
        //Destruction cross-sections
        pt_Scan_Parameters_and_Results->i_min = 14;
        pt_Scan_Parameters_and_Results->i_max = 14;
        //Production cross-sections
        pt_Scan_Parameters_and_Results->k_min = 17;
        pt_Scan_Parameters_and_Results->k_max = 18;
    } else if(pt_Scan_Parameters_and_Results->nuclei=="7Li") {
        pt_Scan_Parameters_and_Results->Y_0 = Y_7Li_0;
        pt_Scan_Parameters_and_Results->Y_min = Y_7Li_Min;
        pt_Scan_Parameters_and_Results->Y_max = Y_7Li_Max;
        //Destruction cross-sections
        pt_Scan_Parameters_and_Results->i_min = 0;
        pt_Scan_Parameters_and_Results->i_max = 6;
        //Production cross-sections
        pt_Scan_Parameters_and_Results->k_min = 0;
        pt_Scan_Parameters_and_Results->k_max = 0;
    } else if(pt_Scan_Parameters_and_Results->nuclei=="7Be") {
        pt_Scan_Parameters_and_Results->Y_0 = Y_7Be_0;
        pt_Scan_Parameters_and_Results->Y_min = Y_7Be_Min;
        pt_Scan_Parameters_and_Results->Y_max = Y_7Be_Max;
        //Destruction cross-sections
        pt_Scan_Parameters_and_Results->i_min = 7;
        pt_Scan_Parameters_and_Results->i_max = 13;
        //Production cross-sections
        pt_Scan_Parameters_and_Results->k_min = 0;
        pt_Scan_Parameters_and_Results->k_max = 0;
    } else {
        cout<<"I couldn't recognize nuclei."<<endl;
        cout << "The string 'nuclei' needs to be a word from this list : '4He', '3He', '2H', '7Li, '7Be'." << endl;
        cout << "Please check. " << endl;

        return;
    }

}
