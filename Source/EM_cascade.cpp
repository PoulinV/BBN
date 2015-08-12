#include "../Include/EM_cascade.h"
#include "../Include/Injected_spectrum.h"
#include "../Include/tools.h"
double  dsigma_compton(double  x, double  z, double g){

	double  dsigma = pi*pow(r_e,2)*m_e/pow(x,2)*(x/g+g/x+pow(m_e/g-m_e/x,2)-2*m_e*(1/g-1/x))*n_e*pow(1+z,3);
	return dsigma;

}
double  dsigma_phph(double  x, double  z,  double g){

	double  dsigma = pow(T_0*(1+z),6)*8*pow(pi,4)*1112*pow(a*r_e,2)*pow(m_e,-6)*pow(x,2)*pow(1-g/x+pow(g/x,2),2)/(63*10125*pi);
	return dsigma;

}


double Function_Integrand_Spectre_Compton(double E_gamma, double E_e, double E_gamma_bar){
	double f;
	double Gamma_e = 4*E_gamma_bar*E_e/(m_e*m_e);
	double q = E_gamma/(Gamma_e*(E_e-E_gamma));

	if(q>=0 && q<=1) {f = 2*q*log(q)+(1+2*q)*(1-q)+pow(Gamma_e*q,2)/(2*(1-Gamma_e*q))*(1-q);
}
	else f=0;
	return f;
}

double Function_Integrand_Spectre_Compton_version_q(double q, double E_e, double E_gamma_bar){
	double f;
	double Gamma_e = 4*E_gamma_bar*E_e/(m_e*m_e);
	double E_gamma=Gamma_e*E_e*q/(1+Gamma_e*q);

	if(q>=0 && q<=1) {f = 2*q*log(q)+(1+2*q)*(1-q)+pow(Gamma_e*q,2)/(2*(1-Gamma_e*q))*(1-q)*pow(E_e-E_gamma,2)*Gamma_e/E_e;
}
	else f=0;
	return f;
}
void  Spectre_electron_compton(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
													 struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
													 struct Structure_Electron_Spectrum * pt_Electron_Spectrum,
													 struct Structure_Gamma_Spectrum * pt_Gamma_Spectrum){

 	double z = pt_Electron_Spectrum->redshift;
	int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
	double E_0 = pt_Particle_Physics_Model->E_0;
	double E_gamma_bb = 2.701*T_0*(1+z);
	double E_x = E_x_0/(1+z), E_c = E_c_0/(1+z);
	double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
	double resultat = 0, Gamma_electron = 0;
	double dE, dE_2, dq, h;
	double E1, E2, E3, f1, f2, f3, F1, F2, F3;
	double E_gamma, E_e;
	double E_gamma_1, E_gamma_2, E_gamma_3;
	double q_min=0.0001, q1, q2, q3;
	double resultat_monochromatique_1, resultat_monochromatique_2;
	dE = (pt_Particle_Physics_Model->E_0 - E_min)/ (double) (Electron_Table_Size-1);
	dq = (1-q_min)/ (double) (n_step-1);
	for(int j =0; j<Electron_Table_Size-1;j++){
		resultat=0;
		Gamma_electron=0;
		dE_2 = (pt_Particle_Physics_Model->E_0 - (E_min+j*dE))/ (double) (n_step-1);
		h = dE_2/2;
		for(int i=0;i<n_step-1;i++){
			if(i==0){
				E1=E_min+i*dE_2+j*dE;
				q1=q_min+i*dq;
				E_e = E1;
				}
			else{
				E1=E3;
				E_gamma_1=E_gamma_3;
				q1=q3;
			}

			E2=E1 + h;
			E3=E2 + h;
			// E_gamma_2=E_gamma_1+h;
			// E_gamma_3=E_gamma_2+h;
			q2=q1+dq/2.;
			q3=q2+dq/2.;
			F1 = Function_Integrand_Spectre_Compton_version_q(q1,E_e, E_gamma_bb);
			F2 = Function_Integrand_Spectre_Compton_version_q(q2,E_e, E_gamma_bb);
			F3 = Function_Integrand_Spectre_Compton_version_q(q3,E_e, E_gamma_bb);
			linearint(pt_Gamma_Spectrum->Gamma_Energy, pt_Gamma_Spectrum->Gamma_Spectrum, pt_Gamma_Spectrum->Gamma_Energy.size(), E1, f1);
			linearint(pt_Gamma_Spectrum->Gamma_Energy, pt_Gamma_Spectrum->Gamma_Spectrum, pt_Gamma_Spectrum->Gamma_Energy.size(), E2, f2);
			linearint(pt_Gamma_Spectrum->Gamma_Energy, pt_Gamma_Spectrum->Gamma_Spectrum, pt_Gamma_Spectrum->Gamma_Energy.size(), E3, f3);
			f1*=(dsigma_compton(E1,z,(E1+m_e-E_e)));
			// cout << "f1 = " << f1<<endl;
			f1/=(gamma_NPC(E1,z)+gamma_compton(E1,z)+gamma_phph(E1,z));
			f2*=(dsigma_compton(E2,z,(E2+m_e-E_e)));
			f2/=(gamma_NPC(E2,z)+gamma_compton(E2,z)+gamma_phph(E2,z));
			f3*=(dsigma_compton(E3,z,(E3+m_e-E_e)));
			f3/=(gamma_NPC(E3,z)+gamma_compton(E3,z)+gamma_phph(E3,z));
			// cout << "f1 après= " << f1<<endl;

			// cout << "E = " << E_min+j*dE<<" E1 = " << E1 << endl;
			resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
			Gamma_electron += h*(F1/3+4.*F2/3+F3/3.);
			if(i==n_step-2 && pt_Spectrum_and_Precision_Parameters->spectrum_choice=="Dirac"){
			}
		}
		cout << " resultat = " << resultat << endl;
		E_e = (E_min+j*dE);
		// if(pt_Spectrum_and_Precision_Parameters->spectrum_choice=="Dirac")resultat_monochromatique_1=pi*pow(r_e,2)*m_e/pow(E_0,2)*n_e*pow(1+z,3)*((E_0+m_e-E_e)/E_0+E_0/(E_0+m_e-E_e)+pow(m_e/(E_0+m_e-E_e)-m_e/E_0,2)-2*m_e*(1/(E_0+m_e-E_e)-1/E_0))/(gamma_NPC(E_0,z)+gamma_compton(E_0,z)+gamma_phph(E_0,z));
		// if(pt_Spectrum_and_Precision_Parameters->spectrum_choice=="Dirac")resultat_monochromatique_2=(dsigma_compton(E_0,z,(E_0+m_e-E_e)))/(gamma_NPC(E_0,z)+gamma_compton(E_0,z)+gamma_phph(E_0,z));
		// cout << " resultat_monochromatique_1 = " << resultat_monochromatique_1 << " resultat_monochromatique_2 = " << resultat_monochromatique_2 << endl;
		// cout << " resultat après = " << resultat << endl;

		pt_Electron_Spectrum->Electron_Energy[j]=E_min+j*dE;
		if(Gamma_electron>0)Gamma_electron*=2*pi*r_e*r_e*m_e*m_e*int_bb/(E_gamma_bb*pt_Electron_Spectrum->Electron_Energy[j]*pt_Electron_Spectrum->Electron_Energy[j]);
		else Gamma_electron=0;
		pt_Electron_Spectrum->Electron_Spectrum[j]=resultat;
		// cout << "E_c = " << E_c <<"E = "  << pt_Electron_Spectrum->Electron_Energy[j]<< "resultat = " << pt_Electron_Spectrum->Electron_Spectrum[j] << " Gamma_electron " << Gamma_electron << endl;
		pt_Electron_Spectrum->Electron_Spectrum[j]/=(Gamma_electron);
		// cout<< "resultat après = " << pt_Electron_Spectrum->Electron_Spectrum[j]<<endl;
	}
	/*for(int j =0; j<Electron_Table_Size;j++){
		resultat=0;
		for(int i=0;i<n_step-1;i++){
			if(i==0){
				E1=E_min+i*dE_2+j*dE;
				E_e = E1;
				q1=q_min+i*dq;

				}
			else{
				E1=E3;
			}

			E2=E1 + h;
			E3=E2 + h;
			if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Electron_Energy, pt_Electron_Spectrum->Electron_Spectrum, pt_Electron_Spectrum->Electron_Energy.size(), E1, f1);
			else f1=0;
			if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Electron_Energy, pt_Electron_Spectrum->Electron_Spectrum, pt_Electron_Spectrum->Electron_Energy.size(), E2, f2);
			else f2=0;
			if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Electron_Energy, pt_Electron_Spectrum->Electron_Spectrum, pt_Electron_Spectrum->Electron_Energy.size(), E3, f3);
			else f3=0;
			F1 = Function_Integrand_Spectre_Compton(E_1+E_gamma_bb-E_e,E1, E_gamma_bb);
			F2 = Function_Integrand_Spectre_Compton(E_2+E_gamma_bb-E_e,E2, E_gamma_bb);
			F3 = Function_Integrand_Spectre_Compton(E_3+E_gamma_bb-E_e,E3, E_gamma_bb);
			f1*=2*F1/(E1*E1);
			f2*=2*F2/(E2*E2);
			f3*=2*F3/(E3*E3);
			// cout << "E = " << E_min+j*dE<<" E1 = " << E1 << endl;
			resultat += h * (f1/3. + 4.*f2/3. + f3/3.);

		}

		pt_Gamma_Spectrum->Gamma_Energy[j]=E_min+j*dE;
		pt_Gamma_Spectrum->Gamma_Spectrum[j]=resultat*2*pi*r_e*r_e*m_e*m_e*int_bb/E_gamma_bb;
		cout << "E = "  << E_min+j*dE<< "resultat = " << 	pt_Gamma_Spectrum->Gamma_Spectrum[j] << endl;

	}
*/


}
void Spectre_gamma_compton(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
													 struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
													 struct Structure_Electron_Spectrum * pt_Electron_Spectrum,
													 struct Structure_Gamma_Spectrum * pt_Gamma_Spectrum){

	double z = pt_Gamma_Spectrum->redshift;
	int n_step = pt_Spectrum_and_Precision_Parameters->n_step;
	double E_gamma_bb = 2.701*T_0*(1+z);
	double E_x = E_x_0/(1+z), E_c = E_c_0/(1+z);
	double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
	double resultat = 0;
	double dE, dE_2, h;
	double E1, E2, E3, f1, f2, f3, F1, F2, F3;
	double E_gamma, E_e;

	dE = (pt_Particle_Physics_Model->E_0 - E_min)/ (double) (Gamma_Table_Size-1);
	for(int j =0; j<Gamma_Table_Size;j++){
		resultat=0;
		dE_2 = (pt_Particle_Physics_Model->E_0 - (E_min+j*dE))/ (double) (n_step-1);
		h = dE_2/2;
		for(int i=0;i<n_step-1;i++){
			if(i==0){
				E1=E_min+i*dE_2+j*dE+m_e;
				E_gamma = E1-m_e;
				}
			else{
				E1=E3;
			}

			E2=E1 + h;
			E3=E2 + h;
			if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Electron_Energy, pt_Electron_Spectrum->Electron_Spectrum, pt_Electron_Spectrum->Electron_Energy.size(), E1, f1);
			else f1=0;
			if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Electron_Energy, pt_Electron_Spectrum->Electron_Spectrum, pt_Electron_Spectrum->Electron_Energy.size(), E2, f2);
			else f2=0;
			if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Electron_Energy, pt_Electron_Spectrum->Electron_Spectrum, pt_Electron_Spectrum->Electron_Energy.size(), E3, f3);
			else f3=0;
			F1 = Function_Integrand_Spectre_Compton(E_gamma,E1, E_gamma_bb);
			F2 = Function_Integrand_Spectre_Compton(E_gamma,E2, E_gamma_bb);
			F3 = Function_Integrand_Spectre_Compton(E_gamma,E3, E_gamma_bb);
			f1*=2*F1/(E1*E1);
			f2*=2*F2/(E2*E2);
			f3*=2*F3/(E3*E3);
			// cout << "E = " << E_min+j*dE<<" E1 = " << E1 << endl;
			resultat += h * (f1/3. + 4.*f2/3. + f3/3.);

		}

		pt_Gamma_Spectrum->Gamma_Energy[j]=E_min+j*dE;
		pt_Gamma_Spectrum->Gamma_Spectrum[j]=resultat*2*pi*r_e*r_e*m_e*m_e*int_bb/E_gamma_bb;
		cout << "E = "  << E_min+j*dE<< "resultat = " << 	pt_Gamma_Spectrum->Gamma_Spectrum[j] << endl;

	}
}
double  gamma_compton(double  x, double  z){

	double X = 2*x/m_e;
	double sigma_cs = 2*pi*pow(r_e,2)/X*((1-4/X-8/pow(X,2))*log(1+X)+1/2+8/X+1/(2*pow(1+X,2)));
	double Gamma = sigma_cs*eta*n_y_0*(1+Y/2)/(1+Y)*pow(1+z,3);
	return Gamma;

}
double  gamma_NPC(double  x, double  z){

	double  k = x/m_e;
	double  rho = (2*k-4)/(k+2+2*pow(2*k,0.5));
	double  sigma_PCN;
	sigma_PCN = a*pow(r_e,2)*(28/9*log(2*k)-218/27) ;
	double  Gamma_2 = sigma_PCN*pow(1+z,3)*eta*n_y_0;
	return Gamma_2;

}

double  gamma_phph(double  x, double  z){

	double  Gamma_3;
	Gamma_3 = 0.1513*pow(a,4)*m_e*pow(x/m_e,3)*pow(T_0*(1+z)/m_e,6);
	return Gamma_3;

}

double Dirac_Spectrum_After_One_Iteration(double  x, double  z, double E_0){


	double Gamma_tot = gamma_NPC(E_0,z)+gamma_compton(E_0,z)+gamma_phph(E_0,z);

	double T = T_0*(1+z);
	double int_BB = 8./63.*pow(pi,4)*pow(T_0*(1+z),6);

	double spectre_gamma_gamma = 1112./(10125*pi)*pow(a*r_e,2)*pow(m_e,-6)*pow(E_0,2)*pow(1-x/E_0+pow(x/E_0,2),2)*int_BB;
	double spectre_compton = pi*r_e*r_e*m_e*pow(E_0,-2)*(x/E_0+E_0/x+pow(m_e/x-m_e/E_0,2)-2*m_e*(1/x-1/E_0))*n_e*pow(1+z,3);
	double f = (spectre_gamma_gamma+spectre_compton)/(Gamma_tot);
	if(x>E_0)f=0;
	return f;

}
void Cascade_Spectrum_Reading_From_File(struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
																				struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum,
																				double z,
																				int iterations){
  ostringstream os;
  string name;
	double tmp_Energy, tmp_Spectrum;
  pt_Cascade_Spectrum->redshift = z;

    os << "Cascade_Spectrum_Folder/Spectrum_m" << pt_Particle_Physics_Model->M_x<<"_z"<< z <<"_" << iterations <<"iterations.dat";
    name = os.str();
    ifstream file(name);
    if(file)cout << "Importing file " << name << " in structure Cascade_Spectrum." << endl;
		else{
			cout << "I couldn't recognize cascade spectrum file. Please check that it is in present in the folder Cascade_Specrum_Folder with proper name : 'Spectrum_mXXX_zXXX_XXXiterations.dat' corresponding to the value of m, z and iterations you are using."<<endl;
			return;
		}
    while(file){
			string line;
	    getline(file, line);
	    // stringstream is(ligne);
	    if(line[0] == '#' or line[0] == '\0') continue;
      file >> tmp_Energy >> tmp_Spectrum ;
			pt_Cascade_Spectrum->Gamma_Energy.push_back(tmp_Energy);
			pt_Cascade_Spectrum->Gamma_Spectrum.push_back(tmp_Spectrum*(gamma_NPC(tmp_Energy,z)+gamma_compton(tmp_Energy,z)+gamma_phph(tmp_Energy,z)));
			  // cout << z << "  " << tmp_Energy << "  " << tmp_Spectrum<<endl;

		}

		file.close();


}

void  Cascade_Spectrum_Calculation(double z,
																	 struct Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
																	 struct Structure_Gamma_Spectrum * pt_Cascade_Spectrum,
																	 struct Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){
	double dE, dE_2, h;
	double E1, E2, E3, E_gamma, Old_spectrum;
	double f1, f2, f3;
	double resultat;
	double E_c = E_c_0/(1+z);
	struct Structure_Gamma_Spectrum Cascade_Spectrum_Old;
	pt_Cascade_Spectrum->redshift = z;
	dE = (pt_Particle_Physics_Model->E_0 - E_min)/ (double) (Gamma_Table_Size-1);
	/*****Initialization****/
	pt_Cascade_Spectrum->Gamma_Energy.resize(Gamma_Table_Size);
	pt_Cascade_Spectrum->Gamma_Spectrum.resize(Gamma_Table_Size);
	Cascade_Spectrum_Old.Gamma_Energy.resize(Gamma_Table_Size);
	Cascade_Spectrum_Old.Gamma_Spectrum.resize(Gamma_Table_Size);
	if(pt_Spectrum_and_Precision_Parameters->spectrum_choice == "universal"){
		for(int i=0;i<Gamma_Table_Size;i++){
			E1=E_min+i*dE;
			pt_Cascade_Spectrum->Gamma_Energy[i]=E1;
			pt_Cascade_Spectrum->Gamma_Spectrum[i]=Universal_Spectrum(E1,z,pt_Particle_Physics_Model->E_0);
		}
	}
	else if(pt_Spectrum_and_Precision_Parameters->spectrum_choice == "Dirac"){
		if(E_c <= pt_Particle_Physics_Model->E_0 ){

			for(int i=0;i<Gamma_Table_Size;i++){
				E1=E_min+i*dE;
				pt_Cascade_Spectrum->Gamma_Energy[i]=E1;
				pt_Cascade_Spectrum->Gamma_Spectrum[i]=Universal_Spectrum(E1,z,pt_Particle_Physics_Model->E_0);
			}
			if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
					if(verbose>1)cout <<" I will now print the spectrum in files." << endl;
					print_spectrum_automatic_names(0, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
			}
		}


		else{

			for(int i=0;i<Gamma_Table_Size;i++){
				E1=E_min+i*dE;
				pt_Cascade_Spectrum->Gamma_Energy[i]=E1;
				pt_Cascade_Spectrum->Gamma_Spectrum[i]=Dirac_Spectrum_After_One_Iteration(E1,z,pt_Particle_Physics_Model->E_0);

			}
			if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
					if(verbose>1)cout <<" I will now print the spectrum in files." << endl;
					print_spectrum_automatic_names(0, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
			}

			for(int k = 0; k<pt_Spectrum_and_Precision_Parameters->iterations;k++){
							if(verbose>1)cout<<"iteration : " << k+1 << endl;
							// for(int l = 0 ; l < Gamma_Table_Size ; l++)
							// {
							// 		Cascade_Spectrum_Old.Gamma_Energy[l] = pt_Cascade_Spectrum->Gamma_Energy[l];
							// 		Cascade_Spectrum_Old.Gamma_Spectrum[l] = pt_Cascade_Spectrum->Gamma_Spectrum[l];
							// }
							for(int j =0; j<Gamma_Table_Size;j++){
								resultat=0;
								dE_2 = (pt_Particle_Physics_Model->E_0 - (E_min+j*dE))/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
								h = dE_2/2;
								for(int i=0;i<pt_Spectrum_and_Precision_Parameters->n_step-1;i++){
									if(i==0){
										E1=E_min+j*dE+i*dE_2;
										E_gamma = E1;
										}
									else{
										E1=E3;
									}

									E2=E1 + h;
									E3=E2 + h;

									linearint(pt_Cascade_Spectrum->Gamma_Energy, pt_Cascade_Spectrum->Gamma_Spectrum, pt_Cascade_Spectrum->Gamma_Energy.size(), E1, f1);
									linearint(pt_Cascade_Spectrum->Gamma_Energy, pt_Cascade_Spectrum->Gamma_Spectrum, pt_Cascade_Spectrum->Gamma_Energy.size(), E2, f2);
									linearint(pt_Cascade_Spectrum->Gamma_Energy, pt_Cascade_Spectrum->Gamma_Spectrum, pt_Cascade_Spectrum->Gamma_Energy.size(), E3, f3);
									// linearint(Cascade_Spectrum_Old.Gamma_Energy, Cascade_Spectrum_Old.Gamma_Spectrum, Cascade_Spectrum_Old.Gamma_Energy.size(), E1, f1);
									// linearint(Cascade_Spectrum_Old.Gamma_Energy, Cascade_Spectrum_Old.Gamma_Spectrum, Cascade_Spectrum_Old.Gamma_Energy.size(), E2, f2);
									// linearint(Cascade_Spectrum_Old.Gamma_Energy, Cascade_Spectrum_Old.Gamma_Spectrum, Cascade_Spectrum_Old.Gamma_Energy.size(), E3, f3);
									f1 *=(dsigma_phph(E1,z,E_gamma)+dsigma_compton(E1,z,E_gamma))/(gamma_NPC(E1,z)+gamma_compton(E1,z)+gamma_phph(E1,z));
									f2 *=(dsigma_phph(E2,z,E_gamma)+dsigma_compton(E2,z,E_gamma))/(gamma_NPC(E2,z)+gamma_compton(E2,z)+gamma_phph(E2,z));
									f3 *=(dsigma_phph(E3,z,E_gamma)+dsigma_compton(E3,z,E_gamma))/(gamma_NPC(E3,z)+gamma_compton(E3,z)+gamma_phph(E3,z));
									// cout << "E = " << E_min+j*dE<<" E1 = " << E1 << endl;
									resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
								}
								pt_Cascade_Spectrum->Gamma_Energy[j]=E_min+j*dE;
								pt_Cascade_Spectrum->Gamma_Spectrum[j]=resultat+Dirac_Spectrum_After_One_Iteration(E_min+j*dE,z,pt_Particle_Physics_Model->E_0);

							}
							if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
								print_spectrum_automatic_names(k+1, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
							}

		}
		}
		if(pt_Spectrum_and_Precision_Parameters->inverse_compton_scattering == "yes"){
			struct Structure_Gamma_Spectrum Inverse_Compton_Spectrum;
			struct Structure_Electron_Spectrum Electron_Spectrum;

			Inverse_Compton_Spectrum.Gamma_Energy.resize(Gamma_Table_Size);
			Inverse_Compton_Spectrum.Gamma_Spectrum.resize(Gamma_Table_Size);
			Electron_Spectrum.Electron_Energy.resize(Electron_Table_Size);
			Electron_Spectrum.Electron_Spectrum.resize(Electron_Table_Size);
			Inverse_Compton_Spectrum.redshift=z;
			Electron_Spectrum.redshift=z;
			Spectre_electron_compton(pt_Particle_Physics_Model,
															 pt_Spectrum_and_Precision_Parameters,
															 &Electron_Spectrum,
															 pt_Cascade_Spectrum);
			Spectre_gamma_compton( pt_Particle_Physics_Model,
															 pt_Spectrum_and_Precision_Parameters,
															 &Electron_Spectrum,
															 &Inverse_Compton_Spectrum);
			for(int i = 0; i<Gamma_Table_Size ; i++){
				pt_Cascade_Spectrum->Gamma_Spectrum[i]+=Inverse_Compton_Spectrum.Gamma_Spectrum[i];
			}
			if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
				print_spectrum_automatic_names(100, pt_Cascade_Spectrum, pt_Particle_Physics_Model);
			}
			if(pt_Spectrum_and_Precision_Parameters->spectrum_mode == "writing"){
				print_spectrum_automatic_names(101, &Inverse_Compton_Spectrum, pt_Particle_Physics_Model);
			}
		}
	}



}
