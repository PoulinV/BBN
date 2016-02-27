#include "bbn/EM_cascade.h"
#include "bbn/injected_spectrum.h"
#include "bbn/structures.h"
#include "bbn/BBN_constraints.h"
#include "bbn/tools.h"
#include "bbn/test_functions.h"

using namespace std;

void print_spectrum_from_function(ostream &file, double (*func)(double,double,double),double z, Structure_Particle_Physics_Model * pt_Particle_Physics_Model,Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){


  double dE = (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table)/pt_Spectrum_and_Precision_Parameters->Energy_Table_Size;
  double E = pt_Spectrum_and_Precision_Parameters->E_min_table;
  int i = 0;
  while(i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size){
    file << E << "   " << (*func)(E,z,pt_Particle_Physics_Model->E_0) << endl;
    i++;
    E+=dE;
  }

}
void print_spectrum(Structure_Output_Options * pt_Output_Options,
                    Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                    Structure_Spectrum * pt_Spectrum,
                    Structure_Particle_Physics_Model * pt_Particle_Physics_Model){
  ostringstream os;
  string name;


  double dE = (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table)/pt_Spectrum->Energy.size();
  double z = pt_Spectrum->redshift;
  double E = pt_Spectrum_and_Precision_Parameters->E_min_table;
  int i = 0;

    if(pt_Output_Options->spectrum_files=="automatic"){
    if(pt_Spectrum->species == "photon"){
      if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "universal"){
          os << "Output/Cascade_Spectrum_Folder/Spectrum_"<<pt_Spectrum->spectrum_name<<"m" << pt_Particle_Physics_Model->M_x<<"_z"<< z <<".dat";
      }
      else{
        if(pt_Spectrum_and_Precision_Parameters->calculation_mode=="iterative")  os << "Output/Cascade_Spectrum_Folder/Spectrum_"<<pt_Spectrum->spectrum_name<<"m" << pt_Particle_Physics_Model->M_x<<"_z"<< z <<"_" << pt_Spectrum_and_Precision_Parameters->number_iterations_photon <<"iterations.dat";
        else if(pt_Spectrum_and_Precision_Parameters->calculation_mode=="triangular")  os << "Output/Cascade_Spectrum_Folder/Spectrum_"<<pt_Spectrum->spectrum_name<<"m" << pt_Particle_Physics_Model->M_x<<"_z"<< z <<"_" << "triangular.dat";
      }

    }
    else if(pt_Spectrum->species == "electron"){
      if(pt_Spectrum_and_Precision_Parameters->calculation_mode=="iterative")  os << "Output/Cascade_Spectrum_Folder/Spectrum_"<<pt_Spectrum->spectrum_name<<"m" << pt_Particle_Physics_Model->M_x<<"_z"<< z <<"_" << pt_Spectrum_and_Precision_Parameters->number_iterations_electron <<"iterations.dat";
      else if(pt_Spectrum_and_Precision_Parameters->calculation_mode=="triangular")  os << "Output/Cascade_Spectrum_Folder/Spectrum_"<<pt_Spectrum->spectrum_name<<"m" << pt_Particle_Physics_Model->M_x<<"_z"<< z <<"_" << "triangular.dat";
    }

    }

    else{
      os << pt_Output_Options->spectrum_files <<pt_Spectrum->species << ".dat";
    }
    name = os.str();

    ofstream file(name);
    cout << "Printing in file " << name <<"."<<  endl;

    if(pt_Spectrum->species=="photon"){
    while(i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size){
      file << pt_Spectrum->Energy[i] << "   " << pt_Spectrum->Spectrum[i] << endl;
      i++;
    }
  }
    if(pt_Spectrum->species=="electron"){
      while(i<pt_Spectrum_and_Precision_Parameters->Energy_Table_Size){
        file << pt_Spectrum->Energy[i] << "   " << pt_Spectrum->Spectrum[i] << endl;
        i++;
      }
    }

}
void print_results_scan(Structure_Output_Options * pt_Output_Options,
                    Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                    Structure_Scan_Parameters_and_Results * pt_Scan_Parameters_and_Results,
                    Structure_Particle_Physics_Model * pt_Particle_Physics_Model){
  ostringstream os;
  string name;


  int i = 0;
  // os << "Output/Result_Scan_Folder/Results_destruc_and_production_"<< pt_Scan_Parameters_and_Results->nuclei << "_m"<<pt_Particle_Physics_Model->M_x<<"MeV.dat";

    if(pt_Output_Options->results_files=="automatic"){
      if(pt_Output_Options->task=="compute_constraints_from_destruction_and_production"){
        os << "Output/Result_Scan_Folder/Results_destruc_and_production_"<< pt_Scan_Parameters_and_Results->nuclei << "_m"<<pt_Particle_Physics_Model->M_x<<"MeV_"<<pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice<<".dat";
      }
      else if(pt_Output_Options->task=="compute_constraints_from_destruction_only"){
        os << "Output/Result_Scan_Folder/Results_destruc_only_"<< pt_Scan_Parameters_and_Results->nuclei << "_m"<<pt_Particle_Physics_Model->M_x<<"MeV_"<<pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice<<".dat";
      }
    }

    else{
      os << pt_Output_Options->results_files << ".dat";
    }
    name = os.str();

    ofstream file(name);
    if(file){
      cout << "Printing in file " << name <<"."<< endl;
      file<< "#Scan parameters : tau in [" <<pt_Scan_Parameters_and_Results->tau_min<<","<<pt_Scan_Parameters_and_Results->tau_max<<"] and zeta in ["<<pt_Scan_Parameters_and_Results->zeta_min<<","<<pt_Scan_Parameters_and_Results->zeta_max<<"]"<<endl;
      file <<"#You have injected photon :" <<pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice <<"and electron :" << pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice << " and used the " <<pt_Spectrum_and_Precision_Parameters->calculation_mode<<" method."<<endl;
      if(pt_Spectrum_and_Precision_Parameters->calculation_mode=="iterative") file<<"#You have asked for "<<pt_Spectrum_and_Precision_Parameters->number_iterations_photon<<" iterations for the photons and " <<pt_Spectrum_and_Precision_Parameters->number_iterations_electron << "iterations for the electrons."<<endl;
        while(i<pt_Scan_Parameters_and_Results->Results_scan_tau_x.size()){
          file << pt_Scan_Parameters_and_Results->Results_scan_tau_x[i] << "   " << pt_Scan_Parameters_and_Results->Results_scan_zeta_x[i] << "   " << pt_Scan_Parameters_and_Results->Results_scan_Abundance[i]<<endl;
          i++;
      }
    }
    else{
      cout<<" I couldn't open file, i'm shutting down." << endl;
      exit(0);
    }


}
void linearint(vector<double> &xa, vector<double> &ya, int n, double x, double &y){
				int i,m,ns=1;
				float x0,x1,y0,y1;
				double dift,dif;
				if(isnan(xa[0])==0){
          dif=sqrt(pow(x-xa[0],2));
          y = ya[0];
          // cout <<" x " << xa[0] << " y = " << y << endl;
        }
        else{
          cout << "careful, you're initial point isn't defined !" << endl;
          return;
        }
				// cout << "dif = " << dif << endl;
				for (i=1;i<n;i++) {
						if ( (dift=sqrt(pow(x-xa[i],2))) < dif) {
						ns=i;
						dif=dift;
						y0=ya[i];
						x0=xa[i];
						y=ya[i];
						}
            // cout << "dift = " << dift << " x " << x << " xa[i] " << xa[i] << endl;

            // cout << "here"<<endl;
						// if(ns!=i)break;
				}
        if(i<n-1){
          x1=xa[ns+1];
          y1=ya[ns+1];
          // if(x1!=x0){
            y=y0+(x-x0)/(x1-x0)*(y1-y0);
          cout << " y = " << y << " x - x0 = " << x-x0 << endl;
        // }
          // else y = y0;
          if(isnan(y)==1){
            cout << "y is a nan ! I automatically put it to 0." << endl;
          y=0;
          }
        }
}

void fill_structure_particle_physics_model(ifstream &file, map_parameters &map_parameters, Structure_Particle_Physics_Model * pt_Particle_Physics_Model){
  int i=0;
  string name = "", value = "";
  string error_name="",error_value="";
  if(file){cout << "(function : fill_structure_particle_physics_model) : Reading input parameters file."<< endl;}
  else{
    cout << "(function : fill_structure_particle_physics_model) : \nI couldn't recognized input parameter file. Please check that it is present in the same folder as the executable."<<endl;
    return;
  }
  while(file){
    string line;
    name = "", value = "";
    error_name="",error_value="";
    getline(file, line);
    i++;
    attribute_name_and_value(line,name,value);
    check_value_and_name_error(name,error_name,value,error_value);
    if(error_value == "yes" || error_name == "yes"){
      cout << "Problem at the line "<<i<<" of your 'input_param.ini' file. I'm shutting down."<<endl;
      exit(0);
    }
    else{
      map_parameters[name]=value;
    }
  }

  pt_Particle_Physics_Model->M_x = atof(map_parameters["m_x"].c_str());
	pt_Particle_Physics_Model->E_0 = atof(map_parameters["m_x"].c_str())/2.;
	pt_Particle_Physics_Model->zeta_x = atof(map_parameters["zeta_x"].c_str());
	pt_Particle_Physics_Model->tau_x = atof(map_parameters["tau_x"].c_str());
	pt_Particle_Physics_Model->T_x	= pow(atof(map_parameters["tau_x"].c_str())*(2*H_r),-0.5)*T_0;
	pt_Particle_Physics_Model->z_x = 	pow(atof(map_parameters["tau_x"].c_str())*(2*H_r),-0.5)-1;
  file.clear();
  file.seekg(0, ios::beg);
}

void fill_structure_scan_parameters_and_results(ifstream &file, map_parameters &map_parameters, Structure_Scan_Parameters_and_Results * pt_Scan_Parameters_and_Results){
  int i=0;
  string name = "", value = "";
  string error_name="",error_value="";
  if(file){cout << "(function : fill_structure_scan_parameters) : Reading input parameters file."<< endl;}
  else{
    cout << "(function : fill_structure_scan_parameters) : \nI couldn't recognized input parameter file. Please check that it is present in the same folder as the executable."<<endl;
    return;
  }
  while(file){
    string line;
    name = "", value = "";
    error_name="",error_value="";
    getline(file, line);
    i++;
    attribute_name_and_value(line,name,value);
    check_value_and_name_error(name,error_name,value,error_value);
    if(error_value == "yes" || error_name == "yes"){
      cout << "Problem at the line "<<i<<" of your 'input_param.ini' file. I'm shutting down."<<endl;
      exit(0);
    }
    else{
      map_parameters[name]=value;
    }
  }
  pt_Scan_Parameters_and_Results->nuclei = map_parameters["nuclei"];
  pt_Scan_Parameters_and_Results->tau_min = atof(map_parameters["tau_min"].c_str());
	pt_Scan_Parameters_and_Results->tau_max = atof(map_parameters["tau_max"].c_str());
	pt_Scan_Parameters_and_Results->tau_step = atof(map_parameters["tau_step"].c_str());
  pt_Scan_Parameters_and_Results->zeta_min = atof(map_parameters["zeta_min"].c_str());
	pt_Scan_Parameters_and_Results->zeta_max = atof(map_parameters["zeta_max"].c_str());
	pt_Scan_Parameters_and_Results->zeta_step = atof(map_parameters["zeta_step"].c_str());
  file.clear();
  file.seekg(0, ios::beg);
}

void fill_structure_spectrum_and_precision_parameters(ifstream &file, map_parameters &map_parameters,
                                                      Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){

  map_spectrum map_spectrum;
  attribute_map_spectrum(map_spectrum,map_parameters);
  int i=0;
  string name = "", value = "";
  string error_name="",error_value="";
  if(file){cout << "(function : fill_structure_spectrum_and_precision_parameters) : Reading input parameters file."<< endl;}
  else{
    cout << "(function : fill_structure_spectrum_and_precision_parameters) : \nI couldn't recognized input parameter file. Please check that it is present in the same folder as the executable."<<endl;
    return;
  }
  while(file){
    string line;
    name = "", value = "";
    error_name="",error_value="";
    getline(file, line);
    i++;
    attribute_name_and_value(line,name,value);
    check_value_and_name_error(name,error_name,value,error_value);
    if(error_value == "yes" || error_name == "yes"){
      cout << "Problem at the line "<<i<<" of your 'input_param.ini' file. I'm shutting down."<<endl;
      exit(0);
    }
    else{
      map_parameters[name]=value;
    }
  }
  pt_Spectrum_and_Precision_Parameters->calculation_mode = map_parameters["calculation_mode"];
	pt_Spectrum_and_Precision_Parameters->number_iterations_photon = atoi(map_parameters["number_iterations_photon"].c_str());
  pt_Spectrum_and_Precision_Parameters->number_iterations_electron = atoi(map_parameters["number_iterations_electron"].c_str());
	pt_Spectrum_and_Precision_Parameters->z_step = atoi(map_parameters["z_step"].c_str());
	pt_Spectrum_and_Precision_Parameters->n_step = atoi(map_parameters["n_step"].c_str());
  pt_Spectrum_and_Precision_Parameters->Energy_Table_Size = atof(map_parameters["Energy_Table_Size"].c_str());
  pt_Spectrum_and_Precision_Parameters->E_min_table = atof(map_parameters["E_min_table"].c_str());
  pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice = map_parameters["photon_spectrum_choice"];
  pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice = map_parameters["electron_spectrum_choice"];
  pt_Spectrum_and_Precision_Parameters->photon_spectrum_file_name = map_parameters["photon_spectrum_file_name"];
  pt_Spectrum_and_Precision_Parameters->electron_spectrum_file_name = map_parameters["electron_spectrum_file_name"];
  pt_Spectrum_and_Precision_Parameters->spectrum_mode = map_parameters["spectrum_mode"];
  pt_Spectrum_and_Precision_Parameters->inverse_compton_scattering = map_parameters["inverse_compton_scattering"];
  pt_Spectrum_and_Precision_Parameters->double_photon_pair_creation = map_parameters["double_photon_pair_creation"];
  pt_Spectrum_and_Precision_Parameters->pair_creation_in_nuclei = map_parameters["pair_creation_in_nuclei"];
  pt_Spectrum_and_Precision_Parameters->compton_scattering = map_parameters["compton_scattering"];
  pt_Spectrum_and_Precision_Parameters->photon_photon_diffusion = map_parameters["photon_photon_diffusion"];
  pt_Spectrum_and_Precision_Parameters->check_energy_conservation = map_parameters["check_energy_conservation"];
  pt_Spectrum_and_Precision_Parameters->integration_method = map_parameters["integration_method"];


  if(map_parameters["calculation_mode"] == "triangular" && map_parameters["photon_spectrum_choice"] == "Dirac"){
    pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum = no_photons_injected;
  }
  else if(map_parameters["calculation_mode"] == "iterative" && map_parameters["photon_spectrum_choice"] == "Dirac"){
    pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum = no_photons_injected;
    // pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum = photon_dirac_spectrum_after_one_iteration;
  }
  else{
    if(map_parameters["photon_spectrum_choice"]=="none")pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum = no_photons_injected;
    else pt_Spectrum_and_Precision_Parameters->Injected_Gamma_Spectrum = map_spectrum[map_parameters["Gamma_Spectrum"]];
  }
  if(map_parameters["calculation_mode"] == "triangular" && map_parameters["electron_spectrum_choice"] == "Dirac"){
    pt_Spectrum_and_Precision_Parameters->Injected_Electron_Spectrum = no_electrons_injected;
  }
  else if(map_parameters["calculation_mode"] == "iterative" && map_parameters["electron_spectrum_choice"] == "Dirac"){
    pt_Spectrum_and_Precision_Parameters->Injected_Electron_Spectrum = no_electrons_injected;
    // pt_Spectrum_and_Precision_Parameters->Injected_Electron_Spectrum = electron_dirac_spectrum_after_one_iteration;
  }
  else{
    if(map_parameters["electron_spectrum_choice"]=="none")pt_Spectrum_and_Precision_Parameters->Injected_Electron_Spectrum = no_electrons_injected;
    else pt_Spectrum_and_Precision_Parameters->Injected_Electron_Spectrum = map_spectrum[map_parameters["Electron_Spectrum"]];
  }

  if(map_parameters["integration_method"] == "simpson_1/3"){
    pt_Spectrum_and_Precision_Parameters->eval_max = 3;
    pt_Spectrum_and_Precision_Parameters->divisor = 6;
    pt_Spectrum_and_Precision_Parameters->weight.push_back(1);
    pt_Spectrum_and_Precision_Parameters->weight.push_back(4);
    pt_Spectrum_and_Precision_Parameters->weight.push_back(1);
  }
  else if(map_parameters["integration_method"] == "simpson_3/8"){
    pt_Spectrum_and_Precision_Parameters->eval_max = 4;
    pt_Spectrum_and_Precision_Parameters->divisor = 8;
    pt_Spectrum_and_Precision_Parameters->weight.push_back(1);
    pt_Spectrum_and_Precision_Parameters->weight.push_back(3);
    pt_Spectrum_and_Precision_Parameters->weight.push_back(3);
    pt_Spectrum_and_Precision_Parameters->weight.push_back(1);
  }
  else if(map_parameters["integration_method"] == "weddle_hardy"){
    pt_Spectrum_and_Precision_Parameters->eval_max = 7;
    pt_Spectrum_and_Precision_Parameters->divisor = 840;
    pt_Spectrum_and_Precision_Parameters->weight.push_back(41);
    pt_Spectrum_and_Precision_Parameters->weight.push_back(216);
    pt_Spectrum_and_Precision_Parameters->weight.push_back(27);
    pt_Spectrum_and_Precision_Parameters->weight.push_back(272);
    pt_Spectrum_and_Precision_Parameters->weight.push_back(27);
    pt_Spectrum_and_Precision_Parameters->weight.push_back(216);
    pt_Spectrum_and_Precision_Parameters->weight.push_back(41);

  }



  file.clear();
  file.seekg(0, ios::beg);
}

void fill_structure_output_options(ifstream &file, map_parameters &map_parameters, Structure_Output_Options * pt_Structure_Output_Options){

  int i=0;
  string name = "", value = "";
  string error_name="",error_value="";
  if(file){cout << "(function : fill_output_options) : Reading input parameters file."<< endl;}
  else{
    cout << "(function : fill_output_options) : \nI couldn't recognized input parameter file. Please check that it is present in the same folder as the executable."<<endl;
    return;
  }
  while(file){
    string line;
    name = "", value = "";
    error_name="",error_value="";
    getline(file, line);
    i++;
    attribute_name_and_value(line,name,value);
    check_value_and_name_error(name,error_name,value,error_value);
    if(error_value == "yes" || error_name == "yes"){
      cout << "Problem at the line "<<i<<" of your 'input_param.ini' file. I'm shutting down."<<endl;
      exit(0);
    }
    else{
      map_parameters[name]=value;
    }
  }
  pt_Structure_Output_Options->results_files = map_parameters["results_files"];
  pt_Structure_Output_Options->spectrum_files = map_parameters["spectrum_files"];
  pt_Structure_Output_Options->interaction_rate_files = map_parameters["interaction_rate_files"];
  pt_Structure_Output_Options->test_files = map_parameters["test_files"];
  pt_Structure_Output_Options->task = map_parameters["task"];
  pt_Structure_Output_Options->task_test = map_parameters["task_test"];
  pt_Structure_Output_Options->EM_cascade_verbose = atoi(map_parameters["EM_cascade_verbose"].c_str());
  pt_Structure_Output_Options->BBN_constraints_verbose = atoi(map_parameters["BBN_constraints_verbose"].c_str());
  pt_Structure_Output_Options->Input_verbose = atoi(map_parameters["Input_verbose"].c_str());
  pt_Structure_Output_Options->Output_verbose = atoi(map_parameters["Output_verbose"].c_str());
  file.clear();
  file.seekg(0, ios::beg);
}

void fill_default_parameters(ifstream &file, map_parameters &map_parameters){

  int i=0;
  string name = "", value = "";
  string error_name="",error_value="";
  string line;
  if(file){cout << "(function : fill_default_parameters) : Reading default parameters file."<< endl;}
  else{
    cout << "(function : fill_default_parameters) : \nI couldn't recognized default parameter file. Please check that it is present in the same folder as the executable."<<endl;
    return;
  }
  while(file){
    name = "", value = "";
    error_name="",error_value="";
    getline(file, line);
    i++;
    attribute_name_and_value(line,name,value);
    check_value_and_name_error(name,error_name,value,error_value);
    if(error_value == "yes" || error_name == "yes"){
      cout << "Problem at the line "<<i<<" of your 'default_param.ini' file. I'm shutting down."<<endl;
      exit(0);
    }
    else{
        if(name!="ignore_line" && value!="ignore_line")map_parameters[name]=value;
    }
  }
  file.clear();
  file.seekg(0, ios::beg);
}
void get_parameter_from_file(ifstream &file, string &parameter){
  string error_parameter;
  check_parameter_error(parameter,error_parameter);
  if(error_parameter=="yes"){
    cout << "I couldn't reckognised parameter " << parameter << " please check it's spelling." << endl;
    return;
  }
  int i=0;
  string name = "", value ="";
  string error_name="",error_value="";
    if(file){cout << "(function : get_parameter_from_file) : Attribute " << parameter << " parameter."<< endl;}
  else{
    cout << "(function : get_parameter_from_file) : \nI couldn't recognized input parameter file. Please check that it is present in the same folder as the executable."<<endl;
    return;
  }

  while(file){
    string line;
    name = "", value = "";
    error_name="",error_value="";
    getline(file, line);
    i++;
    if(name!=parameter){
    attribute_name_and_value(line,name,value);
    check_value_and_name_error(name,error_name,value,error_value);
    }
    if(error_value == "yes" || error_name == "yes"){
      cout << "Problem at the line "<<i<<" of your 'default_param.ini' file. I'm shutting down."<<endl;
      exit(0);
    }

    if(name==parameter){
      parameter = value;
      break;
    }
    if(file.eof() && name!=parameter){
      cout << "I have not found the parameter " << parameter <<" in your input file. I will now attribute the default one." << endl;
      parameter = "default";
    }
  }
  file.clear();
  file.seekg(0, ios::beg);

}
void attribute_name_and_value(const string &line,string &name,string &value){
  int i=0;
  int j=0;
  while(line[i]!='='&&i<=line.size())i++;
  while(line[j]!='#'&&j<=line.size())j++;

  if(i<j){
    for(int k=0;k<i;k++){
      name+=line[k];
    }
    name.erase(remove(name.begin(),name.end(),' '));
    for(int k=i+1;k<j-1;k++){
      value+=line[k];
      value.erase(remove(value.begin(),value.end(),' '));
    }
  }
  else if(i>=j){
    name = "ignore_line";
    value = "ignore_line";
  }

}
void check_parameter_error(const string &parameter, string &error_parameter){
      if(parameter == "calculation_mode" || parameter == "task" || parameter == "task_test" ||parameter == "number_iterations_photon" || parameter == "number_iterations_electron" || parameter == "z_step" || parameter == "n_step"
      || parameter == "photon_spectrum_choice" || parameter == "electron_spectrum_choice" || parameter == "spectrum_mode" || parameter == "inverse_compton_scattering" || parameter == "m_x" || parameter == "tau_x"
      || parameter == "zeta_x" || parameter == "tau_min" || parameter == "tau_max" || parameter == "tau_step" || parameter == "zeta_min" || parameter == "zeta_max" || parameter == "zeta_step" || parameter == "temperature" || parameter == "redshift" || parameter == "results_file" || parameter == "spectrum_files"){
        error_parameter = "no";
      }
      else error_parameter = "yes";
}
void check_value_and_name_error(string &name,string &error_name, string &value,string &error_value){

    char * pEnd;
    // cout << "value =_"<< value << "_ name =_" << name<<"_"<<endl;
    error_name = "no";
    if(name == "calculation_mode"){
      if(value == "iterative" || value == "triangular"){
        error_value = "no";
        error_name =  "no";
      }
      else{
        error_value = "yes";
        cout << "The calculation_mode isn't reckognised, it has to be one among : 'iterative' and 'triangular'." << endl;
      }
    }
    else if(name == "task"){
      if(value == "print_interaction_rate" || value == "compute_cascade_spectrum" || value == "compute_constraints_from_destruction_only" || value == "compute_constraints_from_destruction_and_production"){
        error_value = "no";
        error_name = "no";

      }
      else{
        error_value = "yes";
        cout << "The task isn't reckognised, it has to be one among : 'compute_constraints_from_destruction_and_production', 'compute_constraints_from_destruction_only' and 'compute_cascade_spectrum'." << endl;
      }
    }
    else if(name == "task_test"){
      if(value == "print_interaction_rate" || value == "print_polylog"){
        error_value = "no";
        error_name = "no";

      }
      else{
        error_value = "yes";
        cout << "The task_test isn't reckognised, it has to be one among : 'print_interaction_rate', 'print_polylog'." << endl;
      }
    }
    else if(name == "number_iterations_photon"){
      if(atoi(value.c_str())>20){
        error_value = "yes";
        cout << "The number of iterations for the photon spectrum is too big. Please choose a number <= 20." << endl;
      }
      else {
        error_value = "no";
        error_name =  "no";
      }
    }
    else if(name == "number_iterations_electron"){
      if(atoi(value.c_str())>100){
        error_value = "yes";
        cout << "The number of iterations for the photon spectrum is too big. Please choose a number <= 100." << endl;
      }
      else {
        error_value = "no";
        error_name =  "no";
      }
    }
    else if(name == "z_step"){
      if(atoi(value.c_str())>200){
        error_value = "yes";
        cout << "The number z_step is too big. Please choose a number <= 200." << endl;
      }
      else {
        error_value = "no";
        error_name =  "no";
      }
    }
    else if(name == "n_step"){
      if(atoi(value.c_str())>2000){
        error_value = "yes";
        cout << "The number n_step is too big. Please choose a number <= 2000." << endl;
      }
      else {
        error_value = "no";
        error_name =  "no";
      }
    }
    else if(name == "integration_method"){
      if(value=="simpson_1/3" || value=="simpson_3/8" || value=="weddle_hardy"){
        error_value = "no";
      }
      else {
        cout << "Please choose an integration method between 'simpson_1/3', 'simpson_3/8' and 'weddle_hardy'." << endl;
        error_value = "yes";
      }
    }
    else if(name == "photon_spectrum_choice"){
      check_name_spectrum(value, error_value);
      if(error_value=="yes"){
        cout << "The name of the injected photon spectrum is not reckognised. Please check that it is one among : 'Dirac', 'none', universal' or user specified 'name_of_your_function'."<<endl;
        cout << "The function should have been correctly added in 'injected_spectrum.cpp' and 'injected_spectrum.h' before."<<endl;
      }
    }
    else if(name == "electron_spectrum_choice"){
      check_name_spectrum(value, error_value);
      if(error_value=="yes"){
        cout << "The name of the injected electron spectrum is not reckognised. Please check that it is one among : 'Dirac', 'none', universal' or user specified 'name_of_your_function'."<<endl;
        cout << "The function should have been correctly added in 'injected_spectrum.cpp' and 'injected_spectrum.h' before."<<endl;
      }
    }
    else if(name == "check_energy_conservation"){
      if(value == "yes" || value == "no")error_value="no";
      else{
        cout << "check_energy_conservation : I couldn't understand your choice, please choose 'yes' if you want to take it into account and 'no' otherwise." <<endl;
      }
    }
    else if(name == "spectrum_mode"){
      if(value == "writing" || value == "nothing" || value == "reading"){
        error_value = "no";
        error_name =  "no";
      }
      else {
        error_value = "yes";
        cout << "The parameter 'spectrum_mode' isn't reckognised, please check that it is one among : 'writing', 'reading' and 'nothing'."<<endl;
      }
    }
    else if(name == "inverse_compton_scattering" || name == "double_photon_pair_creation" || name == "pair_creation_in_nuclei" || name == "compton_scattering" || name =="photon_photon_diffusion"){
      if(value == "yes" || value == "no"){
        error_value = "no";
        error_name =  "no";
      }
      else {
        error_value = "yes";
        cout << "The choice for parameter '"<< name <<"' isn't reckognised, please check that it is either 'yes' or 'no'."<<endl;
      }
    }
    else if(name == "Energy_Table_Size"){
      strtod(value.c_str(),&pEnd);
      if(*pEnd != '\0' || pEnd == value.c_str()){
        error_value = "yes";
        cout << "The parameter '" << name << "' is not a double, please check."<<endl;
      }
    }
    else if(name == "E_min_table"){
      strtod(value.c_str(),&pEnd);
      if(*pEnd != '\0' || pEnd == value.c_str()){
        error_value = "yes";
        cout << "The parameter 'E_min_table' is not a double, please check."<<endl;
      }
      if(atof(value.c_str())<1){
        error_value = "yes";
        cout << "The code is not reliable for E < 1 MeV. Please change your minimal energy." << endl;
      }
    }
    else if(name == "m_x"){
      strtod(value.c_str(),&pEnd);
      if(*pEnd != '\0' || pEnd == value.c_str()){
        error_value = "yes";
        cout << "The parameter 'M_x' is not a double, please check."<<endl;
      }
    }
    else if(name == "tau_x"){
      strtod(value.c_str(),&pEnd);
      if(*pEnd != '\0' || pEnd == value.c_str()){
        error_value = "yes";
        cout << "The parameter 'tau_x' is not a double, please check."<<endl;
      }
      if(atof(value.c_str())>1e10 || atof(value.c_str())<1e3){
        error_value = "yes";
        cout << "The parameter 'tau_x' isnt valid. Please choose a value in [1e3,1e10], the code isn't trustable for other lifetimes."<<endl;
      }
    }
    else if(name == "zeta_x"){
      strtod(value.c_str(),&pEnd);
      if(*pEnd != '\0' || pEnd == value.c_str()){
        error_value = "yes";
        cout << "The parameter 'zeta_x' is not a double, please check."<<endl;
      }
    }
    else if(name == "tau_min"){
      strtod(value.c_str(),&pEnd);
      if(*pEnd != '\0' || pEnd == value.c_str()){
        error_value = "yes";
        cout << "The parameter 'tau_min' is not a double, please check."<<endl;
      }
      if(atof(value.c_str())>1e10 || atof(value.c_str())<1e3){
        error_value = "yes";
        cout << "The parameter 'tau_min' isn't valid. Please choose a value in [1e3,1e10], the code isn't trustable for other lifetimes."<<endl;
      }
    }
    else if(name == "tau_max"){
      strtod(value.c_str(),&pEnd);
      if(*pEnd != '\0' || pEnd == value.c_str()){
        error_value = "yes";
        cout << "The parameter 'tau_max' is not a double, please check."<<endl;
      }
      if(atof(value.c_str())>1e12 || atof(value.c_str())<1e3){
        error_value = "yes";
        cout << "The parameter 'tau_max' isn't valid. Please choose a value in [1e3,1e10], the code isn't trustable for other lifetimes."<<endl;
      }
    }
    else if(name == "tau_step"){
      strtod(value.c_str(),&pEnd);
      if(*pEnd != '\0' || pEnd == value.c_str()){
        error_value = "yes";
        cout << "The parameter 'tau_step' is not a double, please check."<<endl;
      }
      else {
        error_value = "no";
      }
    }
    else if(name == "zeta_min"){
      strtod(value.c_str(),&pEnd);
      if(*pEnd != '\0' || pEnd == value.c_str()){
        error_value = "yes";
        cout << "The parameter 'zeta_min' is not a double, please check."<<endl;
      }
      else {
        error_value = "no";
      }
    }
    else if(name == "zeta_max"){
      strtod(value.c_str(),&pEnd);
      if(*pEnd != '\0' || pEnd == value.c_str()){
        error_value = "yes";
        cout << "The parameter 'zeta_max' is not a double, please check."<<endl;
      }
      else {
        error_value = "no";
      }
    }
    else if(name == "zeta_step"){
      strtod(value.c_str(),&pEnd);
      if(*pEnd != '\0' || pEnd == value.c_str()){
        error_value = "yes";
        cout << "The parameter 'zeta_step' is not a double, please check."<<endl;
      }
      else {
        error_value = "no";
      }
    }
    else if(name == "nuclei"){
      if(value == "4He" || value == "3He" || value == "2H" || value == "7Li" || value == "7Be"){
        error_name = "no";
        error_value = "no";
      }
      else {
        error_value = "yes";
      }
    }
    else if(name == "redshift"){
      strtod(value.c_str(),&pEnd);
      if(*pEnd != '\0' || pEnd == value.c_str()){
        error_value = "yes";
        cout << "The parameter 'redshift' is not a double, please check."<<endl;
      }
      else {
        error_value = "no";
      }
    }
    else if(name == "temperature"){
      strtod(value.c_str(),&pEnd);
      if(*pEnd != '\0' || pEnd == value.c_str()){
        error_value = "yes";
        cout << "The parameter 'temperature' is not a double, please check."<<endl;
      }
      else {
        error_value = "no";
      }
    }
    else if(name == "results_files"){
      if(value=="automatic")error_value="no";
      else{
        ofstream file(value);
        if(file){
          error_value="no";
        }
        else{
          cout << "I cannot open file " << value << " please make sure you have created the folder(s) before starting cBBNfast."<<endl;
          error_value = "yes";
        }
      }
    }
    else if(name == "interaction_rate_files"){
      if(value=="automatic")error_value="no";
      else{
        ofstream file(value);
        if(file){
          error_value="no";
        }
        else{
          cout << "I cannot open file " << value << " please make sure you have created the folder(s) before starting cBBNfast."<<endl;
          error_value = "yes";
        }
      }
    }    else if(name == "test_files"){
          if(value=="automatic")error_value="no";
          else{
            ofstream file(value);
            if(file){
              error_value="no";
            }
            else{
              cout << "I cannot open file " << value << " please make sure you have created the folder(s) before starting cBBNfast."<<endl;
              error_value = "yes";
            }
          }
        }
    else if(name == "photon_spectrum_file_name"){
      if(value=="automatic")error_value="no";
      else{
        ifstream file(value);
        if(file){
          error_value="no";
        }
        else{
          cout << "I cannot open file " << value << " please make sure you have created the file before starting cBBNfast."<<endl;
          error_value = "yes";
        }
      }
    }
    else if(name == "electron_spectrum_file_name"){
      if(value=="automatic")error_value="no";
      else{
        ifstream file(value);
        if(file){
          error_value="no";
        }
        else{
          cout << "I cannot open file " << value << " please make sure you have created the file before starting cBBNfast."<<endl;
          error_value = "yes";
        }
      }
    }
    else if(name == "spectrum_files"){
      if(value=="automatic")error_value="no";
      else{
        ofstream file(value);
        if(file){
          error_value="no";
        }
        else{
          cout << "I cannot open file " << value << " please make sure you have created the folder(s) before starting cBBNfast."<<endl;
          error_value = "yes";
        }
      }
    }
    else if(name == "EM_cascade_verbose" || name == "BBN_constraints_verbose" || name == "Input_verbose" || name == "Output_verbose"){
      if(atoi(value.c_str())<0){
        error_value = "yes";
        cout << "The number '" << name <<"' has to be bigger or equal to 0. Please change accordingly." << endl;
      }
      else {
        error_value = "no";
        error_name =  "no";
      }
    }
    else if(name == "ignore_line" || value == "ignore_line"){
      error_value = "no";
      error_name= "no";
    }
    else {
      error_name = "yes";
      cout << " I couldn't reckognised the name " << name << " please check your input file." << endl;
    }
}
double polylog_2(double z, Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){

	double ds;
	double t[pt_Spectrum_and_Precision_Parameters->eval_max], f[pt_Spectrum_and_Precision_Parameters->eval_max];
	double h;
	double resultat=0;
	double t_ini, t_max, factor;
	int y;
  if(z>1)cout << " z bigger than 1 ! " << endl;
  if(fabs(z)<=1){
    factor = 1;
      for(int i = 1 ; i < 1000; i++){
  		resultat += pow(z,i)/(i*i);
  	}
  }
  else {
    if(z<0){
      t_ini=z;
      t_max = 0;
      factor = 1;
    }
    if(z>0){
      t_ini=0;
      t_max = z;
      factor = -1;
    }
    ds = (t_max - t_ini)/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
    y = 0;
    // while(ds>z){
    // 	ds/=10.;
    // 	y++;
    // }
    h = ds/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
    resultat=0;
    // cout << "ds = " << ds << " z = " << z << " y = " << y << endl;
    for(int i=0; i<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1;i++){

          // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
          for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++)
          {
            if(eval == 0){
            if(i==0)	t[eval]=t_ini;
            else t[eval]=t[pt_Spectrum_and_Precision_Parameters->eval_max-1];
            }
            else{
              t[eval]=t[0]+eval*h;
            }
            f[eval]=log(1-t[eval])/t[eval];
            if(isnan(f[eval])==1) f[eval]=0;
            resultat += ds/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
        }
    }
  }



	// cout << " (polylog_2 : )" <<  " z = " << z << " resultat " << factor*resultat << endl;

	return factor*resultat;
}
double expint(const int n, const double x)
{
	static const int MAXIT=100;
	static const double EULER=0.577215664901533,
		EPS=numeric_limits<double>::epsilon(),
		BIG=numeric_limits<double>::max()*EPS;
	int i,ii,nm1=n-1;
	double a,b,c,d,del,fact,h,psi,ans;
	if (n < 0 || x < 0.0 || (x==0.0 && (n==0 || n==1)))
		throw("bad arguments in expint");
	if (n == 0) ans=exp(-x)/x;
	else {
		if (x == 0.0) ans=1.0/nm1;
		else {
			if (x > 1.0) {
				b=x+n;
				c=BIG;
				d=1.0/b;
				h=d;
				for (i=1;i<=MAXIT;i++) {
					a = -i*(nm1+i);
					b += 2.0;
					d=1.0/(a*d+b);
					c=b+a/c;
					del=c*d;
					h *= del;
					if (fabs(del-1.0) <= EPS) {
						ans=h*exp(-x);
						return ans;
					}
				}
				throw("continued fraction failed in expint");
			} else {
				ans = (nm1!=0 ? 1.0/nm1 : -log(x)-EULER);
				fact=1.0;
				for (i=1;i<=MAXIT;i++) {
					fact *= -x/i;
					if (i != nm1) del = -fact/(i-nm1);
					else {
						psi = -EULER;
						for (ii=1;ii<=nm1;ii++) psi += 1.0/ii;
						del=fact*(-log(x)+psi);
					}
					ans += del;
					if (fabs(del) < fabs(ans)*EPS) return ans;
				}
				throw("series failed in expint");
			}
		}
	}
	return ans;
}
double exponential_integral(double x,
														Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters){
	double resultat = 0, h, dt, t[pt_Spectrum_and_Precision_Parameters->eval_max], f[pt_Spectrum_and_Precision_Parameters->eval_max];
	double t_ini = x, t_max = 10000;
	int y = 0;
	dt=(t_max-t_ini)/(double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
	while(dt>t_ini){
		dt/=10.;
		y++;
	}
	h = dt/(pt_Spectrum_and_Precision_Parameters->eval_max-1);
	for(int j=0; j<pow(10,y)*pt_Spectrum_and_Precision_Parameters->n_step-1;j++){

		for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++)
		{
			if(eval == 0){
			if(j==0)	t[eval]=t_ini;
			else t[eval]=t[pt_Spectrum_and_Precision_Parameters->eval_max-1];
			}
			else{
				t[eval]=t[0]+eval*h;
			}

		f[eval]=exp(-t[eval])/t[eval];
		resultat += dt/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];
		// cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
	}
	}
	return resultat;

}
double factorial(int n){
	int factorial;
	for (int i = 0; i <= n; i++){
		if (i == 0)
		factorial = 1;
		else
		factorial *=i;
	}
	return factorial;
}
void check_energy_conservation(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                               Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                               Structure_Spectrum * pt_Gamma_Spectrum,
                               Structure_Spectrum * pt_Electron_Spectrum,
                               double &integrale){
cout << " **** Currently checking energy conservation *** " << endl;
double dE = (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table) / (double) (pt_Spectrum_and_Precision_Parameters->n_step - 1);
double h = dE/6.;
double dE_2, h2;
double E1, E2, E3, E4, E5, E6, E7, f1, f2, f3, f4, f5, f6, f7, g1, g2, g3, g4, g5, g6, g7, E_gamma, E_e, rate_E;
double E[pt_Spectrum_and_Precision_Parameters->eval_max],f[pt_Spectrum_and_Precision_Parameters->eval_max],g[pt_Spectrum_and_Precision_Parameters->eval_max];
double rate_E1,rate_E2,rate_E3,rate_E4,rate_E5,rate_E6,rate_E7;
double resultat_1 = 0, resultat_2 = 0, resultat_3 = 0, resultat_4 = 0, resultat_5 = 0;
double F1, F2, F3, F4, F5, F6, F7;
double z = pt_Gamma_Spectrum->redshift;
double E_c = E_c_0/(1+z);

double E_gamma_bb = 2.701*T_0*(1+z);
double E_cmb_max = 10*E_gamma_bb;
double E_cmb_min = E_gamma_bb/100.;
vector<double> Gamma_Spectrum_Integrated_Over_Kernel;
vector<double> Gamma_Spectrum_Integrated_Over_Kernel_energy;
vector<double> Electron_Spectrum_Integrated_Over_Kernel;
vector<double> Electron_Spectrum_Integrated_Over_Kernel_energy;
double int_bb = 2*pow(T_0*(1+z),3)*1.20205/(pi*pi);
double dE_3 = (pt_Particle_Physics_Model->E_0 - pt_Spectrum_and_Precision_Parameters->E_min_table) / (double) (pt_Spectrum_and_Precision_Parameters->z_step - 1);


if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice=="universal"){

          for(int j=0; j<pt_Spectrum_and_Precision_Parameters->n_step-1;j++){
                            // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
                            for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++)
                            {
                              if(eval == 0){
                              if(j==0)	E[eval]=pt_Spectrum_and_Precision_Parameters->E_min_table;
                              else E[eval]=E[pt_Spectrum_and_Precision_Parameters->eval_max-1];
                              }
                              else{
                                E[eval]=E[0]+eval*h;
                              }

                            if(E[eval]<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E[eval], f[eval]);
                            else f[eval]=0;
                            rate_E = rate_NPC(E[eval],z)+rate_compton(E[eval],z)+rate_gg_scattering(E[eval],z);
                            f[eval]*=rate_E*E[eval];
                            resultat_1 += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];

                          }

                          // if(res_initial==0 && resultat !=0)res_initial = resultat;
                          // if(res_initial!=0 && resultat/res_initial<precision)break;
                          // cout << E1 << " f1 " << f1 << " (exp(E1/T)-1) "  << (exp(E1/T)-1) << " f7 " << f7 << endl;


                          // 	cout << "Egamma = " << E_gamma << " E7 = " << E7 << " resultat = " << resultat<< " j = " << j << endl;
                        }
        integrale = resultat_1;
}
else if(pt_Spectrum_and_Precision_Parameters->calculation_mode == "triangular"){


  for(int j=0; j<pt_Spectrum_and_Precision_Parameters->n_step-1;j++){
        cout << " step j : " << j << " still " << pt_Spectrum_and_Precision_Parameters->n_step-1-j << " to go."<< endl;
                    // cout << "pt_Spectrum_and_Precision_Parameters->eval_max = " << pt_Spectrum_and_Precision_Parameters->eval_max << " h2 " << h2 << endl;
                    for(int eval=0; eval < pt_Spectrum_and_Precision_Parameters->eval_max; eval++)
                    {
                      if(eval == 0){
                      if(j==0)	E[eval]=pt_Spectrum_and_Precision_Parameters->E_min_table;
                      else E[eval]=E[pt_Spectrum_and_Precision_Parameters->eval_max-1];
                      }
                      else{
                        E[eval]=E[0]+eval*h;
                      }

                    if(E[eval]<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E[eval], f[eval]);
                    else f[eval]=0;
                    // rate_E = rate_NPC(E[eval],z)+rate_compton(E[eval],z)+rate_gg_scattering(E[eval],z);
                    // if(pt_Spectrum_and_Precision_Parameters->double_photon_pair_creation=="yes")rate_E+=rate_pair_creation_v2(E[eval],z,pt_Spectrum_and_Precision_Parameters);
                    rate_E = 1;
                    f[eval]*=rate_E*E[eval];
                    resultat_1 += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*f[eval];

                    if(E[eval]<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E[eval], g[eval]);
                    else g[eval]=0;
                    g[eval]*=E[eval];
                    // g[eval]*=E[eval]*(integrator_simpson_rate_inverse_compton_v2(z,E_cmb_min,E_cmb_max,E[eval],pt_Spectrum_and_Precision_Parameters));
                    resultat_2 += dE/pt_Spectrum_and_Precision_Parameters->divisor*pt_Spectrum_and_Precision_Parameters->weight[eval]*g[eval];

                    // cout << "eval " << eval << "E = " << E[eval] << " weight = " << pt_Spectrum_and_Precision_Parameters->weight[eval] << " f[eval] = "<< f[eval] <<" resultat = " << resultat << endl;
                  }

                  // if(res_initial==0 && resultat !=0)res_initial = resultat;
                  // if(res_initial!=0 && resultat/res_initial<precision)break;
                  // cout << E1 << " f1 " << f1 << " (exp(E1/T)-1) "  << (exp(E1/T)-1) << " f7 " << f7 << endl;


                  // 	cout << "Egamma = " << E_gamma << " E7 = " << E7 << " resultat = " << resultat<< " j = " << j << endl;
                }


  integrale = resultat_1;

}


        cout << "The total energy contained in " <<   pt_Gamma_Spectrum->spectrum_name  << "spectrum is " << resultat_1 << " - " << resultat_3 << " = " << resultat_1-resultat_3 << " MeV";
        cout << " and in " << pt_Electron_Spectrum->spectrum_name << "spectrum is " << resultat_2 << " - " << resultat_4 << " = " << resultat_2-resultat_4 << " MeV ";
        cout << "for a total of "<< integrale <<" MeV, you had injected " << pt_Particle_Physics_Model->E_0 << " MeV." << endl;


}





/*********************************FONCTIONS A DEVELOPPER POUR METHODE ITERATIVE*********************************/
// else{
// for(int i=0;i<pt_Spectrum_and_Precision_Parameters->z_step;i++){
//   E_gamma = pt_Spectrum_and_Precision_Parameters->E_min_table+i*dE_3;
//   E_e = E_gamma;
//   resultat_1 = 0;
//   resultat_2 = 0;
//   resultat_3 = 0;
//   resultat_4 = 0;
//   double resultat_5 = 0;
//   dE_2 = (pt_Particle_Physics_Model->E_0 - (E_gamma))/ (double) (pt_Spectrum_and_Precision_Parameters->n_step-1);
//   h2 = dE_2/6.;
//
//         for(int j=0; j<pt_Spectrum_and_Precision_Parameters->n_step-1;j++){
//
//           if(j==0){
//             E1=E_gamma;
//             }
//           else{
//             E1=E7;
//           }
//
//           E2=E1 + h2;
//           E3=E1 + 2*h2;
//           E4=E1 + 3*h2;
//           E5=E1 + 4*h2;
//           E6=E1 + 5*h2;
//           E7=E1 + 6*h2;
//           /********** First : Integrate photon spectrum over gamma->gamma kernel ******/
//           if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E1, f1);
//           else f1=0;
//           if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E2, f2);
//           else f2=0;
//           if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E3, f3);
//           else f3=0;
//           if(E4<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E4, f4);
//           else f4=0;
//           if(E5<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E5, f5);
//           else f5=0;
//           if(E6<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E6, f6);
//           else f6=0;
//           if(E7<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E7, f7);
//           else f7=0;
//
//           F1=f1*(dsigma_phph(E1,z,E_gamma)+dsigma_compton(E1,z,E_gamma));
//           F2=f2*(dsigma_phph(E2,z,E_gamma)+dsigma_compton(E2,z,E_gamma));
//           F3=f3*(dsigma_phph(E3,z,E_gamma)+dsigma_compton(E3,z,E_gamma));
//           F4=f4*(dsigma_phph(E4,z,E_gamma)+dsigma_compton(E4,z,E_gamma));
//           F5=f5*(dsigma_phph(E5,z,E_gamma)+dsigma_compton(E5,z,E_gamma));
//           F6=f6*(dsigma_phph(E6,z,E_gamma)+dsigma_compton(E6,z,E_gamma));
//           F7=f7*(dsigma_phph(E7,z,E_gamma)+dsigma_compton(E7,z,E_gamma));
//           // cout << " E1 = " << E1 << endl;
//
//           resultat_1 += dE_2/840. * (41*F1+216*F2+27*F3+272*F4+27*F5+216*F6+41*F7);
//           if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "Dirac" && j==pt_Spectrum_and_Precision_Parameters->n_step-2 && pt_Spectrum_and_Precision_Parameters->calculation_mode == "iterative"){
//             resultat_1+=(dsigma_phph(pt_Particle_Physics_Model->E_0,z,E_gamma)+dsigma_compton(pt_Particle_Physics_Model->E_0,z,E_gamma))/(rate_NPC(pt_Particle_Physics_Model->E_0,z)+rate_compton(pt_Particle_Physics_Model->E_0,z)+rate_gg_scattering(pt_Particle_Physics_Model->E_0,z));
//           }
//           /******** Second : Integrate electron spectrum over e->gamma kernel ********/
//
//
//           if(E1+m_e<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E1+m_e, g1);
//           else g1=0;
//           if(E2+m_e<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E2+m_e, g2);
//           else g2=0;
//           if(E3+m_e<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E3+m_e, g3);
//           else g3=0;
//           if(E4+m_e<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E4+m_e, g4);
//           else g4=0;
//           if(E5+m_e<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E5+m_e, g5);
//           else g5=0;
//           if(E6+m_e<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E6+m_e, g6);
//           else g6=0;
//           if(E7+m_e<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E7+m_e, g7);
//           else g7=0;
//
//
//           F1=2*g1*Function_Integrand_Spectre_Compton(E_gamma,E1+m_e, E_gamma_bb)/((E1+m_e)*(E1+m_e));
//           F2=2*g2*Function_Integrand_Spectre_Compton(E_gamma,E2+m_e, E_gamma_bb)/((E2+m_e)*(E2+m_e));
//           F3=2*g3*Function_Integrand_Spectre_Compton(E_gamma,E3+m_e, E_gamma_bb)/((E3+m_e)*(E3+m_e));
//           F4=2*g4*Function_Integrand_Spectre_Compton(E_gamma,E4+m_e, E_gamma_bb)/((E4+m_e)*(E4+m_e));
//           F5=2*g5*Function_Integrand_Spectre_Compton(E_gamma,E5+m_e, E_gamma_bb)/((E5+m_e)*(E5+m_e));
//           F6=2*g6*Function_Integrand_Spectre_Compton(E_gamma,E6+m_e, E_gamma_bb)/((E6+m_e)*(E6+m_e));
//           F7=2*g7*Function_Integrand_Spectre_Compton(E_gamma,E7+m_e, E_gamma_bb)/((E7+m_e)*(E7+m_e));
//           // cout << "E7 "<<E7<<" F7  = " << F7 << endl;
//
//           resultat_2 += dE_2/840. * (41*F1+216*F2+27*F3+272*F4+27*F5+216*F6+41*F7);
//           if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice == "Dirac" && j==pt_Spectrum_and_Precision_Parameters->n_step-2 && pt_Spectrum_and_Precision_Parameters->calculation_mode == "iterative"){
//             resultat_2+=2*Function_Integrand_Spectre_Compton(E_gamma,pt_Particle_Physics_Model->E_0, E_gamma_bb)/((pt_Particle_Physics_Model->E_0)*(pt_Particle_Physics_Model->E_0))/Rate_Inverse_Compton(pt_Particle_Physics_Model->E_0,z,pt_Spectrum_and_Precision_Parameters);
//           }
//           /*******Third : Integrate photon spectrum over gamma->e kernel*******/
//           F1=f1*(dsigma_compton(E1,z,(E1+m_e-E_e))+dsigma_NPC(E1+m_e,z,E_e));
//           F2=f2*(dsigma_compton(E2,z,(E2+m_e-E_e))+dsigma_NPC(E2+m_e,z,E_e));
//           F3=f3*(dsigma_compton(E3,z,(E3+m_e-E_e))+dsigma_NPC(E3+m_e,z,E_e));
//           F4=f4*(dsigma_compton(E4,z,(E4+m_e-E_e))+dsigma_NPC(E4+m_e,z,E_e));
//           F5=f5*(dsigma_compton(E5,z,(E5+m_e-E_e))+dsigma_NPC(E5+m_e,z,E_e));
//           F6=f6*(dsigma_compton(E6,z,(E6+m_e-E_e))+dsigma_NPC(E6+m_e,z,E_e));
//           F7=f7*(dsigma_compton(E7,z,(E7+m_e-E_e))+dsigma_NPC(E7+m_e,z,E_e));
//           resultat_3 += dE_2/840. * (41*F1+216*F2+27*F3+272*F4+27*F5+216*F6+41*F7);
//           if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "Dirac" && j==pt_Spectrum_and_Precision_Parameters->n_step-2 && pt_Spectrum_and_Precision_Parameters->calculation_mode == "iterative"){
//             resultat_3+=(dsigma_compton(pt_Particle_Physics_Model->E_0,z,(pt_Particle_Physics_Model->E_0+m_e-E_e))+dsigma_NPC(pt_Particle_Physics_Model->E_0+m_e,z,pt_Particle_Physics_Model->E_0))/(rate_NPC(pt_Particle_Physics_Model->E_0,z)+rate_compton(pt_Particle_Physics_Model->E_0,z)+rate_gg_scattering(pt_Particle_Physics_Model->E_0,z));
//           }
//
//
//           /******** Fourth : Integrate electron spectrum over e->e kernel ********/
//
//           if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E1, g1);
//           else g1=0;
//           if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E2, g2);
//           else g2=0;
//           if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E3, g3);
//           else g3=0;
//           if(E4<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E4, g4);
//           else g4=0;
//           if(E5<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E5, g5);
//           else g5=0;
//           if(E6<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E6, g6);
//           else g6=0;
//           if(E7<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E7, g7);
//           else g7=0;
//
//           F1=g1*Function_Integrand_Spectre_Compton(E1+E_gamma_bb-E_e,E1, E_gamma_bb)/((E1)*(E1));
//           F2=g2*Function_Integrand_Spectre_Compton(E2+E_gamma_bb-E_e,E2, E_gamma_bb)/((E2)*(E2));
//           F3=g3*Function_Integrand_Spectre_Compton(E3+E_gamma_bb-E_e,E3, E_gamma_bb)/((E3)*(E3));
//           F4=g4*Function_Integrand_Spectre_Compton(E4+E_gamma_bb-E_e,E4, E_gamma_bb)/((E4)*(E4));
//           F5=g5*Function_Integrand_Spectre_Compton(E5+E_gamma_bb-E_e,E5, E_gamma_bb)/((E5)*(E5));
//           F6=g6*Function_Integrand_Spectre_Compton(E6+E_gamma_bb-E_e,E6, E_gamma_bb)/((E6)*(E6));
//           F7=g7*Function_Integrand_Spectre_Compton(E7+E_gamma_bb-E_e,E7, E_gamma_bb)/((E7)*(E7));
//
//           resultat_4+= dE_2/840. * (41*F1+216*F2+27*F3+272*F4+27*F5+216*F6+41*F7);
//           if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice == "Dirac" && j==pt_Spectrum_and_Precision_Parameters->n_step-2 && pt_Spectrum_and_Precision_Parameters->calculation_mode == "iterative"){
//             resultat_4+=Function_Integrand_Spectre_Compton(pt_Particle_Physics_Model->E_0+E_gamma_bb-E_e,pt_Particle_Physics_Model->E_0, E_gamma_bb)/((pt_Particle_Physics_Model->E_0)*(pt_Particle_Physics_Model->E_0))/Rate_Inverse_Compton(pt_Particle_Physics_Model->E_0,z,pt_Spectrum_and_Precision_Parameters);
//           }
//
//         }
//
//         resultat_2*=2*pi*r_e*r_e*m_e*m_e*int_bb/E_gamma_bb;
//         resultat_4*=2*pi*r_e*r_e*m_e*m_e*int_bb/E_gamma_bb;
//
//
//         // cout << "Egamma = " << E_gamma << " resultat = " << resultat_1   << "resultat 2 =  "<< resultat_2<< " resultat 3 = " << resultat_3 << "resultat 4 = " << resultat_4<<endl;
//         // if(pt_Spectrum_and_Precision_Parameters->inverse_compton_scattering=="no" && pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice=="none"){
//         //   resultat_2 = 0;
//         //   resultat_3 = 0;
//         //   resultat_4 = 0;
//         // }
//         // if(pt_Spectrum_and_Precision_Parameters->inverse_compton_scattering=="no" && pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice=="none"){
//         //   resultat_1 = 0;
//         //   resultat_2 = 0;
//         //   resultat_3 = 0;
//         // }
//         Gamma_Spectrum_Integrated_Over_Kernel_energy.push_back(E_gamma);
//         Gamma_Spectrum_Integrated_Over_Kernel.push_back(resultat_1+resultat_2);
//         Electron_Spectrum_Integrated_Over_Kernel_energy.push_back(E_e);
//         Electron_Spectrum_Integrated_Over_Kernel.push_back(resultat_3+resultat_4);
//   }
//   resultat_1 = 0;
//   resultat_2 = 0;
//   resultat_3 = 0;
//   resultat_4 = 0;
//         for(int i=0;i<pt_Spectrum_and_Precision_Parameters->n_step-1;i++){
//           if(i==0){
//             E1=pt_Spectrum_and_Precision_Parameters->E_min_table;
//             }
//           else{
//             E1=E7;
//           }
//
//           E2=E1 + h;
//           E3=E2 + h;
//           E4=E3 + h;
//           E5=E4 + h;
//           E6=E5 + h;
//           E7=E6 + h;
//           if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E1, f1);
//           else f1=0;
//           if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E2, f2);
//           else f2=0;
//           if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E3, f3);
//           else f3=0;
//           if(E4<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E4, f4);
//           else f4=0;
//           if(E5<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E5, f5);
//           else f5=0;
//           if(E6<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E6, f6);
//           else f6=0;
//           if(E7<pt_Particle_Physics_Model->E_0)linearint(pt_Gamma_Spectrum->Energy, pt_Gamma_Spectrum->Spectrum, pt_Gamma_Spectrum->Energy.size(), E7, f7);
//           else f7=0;
//
//           f1*=E1*(rate_NPC(E1,z)+rate_compton(E1,z)+rate_gg_scattering(E1,z));
//           f2*=E2*(rate_NPC(E2,z)+rate_compton(E2,z)+rate_gg_scattering(E2,z));
//           f3*=E3*(rate_NPC(E3,z)+rate_compton(E3,z)+rate_gg_scattering(E3,z));
//           f4*=E4*(rate_NPC(E4,z)+rate_compton(E4,z)+rate_gg_scattering(E4,z));
//           f5*=E5*(rate_NPC(E5,z)+rate_compton(E5,z)+rate_gg_scattering(E5,z));
//           f6*=E6*(rate_NPC(E6,z)+rate_compton(E6,z)+rate_gg_scattering(E6,z));
//           f7*=E7*(rate_NPC(E7,z)+rate_compton(E7,z)+rate_gg_scattering(E7,z));
//           resultat_1 += dE/840. * (41*f1+216*f2+27*f3+272*f4+27*f5+216*f6+41*f7);
//           if(pt_Spectrum_and_Precision_Parameters->photon_spectrum_choice == "Dirac" &&  i==pt_Spectrum_and_Precision_Parameters->n_step-2 && pt_Spectrum_and_Precision_Parameters->calculation_mode == "iterative"){
//             resultat_1+=pt_Particle_Physics_Model->E_0;
//           }
//
//           if(E1<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E1, f1);
//           else f1=0;
//           if(E2<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E2, f2);
//           else f2=0;
//           if(E3<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E3, f3);
//           else f3=0;
//           if(E4<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E4, f4);
//           else f4=0;
//           if(E5<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E5, f5);
//           else f5=0;
//           if(E6<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E6, f6);
//           else f6=0;
//           if(E7<pt_Particle_Physics_Model->E_0)linearint(pt_Electron_Spectrum->Energy, pt_Electron_Spectrum->Spectrum, pt_Electron_Spectrum->Energy.size(), E7, f7);
//           else f7=0;
//
//           f1*=E1*(Rate_Inverse_Compton(E1,z,pt_Spectrum_and_Precision_Parameters));
//           f2*=E2*(Rate_Inverse_Compton(E2,z,pt_Spectrum_and_Precision_Parameters));
//           f3*=E3*(Rate_Inverse_Compton(E3,z,pt_Spectrum_and_Precision_Parameters));
//           f4*=E4*(Rate_Inverse_Compton(E4,z,pt_Spectrum_and_Precision_Parameters));
//           f5*=E5*(Rate_Inverse_Compton(E5,z,pt_Spectrum_and_Precision_Parameters));
//           f6*=E6*(Rate_Inverse_Compton(E6,z,pt_Spectrum_and_Precision_Parameters));
//           f7*=E7*(Rate_Inverse_Compton(E7,z,pt_Spectrum_and_Precision_Parameters));
//           resultat_2 += dE/840. * (41*f1+216*f2+27*f3+272*f4+27*f5+216*f6+41*f7);
//           if(pt_Spectrum_and_Precision_Parameters->electron_spectrum_choice == "Dirac" &&  i==pt_Spectrum_and_Precision_Parameters->n_step-2 && pt_Spectrum_and_Precision_Parameters->calculation_mode == "iterative"){
//             resultat_2+=pt_Particle_Physics_Model->E_0;
//           }
//
//           if(E1<pt_Particle_Physics_Model->E_0)linearint(Gamma_Spectrum_Integrated_Over_Kernel_energy, Gamma_Spectrum_Integrated_Over_Kernel, Gamma_Spectrum_Integrated_Over_Kernel_energy.size(), E1, g1);
//           else g1=0;
//           if(E2<pt_Particle_Physics_Model->E_0)linearint(Gamma_Spectrum_Integrated_Over_Kernel_energy, Gamma_Spectrum_Integrated_Over_Kernel, Gamma_Spectrum_Integrated_Over_Kernel_energy.size(), E2, g2);
//           else g2=0;
//           if(E3<pt_Particle_Physics_Model->E_0)linearint(Gamma_Spectrum_Integrated_Over_Kernel_energy, Gamma_Spectrum_Integrated_Over_Kernel, Gamma_Spectrum_Integrated_Over_Kernel_energy.size(), E3, g3);
//           else g3=0;
//           if(E4<pt_Particle_Physics_Model->E_0)linearint(Gamma_Spectrum_Integrated_Over_Kernel_energy, Gamma_Spectrum_Integrated_Over_Kernel, Gamma_Spectrum_Integrated_Over_Kernel_energy.size(), E4, g4);
//           else g4=0;
//           if(E5<pt_Particle_Physics_Model->E_0)linearint(Gamma_Spectrum_Integrated_Over_Kernel_energy, Gamma_Spectrum_Integrated_Over_Kernel, Gamma_Spectrum_Integrated_Over_Kernel_energy.size(), E5, g5);
//           else g5=0;
//           if(E6<pt_Particle_Physics_Model->E_0)linearint(Gamma_Spectrum_Integrated_Over_Kernel_energy, Gamma_Spectrum_Integrated_Over_Kernel, Gamma_Spectrum_Integrated_Over_Kernel_energy.size(), E6, g6);
//           else g6=0;
//           if(E7<pt_Particle_Physics_Model->E_0)linearint(Gamma_Spectrum_Integrated_Over_Kernel_energy, Gamma_Spectrum_Integrated_Over_Kernel, Gamma_Spectrum_Integrated_Over_Kernel_energy.size(), E7, g7);
//           else g7=0;
//
//           g1*=E1;
//           g2*=E2;
//           g3*=E3;
//           g4*=E4;
//           g5*=E5;
//           g6*=E6;
//           g7*=E7;
//
//           resultat_3+=dE/840. * (41*g1+216*g2+27*g3+272*g4+27*g5+216*g6+41*g7);
//
//
//           if(E1<pt_Particle_Physics_Model->E_0)linearint(Electron_Spectrum_Integrated_Over_Kernel_energy, Electron_Spectrum_Integrated_Over_Kernel, Electron_Spectrum_Integrated_Over_Kernel_energy.size(), E1, g1);
//           else g1=0;
//           if(E2<pt_Particle_Physics_Model->E_0)linearint(Electron_Spectrum_Integrated_Over_Kernel_energy, Electron_Spectrum_Integrated_Over_Kernel, Electron_Spectrum_Integrated_Over_Kernel_energy.size(), E2, g2);
//           else g2=0;
//           if(E3<pt_Particle_Physics_Model->E_0)linearint(Electron_Spectrum_Integrated_Over_Kernel_energy, Electron_Spectrum_Integrated_Over_Kernel, Electron_Spectrum_Integrated_Over_Kernel_energy.size(), E3, g3);
//           else g3=0;
//           if(E4<pt_Particle_Physics_Model->E_0)linearint(Electron_Spectrum_Integrated_Over_Kernel_energy, Electron_Spectrum_Integrated_Over_Kernel, Electron_Spectrum_Integrated_Over_Kernel_energy.size(), E4, g4);
//           else g4=0;
//           if(E5<pt_Particle_Physics_Model->E_0)linearint(Electron_Spectrum_Integrated_Over_Kernel_energy, Electron_Spectrum_Integrated_Over_Kernel, Electron_Spectrum_Integrated_Over_Kernel_energy.size(), E5, g5);
//           else g5=0;
//           if(E6<pt_Particle_Physics_Model->E_0)linearint(Electron_Spectrum_Integrated_Over_Kernel_energy, Electron_Spectrum_Integrated_Over_Kernel, Electron_Spectrum_Integrated_Over_Kernel_energy.size(), E6, g6);
//           else g6=0;
//           if(E7<pt_Particle_Physics_Model->E_0)linearint(Electron_Spectrum_Integrated_Over_Kernel_energy, Electron_Spectrum_Integrated_Over_Kernel, Electron_Spectrum_Integrated_Over_Kernel_energy.size(), E7, g7);
//           else g7=0;
//
//           g1*=E1;
//           g2*=E2;
//           g3*=E3;
//           g4*=E4;
//           g5*=E5;
//           g6*=E6;
//           g7*=E7;
//
//           resultat_4+=dE/840. * (41*g1+216*g2+27*g3+272*g4+27*g5+216*g6+41*g7);
//
//           // cout << " E7 = " << E7 << "resultat_1 " << resultat_1 << " resultat_2 = " << resultat_2 << "resultat_3 = " << resultat_3 <<  " i = " << i <<  endl;
//
//         }
//         integrale = (resultat_1-resultat_3)+(resultat_2-resultat_4);
//       }
