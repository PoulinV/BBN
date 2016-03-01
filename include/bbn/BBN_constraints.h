#ifndef __BBN_BBN_CONSTRAINTS_H__
#define __BBN_BBN_CONSTRAINTS_H__

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "bbn/common.h"
#include "bbn/structures.h"

void Spectrum_and_cross_sections_convolution(Structure_Spectrum * pt_Cascade_Spectrum,
                                            Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                                            int i_min,
                                            int i_max,
                                            double &resultat,
                                            double z,
                                            long n_step,
                                            const std::string &spectrum_choice);
static void  Compute_Constraints_from_destruction_only_loop(const int step,
        Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Scan_Parameters_and_Results * pt_Scan_Parameters_and_Results,
        Structure_Output_Options * pt_Output_Options,
        std::vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei,
        std::vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei,
        const int i_min,
        const int i_max);
static void  compute_constraints_from_destruction_and_production_loop(const int step,
        Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
        Structure_Scan_Parameters_and_Results * pt_Scan_Parameters_and_Results,
        Structure_Output_Options * pt_Output_Options,
        std::vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_Nuclei,
        std::vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_Nuclei,
        std::vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Destruction_4He,
        std::vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_Destruction_4He,
        std::vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_redshift_Production_Nuclei,
        std::vector<double> &Cascade_Spectrum_Integrated_Over_Cross_Section_Production_Nuclei,
        const int i_min,
        const int i_max,
        const int j_min,
        const int j_max,
        const int k_min,
        const int k_max);
void Compute_Constraints_from_destruction_only(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                                              Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                                              Structure_Scan_Parameters_and_Results * pt_Scan_Parameters_and_Results,
                                              Structure_Output_Options * pt_Output_Options);

void Compute_constraints_from_destruction_and_production(Structure_Particle_Physics_Model * pt_Particle_Physics_Model,
                                                        Structure_Spectrum_and_Precision_Parameters * pt_Spectrum_and_Precision_Parameters,
                                                        Structure_Scan_Parameters_and_Results * pt_Scan_Parameters_and_Results,
                                                        Structure_Output_Options * pt_Output_Options);


double cross_section(double  x, int i);
void Check_nuclei(const std::string &nuclei,
                  int &i_min,
                  int &i_max,
                  int &k_min,
                  int &k_max,
                  double &Y_min,
                  double &Y_max,
                  double &Y_0);

#endif // __BBN_BBN_CONSTRAINTS_H__
