set logscale
set xrange [2:200]
plot "Output/Cascade_Spectrum_Folder/Spectrum_universal_photon_m2e+06_z425643.dat" u 1:2 w l , \
 "Output/Cascade_Spectrum_Folder/Spectrum_total_photons_m2e+06_z425643_triangular.dat" u 1:2 , \
 "Output/Cascade_Spectrum_Folder/Spectrum_total_electrons_m2e+06_z425643_triangular.dat" u 1:2
