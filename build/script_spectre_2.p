set logscale
set xrange [2:5001]
plot "Output/Cascade_Spectrum_Folder/Spectrum_universal_photon_m500_z500000.dat" u 1:2 w l , \
 "Output/Cascade_Spectrum_Folder/Spectrum_total_photons_m500_z500000_triangular.dat" u 1:2 , \
 "Output/Cascade_Spectrum_Folder/Spectrum_total_electrons_m500_z500000_triangular.dat" u 1:2
