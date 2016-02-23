set logscale
#set xrange [2:1000]
plot "Output/Cascade_Spectrum_Folder/Spectrum_total_photons_m200000_z425643_triangular.dat" u 1:2 w l t 'photons', \
"Output/Cascade_Spectrum_Folder/Spectrum_universal_photon_m200000_z425643.dat" u 1:2 w l t 'universal', \
"Output/Cascade_Spectrum_Folder/Spectrum_total_electrons_m200000_z425643_triangular.dat" u 1:2 w l t 'electrons'
