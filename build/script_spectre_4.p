set logscale x
set xrange [1:1e5]
set yrange [0:40]


plot "Output/Cascade_Spectrum_Folder/Spectrum_universal_photon_m200000_z425643.dat" u 1:(log10($2)) w l t 'universal' lc 'black', \
"test/DP_electron.dat" u 1:(log10($2)) w l t 'DP electrons', \
"test/DP+ICS_electron.dat" u 1:(log10($2)) w l t 'DP+ICS electrons', \
"test/DP+ICS+PP_electron.dat" u 1:(log10($2)) w l t 'DP+ICS+PP electrons', \
"test/DP+ICS+PP+CS_electron.dat" u 1:(log10($2)) w l t 'DP+ICS+PP+CS electrons', \
"Output/Cascade_Spectrum_Folder/Spectrum_total_electrons_m200000_z425643_triangular.dat" u 1:(log10($2)) w l t 'Total electrons', \
"test/DP_photon.dat" u 1:(log10($2)) w l t 'DP photons', \
"test/DP+ICS_photon.dat" u 1:(log10($2)) w l t 'DP+ICS photons', \
"test/DP+ICS+PP_photon.dat" u 1:(log10($2)) w l t 'DP+ICS+PP photons', \
"test/DP+ICS+PP+CS_photon.dat" u 1:(log10($2)) w l t 'DP+ICS+PP+CS photons', \
"Output/Cascade_Spectrum_Folder/Spectrum_total_photons_m200000_z425643_triangular.dat" u 1:(log10($2)) w l t 'Total photons'
