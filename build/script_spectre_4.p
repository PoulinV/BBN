set logscale x
set xrange [1:1e3]
set yrange [0:40]
#f1(x)=log10(1.4e2*(1e-12)*(0.001*x)**(-1.56))
#f2(x)=log10(6.2e-12*(1e-12)*(x*10**-3)**(-5.07)*1e30)

plot "test/spectrum_all_1000MeV_steps500_electron.dat" u 1:(log10($2)) w l t 'total photons', \
"test/spectrum_all_1000MeV_steps500_electron.dat" u 1:(log10($2)) w l t 'total electrons', \
"Output/Cascade_Spectrum_Folder/Spectrum_universal_photon_m2000_z425643.dat" u 1:(log10($2)) w l t 'universal' lc 'black', \
"test/spectrum_DP_1000MeV_steps500_electron.dat" u 1:(log10($2)) w l t 'DP electrons 500', \
"test/spectrum_DP+ICS_1000MeV_steps500_electron.dat" u 1:(log10($2)) w l t 'DP+ICS electrons 500', \
"test/spectrum_DP+ICS_1000MeV_steps500_photon.dat" u 1:(log10($2)) w l t 'DPC+ICS photons 500', \
"test/spectrum_DP+ICS+PP_1000MeV_steps500_electron.dat" u 1:(log10($2)) w l t 'DP+ICS+PP electrons 500', \
"test/spectrum_DP+ICS+PP_1000MeV_steps500_photon.dat" u 1:(log10($2)) w l t 'DP+ICS+PP photons 500', \
"test/spectrum_DP+ICS+PP+CS_1000MeV_steps500_electron.dat" u 1:(log10($2)) w l t 'DP+ICS+PP+CS electrons 500', \
"test/spectrum_DP+ICS+PP+CS_1000MeV_steps500_photon.dat" u 1:(log10($2)) w l t 'DP+ICS+PP+CS photons 500', \
"test/spectrum_DP_1000MeV_steps200_electron.dat" u 1:(log10($2)) w l t 'DP electrons 200', \
"test/spectrum_DP+ICS_1000MeV_steps200_electron.dat" u 1:(log10($2)) w l t 'DP+ICS electrons 200', \
"test/spectrum_DP+ICS_1000MeV_steps200_photon.dat" u 1:(log10($2)) w l t 'DPC+ICS photons 200', \
"test/spectrum_DP+ICS+PP_1000MeV_steps200_electron.dat" u 1:(log10($2)) w l t 'DP+ICS+PP electrons 200', \
"test/spectrum_DP+ICS+PP_1000MeV_steps200_photon.dat" u 1:(log10($2)) w l t 'DP+ICS+PP photons 200', \
"test/spectrum_DP+ICS+PP+CS_1000MeV_steps200_electron.dat" u 1:(log10($2)) w l t 'DP+ICS+PP+CS electrons 200', \
"test/spectrum_DP+ICS+PP+CS_1000MeV_steps200_photon.dat" u 1:(log10($2)) w l t 'DP+ICS+PP+CS photons 200'
#f1(x) w l t 'Kawasaki & Moroi', \
#f2(x) w l
