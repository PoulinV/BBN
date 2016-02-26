set logscale x
set xrange [1:1e3]
set yrange [0:40]
#f1(x)=log10(1.4e2*(1e-12)*(0.001*x)**(-1.56))
#f2(x)=log10(6.2e-12*(1e-12)*(x*10**-3)**(-5.07)*1e30)

plot "test/spectrum_ICS_1000MeV_steps500_electron.dat" u 1:(log10($2)) w l t 'ICS v1 500', \
"test/spectrum_ICS_1000MeV_steps200_electron.dat" u 1:(log10($2)) w l t 'ICS v1 200', \
"test/spectrum_ICSv2_1000MeV_steps200_electron.dat" u 1:(log10($2)) w l t 'ICS v2 200', \
"test/spectrum_ICSv2_1000MeV_steps500_electron.dat" u 1:(log10($2)) w l t 'ICS v2 500', \
"test/spectrum_ICSv3_1000MeV_steps500_electron.dat" u 1:(log10($2)) w l t 'ICS v3 500'
#f1(x) w l t 'Kawasaki & Moroi', \
#f2(x) w l
