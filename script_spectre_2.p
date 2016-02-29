set logscale x

plot 'test/DP_v2_electron.dat' u 1:(log10($2)) w l t 'DP electron' , \
'test/DP+ICS_v2_electron.dat' u 1:(log10($2)) w l t 'DP+ICS electron', \
'test/DP+ICS_v2_photon.dat' u 1:(log10($2)) w l t 'DP+ICS photon'
