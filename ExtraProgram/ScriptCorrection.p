set   autoscale                        # scale axes automatically

unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically

set xlabel "E"
set ylabel "f"
set logscale

plot "../OutputExtra/Correction4iterations_70MeV.dat"  using 1:2 title 'correction 70MeV', \
"../OutputExtra/Correction4iterations_70MeV_T10000.dat"  using 1:2 title 'correction 70MeV T10000', \
"../OutputExtra/Correction4iterations_30MeV.dat"  using 1:2 title 'correction 30MeV', \
"../OutputExtra/Correction4iterations_20MeV.dat"  using 1:2 title 'correction 20MeV', \
"../OutputExtra/Correction4iterations_4MeV.dat"  using 1:2 title 'correction 4MeV', \
