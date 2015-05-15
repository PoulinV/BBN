      set   autoscale                        # scale axes automatically

      unset label                            # remove any previous labels
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically

      set xlabel "E"
      set ylabel "f"
      set logscale

   plot "../OutputExtra/SpectreTotal_70MeV.dat"  using 1:2 title 'spectre une iterations' with lines, \
"../OutputExtra/SpectreTotal_70MeV.dat"  using 1:3 title 'spectre deux iterations' with lines, \
"../OutputExtra/SpectreTotal_70MeV.dat"  using 1:4 title 'spectre trois iterations' with lines, \
"../OutputExtra/SpectreTotal_70MeV.dat"  using 1:5 title 'spectre quatre iterations' with lines, \
"../OutputExtra/SpectreTotal_70MeV.dat"  using 1:6 title 'spectre cinq iterations' with lines, \
"../OutputExtra/SpectreTotal_70MeV.dat"  using 1:7 title 'spectre six iterations' with lines, \
"../OutputExtra/SpectreTotal_70MeV.dat"  using 1:8 title 'spectre sept iterations' with lines, \
"../OutputExtra/SpectreTotal_70MeV.dat"  using 1:10 title 'spectre compton' with lines, \
"../OutputExtra/SpectreTotal_70MeV.dat"  using 1:9 title 'spectre total' with lines, \
