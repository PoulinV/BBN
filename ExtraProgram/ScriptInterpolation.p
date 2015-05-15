      set   autoscale                        # scale axes automatically

      unset label                            # remove any previous labels
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically

      set xlabel "z"
      set ylabel "Integrale"
      set logscale
      plot "../OutputExtra/Spectre4Iterations.dat"  using 1:2 title 'interpolation' , \
      "../OutputExtra/results_interpolations_spectre.dat"  using 1:2 title 'results' , \
