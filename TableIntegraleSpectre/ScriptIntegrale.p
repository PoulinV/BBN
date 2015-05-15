      set   autoscale                        # scale axes automatically

      unset label                            # remove any previous labels
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically

      set xlabel "z"
      set ylabel "Integrale"
      set logscale
      plot "ExtraOutput/Interpolation.dat"  using 1:4 title 'interpolation' , \
      "TableIntegraleSpectre/results_trois_iterations_Produc2H_70MeV.dat"  using 1:2 title 'results' w l, \
