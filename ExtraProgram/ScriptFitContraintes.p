      set   autoscale                        # scale axes automatically

      unset label                            # remove any previous labels
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically

      set xlabel "z"
      set ylabel "f(z)"
      set autoscale
      unset logscale;
	  f1(x)=a*x+b
	  f2(x)=c*x+d
	  f3(x)=e*x+f
	  f4(x)=g*x+h
	  f5(x)=i*x+j
	  f6(x)=k*x+l
    f7(x)=m*x+n

    fit [x = 7.2:7.3] f7(x) "results_Reinj_Monochromatique_Destruc2H.dat" using 1:2 via m,n
	  fit [x = 7:7.2] f6(x) "results_Reinj_Monochromatique_Destruc2H.dat" using 1:2 via k,l
	  fit [x = 6.8:7] f5(x) "results_Reinj_Monochromatique_Destruc2H.dat" using 1:2 via i,j
    fit [x = 6.4:6.7] f4(x) "results_Reinj_Monochromatique_Destruc2H.dat" using 1:2 via g,h
    fit [x =5.3:6] f3(x) "results_Reinj_Monochromatique_Destruc2H.dat" using 1:2 via e,f
    fit [x = 5.:5.3] f2(x) "results_Reinj_Monochromatique_Destruc2H.dat" using 1:2 via c,d
    fit [x = 4.2:5.] f1(x) "results_Reinj_Monochromatique_Destruc2H.dat" using 1:2 via a,b
    plot "results_Reinj_Monochromatique_Destruc2H.dat" using 1:2 notitle with lines , \
    f1(x) ,\
    f2(x) ,\
    f3(x) ,\
    f4(x) ,\
    f5(x) ,\
    f6(x) ,\
    f7(x)
