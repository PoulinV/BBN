set logscale
plot "test/interaction_rate_z425643.dat" u 1:2 w l t 'yp->ee', \
"test/interaction_rate_z425643.dat" u 1:3 w l t 'ye->ye', \
"test/interaction_rate_z425643.dat" u 1:4 w l t 'yy->yy', \
"test/interaction_rate_z425643.dat" u 1:5 w l t 'yy->ee', \
"test/interaction_rate_z425643.dat" u 1:6 w l t 'tot', \
"test/interaction_rate_z425643.dat" u 1:7 w l t 'ey->ye', \
"test/interaction_rate_z425643.dat" u 1:8 w l t 'ey->ye v2'
