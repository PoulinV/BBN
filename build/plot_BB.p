T = 0.0001
set xrange [0.00001:1]
set logscale
f(x) = x*x/(exp(x/T)-1)
plot f(x) title 'BB' w l
