set yrange [*:*]
set xrange [*:*]
set xlabel "step"
set ylabel "norm"

p "norm.dat" u 1:2 w l ti "norm_psi" 

pause -1


