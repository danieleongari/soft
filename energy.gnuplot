set yrange [*:*]
set xrange [*:*]
set xlabel "step"
set ylabel "energies"

  p "energy.dat" u 1:2 w l ti "E_kin" 
rep "energy.dat" u 1:3 w l ti "E_pot"
rep "energy.dat" u 1:4 w l ti "E_tot"

pause -1
