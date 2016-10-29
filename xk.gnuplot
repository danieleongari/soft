set xlabel "Steps"

  p "norm.dat" u 1:($3-25) w l lt 1 ti "xave-xeq"
rep "norm.dat" u 1:($4*5000) w l lt 3 ti "kave/5000"

pause -1
