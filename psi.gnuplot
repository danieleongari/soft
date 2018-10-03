set yrange [-1:1]
set xlabel "x"

do for [IDX=0:999]{ 

plot  "potential.dat"     u 2:3 w l lt rgb "green" lw 3 notitle,\
      "psi.dat" index IDX u 2:3 w l lt rgb "blue"  lw 1 ti "Re",\
      "psi.dat" index IDX u 2:4 w l lt rgb "red"   lw 1 ti "Im",\
      "psi.dat" index IDX u 2:5 w l lt rgb "black" lw 1 ti "sqr"

pause -1

}
