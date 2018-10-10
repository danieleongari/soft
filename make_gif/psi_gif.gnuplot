set yrange [-0.7:0.7]
set xlabel "x"

set term png size 600,300

do for [IDX=0:99]{ 
plot  "potential.dat"     u ($2):($3/10) w l lt rgb "green" lw 3 notitle,\
      "psi.dat" index IDX u 2:3          w l lt rgb "blue"  lw 1 ti "Re",\
      "psi.dat" index IDX u 2:4          w l lt rgb "red"   lw 1 ti "Im",\
      "psi.dat" index IDX u 2:5          w l lt rgb "black" lw 1 ti "sqr"

#pause -1

}
