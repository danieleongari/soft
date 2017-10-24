set xrange [-10:+10]
set yrange [-0.3:+0.3]
set xlabel "k"

do for [IDX=0:999]{ 

plot  "psip.dat" index IDX u 2:3 w l  ls 3 ti "Re",\
      "psip.dat" index IDX u 2:4 w l  ls 4 ti "Im",\
      "psip.dat" index IDX u 2:5 w l lc "k" ti "sqr"

pause -1

}
