#set yrange [-0.3:0.3]
set xlabel "sx"

do for [IDX=0:999]{ 

plot  "psi.dat" index IDX u 1:5 w l  ls 3 ti "Re",\
      "psi.dat" index IDX u 1:6 w l  ls 4 ti "Im",\
      "psisq.dat" index IDX u 1:4 w l lc "k" ti "sqr"

pause 0.1

}
