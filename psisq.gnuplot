set yrange [0:1]
set xlabel "x"

do for [IDX=1:999]{ 

plot "psisq.dat" index IDX u 2:3 w l  ls 3 notitle, "potential.dat" u 2:3 w l ls 1 notitle
#pause 0.1
}
