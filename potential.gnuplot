#set yrange [*:2]
#set xrange [5:40]
set xlabel "X"
set ylabel "E"

p "potential.dat" u 2:3 w l lt rgb "green" lw 3 t "V(x)" 


pause -1
