unset key
set xrange [0:2000]
set yrange [700:0]
plot "goose_positions.txt" using 3:2 with lines lt rgb "black", "winterbreedingposition.txt" using 2:1 with points pt 3 ps 5
set yrange [0:1000]
unset xrange
plot "goose_positions.txt" using 1:4 with lines
