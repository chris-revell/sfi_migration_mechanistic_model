unset key
set xrange [0:2000]
set yrange [700:0]
plot "output_data/goose_positions.txt" using 3:2 with lines lt rgb "black", "output_data/winterbreedingposition.txt" using 2:1 with points pt 3 ps 5
