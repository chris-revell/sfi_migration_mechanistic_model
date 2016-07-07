set xrange [0:2000]
set yrange [700:0]
unset key
plot "NDVI_data/goose_positions.txt" using 3:2 with lines lt rgb "black", "NDVI_data/winterbreedingposition.txt" using 2:1 with points pt 3 ps 5
