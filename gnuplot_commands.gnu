set terminal png
unset key
set output output_data/20167171730/path.png
set xrange [0:2000]
set yrange [700:0]
plot "goose_positions.txt" using 3:2 with lines lt rgb "black", "winterbreedingposition.txt" using 2:1 with points pt 3 ps 5
set yrange [0:1000]
unset xrange
set output output_data/20167171730/distance.png
plot "goose_positions.txt" using 1:4 with lines
