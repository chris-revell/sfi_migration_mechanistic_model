#Script to run seabird simulation separately on each core of a machine.

initial_lat=-54.2
initial_lon=-36.3
a=0.1
kT=0.1
start_month=1
end_month=3

num_cores=$(getconf _NPROCESSORS_ONLN)
echo "Number of cores: "$num_cores
i=1
while ((i<=num_cores)); do
  python3 seabird_sim.py $initial_lat $initial_lon $a $kT $start_month $end_month $i >> run$i_out.txt &
  let i++
done
echo "Running "$i" simulations"
wait
