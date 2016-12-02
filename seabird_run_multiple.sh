#Script to run seabird simulation separately on each core of a machine.

initial_lat=
initial_lon=
a=
kT=
start_month=
end_month=

num_cores=$(getconf _NPROCESSORS_ONLN)
echo "Number of cores: "$num_cores
i=1
while ((i<=num_cores)); do
  python3 seabird_sim.py initial_lat initial_lon a kT start_month end_month $i >> run$i.txt &&
done
