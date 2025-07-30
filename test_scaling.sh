#!/bin/bash

output="scaling_results_parallel.txt"


N_values=(5000 10000 15000 20000 30000 40000 )
#25000 30000 35000 40000 45000 50000 75000 100000
for N in "${N_values[@]}"
do
    {
        runtime=$(/usr/bin/time -f "%e" ./app/fmm "$N" 2>&1 >/dev/null)
        echo "$N $runtime" >> "$output"
    } &
done

wait
echo "Appended results to $output"
