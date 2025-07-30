#!/bin/bash

output="scaling_results_parallel.txt"

# Add header only if file does not exist
if [ ! -f "$output" ]; then
    echo "N runtime_seconds" > "$output"
fi

N_values=(5000 8000 10000 15000 )
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
