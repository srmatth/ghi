#!/bin/bash

# Loop to run the R script 100 times, passing the argument and saving output to a file
for i in {1..50}
do
    Rscript ghi_sim_bootstrapped_min_resid.R "$i" > "out_${i}.txt" 2>&1 &  # Pass {i} as argument and save output to out_{i}.txt
done
