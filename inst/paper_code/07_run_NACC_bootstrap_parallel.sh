#!/bin/bash

for i in {1..50}
do
  Rscript NACC_bootstrap.R "$i" > "out_${i}.txt" 2>&1 &  # Pass {i} as argument and save output to out_{i}.txt
done
