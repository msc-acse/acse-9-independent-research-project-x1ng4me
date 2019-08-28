#!/bin/bash
# Xianzheng Li, ICL
# GitHub: x1ng4me

# Set how many threads should be used in a test
export OMP_NUM_THREADS=6
# Set OMP_PINNING to bind all threads on one socket
export OMP_PROC_BIND=true
echo $OMP_NUM_THREADS

#Run iSALE, output with a txt file
echo "Running iSALE!"
./iSALE2D -i asteroid_6.inp >> vtune_output_6.txt
echo "Finished!"
