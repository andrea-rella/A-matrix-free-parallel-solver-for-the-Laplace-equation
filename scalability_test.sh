#!/bin/bash

# Test scalability of the Laplace solver

# Ensure the script exits if any command fails
set -e

# Compile the program with the _SCALABILITY environment variable set
make clean
_SCALABILITY=1 make

# Specify the path to the input file
INPUT_FILE=data1.json

# Define the number of processors and grid sizes for testing
processor_counts=(1 2 4 8)
cell_sizes="(15,31,63,127,255)"

# Output file for results
output_file="scalability_results/scalability_results.txt"
mkdir -p scalability_results
echo "Scalability Test Results" > $output_file
echo "========================" >> $output_file

# Run the tests and capture the output
for np in "${processor_counts[@]}"; do  
        echo "Running test with ${np} processors"
        echo "-------------------------------------" >> $output_file
        mpirun -np ${np} ./laplace_solver ${INPUT_FILE} ${cell_sizes} >> $output_file
        echo "" >> $output_file    
done

# Summarize results
echo "Scalability test completed. Results are stored in $output_file."