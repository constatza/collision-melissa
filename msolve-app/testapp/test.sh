#!/bin/bash

# Define the number of timesteps
timesteps=10

# Define the number of numbers to print
numbers_to_print=15402

# Define the wait time in seconds
wait_time=5

# Loop through each timestep
for (( timestep=1; timestep<=$timesteps; timestep++ )); do
    echo "Timestep $timestep:"

    # Loop to print the specified number of numbers
    for (( i=1; i<=$numbers_to_print; i++ )); do
        echo $i
    done

    # Wait for the specified time before proceeding to the next timestep
    sleep $wait_time
done

