#!/bin/bash

usage() {
  echo "Usage: $0 -t 'time' -n 'num_processes'"
  exit 1
}

# Default number of processes if not specified
num_processes=5

while getopts "t:n:" opt; do
  case ${opt} in
    t )
      time_arg=$OPTARG
      ;;
    n )
      num_processes=$OPTARG
      ;;
    \? )
      usage
      ;;
  esac
done

if [ -z "$time_arg" ]; then
  usage
fi

# Try to find the Julia executable and update the PATH if necessary
julia_exec=$(which julia)
if [ -z "$julia_exec" ]; then
  echo "Julia executable not found in PATH. Please ensure Julia is installed and available in your PATH."
  exit 1
fi

echo "Using Julia at: $julia_exec"

# Update Julia packages
$julia_exec --project -e 'using Pkg; Pkg.update()'

# Loop to run multiple processes
for (( i=1; i<=num_processes; i++ ))
do
  $julia_exec --project -e "include(\"experiments/CalabiYau/char7/RunCalabiYau.jl\"); run_experiment(time = $time_arg)" &
done

# Wait for all background processes to finish
wait
