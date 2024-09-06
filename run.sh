#!/bin/bash

usage() {
  echo "Usage: $0 -t 'time' -n 'num_processes'"
  exit 1
}

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

julia_exec=$(which julia)
if [ -z "$julia_exec" ]; then
  echo "Julia executable not found in PATH. Please ensure Julia is installed and available in your PATH."
  exit 1
fi

echo "Using Julia at: $julia_exec"

$julia_exec --project -e 'using Pkg; Pkg.update()'

for (( i=1; i<=num_processes; i++ ))
do
  $julia_exec --project -e "include(\"experiments/CalabiYau/char7/RunCalabiYau.jl\"); run_experiment(time = $time_arg)" &
done

wait
