#!/bin/bash

usage() {
  echo "Usage: $0 -t 'time'"
  exit 1
}

while getopts "t:" opt; do
  case ${opt} in
    t )
      time_arg=$OPTARG
      ;;
    \? )
      usage
      ;;
  esac
done

if [ -z "$time_arg" ]; then
  usage
fi

julia --project -e 'using Pkg; Pkg.update()'
julia --project -e "include(\"experiments/CalabiYau/RunCalabiYau.jl\"); run_experiment(time = $time_arg)" &
julia --project -e "include(\"experiments/CalabiYau/RunCalabiYau.jl\"); run_experiment(time = $time_arg)" &
julia --project -e "include(\"experiments/CalabiYau/RunCalabiYau.jl\"); run_experiment(time = $time_arg)" &
julia --project -e "include(\"experiments/CalabiYau/RunCalabiYau.jl\"); run_experiment(time = $time_arg)" &
