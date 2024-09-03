#!/bin/bash

usage() {
  echo "Usage: $0 -t 'time'"
  exit 1
}


julia --project -e 'using Pkg; Pkg.update()'
julia --project -e "include(\"experiments/CalabiYau/char7/CPURunCalabiYau.jl\")" &
