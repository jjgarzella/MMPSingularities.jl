#!/usr/bin/env bash

# Check that an argument was provided
if [ -z "$1" ]; then
    echo "Usage: $0 {3|5|7|11|13}"
    exit 1
fi

VALUE="$1"

# Map input value to thread count
case "$VALUE" in
    3|5|7)
        THREADS=5
        ;;
    11|13)
        THREADS=2
        ;;
    *)
        echo "Invalid value: $VALUE"
        echo "Allowed values: 3, 5, 7, 11, 13"
        exit 1
        ;;
esac

# Set up project
julia --project=cluster/ -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'

# Run the actual job with selected thread count
julia --threads "$THREADS" --project=cluster/ -e "include(\"cluster/run.jl\"); run(4, $1)" &
