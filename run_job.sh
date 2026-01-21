#!/bin/bash
# When running with a volume mount, the Manifest.toml inside the container is hidden.
# We need to re-instantiate the environment to ensure Julia knows where to find the pre-installed packages.
julia --project=cluster -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'

# Run the actual job
exec julia --threads 4 --project=cluster/ -e 'include("cluster/run.jl")'
