ARG IMAGE=nvidia/cuda:12.1.1-devel-ubuntu20.04
FROM $IMAGE

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get install --yes --no-install-recommends \
                    curl ca-certificates vim git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# julia
ARG JULIA_RELEASE=1.12
ARG JULIA_VERSION=1.12.4
RUN curl -s -L https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_RELEASE}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz | \
    tar -C /usr/local -x -z --strip-components=1 -f -

# Setup depot with open permissions
RUN mkdir -m 777 /data
ENV JULIA_DEPOT_PATH=/data
ENV JULIA_HISTORY=/data/logs/repl_history.jl

# Setup workdir
WORKDIR /app
