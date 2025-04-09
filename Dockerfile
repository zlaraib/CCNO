FROM nvidia/cuda:12.1.1-devel-ubuntu22.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y python3 python3-pip gfortran build-essential libhdf5-openmpi-dev \
                       openmpi-bin pkg-config libopenmpi-dev openmpi-bin libblas-dev \
                       liblapack-dev libpnetcdf-dev git python-is-python3 wget && \
    pip3 install numpy matplotlib h5py sympy scipy &&\
    rm -rf /var/lib/apt/lists/*
ENV USER=jenkins
ENV LOGNAME=jenkins
ENV JULIA_VERSION=1.11.4
RUN wget -q https://julialang-s3.julialang.org/bin/linux/x64/1.11/julia-${JULIA_VERSION}-linux-x86_64.tar.gz \
    && tar -xzf julia-${JULIA_VERSION}-linux-x86_64.tar.gz -C /opt \
    && mv /opt/julia-${JULIA_VERSION} /opt/julia \
    && ln -s /opt/julia/bin/julia /usr/local/bin/julia \
    && rm julia-${JULIA_VERSION}-linux-x86_64.tar.gz

# create directory where julia downloads and stores dependencies
ENV JULIA_DEPOT_PATH=/opt/julia
RUN mkdir -p $JULIA_DEPOT_PATH

# install some dependencies to hopefully make the jenkins process faster
RUN julia -e 'using Pkg; Pkg.add(["ITensors", "ITensorMPS", "Plots", "Measures", "LinearAlgebra", "ITensorTDVP", "HDF5"]); Pkg.precompile()'
RUN chmod -R 777 $JULIA_DEPOT_PATH