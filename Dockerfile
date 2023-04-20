FROM nvidia/cuda:11.6.2-devel-ubuntu20.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y python3 python3-pip gfortran build-essential libhdf5-openmpi-dev openmpi-bin pkg-config libopenmpi-dev openmpi-bin libblas-dev liblapack-dev libpnetcdf-dev git python-is-python3 wget

# Download and install Julia 1.8.5
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.5-linux-x86_64.tar.gz && \
    tar zxvf julia-1.8.5-linux-x86_64.tar.gz && \
    rm julia-1.8.5-linux-x86_64.tar.gz && \
    mv julia-1.8.5 /opt/julia 
# RUN apt-get install -y python3 python3-pip gfortran build-essential libhdf5-openmpi-dev openmpi-bin pkg-config libopenmpi-dev openmpi-bin libblas-dev liblapack-dev libpnetcdf-dev git python-is-python3 julia wget
RUN pip3 install numpy matplotlib h5py sympy scipy
ENV USER=jenkins
ENV LOGNAME=jenkins
ENV JULIA_DEPOT_PATH=/opt/julia:$JULIA_DEPOT_PATH
ENV JULIA_LOAD_PATH=/opt/julia:$JULIA_LOAD_PATH
RUN julia -e 'import Pkg; Pkg.add("ITensors"); Pkg.precompile()'
RUN chmod a+rX -R /opt/julia
