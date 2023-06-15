FROM nvidia/cuda:12.1.1-devel-ubuntu22.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y python3 python3-pip gfortran build-essential libhdf5-openmpi-dev openmpi-bin pkg-config libopenmpi-dev openmpi-bin libblas-dev liblapack-dev libpnetcdf-dev git python-is-python3 wget
RUN pip3 install numpy matplotlib h5py sympy scipy
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.7-linux-x86_64.tar.gz
RUN tar -zxvf julia-1.6.7-linux-x86_64.tar.gz
ENV PATH=/julia-1.6.7/bin:$PATH
ENV USER=jenkins
ENV LOGNAME=jenkins
ENV JULIA_DEPOT_PATH=/opt/julia:$JULIA_DEPOT_PATH
ENV JULIA_LOAD_PATH=/opt/julia:$JULIA_LOAD_PATH
RUN julia -e 'import Pkg; Pkg.add("ITensors"); Pkg.add("Plots"); Pkg.add("Measures"); Pkg.precompile()' 
RUN chmod a+rX -R /opt/julia