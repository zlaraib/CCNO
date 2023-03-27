FROM nvidia/cuda:11.6.2-base-ubuntu20.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y python3 python3-pip gfortran build-essential libhdf5-openmpi-dev openmpi-bin pkg-config libopenmpi-dev openmpi-bin libblas-dev liblapack-dev libpnetcdf-dev git python-is-python3 julia
RUN pip3 install numpy matplotlib h5py sympy scipy
ENV USER=jenkins
ENV LOGNAME=jenkins
