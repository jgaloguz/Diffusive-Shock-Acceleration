#!/bin/bash

# configure code with autotools
autoreconf
automake --add-missing
./configure CXXFLAGS="-Ofast" --with-mpi=$1 --with-execution=PARALLEL --with-trajectory=PARKER --with-time_flow=$2 --with-rkmethod=0 --with-server=SELF

# make folder(s) to output results
mkdir runs_planar_shock/dsa_results
mkdir runs_spherical_shock/dsa_results