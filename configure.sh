#!/bin/bash

# configure code with autotools
autoreconf
automake --add-missing
if [ "$3" == "planar" ]; then
   ./configure CXXFLAGS="-Ofast" --with-mpi=$1 --with-execution=PARALLEL --with-trajectory=PARKER --with-time_flow=$2 --with-rkmethod=0 --with-server=SELF
   mkdir -p runs_planar_shock/dsa_results
   echo "Code configured for:"
   echo "    forward planar shock simulations;"
   echo "    backward planar shock simulations;"
   echo "    backward spherical simulations without a source."
elif [ "$3" == "spherical" ]; then
   ./configure CXXFLAGS="-Ofast" --with-mpi=$1 --with-execution=PARALLEL --with-trajectory=PARKER_SOURCE --with-time_flow=$2 --with-rkmethod=0 --with-server=SELF
   mkdir -p runs_spherical_shock/dsa_results
   echo "Code configured for:"
   echo "    forward spherical shock simulations with a source;"
   echo "    backward spherical shock simulations with a source."
else
   echo "Shock geometry parameter unrecognized."
   echo "Valid options are: 'planar' or 'spherical'."
fi
