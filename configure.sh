#!/bin/bash

# Check arguments
if [ "$2" == "forward" ]; then
   TIME_FLOW="FORWARD"
elif [ "$2" == "backward" ]; then
   TIME_FLOW="BACKWARD"
else
   echo ""
   echo "****************************************"
   echo "ERROR: Invalid time flow parameter."
   echo "Valid options are: 'forward' or 'backward'."
   echo "****************************************"
   echo ""
   exit 1
fi

if [ "$3" == "planar" ]; then
   TRAJECTORY_TYPE="PARKER"
   mkdir -p runs_planar_shock/dsa_results
   mkdir -p runs_planar_shock_acc/dsa_results
elif [ "$3" == "spherical" ]; then
   TRAJECTORY_TYPE="PARKER_SOURCE"
   mkdir -p runs_spherical_shock/dsa_results
else
   echo ""
   echo "****************************************"
   echo "ERROR: Invalid shock geometry parameter."
   echo "Valid options are: 'planar' or 'spherical'."
   echo "****************************************"
   echo ""
   exit 2
fi

# configure code with autotools
autoreconf
automake --add-missing
./configure CXXFLAGS="-Ofast" --with-mpi=${1} --with-execution=PARALLEL --with-trajectory=${TRAJECTORY_TYPE} --with-time_flow=${TIME_FLOW} --with-rkmethod=0 --with-server=SELF

# Output status if successful
if [ $? == 0 ]; then
   echo ""
   echo "****************************************"
   echo "Code configured for:"
   echo "    ${2} ${3} shock simulations."
   echo "****************************************"
   echo ""
fi