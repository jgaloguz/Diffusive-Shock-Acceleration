#!/bin/bash

autoreconf
automake --add-missing
./configure CXXFLAGS="-Ofast" --enable-silo --with-mpi=$1 --with-execution=PARALLEL --with-trajectory=PARKER --with-time_flow=$2 --with-rkmethod=0 --with-server=SELF