# Diffusive Shock Acceleration

This is a specialization of the SPECTRUM software applied to modeling diffusive shock acceleration. Specifically, the codes in this repository simulate the acceleration of test-particles in the presence of a planar or spherical shock for a variety of diffusion profiles using stochastic methods, analyze the results, and generate useful figures for a future scientific publication. In addition, both forward and backward time-flow approaches are used to solve the relevant equations so their results can be compared and their relative advantages or disadvantes explored.

## Using the code

Before first use, the code must be configured with autotools. After cloning this repository, execute the configure script in the working directory
```
git clone https://github.com/jgaloguz/Diffusive-Shock-Acceleration
cd Diffusive-Shock-Acceleration
./configure.sh <mpi-option> <time-flow>
where `<mpi-option>` is either `openmpi` or `mpich`, whichever is installed in your system, and `<time-flow>` is either `FORWARD` or `BACKWARD`, depending on which algorithm you want to use to solve the relevant transport equations. You may have to change the permissions of `configure.sh` before you can execute it. You will know the configuration stage ran successfully if a `config.h` file was generated in the working directory.
```


## Important note

**This is NOT the official SPECTRUM repository.** For information about SPECTRUM, go to https://github.com/vflorins/SPECTRUM.
