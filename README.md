# Diffusive Shock Acceleration

This is a specialization of the SPECTRUM software applied to modeling diffusive shock acceleration. Specifically, the codes in this repository simulate the acceleration of test-particles in the presence of a planar or spherical shock for a variety of diffusion profiles using stochastic methods, analyze the results, and generate useful figures for a future scientific publication. In addition, both forward and backward time-flow approaches are used to solve the relevant equations so their results can be compared and their relative advantages or disadvantes explored.

## Using the code

Before first use, the code must be configured with autotools. After cloning this repository, execute the configure script in the working directory
```
git clone https://github.com/jgaloguz/Diffusive-Shock-Acceleration
cd Diffusive-Shock-Acceleration
./configure.sh <mpi-option> <time-flow>
```
where `<mpi-option>` is either `openmpi` or `mpich`, whichever is installed in your system, and `<time-flow>` is either `FORWARD` or `BACKWARD`, depending on which algorithm you want to use to solve the relevant transport equations. You may have to change the permissions of `configure.sh` before you can execute it. You will know the configuration stage ran successfully if a `config.h` file was generated in the working directory.

**Planar Shock**

To compile and run code, first navigate to the `runs` subdirectory
```
cd runs
```
To obtain analytic results in the case of a planar shock with a diffusion coefficient that is proportional to the square of the flow speed, compile and run the code using
```
make dsa_analytic
./dsa_analytic
```

For forward-in-time simulations, configure the code using the `FORWARD` option.
Compile and run the code using
```
make dsa_forward
mpirun -np <N> dsa_forward <number-of-trajectories> <batch-size>
```
where `<N>` is the number of processors you want to use for running the code, `<number-of-trajectories>` is the total number of trajectories to average over during this run, and `<batch-size>` is the number of trajectories assigned per batch as each processor becomes available to perform work.
The results can be post-processed with the command
```
make dsa_forward_postprocess
./dsa_forward_postprocess <number-of-trajectories>
```

Similarly, for backward-in-time simulations, configure the code using the `BACKWARD` option.
Compile and run the code using
```
make dsa_backward
mpirun -np <N> dsa_backward <number-of-trajectories> <batch-size>
```
where the options mean the same as in the forward-in-time runs.
The results can be post-processed with the command
```
make dsa_backward_postprocess
./dsa_backward_postprocess
```
where no options are needed in this case.

The `runs/params.dat` file contains the parameters that control the simulation execution and output.
The first parameter controls maximum spatial displacement per step for the trajectories away from the shock.
The second parameter controls the width of the shock encountered by the particles.
The third parameter controls maximum spatial displacement per step for the trajectories near the shock.
The fourth parameter controls the initial spatial location of the particles for the backward runs.

All results are stored in the `runs/dsa_results` folder.
They can be visualized by running a Python script
```
python dsa_plots.py <which-time-flow>
```
where `<which-time-flow>` is either `forward` or `backward`, depending on which result you want to compare with analytic results.

## Important note

**This is NOT the official SPECTRUM repository.** For information about SPECTRUM, go to https://github.com/vflorins/SPECTRUM.
