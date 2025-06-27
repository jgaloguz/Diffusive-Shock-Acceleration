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

After a successful configuration, to compile and run code, first navigate to either the `runs_planar_shock` subdirectory for simulations using a planar shock
```
cd runs_planar_shock
```
or the `runs_spherical_shock` subdirectory for simulations using a spherical shock
```
cd runs_spherical_shock
```
In either folder, the code is run the same way.
Here are some notes on the differences between the codes.

 - In the case of a planar shock, the diffusion coefficient is proportional to the square of the flow speed everywhere. In the spherical shock, this is only true within the shock and the downstream region, with the upstream having a diffusion coefficient that is proportional to the radial coordinate.
 - The planar shock has a constant flow speed on both sides with the downstream speed equaling the upstream speed divided by the shock strength, with a hyperbolic tangent transition between them. In the spherical shock, the radial speed is constant upstream and decreases with the square of the radial coordinate downstream, with a hyperbolic tangent transition as well.

**Analytic Solutions**

To obtain analytic results, compile and run the code using
```
make dsa_analytic
./dsa_analytic
```

**Forward-in-time Simulations**

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

**Backward-in-time Simulations**

Similarly, for backward-in-time simulations, configure the code using the `BACKWARD` option.
Compile and run the code using
```
make dsa_backward
mpirun -np <N> dsa_backward <number-of-trajectories> <batch-size> <time-index>
```
where the first two options mean the same as in the forward-in-time runs and the third is the index of the time within the common time array at which to initialize the pseudo-particles.
The results can be post-processed with the command
```
make dsa_backward_postprocess
./dsa_backward_postprocess
```
where no options are needed in this case.

The `params.dat` file contains the parameters that control the simulation execution and output.
The first parameter equals the maximum spatial displacement per step for the trajectories away from the shock in units of au.
The second parameter equals the width of the shock in units of au.
The third parameter equals the maximum spatial displacement per step for the trajectories near the shock as a fraction of the first parameter.
The fourth parameter equals the initial spatial location of the pseudo-particles for the backward runs as a factor of `z_shock` or `r_shock`, which is the location of the bin center closest to the shock front on the downstream region.

**Plotting Results**

All results are stored in the `runs/dsa_results` folder.
They can be visualized by running a Python script
```
python dsa_plots.py <which-time-flow>
```
where `<which-time-flow>` is either `forward` or `backward`, depending on which result you want to compare with analytic results.

## Important note

**This is NOT the official SPECTRUM repository.** For information about SPECTRUM, go to https://github.com/vflorins/SPECTRUM.
