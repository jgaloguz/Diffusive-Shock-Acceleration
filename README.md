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
This command will post-process the results for all time, spatial, and momentum bins, since a single forward-in-time simulation, provided high enough statistics, will approximate the momentum spectrum at all times and all spatial locations.

**Backward-in-time Simulations**

Similarly, for backward-in-time simulations, configure the code using the `BACKWARD` option.
Compile and run the code using
```
make dsa_backward
mpirun -np <N> dsa_backward <number-of-trajectories> <batch-size> <time-index>
```
where `<N>`, `<number-of-trajectories>`, and `<batch-size>` mean the same as in the forward-in-time runs while `<time-index>` is the index of the time within the common time array at which to initialize the pseudo-particles.
The results can be post-processed with the command
```
make dsa_backward_postprocess
./dsa_backward_postprocess <time-index>
```
This command will only post-process the results for the time corresponding to `<time-index>` and the spatial location given by input *Parameter 4* (see list of parameters below).
This is because each backward-in-time run will only compute the spectrum at a single time and location.

**Parameters File**
The `params.dat` file contains the parameters that control the simulation execution and output.
Below is a list of all of the parameters in the order that they should be specified within the file (separated by white space/return) and their meaning/usage.
The units are mentioned within parentheses.

 - *Parameter 1*: maximum spatial displacement per step for the trajectories away from the shock (au).
 - *Parameter 2*: width of the shock (au).
 - *Parameter 3*: maximum spatial displacement per step for the trajectories near the shock given as a fraction of the shock width (unitless).
 - *Parameter 4*: spatial/radial location where to compute/plot spectrum (au).
 - *Parameter 5*: injection momentum (MeV).
 - *Parameter 6*: shock strength, which should be in the range (1,4] (unitless).
 - *Parameter 7*: upstream flow speed near the shock (cm / s).
 - *Parameter 8*: upstream diffusion coefficient near the shock (cm^2 / s).
 - *Parameter 9*: injection rate (1 / cm^2 / s)
 - *Parameter 10*: lower bound of momentum range (MeV). **Used only in spherical shock simulations**. For the planar shock problem, the injection momentum is used as the lower bound of the momentum range.
 - *Parameter 11*: upper bound of momentum range (MeV).
 - *Parameter 12*: Boolean flag (0 or 1) to indicate whether bins should be logarithmic (1) or linear (0). **Used only in spherical shock simulations**. For the planar shock problem, the spatial bins are always linear.
 - *Parameter 13*: lower bound of spatial/radius range (au).
 - *Parameter 14*: upper bound of spatial/radius range (au).
 - *Parameter 15*: lower bound of temporal range (day).
 - *Parameter 16*: upper bound of temporal range (day).
 - *Parameter 17*: radius of spherical shock (au). **Used only in spherical shock simulations**. For the planar shock problem, the shock is placed at x = 0.

 Note that the default binning resolutions of the momentum, spatial, and temporal ranges are 100, 100, and 5 bins respectively.
 These do not affect the physics or the execution time, and should be the same across all runs for a fair comparison.
 They can be manually changed in the `dsa_common.hh` file prior to compilation.

**Plotting Results**

All results are stored in the `runs/dsa_results` folder.
They can be visualized by running a Python script
```
python dsa_plots.py <Nt1> <Nt2> <which-variables> <which-time-flow>
```
where `<Nt1>` and `<Nt2>` are the lower and upper time indices specifying a subset of the results to be plotted, `<which-variables>` can be `none`, `pos`, `mom`, or `both`, and `<which-time-flow>` is either `forward` or `backward`.
In particular, setting `<Nt1> = 0` and `<Nt2> = 5` will plot all the available results for the default time resolution of 5 bins.
Using `<which-variables> = none` will only plot the analytic results (which are always plotted).
Using `<which-variables> = pos` will add the number density (spectrum integrated over momentum) vs space plots in the top panel, `<which-variables> = mom` will add the spectrum vs momentum plots at the location indicated by *Parameter 4* in the bottom panel, and `<which-variables> = both` will add both.
Finally, setting `<which-time-flow> = forward` will plot the results from the forward-in-time runs, while setting `<which-time-flow> = backward` will plot the results from the backward-in-time runs.

## Important note

**This is NOT the official SPECTRUM repository.** For information about SPECTRUM, go to https://github.com/vflorins/SPECTRUM.
