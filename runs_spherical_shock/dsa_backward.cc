#include "src/simulation.hh"
#include "src/distribution_other.hh"
#include "src/background_solarwind_termshock.hh"
#include "src/diffusion_other.hh"
#include "src/boundary_time.hh"
#include "src/boundary_momentum.hh"
#include "src/initial_time.hh"
#include "src/initial_space.hh"
#include "src/initial_momentum.hh"
#include "dsa_common.hh"
#include <iostream>
#include <iomanip>

using namespace Spectrum;

int main(int argc, char** argv)
{
   int i;
   DataContainer container;
   ReadParams();
   DefineArrays();

//--------------------------------------------------------------------------------------------------
// Create a simulation object
//--------------------------------------------------------------------------------------------------

   std::unique_ptr<SimulationWorker> simulation;
   simulation = CreateSimulation(argc, argv);

//--------------------------------------------------------------------------------------------------
// Particle type
//--------------------------------------------------------------------------------------------------

   simulation->SetSpecie(specie);

//--------------------------------------------------------------------------------------------------
// Background
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Initial time
   container.Insert(0.0);

// Origin
   container.Insert(gv_zeros);

// Upstream velocity
   GeoVector u0(U_up, 0.0, 0.0);
   container.Insert(u0);

// Upstream magnetic field
   double RS = 6.957e10 / unit_length_fluid;
   double r_ref = 3.0 * RS;
   double BmagE = 5.0e-5 / unit_magnetic_fluid;
   double Bmag_ref = BmagE * Sqr(one_au / r_ref);
   GeoVector B0(Bmag_ref, 0.0, 0.0);
   container.Insert(B0);

// Maximum displacement
   container.Insert(dmax);

// Solar rotation vector
   GeoVector Omega(0.0, 0.0, 0.0);
   container.Insert(Omega);

// Reference equatorial distance
   container.Insert(r_ref);

// dmax fraction
   container.Insert(dmax_fraction);

// Spherical shock radius
   container.Insert(R_sh); 

// Spherical shock width
   container.Insert(w_sh);
   
// Spherical shock strength
   container.Insert(s);

   simulation->AddBackground(BackgroundSolarWindTermShock(), container);

//--------------------------------------------------------------------------------------------------
// Time initial condition
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Time to obtain solution
   int t_idx = 0;
   if(argc > 3) {
      t_idx = atoi(argv[3]);
      if(t_idx >= Nt) t_idx = Nt - 1;
   };
   container.Insert(t_arr[t_idx]);

   simulation->AddInitial(InitialTimeFixed(), container);

//--------------------------------------------------------------------------------------------------
// Spatial initial condition
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Point to obtain solution
   GeoVector init_pos(r_spectrum, 0.0, 0.0);
   container.Insert(init_pos);

   simulation->AddInitial(InitialSpaceFixed(), container);

//--------------------------------------------------------------------------------------------------
// Momentum initial condition
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Lower bound for momentum
   container.Insert(p0);

// Upper bound for momentum
   container.Insert(pf);

// Log bias
   bool log_bias = true;
   container.Insert(log_bias);

   simulation->AddInitial(InitialMomentumThickShell(), container);

//--------------------------------------------------------------------------------------------------
// Time limit
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Maximum crossings
   int max_crossings_time = 1;
   container.Insert(max_crossings_time);

// Actions vector
   std::vector<int> actions_time;
   actions_time.push_back(-1);
   actions_time.push_back(-1);
   actions_time.push_back(1);
   container.Insert(actions_time);

// Initial time at which to stop integrating
   container.Insert(0.0);

   simulation->AddBoundary(BoundaryTimeExpire(), container);

//--------------------------------------------------------------------------------------------------
// Injection boundary
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Maximum crossings
   int max_crossings_mom = 1;
   container.Insert(max_crossings_mom);

// Actions vector
   std::vector<int> actions_mom;
   actions_mom.push_back(0);
   actions_mom.push_back(0);
   actions_mom.push_back(0);
   container.Insert(actions_mom);

// Injection momentum
   container.Insert(p_inj);

// Center of spherical shells
   container.Insert(gv_zeros);

// Radius of inner shell
   double inner_shell_radius = R_sh;
   container.Insert(inner_shell_radius);

// Radius of outer shell
   double outer_shell_radius = R_sh + w_sh;
   container.Insert(outer_shell_radius);

   simulation->AddBoundary(BoundaryMomentumInjectRestrictShell(), container);

//--------------------------------------------------------------------------------------------------
// Diffusion model
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Reference diffusion coefficient
   container.Insert(kappa_up);

// Normalization of kinetic energy
   double T0 = one_MeV;
   container.Insert(T0);

// Normalization of radius
   container.Insert(R_sh);

// Power of kinetic energy dependance
   double power_law_T = 0.0;
   container.Insert(power_law_T);

// Power of radial dependance
   double power_law_r = 1.0;
   container.Insert(power_law_r);

// Ratio of perpendicular to parallel diffusion
   double kap_rat = 0.0;
   container.Insert(kap_rat);

// Downstream dependance index
   int stream_dep_idx = 1;
   container.Insert(stream_dep_idx);

// Upstream flow at the start of shock
   container.Insert(U_up);

// Shock position
   container.Insert(R_sh);

// Shock width
   container.Insert(w_sh);

   simulation->AddDiffusion(DiffusionKineticEnergyRadialDistancePowerLaw(), container);

//--------------------------------------------------------------------------------------------------
// Distribution 1 (time)
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Number of bins
   MultiIndex n_bins1(100, 0, 0);
   container.Insert(n_bins1);
   
// Smallest value
   GeoVector minval1(0.0, 0.0, 0.0);
   container.Insert(minval1);

// Largest value
   GeoVector maxval1(tf, 0.0, 0.0);
   container.Insert(maxval1);

// Linear or logarithmic bins
   MultiIndex log_bins1(0, 0, 0);
   container.Insert(log_bins1);

// Add outlying events to the end bins
   MultiIndex bin_outside1(0, 0, 0);
   container.Insert(bin_outside1);

// Physical units of the distro
   double unit_distro1 = 1.0;
   container.Insert(unit_distro1);

// Physical units of the bin value
   GeoVector unit_val1 = {unit_time_fluid, 1.0, 1.0};
   container.Insert(unit_val1);

// Don't keep records
   bool keep_records1 = false;
   container.Insert(keep_records1);

// Value for the "hot" condition
   double val_hot1 = 1.0;
   container.Insert(val_hot1);

// Value for the "cold" condition
   double val_cold1 = 0.0;
   container.Insert(val_cold1);

// Coordinates to use (initial or final)
   int val_time1 = 1;
   container.Insert(val_time1);

   simulation->AddDistribution(DistributionTimeUniform(), container);

//--------------------------------------------------------------------------------------------------
// Distribution 2 (position)
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Number of bins
   MultiIndex n_bins2(Nr, 0, 0);
   container.Insert(n_bins2);

// Smallest value
   GeoVector minval2(R_sh, 0.0, 0.0);
   container.Insert(minval2);

// Largest value
   GeoVector maxval2(R_sh + w_sh, 0.0, 0.0);
   container.Insert(maxval2);

// Linear or logarithmic bins
   MultiIndex log_bins2(0, 0, 0);
   container.Insert(log_bins2);

// Add outlying events to the end bins
   MultiIndex bin_outside2(0, 0, 0);
   container.Insert(bin_outside2);

// Physical units of the distro
   double unit_distro2 = 1.0;
   container.Insert(unit_distro2);

// Physical units of the bin values
   GeoVector unit_val2 = {unit_length_fluid, 1.0, 1.0};
   container.Insert(unit_val2);

// Don't keep records
   bool keep_records2 = false;
   container.Insert(keep_records2);

// Constant value for the "hot" condition
   double val_hot2 = 1.0;
   container.Insert(val_hot2);

// Constant value for the "cold" condition
   double val_cold2 = 0.0;
   container.Insert(val_cold2);

// Which coordinates to use for value: 0 initial, 1 final
   int val_time2 = 1;
   container.Insert(val_time2);

// Which coordinate representation to use for value: 0 "native coordinates", 1 locally spherical with B || z
   int val_coord2 = 0;
   container.Insert(val_coord2);

   simulation->AddDistribution(DistributionPositionUniform(), container);

//--------------------------------------------------------------------------------------------------
// Distribution 3 (momentum)
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Number of bins
   MultiIndex n_bins3(Np, 0, 0);
   container.Insert(n_bins3);

// Smallest value
   GeoVector minval3(p0, 0.0, 0.0);
   container.Insert(minval3);

// Largest value
   GeoVector maxval3(pf, 0.0, 0.0);
   container.Insert(maxval3);

// Linear or logarithmic bins
   MultiIndex log_bins3(1, 0, 0);
   container.Insert(log_bins3);

// Add outlying events to the end bins
   MultiIndex bin_outside3(0, 0, 0);
   container.Insert(bin_outside3);

// Physical units of the distro variable
   double unit_distro3 = 1.0;
   container.Insert(unit_distro3);

// Physical units of the bin variable which is momentum here. This is for x axis.
   GeoVector unit_val3 = {unit_momentum_particle, 1.0, 1.0};
   container.Insert(unit_val3);

// Don't keep records
   bool keep_records3 = false;
   container.Insert(keep_records3);

// Constant value for the "hot" condition
   double val_hot3 = 1.0;
   container.Insert(val_hot3);

// Constant value for the "cold" condition
   double val_cold3 = 0.0;
   container.Insert(val_cold3);

// Which coordinates to use for value: 0 initial, 1 final
   int val_time3 = 0;
   container.Insert(val_time3);

// Which coordinate representation to use for value: 0 "native coordinates", 1 locally spherical with B || z
   int val_coord3 = 0;
   container.Insert(val_coord3);

   simulation->AddDistribution(DistributionMomentumUniform(), container);

//--------------------------------------------------------------------------------------------------
// Run the simulation
//--------------------------------------------------------------------------------------------------

   int n_traj;
   int batch_size;

   batch_size = n_traj = 1;
   if(argc > 1) n_traj = atoi(argv[1]);
   if(argc > 2) batch_size = atoi(argv[2]);

   std::string simulation_files_prefix = "dsa_results/dsa_backward_";
   simulation->DistroFileName(simulation_files_prefix);
   simulation->SetTasks(n_traj, batch_size);
   simulation->MainLoop();
   simulation->PrintDistro1D(0, 0, simulation_files_prefix + "time_" + std::to_string(t_idx) + ".dat", false);
   simulation->PrintDistro1D(1, 0, simulation_files_prefix + "pos_" + std::to_string(t_idx) + ".dat", false);
   simulation->PrintDistro1D(2, 0, simulation_files_prefix + "mom_" + std::to_string(t_idx) + ".dat", false);

   return 0;
};


