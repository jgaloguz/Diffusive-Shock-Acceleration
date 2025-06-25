#include "src/simulation.hh"
#include "src/distribution_other.hh"
#include "src/background_smooth_shock.hh"
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
   DefineArrays();
   ReadParams();

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
   double B_up = 5.0e-7 / unit_magnetic_fluid;
   GeoVector B0(B_up, 0.0, 0.0);
   container.Insert(B0);

// Maximum displacement
   double one_au = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   double dmax = params[0] * one_au;
   container.Insert(dmax);

// Shock starting position
   container.Insert(gv_zeros);

// Shock normal
   GeoVector n_shock (-1.0, 0.0, 0.0);
   container.Insert(n_shock);

// Shock velocity
   double v_shock = 0.0;
   container.Insert(v_shock);

// Downstream velocity
   GeoVector u1 (U_dn, 0.0, 0.0);
   container. Insert(u1);

// Downstream magnetic field
   double B_dn = B_up * s;   
   container. Insert(B_dn);

// Shock width
   double width_shock = params[1] * one_au;
   container.Insert(width_shock);

// dmax fraction
   double dmax_fraction = params[2];
   container.Insert(dmax_fraction);

   simulation->AddBackground(BackgroundSmoothShock(), container);

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
   GeoVector init_pos(params[3] * z_shock, 0.0, 0.0);
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
   container.Insert(p0);

   simulation->AddBoundary(BoundaryMomentumInject(), container);

//--------------------------------------------------------------------------------------------------
// Diffusion model
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Reference diffusion coefficient
   container.Insert(kappa_up);

// Normalization of bulk velocity
   container.Insert(U_up);

// Power of bulk velocity dependance
   double power_law_U = 2.0;
   container.Insert(power_law_U);

// Normalization of particle momentum
   double p_0 =  0.1 * mass[specie] * c_code;
   container.Insert(p_0);

// Power of particle momentum dependance
   double power_law_p = 0.0;
   container.Insert(power_law_p);

// Ratio of perpendicular to parallel diffusion
   double kap_rat = 0.0;
   container.Insert(kap_rat);

   simulation->AddDiffusion(DiffusionFlowMomentumPowerLaw(), container);

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
   MultiIndex n_bins2(Nz, 0, 0);
   container.Insert(n_bins2);

// Smallest value
   GeoVector minval2(-width_shock, 0.0, 0.0);
   container.Insert(minval2);

// Largest value
   GeoVector maxval2(width_shock, 0.0, 0.0);
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


