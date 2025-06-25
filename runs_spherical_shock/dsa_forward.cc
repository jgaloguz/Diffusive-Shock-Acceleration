#include "src/simulation.hh"
#include "src/distribution_other.hh"
#include "src/background_solarwind_termshock.hh"
#include "src/diffusion_other.hh"
#include "src/boundary_time.hh"
#include "src/boundary_space.hh"
#include "src/initial_time.hh"
#include "src/initial_space.hh"
#include "src/initial_momentum.hh"
#include "dsa_common.hh"
#include <iostream>
#include <iomanip>

using namespace Spectrum;

int main(int argc, char** argv)
{
   int i, j;
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
   double RS = 6.957e10 / unit_length_fluid;
   double r_ref = 3.0 * RS;
   double BmagE = 5.0e-5 / unit_magnetic_fluid;
   double Bmag_ref = BmagE * Sqr((GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid) / r_ref);
   GeoVector B0(Bmag_ref, 0.0, 0.0);
   container.Insert(B0);

// Maximum displacement
   double dmax = params[0] * one_au;
   container.Insert(dmax);

// Solar rotation vector
   GeoVector Omega(0.0, 0.0, 0.0);
   container.Insert(Omega);

// Reference equatorial distance
   container.Insert(r_ref);

// dmax fraction
   double dmax_fraction = params[2];
   container.Insert(dmax_fraction);

// Spherical shock radius
   container.Insert(R_sh); 

// Spherical shock width
   double w_sh = params[1] * one_au;
   container.Insert(w_sh);
   
// Spherical shock strength
   container.Insert(s);

   simulation->AddBackground(BackgroundSolarWindTermShock(), container);

//--------------------------------------------------------------------------------------------------
// Time initial condition
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Initial time
   container.Insert(0.0);

// Final time
   container.Insert(tf);

// Number of subintervals (0 for random points)
   int n_intervals = 0;
   container.Insert(n_intervals);

   simulation->AddInitial(InitialTimeInterval(), container);

//--------------------------------------------------------------------------------------------------
// Spatial initial condition
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Lower bound for position
   GeoVector init_pos(R_sh + 0.5 * w_sh, 0.0, 0.0);
   container.Insert(init_pos);

   simulation->AddInitial(InitialSpaceFixed(), container);

//--------------------------------------------------------------------------------------------------
// Momentum initial condition
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Injection momentum
   container.Insert(p0);

   simulation->AddInitial(InitialMomentumShell(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Spatial Inner boundary
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   int max_crossings_Sun = 1;
   container.Insert(max_crossings_Sun);

// Action Vector
   std::vector<int> actions_Sun;
   actions_Sun.push_back(-1);
   for (i = 0; i < Nt; i++) actions_Sun.push_back(-1);
   container.Insert(actions_Sun);

// Origin
   container.Insert(gv_zeros);

// Radius
   double inner_boundary = 1.0 * one_au;
   container.Insert(inner_boundary);

   simulation->AddBoundary(BoundarySphereAbsorb(), container);

//--------------------------------------------------------------------------------------------------
// Time check-points
//--------------------------------------------------------------------------------------------------

   for (j = 0; j < Nt-1; j++) {
      container.Clear();

// Maximum crossings
      int max_crossings_time = -1;
      container.Insert(max_crossings_time);

// Actions vector
      std::vector<int> actions_time;
// Time distro
      actions_time.push_back(-1);
// Position-momentum distros
      for (i = 0; i < Nt-1; i++) actions_time.push_back((i == j ? 0 : -1));
      actions_time.push_back(-1);
      container.Insert(actions_time);

// Temporal check-points
      container.Insert(t_arr[j]);

      simulation->AddBoundary(BoundaryTimePass(), container);
   };

//--------------------------------------------------------------------------------------------------
// Time limit
//--------------------------------------------------------------------------------------------------

   container.Clear();

// Maximum crossings
   int max_crossings_time = 1;
   container.Insert(max_crossings_time);

// Actions vector
   std::vector<int> actions_time_lim;
// Time distro
   actions_time_lim.push_back(0);
// Position-momentum distros
   for (i = 0; i < Nt-1; i++) actions_time_lim.push_back(-1);
   actions_time_lim.push_back(0);
   container.Insert(actions_time_lim);

// Final time at which to stop integrating
   container.Insert(tf);

   simulation->AddBoundary(BoundaryTimeExpire(), container);

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
   int val_time1 = 0;
   container.Insert(val_time1);

   simulation->AddDistribution(DistributionTimeUniform(), container);

//--------------------------------------------------------------------------------------------------
// Distribution 2 (position-momentum)
//--------------------------------------------------------------------------------------------------

   for (i = 0; i < Nt; i++) {
      container.Clear();

// Number of bins
      MultiIndex n_bins2(Nr, Np, 0);
      container.Insert(n_bins2);

// Smallest value
      GeoVector minval2(r0, p0, 0.0);
      container.Insert(minval2);

// Largest value
      GeoVector maxval2(rf, pf, 0.0);
      container.Insert(maxval2);

// Linear or logarithmic bins
      MultiIndex log_bins2(0, 1, 0);
      container.Insert(log_bins2);

// Add outlying events to the end bins
      MultiIndex bin_outside2(0, 0, 0);
      container.Insert(bin_outside2);

// Physical units of the distro
      double unit_distro2 = 1.0;
      container.Insert(unit_distro2);

// Physical units of the bin values
      GeoVector unit_val2 = {unit_length_fluid, unit_momentum_particle, 1.0};
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

// Which component of position to use
      int pos_coord2 = 0;
      container.Insert(pos_coord2);

// Which component of momentum to use
      int mom_coord2 = 0;
      container.Insert(mom_coord2);

      simulation->AddDistribution(DistributionPositionMomentumUniform(), container);
   };

//--------------------------------------------------------------------------------------------------
// Run the simulation
//--------------------------------------------------------------------------------------------------

   int n_traj;
   int batch_size;

   batch_size = n_traj = 1;
   if(argc > 1) n_traj = atoi(argv[1]);
   if(argc > 2) batch_size = atoi(argv[2]);

   std::string simulation_files_prefix = "dsa_results/dsa_forward_";
   simulation->DistroFileName(simulation_files_prefix);
   simulation->SetTasks(n_traj, batch_size);
   simulation->MainLoop();
   simulation->PrintDistro1D(0, 0, simulation_files_prefix + "time.dat", false);
   for (i = 0; i < Nt; i++) {
// Number density vs position
      simulation->PrintDistro1D(i + 1, 0, simulation_files_prefix + "pos_" + std::to_string(i) + ".dat", false);
// Spectrum vs momentum
      simulation->PrintDistro2D(i + 1, 0, 1, simulation_files_prefix + "mom_" + std::to_string(i) + ".dat", false);
   };

   double divK = (kappa_up - kappa_dn) / w_sh;
   if (MPI_Config::is_master) {
      std::cout << "divK = " << divK << std::endl;
      std::cout << "width = " << w_sh << std::endl;
      std::cout << "K / U (up) = " << kappa_up / U_up << std::endl;
      std::cout << "K / U (dn) = " << kappa_dn / U_dn << std::endl;
      std::cout << "K / U^2 = " << 0.25 * tau / one_day << std::endl;
      std::cout << "min dt_adv = " << w_sh / (U_up + divK) << std::endl;
      std::cout << "max dt_dif = " << Sqr(w_sh) / kappa_dn << std::endl;
   };

   return 0;
};


