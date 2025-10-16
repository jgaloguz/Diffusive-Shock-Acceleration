#include "src/background_smooth_shock.hh"
#include "src/diffusion_other.hh"
#include "common/physics.hh"
#include "common/random.hh"
#include "common/spatial_data.hh"
#include "dsa_common.hh"
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace Spectrum;

int main(int argc, char** argv)
{
// Initialize the MPI environment
   MPI_Init(&argc, &argv);
   int comm_size, comm_rank;
   MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
   MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

// Header message
   if (comm_rank == 0) {
      std::cerr << "Number of CPUs: " << comm_size << std::endl;
      std::cerr << "Number of trajectories per CPU: " << n_traj << std::endl;
      std::cerr << std::endl;
   };

   int i, j, k, counter, n_traj_10;
   double t, dt, Kpara, lnw, alpha_3;
   double distro, distro_out;
   GeoVector x, p, divK;
   SpatialData spdata;
   std::string outfilename;

   DataContainer container;
   BackgroundSmoothShock background;
   DiffusionFlowMomentumPowerLaw diffusion;
   RNG rng(time(NULL) + comm_rank);

   ReadParams();
   DefineArrays();

   alpha_3 = n_thrs * log(n_chld) / log(pf / p0) / 3.0;
   spdata._mask = BACKGROUND_U | BACKGROUND_B | BACKGROUND_gradU | BACKGROUND_gradB;

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
   container.Insert(w_sh);

// dmax fraction
   container.Insert(dmax_fraction);

   background.SetupObject(container);

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

   diffusion.SetupObject(container);

//--------------------------------------------------------------------------------------------------
// Trajectory Loop
//--------------------------------------------------------------------------------------------------

// Load command line arguments
   double x_init = 0.0;
   if (argc > 1) x_init = atof(argv[1]);
   double t_final = t_arr[0];
   if (argc > 2) t_final = t_arr[atoi(argv[2])];

// Initialize variables
   x = gv_zeros;
   p = gv_zeros;
   p[0] = p0;
   divK = gv_zeros;
   distro = 0.0;
   counter = n_traj;
   n_traj_10 = n_traj / 10;
   while (counter > 0) {
      if (comm_rank == 0) {
         if (counter % n_traj_10 == 0) std::cerr << counter << std::endl;
      };
// Initialize particle
      t = 0.0;
      x[0] = x_init;
      lnw = 0.0;

// Time loop
      while (t < t_final) {
         background.GetFields(t, x, p, spdata);

// Compute Kpara and grad(Kpara) and assemble diffusion tensor
         Kpara = diffusion.GetComponent(1, t, x, p, spdata);
         divK[0] = diffusion.GetDirectionalDerivative(0);

// Take step
         dt = fmin(spdata.dmax / (spdata.Uvec + divK).Norm(), Sqr(spdata.dmax) / Kpara);
         t += dt;
         x[0] += (spdata.Uvec[0] + divK[0]) * dt + sqrt(2.0 * Kpara * dt) * rng.GetNormal();
         lnw += alpha_3 * spdata.divU() * dt;
      };
      counter--;

// Bin particle
      distro += exp(-lnw);
   };
   distro /= n_traj;

// Share results
   MPI_Reduce(&distro, &distro_out, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

// Output results
   x[0] = x_init;
   background.GetFields(0.0, x, p, spdata);
   Kpara = diffusion.GetComponent(1, 0.0, x, p, spdata);
   std::cout << std::setprecision(6);
   if (comm_rank == 0) {
      std::cout << std::setw(18) << x_init
                << std::setw(18) << distro_out / comm_size
                << std::setw(18) << Kpara
                << std::endl;
   };

// Finalize the MPI environment.
   MPI_Finalize();

   return 0;
};