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

constexpr int n_traj_10 = n_traj / 10;                // Number of trajectories divided by 10
constexpr double child_lnw = log(1.0 / n_chld);       // Log of inverse of childs per split

// Importance sampling auxiliary function and its derivatives
inline double A(double x)
{
   return pow(cosh(x / w_sh), A0 * w_sh);
};
inline double dlnAdx(double x)
{
   return A0 * tanh(x / w_sh);
};
inline double d2lnAdx2(double x)
{
   return Sqr(dlnAdx(x,w_sh)) * (1.0 - 1.0 / A0 / w_sh) + (A0 / w_sh);
};

int main(int argc, char** argv)
{
// Initialize the MPI environment
   MPI_Init(&argc, &argv);
   int comm_size, comm_rank;
   MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
   MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

   int i, j, k, counter;
   double t, dt, Kpara, lnw, dA;
   GeoVector x, p, divK;
   SpatialData spdata;
   std::string outfilename;

   DataContainer container;
   BackgroundSmoothShock background;
   DiffusionFlowMomentumPowerLaw diffusion;
   RNG rng(time(NULL) + comm_rank);

   bool child = false;
   int i_splt = -1, p_level, p_level_new, n_splits = 0, n_splits_out = 0;
   double t_splt[max_splt] = {0.0};
   double x_splt[max_splt] = {0.0};
   double p_splt[max_splt] = {0.0};
   double lnw_splt[max_splt] = {0.0};
   double p_thrs[n_thrs] = {0.0};

   ReadParams();
   DefineArrays();

   double dlogp_splt = (logpf - logp0) / n_thrs;
   for (i = 0; i < n_thrs; i++) p_thrs[i] = pow(10.0, logp0 + (i + 0.5) * dlogp_splt);

   j = LocateInArray(0, Nz-1, z_arr, z_spectrum, false);
   double z1 = z_arr[j];
   double z2 = z_arr[j+1];
   double distro[Np] = {0.0};
   double distro_out[Np] = {0.0};
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

   x = gv_zeros;
   p = gv_zeros;
   divK = gv_zeros;
   counter = n_traj;
   while (counter > 0) {
      if (comm_rank == 0) {
         if (counter % n_traj_10 == 0) std::cerr << counter << std::endl;
      };
// Initialize particle
      if (child) {
         t = t_splt[i_splt];
         x[0] = x_splt[i_splt];
         p[0] = p_splt[i_splt];
         p_level = LocateInArray(0, n_thrs-1, p_thrs, p[0], true);
         lnw = lnw_splt[i_splt];
         i_splt--;
         if (i_splt < 0) child = false;
      }
      else {
         t = rng.GetUniform() * tf;
         x[0] = 0.0;
         p[0] = p0;
         p_level = -1;
         lnw = 0.0;
      };

// Time loop
      while (t < tf) {
         background.GetFields(t, x, p, spdata);

// Compute Kpara and grad(Kpara) and assemble diffusion tensor
         Kpara = diffusion.GetComponent(1, t, x, p, spdata);
         divK[0] = diffusion.GetDirectionalDerivative(0);

// Take step
         dt = fmin(spdata.dmax / (spdata.Uvec + divK).Norm(), Sqr(spdata.dmax) / Kpara);
         dt = 0.5 * fmin(dt, 3.0 * 0.01 / fabs(spdata.divU()));
         t += dt;
         x[0] += (spdata.Uvec[0] + divK[0]) * dt + sqrt(2.0 * Kpara * dt) * rng.GetNormal();
         p[0] -= p[0] * spdata.divU() * dt / 3.0;
#ifdef ENABLE_IMPORTANCE
         dA = dlnAdx(x[0]);
         x[0] -= 2.0 * Kpara * dA * dt;
         lnw -= ((spdata.Uvec[0] + divK[0]) * dA + Kpara * d2lnAdx2(x[0])) * dt;
#endif

#ifdef ENABLE_SPLITTING
// Check momentum splitting threshold crossing
         p_level_new = LocateInArray(0, n_thrs-1, p_thrs, p[0], true);
         if (p_level_new > p_level) {
            n_splits++;
            child = true;
            p_level = p_level_new;
            lnw += child_lnw;
            for (i = 1; i < n_chld; i++) {
               i_splt++;
               t_splt[i_splt] = t;
               x_splt[i_splt] = x[0];
               p_splt[i_splt] = p[0];
               lnw_splt[i_splt] = lnw;
               counter++;
            };
         };
#endif
      };
      counter--;

// Bin particle
      if (z1 < x[0] && x[0] < z2) {
         k = LocateInArray(0, Np-1, p_arr, p[0], false);
         if (k + 1) distro[k] += exp(lnw);
      };
   };

// Share results
   MPI_Reduce(distro, distro_out, 100, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&n_splits, &n_splits_out, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

// Output results
   if (comm_rank == 0) {
      std::cerr << std::endl << "Number of splits = " << n_splits_out << std::endl;
      outfilename = "dsa_results/dsa_forward_mom_" + std::to_string(Nt-1) + "_pp";
#ifdef ENABLE_SPLITTING
      outfilename += "_split.dat";
      std::cerr << std::endl << "SPLITTING" << std::endl;
#else
#ifdef ENABLE_IMPORTANCE
      outfilename += "_imps.dat";
      std::cerr << std::endl << "IMPORTANCE" << std::endl;
#else
      outfilename += "_baseline.dat";
      std::cerr << std::endl << "BASELINE" << std::endl;
#endif
#endif
      std::ofstream output_sda_file(outfilename);
      for(k = 0; k < Np; k++) {
         output_sda_file << std::setw(20) << EnrKin(p_arr[k], specie) / one_MeV
                         << std::setw(20) << distro_out[k] / (dz * dp_arr[k]) * Q * tf / (n_traj * comm_size)
                         << std::endl;
      };
      output_sda_file.close();
   };

// Finalize the MPI environment.
   MPI_Finalize();

   return 0;
};