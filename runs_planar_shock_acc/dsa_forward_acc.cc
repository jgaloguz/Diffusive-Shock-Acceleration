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
   return Sqr(dlnAdx(x)) * (1.0 - 1.0 / A0 / w_sh) + (A0 / w_sh);
};

// Function to accumulate path density
void ConditionalPathDensity(double **pd, double bin_time, std::vector<double> t_hist,
                            std::vector<double> x_hist, std::vector<double> p_hist)
{
   int i = 0, j, k;

// Bin path
   i = LocateInArray(0, t_hist.size()-1, t_hist.data(), bin_time, false);
   if (i >= 0) {
      j = LocateInArray(0, Nz, z_arr_edges, x_hist[i], false);
      k = LocateInArray(0, Np, p_arr_edges, p_hist[i], false);
      if (j >= 0 && k >= 0) pd[j][k] += 1.0;
   };
};

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

   int i, j, k, counter, idx;
   double t, dt, Kpara, lnw, dA;
   double **pd0, **pd1, **pd2, **pd3;
   double **pd0_out, **pd1_out, **pd2_out, **pd3_out;
   GeoVector x, p, divK;
   SpatialData spdata;
   std::string outfilename1;

   DataContainer container;
   BackgroundSmoothShock background;
   DiffusionFlowMomentumPowerLaw diffusion;
   RNG rng(time(NULL) + comm_rank);

   bool child = false;
   int i_splt = -1, p_level, p_level_new, n_splits = 0, n_splits_out = 0;
   std::vector<int> idx_splt;
   std::vector<double> t_hist;
   std::vector<double> x_hist;
   std::vector<double> p_hist;
   std::vector<double> lnw_splt;
   double p_thrs[n_thrs] = {0.0};

// Allocate memory
   pd0 = Create2D<double>(Nz, Np);
   pd1 = Create2D<double>(Nz, Np);
   pd2 = Create2D<double>(Nz, Np);
   pd3 = Create2D<double>(Nz, Np);
   pd0_out = Create2D<double>(Nz, Np);
   pd1_out = Create2D<double>(Nz, Np);
   pd2_out = Create2D<double>(Nz, Np);
   pd3_out = Create2D<double>(Nz, Np);
   for (j = 0; j < Nz; j++) {
      for (k = 0; k < Np; k++) {
         pd0[j][k] = 0.0;
         pd1[j][k] = 0.0;
         pd2[j][k] = 0.0;
         pd3[j][k] = 0.0;
         pd0_out[j][k] = 0.0;
         pd1_out[j][k] = 0.0;
         pd2_out[j][k] = 0.0;
         pd3_out[j][k] = 0.0;
      };
   };

   ReadParams();
   DefineArrays();

   double dlogp_splt = (logpf - logp0) / n_thrs;
   for (i = 0; i < n_thrs; i++) p_thrs[i] = pow(10.0, logp0 + (i + 0.5) * dlogp_splt);

   j = LocateInArray(0, Nz, z_arr_edges, z_spectrum, false);
   double z1 = z_arr_edges[j];
   double z2 = z_arr_edges[j+1];
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
// State variables
         idx = idx_splt.back();
         t = t_hist[idx];
         x[0] = x_hist[idx];
         p[0] = p_hist[idx];
         p_level = LocateInArray(0, n_thrs-1, p_thrs, p[0], true);
         lnw = lnw_splt.back();
         i_splt--;
         if (i_splt < 0) child = false;
// History arrays
         idx_splt.pop_back();
         t_hist.resize(idx+1);
         x_hist.resize(idx+1);
         p_hist.resize(idx+1);
         lnw_splt.pop_back();
      }
      else {
// State variables
         t = rng.GetUniform() * tf;
         x[0] = 0.0;
         p[0] = p0;
         p_level = -1;
         lnw = 0.0;
// History arrays
         t_hist.clear();
         x_hist.clear();
         p_hist.clear();
         idx = 0;
         t_hist.push_back(t);
         x_hist.push_back(x[0]);
         p_hist.push_back(p[0]);
      };

// Time loop
      while (t < t_arr[3]) {
         background.GetFields(t, x, p, spdata);

// Compute Kpara and grad(Kpara) and assemble diffusion tensor
         Kpara = diffusion.GetComponent(1, t, x, p, spdata);
         divK[0] = diffusion.GetDirectionalDerivative(0);

// Take step and update state variables
         dt = fmin(spdata.dmax / (spdata.Uvec + divK).Norm(), Sqr(spdata.dmax) / Kpara);
         dt = 0.5 * fmin(dt, 3.0 * 0.01 / fabs(spdata.divU()));
         idx++;
         t += dt;
         x[0] += (spdata.Uvec[0] + divK[0]) * dt + sqrt(2.0 * Kpara * dt) * rng.GetNormal();
         p[0] -= p[0] * spdata.divU() * dt / 3.0;
#if defined(ENABLE_IMPORTANCE)
         dA = dlnAdx(x[0]);
         x[0] -= 2.0 * Kpara * dA * dt;
         lnw -= ((spdata.Uvec[0] + divK[0]) * dA + Kpara * d2lnAdx2(x[0])) * dt;
#elif defined(ENABLE_SPLITTING)
// Check momentum splitting threshold crossing
         p_level_new = LocateInArray(0, n_thrs-1, p_thrs, p[0], true);
         if (p_level_new > p_level) {
            n_splits++;
            child = true;
            p_level = p_level_new;
            lnw += child_lnw;
            for (i = 1; i < n_chld; i++) {
               i_splt++;
               idx_splt.push_back(idx);
               lnw_splt.push_back(lnw);
               counter++;
            };
         };
#endif
// Update history arrays
         t_hist.push_back(t);
         x_hist.push_back(x[0]);
         p_hist.push_back(p[0]);
      };
      counter--;
      ConditionalPathDensity(pd0, t_arr[0], t_hist, x_hist, p_hist);
      ConditionalPathDensity(pd1, t_arr[1], t_hist, x_hist, p_hist);
      ConditionalPathDensity(pd2, t_arr[2], t_hist, x_hist, p_hist);
      ConditionalPathDensity(pd3, t_arr[3], t_hist, x_hist, p_hist);

// Bin particle
      if (z1 < x[0] && x[0] < z2) {
         k = LocateInArray(0, Np, p_arr_edges, p[0], false);
         if (k + 1) distro[k] += exp(lnw);
      };
   };

// Share results
   MPI_Reduce(pd0[0], pd0_out[0], Nz*Np, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(pd1[0], pd1_out[0], Nz*Np, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(pd2[0], pd2_out[0], Nz*Np, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(pd3[0], pd3_out[0], Nz*Np, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(distro, distro_out, Np, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&n_splits, &n_splits_out, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

// Output results
   if (comm_rank == 0) {
      std::cerr << std::endl << "Number of splits = " << n_splits_out << std::endl;
      outfilename1 = "dsa_results/dsa_forward_path_dens_pp";
#if defined(ENABLE_IMPORTANCE)
      outfilename1 += "_imps.dat";
      std::cerr << std::endl << "IMPORTANCE" << std::endl;
#elif defined(ENABLE_SPLITTING)
      outfilename1 += "_split.dat";
      std::cerr << std::endl << "SPLITTING" << std::endl;
#else
      outfilename1 += "_baseline.dat";
      std::cerr << std::endl << "BASELINE" << std::endl;
#endif
// Path densities
      std::ofstream output_dsa_file1(outfilename1);
      for(j = 0; j < Nz; j++) {
         for(k = 0; k < Np; k++) {
            output_dsa_file1 << std::setw(20) 
                             << pd0_out[j][k] * Q * tf / (dz * dp_arr[k]) / (n_traj * comm_size);
         };
         output_dsa_file1 << std::endl;
      };
      for(j = 0; j < Nz; j++) {
         for(k = 0; k < Np; k++) {
            output_dsa_file1 << std::setw(20) 
                             << pd1_out[j][k] * Q * tf / (dz * dp_arr[k]) / (n_traj * comm_size);
         };
         output_dsa_file1 << std::endl;
      };
      for(j = 0; j < Nz; j++) {
         for(k = 0; k < Np; k++) {
            output_dsa_file1 << std::setw(20) 
                             << pd2_out[j][k] * Q * tf / (dz * dp_arr[k]) / (n_traj * comm_size);
         };
         output_dsa_file1 << std::endl;
      };
      for(j = 0; j < Nz; j++) {
         for(k = 0; k < Np; k++) {
            output_dsa_file1 << std::setw(20) 
                             << pd3_out[j][k] * Q * tf / (dz * dp_arr[k]) / (n_traj * comm_size);
         };
         output_dsa_file1 << std::endl;
      };
      output_dsa_file1.close();
   };

// Free memory
   Delete2D(pd0);
   Delete2D(pd1);
   Delete2D(pd2);
   Delete2D(pd3);

// Finalize the MPI environment.
   MPI_Finalize();

   return 0;
};