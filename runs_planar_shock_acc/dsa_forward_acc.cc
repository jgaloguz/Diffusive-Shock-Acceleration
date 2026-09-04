#include "src/background_smooth_shock.hh"
#include "src/diffusion_other.hh"
#include "common/physics.hh"
#include "common/random.hh"
#include "common/spatial_data.hh"
#include "dsa_common.hh"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <chrono>

using namespace Spectrum;

// Accumulate path density
void CompCondPathDens(double **pd, double bin_time, std::vector<double> t_hist,
                      std::vector<double> x_hist, std::vector<double> p_hist)
{
   int i, j, k;
// Bin path
   i = LocateInArray(0, t_hist.size()-1, t_hist.data(), bin_time, false);
   if (i >= 0) {
      j = (x_hist[i] - z0) / dz;
      k = (log10(p_hist[i]) - logp0) / dlogp;
      if (k < 0) k = 0;
      if (0 <= j && j < Nz && k < Np) {
#if defined(LIKELIHOOD_TEST)
         pd[j][k] += pow(p_hist.back() / p0, alpha);
#else
         pd[j][k] += 1.0;
#endif
      };
   };
};

// Output path density
void OutputCondPathDens(std::ofstream& output_file, double **pd, double bin_time, int traj_total)
{
   int j, k;
   for(j = 0; j < Nz; j++) {
      for(k = 0; k < Np; k++) {
         output_file << std::setw(20) 
                     << pd[j][k] / (dz * dp_arr[k])
                      * Qtf / traj_total;
      };
      output_file << std::endl;
   };
};

int main(int argc, char** argv)
{
// Initialize the MPI environment
   MPI_Init(&argc, &argv);
   int comm_size, comm_rank;
   MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
   MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

   int n_traj = 1000;
   if (argc > 1) n_traj = atoi(argv[1]);
   int n_traj_10 = n_traj / 10;

// Header message
   if (comm_rank == 0) {
      std::cerr << "Number of CPUs: " << comm_size << std::endl;
      std::cerr << "Number of trajectories per CPU: " << n_traj << std::endl;
#if defined(ENABLE_IMP_SAMP) && (ENABLE_IMP_SAMP == 0)
      std::cerr << "Reading likelihood from file:"
                << "dsa_results/feynman_kac/dsa_likelihood_alpha=" + alpha_string + ".dat"
                << std::endl;
#endif
      std::cerr << std::endl;
   };

// Declare variables, arrays, objects, etc
   int i, j, k, counter, idx;
   long n_steps_traj;
   long n_traj_total;
   double t, dt, Kpara, lnw, bin_wgt;
   double dAdt, dAdx, d2Adx2, dAdp;
   double exec_time, exec_time_out;
   double var_time, var_time_out;
   double **pd0, **pd1, **pd2, **pd3;
   double **pd0_out, **pd1_out, **pd2_out, **pd3_out;
   GeoVector x, p, divK;
   SpatialData spdata;
   std::string outfilename1, outfilename2;

   DataContainer container;
   BackgroundSmoothShock background;
   DiffusionFlowMomentumPowerLaw diffusion;

   bool child = false;
   int i_splt = -1, p_level, p_level_new;
   int num_child;
   double n_splits, n_splits_out;
   double var_splits, var_splits_out;
   double n_steps, n_steps_out;
   double var_steps, var_steps_out;
   std::vector<int> idx_splt;
   std::vector<double> t_hist;
   std::vector<double> x_hist;
   std::vector<double> p_hist;
   std::vector<double> lnw_splt;

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

// Initialize parameters and random number generator
   ReadParams();
   DefineArrays();
   ImportanceSamplingSetup();
   RNG rng(time(NULL) + comm_rank);
   spdata._mask = BACKGROUND_U | BACKGROUND_B | BACKGROUND_gradU | BACKGROUND_gradB;

#if ENABLE_SPLITTING > 0
// Initialize splitting thresholds and child probability
   double p_thrs[n_thrs] = {0.0};
   double dlogp_splt = (logpf - logp0) / n_thrs;
   double p_chld = pow(pf / p0, alpha / n_thrs);
   for (i = 0; i < n_thrs; i++) p_thrs[i] = pow(10.0, logp0 + (i + 0.5) * dlogp_splt);
#endif

// Initialize distribution arrays
   j = LocateInArray(0, Nz, z_arr_edges, z_spectrum, false);
   double z1 = z_arr_edges[j];
   double z2 = z_arr_edges[j+1];
   double distro[Np] = {0.0};
   double distro_out[Np] = {0.0};
   double distro_ref;
   double distro_err[Np] = {0.0};
   double distro_err_out[Np] = {0.0};

// Compute phase-space integral of interest for e^2 metric
   double H1;
   double H[Np] = {0.0};
   double Hh[Np] = {0.0};
   double Hh_out[Np] = {0.0};
   double e2[Np] = {0.0};
   double e2_out[Np] = {0.0};
   for (i = 0; i < Np; i++) H[i] = SpaceMomentumIntegral(z1, z2, 10, p_arr_edges[i], p_arr_edges[i+1], 100, tf);
   for (i = Np-1; i > 0; i--) H[i-1] += H[i];

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

// Power of bulk velocity dependence
   double power_law_U = 2.0;
   container.Insert(power_law_U);

// Normalization of particle momentum
   double p_0 =  0.1 * mass[specie] * c_code;
   container.Insert(p_0);

// Power of particle momentum dependence
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
   n_splits = 0.0;
   n_steps = 0.0;
   std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now(); // Start timer
   while (counter > 0) {
      if (comm_rank == 0) {
         if (counter % n_traj_10 == 0) std::cerr << counter << std::endl;
      };
// Initialize particle
      n_steps_traj = 0;
      if (child) {
// State variables
         idx = idx_splt.back();
         t = t_hist[idx];
         x[0] = x_hist[idx];
         p[0] = p_hist[idx];
#if ENABLE_SPLITTING > 0
         p_level = LocateInArray(0, n_thrs-1, p_thrs, p[0], true);
#endif
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
         t = init_time(rng.GetUniform());
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
      while (t < tf) {
         background.GetFields(t, x, p, spdata);

// Compute Kpara and grad(Kpara) and assemble diffusion tensor
         Kpara = diffusion.GetComponent(1, t, x, p, spdata);
         divK[0] = diffusion.GetDirectionalDerivative(0);

// Take step and update state variables
         dt = fmin(Sqr(spdata.dmax) / Kpara, 3.0 * 0.01 / fabs(spdata.divU()));
#if defined(ENABLE_IMP_SAMP)
         dAdt = dlnAdt(t, x[0]);
         dAdx = dlnAdx(t, x[0]);
         d2Adx2 = d2lnAdx2(t, x[0]);
         dAdp = dlnAdp(p[0]);
         dt = 0.5 * fmin(dt, spdata.dmax / fabs(spdata.Uvec[0] + divK[0] + 2.0 * Kpara * dAdx));
         x[0] += 2.0 * Kpara * dAdx * dt;
         lnw += ( dAdt + (spdata.Uvec[0] + divK[0] + Kpara * dAdx) * dAdx
                + Kpara * d2Adx2 - (p[0] / 3.0) * spdata.divU() * dAdp ) * dt;
#else
         dt = 0.5 * fmin(dt, spdata.dmax / fabs(spdata.Uvec[0] + divK[0]));
#endif
         idx++;
         t += dt;
         x[0] += (spdata.Uvec[0] + divK[0]) * dt + sqrt(2.0 * Kpara * dt) * rng.GetNormal();
         p[0] -= p[0] * spdata.divU() * dt / 3.0;
#if defined(ENABLE_SPLITTING)
#if ENABLE_SPLITTING == 0
// Check momentum splitting threshold
         num_child = pow(1.0 - spdata.divU() * dt / 3.0, alpha) + rng.GetUniform();
#else
// Check momentum splitting threshold crossing
         p_level_new = LocateInArray(0, n_thrs-1, p_thrs, p[0], true);
         num_child = 1;
         for (i = p_level; i < p_level_new; i++) {
            num_child *= static_cast<int>(p_chld + rng.GetUniform());
         };
         p_level = p_level_new;
#endif
// Check if particle splits based number of children (computed from splitting probability)
         if (num_child > 1) {
            child = true;
            lnw += log(1.0 / num_child);
            for (i = 1; i < num_child; i++) {
               n_splits += 1.0;
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
         n_steps_traj += 1;
      };
      counter--;
// Bin conditional path densities
      CompCondPathDens(pd0, t_arr[0], t_hist, x_hist, p_hist);
      CompCondPathDens(pd1, t_arr[1], t_hist, x_hist, p_hist);
      CompCondPathDens(pd2, t_arr[2], t_hist, x_hist, p_hist);
      CompCondPathDens(pd3, t_arr[3], t_hist, x_hist, p_hist);

// Bin particle and compute e^2 metric
      k = (log10(p[0]) - logp0) / dlogp;
      if (z1 < x[0] && x[0] < z2 && k < Np) {
         distro[k] += exp(lnw);
         H1 = Qtf * exp(lnw) / A(tf, x[0], p[0]);
      }
      else H1 = 0.0;
      for (i = 0; i < fmin(k+1,Np); i++) {
         Hh[i] += H1;
         e2[i] += Sqr((H1 - H[i]) / H[i]);
      };
      for (i = k+1; i < Np; i++) {
         Hh[i] += 0.0;
         e2[i] += 1.0;
      };
// Tally total number of steps
      n_steps += n_steps_traj;
   };
   std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now(); // End timer
   exec_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

// Share results
   MPI_Reduce(pd0[0], pd0_out[0], Nz*Np, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(pd1[0], pd1_out[0], Nz*Np, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(pd2[0], pd2_out[0], Nz*Np, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(pd3[0], pd3_out[0], Nz*Np, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(distro, distro_out, Np, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

// Share execution time statistics
   MPI_Reduce(&exec_time, &exec_time_out, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   exec_time_out /= comm_size;
   MPI_Bcast(&exec_time_out, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   var_time = Sqr(exec_time - exec_time_out);
   MPI_Reduce(&var_time, &var_time_out, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   var_time_out /= comm_size;

#if defined(ENABLE_SPLITTING)
// Share splitting statistics
   MPI_Reduce(&n_splits, &n_splits_out, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   n_splits_out /= comm_size;
   MPI_Bcast(&n_splits_out, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   var_splits = Sqr(n_splits - n_splits_out);
   MPI_Reduce(&var_splits, &var_splits_out, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   var_splits_out /= comm_size;
#endif

// Share number of steps statistics
   MPI_Reduce(&n_steps, &n_steps_out, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   n_steps_out /= comm_size;
   MPI_Bcast(&n_steps_out, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   var_steps = Sqr(n_steps - n_steps_out);
   MPI_Reduce(&var_steps, &var_steps_out, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   var_steps_out /= comm_size;

// Share result variance
   for(k = 0; k < Np; k++) {
      distro[k] *= Qtf / A(tf, z_spectrum, p_arr[k]) / (dz * dp_arr[k]) / n_traj;
      distro_ref = N12(z_spectrum, p_arr[k], tf) * M_4PI * Sqr(p_arr[k]);
      distro_err[k] = Sqr((distro[k] - distro_ref) / distro_ref);
   };
   MPI_Reduce(distro_err, distro_err_out, Np, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   for(k = 0; k < Np; k++) distro_err_out[k] /= comm_size;

// Share e^2 statistic
   n_traj_total = n_traj * comm_size;
   MPI_Reduce(Hh, Hh_out, Np, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   for(k = 0; k < Np; k++) Hh_out[k] /= n_traj_total;
   MPI_Reduce(e2, e2_out, Np, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   for(k = 0; k < Np; k++) e2_out[k] /= n_traj_total;

// Output results
   if (comm_rank == 0) {
// Elapsed time
      std::cout << "Execution time per seed particle = " << exec_time_out / n_traj << " +- " << sqrt(var_time_out) / n_traj << " [ms]" << std::endl;

      outfilename1 = "dsa_results/dsa_forward_path_dens_pp";
      outfilename2 = "dsa_results/dsa_forward_mom_pp";
#if defined(ENABLE_SPLITTING)
      outfilename1 += "_split_" + std::to_string(ENABLE_SPLITTING * n_thrs) + ".dat";
      outfilename2 += "_split_" + std::to_string(ENABLE_SPLITTING * n_thrs) + ".dat";
      std::cerr << std::endl << "PART. SPLIT. " << ENABLE_SPLITTING << std::endl;
      std::cout << "Splits per seed particle = " << n_splits_out / n_traj << " +- " << sqrt(var_splits_out) / n_traj << std::endl;
#elif defined(ENABLE_IMP_SAMP)
      outfilename1 += "_imps_" + std::to_string(ENABLE_IMP_SAMP) + ".dat";
      outfilename2 += "_imps_" + std::to_string(ENABLE_IMP_SAMP) + ".dat";
      std::cerr << std::endl << "IMP. SAMP. " << ENABLE_IMP_SAMP << std::endl;
      std::cout << "Importance sampling factor @ x = 0: " << A(tf, 0.0, p0) << std::endl;
#elif defined(LIKELIHOOD_TEST)
      outfilename1 += "_ltest.dat";
      outfilename2 += "_ltest.dat";
      std::cerr << std::endl << "LIKE. TEST" << std::endl;
#else
      outfilename1 += "_base.dat";
      outfilename2 += "_base.dat";
      std::cerr << std::endl << "BASE METHOD" << std::endl;
#endif
// Number of steps
      std::cout << "Steps per seed particle = " << n_steps_out / n_traj << " +- " << sqrt(var_steps_out) / n_traj << std::endl;

// Path densities
      std::ofstream output_dsa_file1(outfilename1);
      OutputCondPathDens(output_dsa_file1, pd0_out, t_arr[0], n_traj_total);
      OutputCondPathDens(output_dsa_file1, pd1_out, t_arr[1], n_traj_total);
      OutputCondPathDens(output_dsa_file1, pd2_out, t_arr[2], n_traj_total);
      OutputCondPathDens(output_dsa_file1, pd3_out, t_arr[3], n_traj_total);
      output_dsa_file1.close();

// Spectrum near shock and related metrics
      std::ofstream output_dsa_file2(outfilename2);
      for(k = 0; k < Np; k++) {
         distro_out[k] *= Qtf / A(tf, z_spectrum, p_arr[k]) / (dz * dp_arr[k]) / n_traj_total;
         output_dsa_file2 << std::setw(20) << EnrKin(p_arr[k], specie) / one_MeV
                          << std::setw(20) << distro_out[k]
                          << std::setw(20) << N12(z_spectrum, p_arr[k], tf) * M_4PI * Sqr(p_arr[k])
                          << std::setw(20) << sqrt(distro_err_out[k])
                          << std::setw(20) << H[k]
                          << std::setw(20) << Hh_out[k]
                          << std::setw(20) << e2_out[k]
                          << std::endl;
      };
      output_dsa_file2.close();
   };

// Free memory
   Delete2D(pd0);
   Delete2D(pd1);
   Delete2D(pd2);
   Delete2D(pd3);
#if defined(ENABLE_IMP_SAMP)
   ImportanceSamplingFree();
#endif

// Finalize the MPI environment.
   MPI_Finalize();

   return 0;
};