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

#define CFL 0.5                                       // CFL factor
#define eps 1.0e-8                                    // tolerance for Jacobi method

// Function to check the (1-norm) difference between two matrices
double norm1(double **A, double **B)
{
   double D = 0.0;
   for (int i = 1; i < Nz+1; i++) {
      for (int j = 1; j < Np+1; j++) {
         D += fabs(A[i][j] - B[i][j]);
      };
   };
   return D;
};

int main(int argc, char** argv)
{
   int i, j, iter, time_index;
   double t, dt, dt_dx, dt_dx2, dt_dlnp, IC;
   double *u, *k, *du_dx, *dk_dx;
   double **h_new, **h_old, **RHS_old, RHS_new;
   GeoVector x, p;
   SpatialData spdata;
   std::string outfilename = "dsa_results/likelihood_expectation.dat";

   DataContainer container;
   BackgroundSmoothShock background;
   DiffusionFlowMomentumPowerLaw diffusion;

   ReadParams();
   DefineArrays();

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
// Iteratively solve linear system of equations
//--------------------------------------------------------------------------------------------------

// Allocate memory
   h_new = Create2D<double>(Nz+2, Np+2);
   h_old = Create2D<double>(Nz+2, Np+2);
   RHS_old = Create2D<double>(Nz+2, Np+2);
   u = new double[Nz];
   k = new double[Nz];
   du_dx = new double[Nz];
   dk_dx = new double[Nz];

   dt = CFL * Sqr(dz) / kappa_up;
   dt_dx = dt / dz / 4.0;
   dt_dx2 = dt / Sqr(dz) / 2.0;
   dt_dlnp = dt / dlogp / 6.0;

// Precompute values
   t = 0.0;
   x = gv_zeros;
   p = gv_zeros;
   for (i = 0; i < Nz; i++) {
      x[0] = z_arr[i];
      background.GetFields(t, x, p, spdata);
      u[i] = spdata.Uvec[0];
      k[i] = diffusion.GetComponent(1, t, x, p, spdata);
      du_dx[i] = spdata.gradUvec[0][0];
      dk_dx[i] = diffusion.GetDirectionalDerivative(0);
   };

// Initialize system: Domain
   for (j = 1; j < Np+1; j++) {
      IC = pow(n_chld, floor(n_thrs * (log10(p_arr[j-1]) - logp0) / (logpf - logp0)));
      for (i = 1; i < Nz+1; i++) {
         h_old[i][j] = IC;
         h_new[i][j] = IC;
      };
   };
// Initialize system: Edges
   for (i = 1; i < Nz+1; i++) {
      h_old[i][0] = 0.0;
      h_old[i][Np+1] = 0.0;
      h_new[i][0] = 0.0;
      h_new[i][Np+1] = 0.0;
   };
   for (j = 1; j < Np+1; j++) {
      h_old[0][j] = 0.0;
      h_old[Nz+1][j] = 0.0;
      h_new[0][j] = 0.0;
      h_new[Nz+1][j] = 0.0;
   };
// Initialize system: Corners
   h_old[0][0] = 0.0;
   h_old[0][Np+1] = 0.0;
   h_old[Nz+1][0] = 0.0;
   h_old[Nz+1][Np+1] = 0.0;
   h_new[0][0] = 0.0;
   h_new[0][Np+1] = 0.0;
   h_new[Nz+1][0] = 0.0;
   h_new[Nz+1][Np+1] = 0.0;

// Iterate through time
   time_index = 3;
   std::cout << "t_f = " << t_arr[time_index] << std::endl;
   while (t < t_arr[time_index]) {
      std::cout << "\rt = " << t;
      for (i = 1; i < Nz+1; i++) {
         for (j = 1; j < Np+1; j++) {
            RHS_old[i][j] = h_old[i][j]
                          - (u[i-1] - dk_dx[i-1]) * dt_dx * (h_old[i+1][j] - h_old[i-1][j])
                          + k[i-1] * dt_dx2 * (h_old[i+1][j] - 2.0 * h_old[i][j] + h_old[i-1][j])
                          + du_dx[i-1] * (h_old[i][j+1] - h_old[i][j-1]);
         };
      };
// Perform Jacobi iterations
      iter = 0;
      do {
         for (i = 1; i < Nz+1; i++) memcpy(h_old[i]+1, h_new[i]+1, Np * sizeof(double));
         for (i = 1; i < Nz+1; i++) {
            for (j = 1; j < Np+1; j++) {
               RHS_new = (u[i-1] - dk_dx[i-1]) * dt_dx * (h_old[i+1][j] - h_old[i-1][j])
                       - k[i-1] * dt_dx2 * (h_old[i+1][j] + h_old[i-1][j])
                       - du_dx[i-1] * (h_old[i][j+1] - h_old[i][j-1]);
               h_new[i][j] = (RHS_old[i][j] - RHS_new) / (1.0 - 2.0 * k[i-1] * dt_dx2);
            };
         };
         iter++;
      } while (norm1(h_old, h_new) > eps);
      for (i = 1; i < Nz+1; i++) memcpy(h_old[i]+1, h_new[i]+1, Np * sizeof(double));
      t += dt;
   };
   std::cout << "\rt = " << t_arr[time_index] << std::endl;

// Output to data
   std::ofstream output_sda_file(outfilename);
   for (i = 1; i < Nz+1; i++) {
      for (j = 1; j < Np+1; j++) {
         output_sda_file << std::setw(20) << h_new[i][j];
      };
      output_sda_file << std::endl;
   };
   output_sda_file.close();

// Clean-up
   Delete2D(h_new);
   Delete2D(h_old);
   Delete2D(RHS_old);
   delete[] u;
   delete[] k;
   delete[] du_dx;
   delete[] dk_dx;

   return 0;
};