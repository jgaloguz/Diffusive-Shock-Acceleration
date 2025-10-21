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

// Numerical solver parameters
#define Nx_lft 400
#define Nx_ctr 200
#define Nx_rgt 100
const int Nx_FD = Nx_lft + Nx_ctr + Nx_rgt + 1;
const double X1_FD =-10.0;
const double X2_FD = -0.01;
const double X3_FD = 0.01;
const double X4_FD = 1.0;
const double Tf_FD = 10.0 * one_day;
const double dlogx_lft_FD = (log(-X2_FD) - log(-X1_FD)) / Nx_lft;
const double dx_ctr = (X3_FD - X2_FD) / Nx_ctr;
const double dlogx_rgt_FD = (log(X4_FD) - log(X3_FD)) / Nx_rgt;
const double dt_FD = 1.0e-8;
const double dt_out_FD = Tf_FD / 100.0;

// Lax-Wendroff
void LaxWendroff(double *fnew, double *fold, double *df, double *Kdf,
                 double *cf_lft, double *cf_ctr, double *cf_rgt,
                 double *u, double *kap, double *du, double alpha_3)
{
   int i;
// Pre-compute first derivative of f
   for (i = 1; i < Nx_FD-1; i++) {
      df[i] = cf_lft[i-1] * fold[i-1] + cf_ctr[i-1] * fold[i] + cf_rgt[i-1] * fold[i+1];
      Kdf[i] = kap[i] * df[i];
   };
   for (i = 2; i < Nx_FD-2; i++) {
// Advection
      fnew[i] += dt_FD * (u[i] * df[i]);
// Diffusion
      fnew[i] += dt_FD * (cf_lft[i-1] * Kdf[i-1] + cf_ctr[i-1] * Kdf[i] + cf_rgt[i-1] * Kdf[i+1]);
// Source
      fnew[i] -= dt_FD * (du[i] * fold[i] * alpha_3);
   };
};


int main(int argc, char** argv)
{
   int i, j;
   int Nt_out;
   double t, t_out;
   double *X, *dX;
   double *Cf_lft, *Cf_ctr, *Cf_rgt;
   double *U, *K, *dU;
   double *fold, *fnew;
   double *df, *Kdf;
   GeoVector pos = gv_zeros, mom = gv_zeros;
   SpatialData spdata;
   std::ofstream params, coeffs, solution;

   DataContainer container;
   BackgroundSmoothShock background;
   DiffusionFlowMomentumPowerLaw diffusion;

   ReadParams();
   DefineArrays();
   double alpha_3 = n_thrs * log(n_chld) / log(pf / p0) / 3.0;

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
// 1D (space) advection-diffusion-source equation solver with Lax-Wendroff method
//--------------------------------------------------------------------------------------------------

// Allocate memory
   X = new double[Nx_FD];
   dX = new double[Nx_FD-1];
   Cf_lft = new double[Nx_FD-2];
   Cf_ctr = new double[Nx_FD-2];
   Cf_rgt = new double[Nx_FD-2];
   U = new double[Nx_FD];
   K = new double[Nx_FD];
   dU = new double[Nx_FD];
   fold = new double[Nx_FD];
   fnew = new double[Nx_FD];
   df = new double[Nx_FD];
   Kdf = new double[Nx_FD];

// Define grid, differences, and reconstruction coefficients
   for (i = 0; i < Nx_lft; i++) X[i] = -exp(log(-X2_FD) - (Nx_lft - i) * dlogx_lft_FD);
   for (i = 0; i < Nx_ctr; i++) X[i+Nx_lft] = X2_FD + i * dx_ctr;
   for (i = 0; i <= Nx_rgt; i++) X[i+Nx_lft+Nx_ctr] = exp(log(X3_FD) + i * dlogx_rgt_FD);
   for (i = 0; i < Nx_FD-1; i++) dX[i] = X[i+1] - X[i];
   for (i = 0; i < Nx_FD-2; i++) {
      Cf_lft[i] = 1.0/(dX[i+1]+dX[i]) - 1.0/dX[i];
      Cf_ctr[i] = 1.0/dX[i] - 1.0/dX[i+1];
      Cf_rgt[i] = 1.0/dX[i+1] - 1.0/(dX[i+1]+dX[i]);
   };

// Set coefficients (and print to file) and initial condition
   coeffs.open("dsa_results/likelihood_coeffs_FD.dat");
   coeffs << std::setprecision(8);
   for (i = 0; i < Nx_FD; i++) {
      pos[0] = X[i];
      background.GetFields(t, pos, mom, spdata);
      U[i] = spdata.Uvec[0];
      K[i] = diffusion.GetComponent(1, t, pos, mom, spdata);
      dU[i] = spdata.gradUvec[0][0];
      coeffs << std::setw(16) << pos[0]
             << std::setw(16) << U[i]
             << std::setw(16) << K[i]
             << std::setw(16) << dU[i]
             << std::endl;
      U[i] += diffusion.GetDirectionalDerivative(0);
      fnew[i] = 1.0;
      df[i] = 0.0;
      Kdf[i] = 0.0;
   };
   coeffs.close();

// Iterate
   solution.open("dsa_results/likelihood_solution_FD.dat");
   solution << std::setprecision(8);
   t = Tf_FD;
   Nt_out = 0;
   t_out = Tf_FD;
   while (t > 0.0) {
// Print solution
      if (t <= t_out) {
         for (i = 0; i < Nx_FD; i++) solution << std::setw(16) << fnew[i];
         solution << std::endl;
         Nt_out++;
         t_out -= dt_out_FD;
      };
// Take a step
      memcpy(fold, fnew, Nx_FD * sizeof(double));
      LaxWendroff(fnew, fold, df, Kdf,
                  Cf_lft, Cf_ctr, Cf_rgt,
                  U, K, dU, alpha_3);
// Increase time
      t -= dt_FD;
      std::cerr << "\rt = " << t << std::endl;
   };
   for (i = 0; i < Nx_FD; i++) solution << std::setw(16) << fnew[i];
   solution << std::endl;
   Nt_out++;
   std::cerr << std::endl;
   solution.close();

// Print parameters to file
   params.open("dsa_results/likelihood_params_FD.dat");
   params << std::setprecision(8);
   params << std::setw(16) << Nx_FD
          << std::setw(16) << Nx_lft
          << std::setw(16) << Nx_ctr
          << std::setw(16) << Nx_rgt
          << std::setw(16) << X1_FD
          << std::setw(16) << X2_FD
          << std::setw(16) << X3_FD
          << std::setw(16) << X4_FD
          << std::setw(16) << dt_FD
          << std::setw(16) << Tf_FD
          << std::setw(16) << Nt_out
          << std::setw(16) << 3.0 * alpha_3
          << std::endl;
   params.close();

// De-allocate memory
   delete[] X;
   delete[] dX;
   delete[] Cf_lft;
   delete[] Cf_ctr;
   delete[] Cf_rgt;
   delete[] U;
   delete[] K;
   delete[] dU;
   delete[] fold;
   delete[] fnew;
   delete[] df;
   delete[] Kdf;

   return 0;
};