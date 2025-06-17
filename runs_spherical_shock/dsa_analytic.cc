#include "common/physics.hh"
#include "dsa_common.hh"
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace Spectrum;

// Phase-space density upstream
inline double N1(double z, double p, double t)
{
   double a = pow(p0 / p, beta);
   double b = 2.0 * z / (tau * U_up);
   double c = sqrt(t / tau);
   double d = 0.5 * beta / c * log(p / p0) - z / (sqrt(t * tau) * U_up);
   return amp * sqrt(Cube(p0 / p)) * exp(U_up * z / (2.0 * kappa_up)) * (exp(b) * erfc(d - c) * a + exp(-b) * erfc(d + c) / a);
};

// Phase-space density downstream
inline double N2(double z, double p, double t)
{
   double a = pow(p0 / p, beta);
   double b = 2.0 * z / (tau * U_dn);
   double c = sqrt(t / tau);
   double d = 0.5 * beta / c * log(p / p0) + z / (sqrt(t * tau) * U_dn);
   return amp * sqrt(Cube(p0 / p)) * exp(U_dn * z / (2.0 * kappa_dn)) * (exp(-b) * erfc(d - c) * a + exp(b) * erfc(d + c) / a);
};

// Integral in momentum
double MomentumIntegral(double z, double t)
{
   int i;
   double S = 0.0;
   if (z < 0.0) {
      for (i = 0; i < Np; i++) S += N1(z, p_arr[i], t) * Sqr(p_arr[i]) * dp_arr[i];
   }
   else {
      for (i = 0; i < Np; i++) S += N2(z, p_arr[i], t) * Sqr(p_arr[i]) * dp_arr[i];
   };
   return S * M_4PI;
};

int main(int argc, char** argv)
{
   int i,j;
   std::ofstream dsa_analytic_file;

// Define initialize momentum, position, and time arrays
   DefineArrays();
   ReadParams();
   std::cout << "z_shock = " << z_shock << std::endl;
   std::cout << "t_acc = " << tau << std::endl;

// Loop over times
   for (i = 0; i < Nt; i++) {
// Number density vs position (integrate over momentum)
      dsa_analytic_file.open("dsa_results/dsa_analytic_pos_" + std::to_string(i) + ".dat");
// Loop over spatial positions
      for (j = 0; j < Nz; j++) {
// Output to file
         dsa_analytic_file << std::setw(16) << z_arr[j]
                           << std::setw(16) << MomentumIntegral(z_arr[j], t_arr[i])
                           << std::endl;
      };
      dsa_analytic_file.close();
// Spectrum vs momentum near shock
      dsa_analytic_file.open("dsa_results/dsa_analytic_mom_" + std::to_string(i) + ".dat");
// Loop over momentum
      for (j = 0; j < Np; j++) {
// Output to file
         dsa_analytic_file << std::setw(16) << p_arr[j]
                           << std::setw(16) << N2(params[3] * z_shock, p_arr[j], t_arr[i]) * M_4PI * Sqr(p_arr[j])
                           << std::endl;
      };
      dsa_analytic_file.close();
   };

// Output times
   dsa_analytic_file.open("dsa_results/dsa_analytic_time.dat");
   for (i = 0; i < Nt; i++) dsa_analytic_file << std::setw(16) << t_arr[i] << std::endl;
   dsa_analytic_file.close();

   return 0;
};
