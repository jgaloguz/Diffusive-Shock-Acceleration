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

// Phase-space density upstream AND downstream
inline double N12(double z, double p, double t)
{
   if (z < 0.0) return N1(z, p, t);
   else return N2(z, p, t);
};

// Integral in momentum
inline double MomentumIntegral(double z, double t)
{
   int i;
   double S = 0.0;
   for (i = 0; i < Np; i++) S += N12(z, p_arr[i], t) * Sqr(p_arr[i]) * dp_arr[i];
   return S * M_4PI;
};

int main(int argc, char** argv)
{
   int i,j;
   std::ofstream dsa_analytic_file;

// Define initialize momentum, position, and time arrays
   ReadParams();
   DefineArrays();
   std::cout << "z_diff = " << kappa_up / U_up / one_au << " au" << std::endl;
   std::cout << "z_spectrum = " << z_spectrum / one_au << " au" << std::endl;
   std::cout << "w_shock = " << w_sh / one_au << " au" << std::endl;
   std::cout << "t_acc = " << tau / one_day << " days" << std::endl;
   std::cout << "t_final = " << t_arr[Nt-1] / one_day << " days" << std::endl;
   std::cout << "min dt_adv = " << w_sh / (U_up + (kappa_up - kappa_dn) / w_sh) / one_day << " days" << std::endl;
   std::cout << "max dt_dif = " << Sqr(w_sh) / kappa_dn << " days" << std::endl;

// Loop over times
   for (i = 0; i < Nt; i++) {
// Number density vs position (integrate over momentum)
      dsa_analytic_file.open("dsa_results/dsa_analytic_pos_" + std::to_string(i) + ".dat");
// Loop over spatial positions
      for (j = 0; j < Nz; j++) {
// Output to file
         dsa_analytic_file << std::setw(16) << z_arr[j] / one_au
                           << std::setw(16) << MomentumIntegral(z_arr[j], t_arr[i])
                           << std::endl;
      };
      dsa_analytic_file.close();
// Spectrum vs momentum near shock
      dsa_analytic_file.open("dsa_results/dsa_analytic_mom_" + std::to_string(i) + ".dat");
// Loop over momentum
      for (j = 0; j < Np; j++) {
// Output to file
         dsa_analytic_file << std::setw(16) << EnrKin(p_arr[j], specie) / one_MeV
                           << std::setw(16) << N12(z_spectrum, p_arr[j], t_arr[i]) * M_4PI * Sqr(p_arr[j])
                           << std::endl;
      };
      dsa_analytic_file.close();
   };

// Output times
   dsa_analytic_file.open("dsa_results/dsa_analytic_time.dat");
   for (i = 0; i < Nt; i++) dsa_analytic_file << std::setw(16) << t_arr[i] / one_day << std::endl;
   dsa_analytic_file.close();

   return 0;
};
