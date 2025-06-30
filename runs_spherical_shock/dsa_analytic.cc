#include "common/physics.hh"
#include "dsa_common.hh"
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace Spectrum;

// Phase-space density upstream AND downstream
double N12(double r, double p)
{
   double p1, p2, r1, r2, arg_p, arg_r, cef1, cef2, f;

   p1 = pow(p / p_inj, 3.0 * (Sqr(xi) - Sqr(zeta1)) / (2.0 * eta_up));
   p2 = pow(p / p_inj, 3.0 * (Sqr(xi) - Sqr(zeta2)) / (2.0 * eta_up));

   if (r < R_sh) {
      r1 = pow(r / R_sh, xi - zeta1);
      r2 = pow(r / R_sh, xi - zeta2);
      arg_r = sqrt(eta_up / (6.0 * log(p_inj / p))) * log(R_sh / r);
   }
   else {
      r1 = (1.0 - exp(-eta_dn * R_sh / r)) / (1.0 - exp(-eta_dn));
      r2 = r1;
      arg_r = 0.0;
   };
   arg_p = sqrt(3.0 / (2.0 * eta_up) * log(p_inj / p));

   if (p < p_inj) {
      cef1 = 1.0 - erf(zeta1 * arg_p + arg_r);
      cef2 = 1.0 + erf(zeta2 * arg_p + arg_r);
      f = -(zeta1 * r1 * p1 * cef1 + zeta2 * r2 * p2 * cef2);
   }
   else {
      f = -2.0 * zeta2 * r2 * p2;
   };

   return amp * f;
};

// Integral in momentum
double MomentumIntegral(double r)
{
   int i;
   double S = 0.0;
   for (i = 0; i < Np; i++) S += N12(r, p_arr[i]) * Sqr(p_arr[i]) * dp_arr[i];
   return S * M_4PI;
};

int main(int argc, char** argv)
{
   int i = Nt-1, j;
   std::ofstream dsa_analytic_file;

// Define initialize momentum, position, and time arrays
   ReadParams();
   DefineArrays();
   std::cout << "r_diff = " << (kappa_up / U_up) / one_au << " au" << std::endl;
   std::cout << "r_spectrum = " << r_spectrum / one_au << " au" << std::endl;
   std::cout << "t_acc = " << tau / one_day << " days" << std::endl;
   std::cout << "t_final = " << t_arr[Nt-1] / one_day << " days" << std::endl;

// Number density vs position (integrate over momentum)
   dsa_analytic_file.open("dsa_results/dsa_analytic_pos_" + std::to_string(i) + ".dat");
// Loop over spatial positions
   for (j = 0; j < Nr; j++) {
// Output to file
      dsa_analytic_file << std::setw(16) << r_arr[j] / one_au
                        << std::setw(16) << MomentumIntegral(r_arr[j])
                        << std::endl;
   };
   dsa_analytic_file.close();
// Spectrum vs momentum near shock
   dsa_analytic_file.open("dsa_results/dsa_analytic_mom_" + std::to_string(i) + ".dat");
// Loop over momentum
   for (j = 0; j < Np; j++) {
// Output to file
      dsa_analytic_file << std::setw(16) << EnrKin(p_arr[j], specie) / one_MeV
                        << std::setw(16) << N12(r_spectrum, p_arr[j]) * M_4PI * Sqr(p_arr[j])
                        << std::endl;
   };
   dsa_analytic_file.close();

// Output times
   dsa_analytic_file.open("dsa_results/dsa_analytic_time.dat");
   for (j = 0; j < Nt; j++) dsa_analytic_file << std::setw(16) << t_arr[j] / one_day << std::endl;
   dsa_analytic_file.close();

   return 0;
};
