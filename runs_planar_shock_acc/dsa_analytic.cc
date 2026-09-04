#include "common/physics.hh"
#include "dsa_common.hh"
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace Spectrum;

int main(int argc, char** argv)
{
   int i, j, k;
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

   std::cout << std::endl;
   std::cout << "========================================" << std::endl;
   std::cout << "Dimensionless parameters:" << std::endl;
   std::cout << "u1 = " << U_up << std::endl;
   std::cout << "u2 = " << U_dn << std::endl;
   std::cout << "kappa1 = " << kappa_up << std::endl;
   std::cout << "Lsh = " << w_sh << std::endl;
   std::cout << "p0 = " << p0 << std::endl;
   std::cout << "alpha = " << alpha << std::endl;
   std::cout << "alpha-r_max slopes:" << std::endl;
   std::cout << "    (n_s = 1.25) " << log10(1.25) / log10(pf / p0) << std::endl;
   std::cout << "    (n_s = 1.50) " << log10(1.5) / log10(pf / p0) << std::endl;
   std::cout << "    (n_s = 1.75) " << log10(1.75) / log10(pf / p0) << std::endl;
   std::cout << "    (n_s = 2.00) " << log10(2.0) / log10(pf / p0) << std::endl;
   std::cout << "t1 = " << t_arr[0] << std::endl;
   std::cout << "t2 = " << t_arr[1] << std::endl;
   std::cout << "t3 = " << t_arr[2] << std::endl;
   std::cout << "t4 = " << t_arr[3] << std::endl;
   std::cout << "========================================" << std::endl;

// Loop over times
   dsa_analytic_file.open("dsa_results/dsa_forward_path_dens_pp_analytic.dat");
   for (i = 0; i < Nt; i++) {
// Loop over 2D grid
      for(j = 0; j < Nz; j++) {
         for(k = 0; k < Np; k++) {
            dsa_analytic_file << std::setw(16) << N12(z_arr[j], p_arr[k], t_arr[i]) * M_4PI * Sqr(p_arr[k]);
         };
         dsa_analytic_file << std::endl;
      };
   };
   dsa_analytic_file.close();

// Output times
   dsa_analytic_file.open("dsa_results/dsa_analytic_time.dat");
   for (i = 0; i < Nt; i++) dsa_analytic_file << std::setw(16) << t_arr[i] / one_day << std::endl;
   dsa_analytic_file.close();
// Output spatial grid
   dsa_analytic_file.open("dsa_results/dsa_analytic_pos.dat");
   for (j = 0; j < Nz; j++) dsa_analytic_file << std::setw(16) << z_arr[j] / one_au << std::endl;
   dsa_analytic_file.close();
// Output momentum grid
   dsa_analytic_file.open("dsa_results/dsa_analytic_mom.dat");
   for (k = 0; k < Np; k++) dsa_analytic_file << std::setw(16) << p_arr[k] / p0 << std::endl;
   dsa_analytic_file.close();
   dsa_analytic_file.open("dsa_results/dsa_analytic_dmom.dat");
   for (k = 0; k < Np; k++) dsa_analytic_file << std::setw(16) << dp_arr[k] / p0 << std::endl;
   dsa_analytic_file.close();
// Output kinetic energy grid
   dsa_analytic_file.open("dsa_results/dsa_analytic_enr.dat");
   for (k = 0; k < Np; k++) dsa_analytic_file << std::setw(16) << EnrKin(p_arr[k], specie) / one_MeV << std::endl;
   dsa_analytic_file.close();

   return 0;
};
