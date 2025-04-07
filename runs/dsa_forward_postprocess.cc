#include "common/physics.hh"
#include "dsa_common.hh"
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace Spectrum;

int main(int argc, char** argv)
{
   int i, j, k;
   int n_traj = 1, sum_c[Nz][Np];
   if(argc > 1) n_traj = atoi(argv[1]);
   double coord1[Nz][Np], coord2[Nz][Np], distro[Nz][Np], sum_w[Nz][Np];
   std::string infilename;
   std::string outfilename;
   std::string line;
   DefineArrays();
   ReadParams();

   for (i = 0; i < Nt; i++) {
// Number density vs position
      infilename = "dsa_results/dsa_forward_pos_" + std::to_string(i) + ".dat";
      outfilename = "dsa_results/dsa_forward_pos_" + std::to_string(i) + "_pp.dat";

// Open input analytic distro file
      std::ifstream input_sda_file(infilename);

// Read first two lines of distro file
      std::getline(input_sda_file, line);
      std::getline(input_sda_file, line);

// Read data
      for(j = 0; j < Nz; j++) {
         input_sda_file >> coord1[j][0];
         input_sda_file >> distro[j][0];
         input_sda_file >> sum_w[j][0];
         input_sda_file >> sum_c[j][0];
      };

// Close input cartesian distro file
      input_sda_file.close();

// Open output distro file
      std::ofstream output_sda_file(outfilename);

// Output data
      output_sda_file << std::setprecision(8);
      for(j = 0; j < Nz; j++) {
         output_sda_file << std::setw(20) << coord1[j][0]
                         << std::setw(20) << sum_c[j][0] / dz * Q * tf / n_traj
                         << std::endl;
      };

// Close output distro file
      output_sda_file.close();
   };

   for (i = 0; i < Nt; i++) {
// Spectrum vs momentum
      infilename = "dsa_results/dsa_forward_mom_" + std::to_string(i) + ".dat";
      outfilename = "dsa_results/dsa_forward_mom_" + std::to_string(i) + "_pp.dat";

// Open input analytic distro file
      std::ifstream input_sda_file(infilename);

// Read first two lines of distro file
      std::getline(input_sda_file, line);
      std::getline(input_sda_file, line);

// Read data
      for(j = 0; j < Nz; j++) {
         for(k = 0; k < Np; k++) {
            input_sda_file >> coord1[j][k];
            input_sda_file >> coord2[j][k];
            input_sda_file >> distro[j][k];
            input_sda_file >> sum_w[j][k];
            input_sda_file >> sum_c[j][k];
         };
      };

// Close input cartesian distro file
      input_sda_file.close();

// Open output distro file
      std::ofstream output_sda_file(outfilename);

// Output data
      output_sda_file << std::setprecision(8);
      j = (params[3] * z_shock - z0) / dz;
      for(k = 0; k < Np; k++) {
         output_sda_file << std::setw(20) << coord2[j][k]
                         << std::setw(20) << sum_c[j][k] / (dz * dp_arr[k]) * Q * tf / n_traj
                         << std::endl;
      };

// Close output distro file
      output_sda_file.close();
   };

   std::cout << "Post-processed position distribution files outputed." << std::endl;
   return 0;
};
