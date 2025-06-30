#include "common/physics.hh"
#include "dsa_common.hh"
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace Spectrum;

int main(int argc, char** argv)
{
   int i, j, k; 
   std::string infilename;
   std::string outfilename;
   std::string line;
   int sum_c[Np];
   double coord[Np], distro[Np], sum_w[Np];
   double S;
   ReadParams();
   DefineArrays();

   for (i = 0; i < Nt; i++) {
// Spectrum vs momentum
      infilename = "dsa_results/dsa_backward_mom_" + std::to_string(i) + ".dat";
      outfilename = "dsa_results/dsa_backward_mom_" + std::to_string(i) + "_pp.dat";

// Open input analytic distro file
      std::ifstream input_sda_file(infilename);

// Read first two lines of distro file
      std::getline(input_sda_file, line);
      std::getline(input_sda_file, line);

// Read data
      for(j = 0; j < Np; j++) {
         input_sda_file >> coord[j];
         input_sda_file >> distro[j];
         input_sda_file >> sum_w[j];
         input_sda_file >> sum_c[j];
      };

// Close input cartesian distro file
      input_sda_file.close();

// Open output distro file
      std::ofstream output_sda_file(outfilename);

// Output data
      output_sda_file << std::setprecision(8);
      for(j = 0; j < Np; j++) {
         output_sda_file << std::setw(20) << EnrKin(coord[j], specie) / one_MeV
                         << std::setw(20) << amp * M_8PI * distro[j] * Sqr(coord[j])
                         << std::endl;
      };

// Close output distro file
      output_sda_file.close();

// Integrate spectrum to obtain number density
      S = 0.0;
      for(j = 0; j < Np; j++) S += distro[j] * Sqr(coord[j]) * dp_arr[j];
      std::cout << S * amp * M_8PI << std::endl;
   };

   return 0;
};
