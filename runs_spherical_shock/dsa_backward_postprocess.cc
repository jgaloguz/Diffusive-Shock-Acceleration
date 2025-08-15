#include "common/physics.hh"
#include "dsa_common.hh"
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace Spectrum;

int main(int argc, char** argv)
{
   int i = Nt-1, j, k; 
   if(argc > 1) i = atoi(argv[1]);
   std::string infilename;
   std::string outfilename;
   std::string line;
   int sum_c[Np];
   double coord[Np], distro[Np], sum_w[Np];
   double S;
   ReadParams();
   DefineArrays();
   double factor = (sqrt(lambda) / s) * (U_up / DeltaU);

// Spectrum vs momentum
   infilename = "dsa_results/dsa_backward_mom_" + std::to_string(i) + ".dat";
   outfilename = "dsa_results/dsa_backward_mom_" + std::to_string(i) + "_pp.dat";

// Open input analytic distro file
   std::ifstream input_dsa_file(infilename);

// Read first two lines of distro file
   std::getline(input_dsa_file, line);
   std::getline(input_dsa_file, line);

// Read data
   for(j = 0; j < Np; j++) {
      input_dsa_file >> coord[j];
      input_dsa_file >> distro[j];
      input_dsa_file >> sum_w[j];
      input_dsa_file >> sum_c[j];
   };

// Close input cartesian distro file
   input_dsa_file.close();

// Open output distro file
   std::ofstream output_dsa_file(outfilename);

// Output data
   output_dsa_file << std::setprecision(8);
   for(j = 0; j < Np; j++) {
      output_dsa_file << std::setw(20) << EnrKin(coord[j], specie) / one_MeV
                      << std::setw(20) << amp * factor * M_8PI * distro[j] * Sqr(coord[j])
                      << std::endl;
   };

// Close output distro file
   output_dsa_file.close();

// Number density vs position
   outfilename = "dsa_results/dsa_backward_pos_" + std::to_string(i) + "_pp.dat";

// Integrate spectrum to obtain number density
   S = 0.0;
   for(j = 0; j < Np; j++) S += distro[j] * Sqr(coord[j]) * dp_arr[j];

// Open output number density file
   output_dsa_file.open(outfilename);

   output_dsa_file << std::setprecision(8);
   output_dsa_file << std::setw(20) << r_spectrum
                   << std::setw(20) << amp * factor * M_8PI * S 
                   << std::endl;

// Close output number density file
   output_dsa_file.close();

   std::cout << "Post-processed distribution files outputed." << std::endl;
   return 0;
};
