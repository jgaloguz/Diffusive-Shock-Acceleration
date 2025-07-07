#include "src/server_config.hh"
#include "src/background_solarwind_termshock.hh"
#include "dsa_common.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>

using namespace Spectrum;

int main(int argc, char** argv)
{
   int i;
   BackgroundSolarWindTermShock background;
   std::ofstream plot_file;
   SpatialData spdata;
   DataContainer container;
   GeoVector pos, mom;

   ReadParams();
   DefineArrays();
   spdata._mask = BACKGROUND_U | BACKGROUND_B | BACKGROUND_gradU;

   container.Clear();

// Initial time
   container.Insert(0.0);

// Origin
   container.Insert(gv_zeros);

// Upstream velocity
   GeoVector u0(U_up, 0.0, 0.0);
   container.Insert(u0);

// Upstream magnetic field
   double RS = 6.957e10 / unit_length_fluid;
   double r_ref = 3.0 * RS;
   double BmagE = 5.0e-5 / unit_magnetic_fluid;
   double Bmag_ref = BmagE * Sqr(one_au / r_ref);
   GeoVector B0(Bmag_ref, 0.0, 0.0);
   container.Insert(B0);

// Maximum displacement
   container.Insert(dmax);

// Solar rotation vector
   GeoVector Omega(0.0, 0.0, 0.0);
   container.Insert(Omega);

// Reference equatorial distance
   container.Insert(r_ref);

// dmax fraction
   container.Insert(dmax_fraction);

// Spherical shock radius
   container.Insert(R_sh); 

// Spherical shock width
   container.Insert(w_sh);
   
// Spherical shock strength
   container.Insert(s);

   background.SetupObject(container);

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Parameters for plotting
   int pts = Nr;
   if(argc > 1) pts = atoi(argv[1]);
   double dr_plot = (rf - r0) / (pts-1);
   pos[0] = r0;
   pos[1] = 0.0;
   pos[2] = 0.0;
   mom[0] = p_inj;

// Plot
   plot_file.open("dsa_results/solarwind.dat");
   for (i = 0; i < pts; i++) {
      background.GetFields(0.0, pos, mom, spdata);
      plot_file << std::setw(18) << pos[0]
                << std::setw(18) << spdata.Uvec.Norm() * unit_velocity_fluid
                << std::setw(18) << spdata.Bmag * unit_magnetic_fluid
                << std::setw(18) << spdata.divU() * unit_velocity_fluid / unit_length_fluid
                << std::endl;
      pos[0] += dr_plot;
   };
   plot_file.close();

   return 0;
};


