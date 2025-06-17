#include <fstream>

#ifndef DSA_COMMON_HH
#define DSA_COMMON_HH

using namespace Spectrum;

// Constants
const int specie = Specie::proton;
const double one_MeV = SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle;
const double one_au = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
const double one_day = 24.0 * 60.0 * 60.0 / unit_time_fluid;

// Shock parameters
const double p0 = Mom(1.0 * one_MeV, specie);
const double s = 4.0;
const double U_up = 4.0e7 / unit_velocity_fluid;
const double U_dn = U_up / s;
const double DeltaU = U_up - U_dn;
const double kappa_up = 1.0e20 / unit_diffusion_fluid;
const double kappa_dn = kappa_up * Sqr(U_dn / U_up);
const double beta = 1.5 * (s + 1.0) / (s - 1.0);
const double tau = 4.0 * kappa_up / Sqr(U_up);
const double Q = 1.0;
const double amp = 3.0 * Q / (2.0 * M_4PI * DeltaU * Cube(p0));

// Numerical integration parameters
const int Np = 100;
const double pf = Mom(100.0 * one_MeV, specie);
const double logp0 = log10(p0);
const double logpf = log10(pf);
const double dlogp = (logpf - logp0) / Np;
const int Nz = 100;
const double z0 = -1.0 * one_au;
const double zf = 4.0 * one_au;
const double dz = (zf - z0) / Nz;
const int Nt = 5;
const double t0 = 1.0 * one_day;
const double tf = 100.0 * one_day;
const double logt0 = log10(t0);
const double logtf = log10(tf);
const double dlogt = (logtf - logt0) / (Nt - 1);

double p_arr[Np];
double dp_arr[Np];
double z_arr[Nz];
double z_shock;
double t_arr[Nt];

const int N_params = 4;
double params[N_params];
std::ifstream params_file;

// Define initialize momentum, position, and time arrays
void DefineArrays(void)
{
   int i;
// Momentum
   for (i = 0; i < Np; i++) {
      p_arr[i] = pow(10.0, logp0 + (i+0.5) * dlogp);
      dp_arr[i] = pow(10.0, logp0 + (i+1) * dlogp) - pow(10.0, logp0 + i * dlogp);
   };
// Position
   for (i = 0; i < Nz; i++) z_arr[i] = z0 + (i+0.5) * dz;
   for (i = Nz-1; i >= 0; i--) {
      if (z_arr[i] <= 0.0) break;
      else z_shock = z_arr[i];
   };
// Time
   for (i = 0; i < Nt; i++) t_arr[i] = pow(10.0, logt0 + i * dlogt);
};

// Read parameters from file
void ReadParams(void)
{
   params_file.open("params.dat");
   for (int i = 0; i < N_params; i++) params_file >> params[i];
   params_file.close();
}

#endif
