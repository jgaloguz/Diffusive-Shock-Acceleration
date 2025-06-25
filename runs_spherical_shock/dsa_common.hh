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
const double p_inj = Mom(1.0 * one_MeV, specie);
const double s = 4.0;
const double R_sh = 80.0 * one_au;
const double U_up = 4.0e7 / unit_velocity_fluid;
const double U_dn = U_up / s;
const double DeltaU = U_up - U_dn;
const double kappa_up = 1.0e21 / unit_diffusion_fluid;
const double kappa_dn = kappa_up * Sqr(U_dn / U_up);
const double eta_up = R_sh * U_up / kappa_up;
const double eta_dn = R_sh * U_dn / kappa_dn;
const double xi = 0.5 * eta_up - 1.0;
const double lambda = Sqr(s + (s - 1.0) * xi) + (s - 1.0) * 2.0 * eta_up * exp(-eta_dn) / (1.0 - exp(-eta_dn));
const double zeta1 = (sqrt(lambda) - s) / (s - 1.0);
const double zeta2 = -(sqrt(lambda) + s) / (s - 1.0);
const double tau = 4.0 * kappa_up / Sqr(U_up);
const double Q = 1.0;
const double amp = 3.0 * s / (8.0 * M_PI * sqrt(lambda) * Cube(p_inj) * U_up);

// const double beta = 1.5 * (s + 1.0) / (s - 1.0);
// const double amp = 3.0 * Q / (2.0 * M_4PI * DeltaU * Cube(p0));

// Numerical integration parameters
const int Np = 100;
const double p0 = Mom(0.01 * one_MeV, specie);
const double pf = Mom(100.0 * one_MeV, specie);
const double logp0 = log10(p0);
const double logpf = log10(pf);
const double dlogp = (logpf - logp0) / Np;
const int Nr = 100;
const double r0 = R_sh - 20.0 * one_au;
const double rf = R_sh + 40.0 * one_au;
const double dr = (rf - r0) / Nr;
const int Nt = 5;
const double t0 = 1.0 * one_day;
const double tf = 100.0 * one_day;
const double logt0 = log10(t0);
const double logtf = log10(tf);
const double dlogt = (logtf - logt0) / (Nt - 1);

double p_arr[Np];
double dp_arr[Np];
double r_arr[Nr];
double r_shock;
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
   for (i = 0; i < Nr; i++) r_arr[i] = r0 + (i+0.5) * dr;
   for (i = Nr-1; i >= 0; i--) {
      if (r_arr[i] <= R_sh) break;
      else r_shock = r_arr[i];
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
