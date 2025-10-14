#include <fstream>

#ifndef DSA_COMMON_HH
#define DSA_COMMON_HH

using namespace Spectrum;

#define ENABLE_SPLITTING                              // Flag to enable particle splitting
#define ENABLE_IMPORTANCE                             // Flag to enable importance sampling
#define n_traj 1000                                   // Number of trajectories per process
#define n_chld 2                                      // Number of child particles per split
#define n_thrs 10                                     // Number of momentum thresholds to split
#define max_splt 100000                               // Maximum number of splits per trajectory
#define A0 0.1                                        // Importance sampling constant

// Constants
const int specie = Specie::proton;
const double one_MeV = SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle;
const double one_au = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
const double one_day = 24.0 * 60.0 * 60.0 / unit_time_fluid;

// Shock parameters
double p0;                    // Injection momentum. Also lower limit of momentum range
double s;                     // Shock strength
double U_up;                  // Upstream flow speed
double U_dn;                  // Downstream flow speed
double DeltaU;                // Difference between up and downstream flow speed
double kappa_up;              // Diffusion coefficient upstream
double kappa_dn;              // Diffusion coefficient downstream
double beta;                  // Parameter of interest (3/2)*(s+1)/(s-1)
double tau;                   // Acceleration time
double Q;                     // Injection rate
double amp;                   // Scaling factor for analytic solution

// Numerical simulation parameters
double dmax;                  // Maximum trajectory displacement away from shock
double w_sh;                  // Shock width
double dmax_fraction;         // Ratio of shock width to impose as maximum trajectory displacement near shock
const int Np = 100;           // Number of momentum bins
double pf;                    // Upper limit in momentum range
double logp0;                 // Logarithm of p0
double logpf;                 // Logarithm of pf
double dlogp;                 // Difference between logarithms of p0 and pf
const int Nz = 100;           // Number of spatial bins
double z0;                    // Lower limit of spatial range
double zf;                    // Upper limit of spatial range
double dz;                    // Spatial bin size
const int Nt = 5;             // Number of time bins
double t0;                    // Lower limit of time
double tf;                    // Upper limit of time
double logt0;                 // Logarithm of t0
double logtf;                 // Logarithm of tf
double dlogt;                 // Difference between logarithms of t0 and tf

double p_arr[Np];             // Momentum bin centers
double dp_arr[Np];            // Momentum bin sizes
double z_arr[Nz];             // Spatial bin centers
double z_spectrum;            // Spatial location where to plot spectrum
double t_arr[Nt];             // Time bin centers

const int N_params = 17;      // Number of parameters to read from file
double params[N_params];      // Array of parameters
std::ifstream params_file;    // Parameter file

// Read parameters from file
void ReadParams(void)
{
   params_file.open("params.dat");
   for (int i = 0; i < N_params; i++) params_file >> params[i];
   params_file.close();

// Unpack parameters
   dmax = params[0] * one_au;
   w_sh = params[1] * one_au;
   dmax_fraction = params[2];
   z_spectrum = params[3] * one_au;
   p0 = Mom(params[4] * one_MeV, specie);
   s = params[5];
   U_up = params[7] / unit_velocity_fluid;
   U_dn = U_up / s;
   DeltaU = U_up - U_dn;
   kappa_up = params[8] / unit_diffusion_fluid;
   kappa_dn = kappa_up * Sqr(U_dn / U_up);
   beta = 1.5 * (s + 1.0) / (s - 1.0);
   tau = 4.0 * kappa_up / Sqr(U_up);
   Q = params[9];
   amp = 3.0 * Q / (M_8PI * DeltaU * Cube(p0));
   pf = Mom(params[11] * one_MeV, specie);
   logp0 = log10(p0);
   logpf = log10(pf);
   dlogp = (logpf - logp0) / Np;
   z0 = params[13] * one_au;
   zf = params[14] * one_au;
   dz = (zf - z0) / Nz;
   t0 = params[15] * one_day;
   tf = params[16] * one_day;
   logt0 = log10(t0);
   logtf = log10(tf);
   dlogt = (logtf - logt0) / (Nt - 1);
};

// Define initialize momentum, position, and time arrays
void DefineArrays(void)
{
   int i;
// Momentum
   for (i = 0; i < Np; i++) {
      p_arr[i] = pow(10.0, logp0 + (i + 0.5) * dlogp);
      dp_arr[i] = pow(10.0, logp0 + (i + 1) * dlogp) - pow(10.0, logp0 + i * dlogp);
   };
// Position
   for (i = 0; i < Nz; i++) z_arr[i] = z0 + (i + 0.5) * dz;
// Time
   for (i = 0; i < Nt; i++) t_arr[i] = pow(10.0, logt0 + i * dlogt);
};

#endif
