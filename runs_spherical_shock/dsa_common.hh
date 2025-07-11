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
double p_inj;                    // Injection momentum
double s;                        // Shock strength
double R_sh;                     // Radius of shock surface
double U_up;                     // Upstream flow speed
double U_dn;                     // Downstream flow speed
double DeltaU;                   // Difference between up and downstream flow speed
double kappa_up;                 // Diffusion coefficient upstream
double kappa_dn;                 // Diffusion coefficient downstream
double eta_up;                   // Parameter of interest upstream U*R/kappa
double eta_dn;                   // Parameter of interest downstream U*R/kappa
double xi;                       // Parameter of interest upstream 0.5*eta - 1
double lambda;                   // Parameter of interest (see manuscript)
double zeta1;                    // Parameter of interest (sqrt(lambda) - s)/(s-1)
double zeta2;                    // Parameter of interest -(sqrt(lambda) + s)/(s-1)
double tau;                      // Acceleration time
double Q;                        // Injection rate
double amp;                      // Scaling factor for analytic solution

// Numerical simulation parameters
double dmax;                     // Maximum trajectory displacement away from shock
double w_sh;                     // Shock width
double dmax_fraction;            // Ratio of shock width to impose as maximum trajectory displacement near shock
const int Np = 100;              // Number of momentum bins
double p0;                       // Lower limit in momentum range
double pf;                       // Upper limit in momentum range
double logp0;                    // Logarithm of p0
double logpf;                    // Logarithm of pf
double dlogp;                    // Difference between logarithms of p0 and pf
const int Nr = 100;              // Number of radial bins
int log_rbins;                   // Boolean flag to indicate linear or logarithmic radial bins
double r0;                       // Lower limit of radial range
double rf;                       // Upper limit of radial range
double dr;                       // Radial bin size
const int Nt = 5;                // Number of time bins
double t0;                       // Lower limit of time
double tf;                       // Upper limit of time
double logt0;                    // Logarithm of t0
double logtf;                    // Logarithm of tf
double dlogt;                    // Difference between logarithms of t0 and tf

double p_arr[Np];                // Momentum bin centers
double dp_arr[Np];               // Momentum bin sizes
double r_arr[Nr];                // Radial bin centers
double dr_arr[Nr];               // Radial bin sizes
double r_spectrum;               // Radial location where to plot spectrum
double t_arr[Nt];                // Time bin centers

const int N_params = 17;          // Number of parameters to read from file
double params[N_params];         // Array of parameters
std::ifstream params_file;       // Parameter file

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
   r_spectrum = params[3] * one_au;
   p_inj = Mom(params[4] * one_MeV, specie);
   s = params[5];
   R_sh = params[6] * one_au;
   U_up = params[7] / unit_velocity_fluid;
   U_dn = U_up / s;
   DeltaU = U_up - U_dn;
   kappa_up = params[8] / unit_diffusion_fluid;
   kappa_dn = kappa_up * Sqr(U_dn / U_up);
   eta_up = R_sh * U_up / kappa_up;
   eta_dn = R_sh * U_dn / kappa_dn;
   xi = 0.5 * eta_up - 1.0;
   lambda = Sqr(s + (s - 1.0) * xi) + (s - 1.0) * 2.0 * eta_up * exp(-eta_dn) / (1.0 - exp(-eta_dn));
   zeta1 = (sqrt(lambda) - s) / (s - 1.0);
   zeta2 = -(sqrt(lambda) + s) / (s - 1.0);
   tau = 4.0 * kappa_up / Sqr(U_up);
   Q = params[9];
   amp = 3.0 * Q * s / (2.0 * M_4PI * sqrt(lambda) * Cube(p_inj) * U_up);
   p0 = Mom(params[10] * one_MeV, specie);
   pf = Mom(params[11] * one_MeV, specie);
   logp0 = log10(p0);
   logpf = log10(pf);
   dlogp = (logpf - logp0) / Np;
   log_rbins = params[12];
   r0 = params[13] * one_au;
   rf = params[14] * one_au;
   if (log_rbins) dr = (log10(rf) - log10(r0)) / Nr;
   else dr = (rf - r0) / Nr;
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
      p_arr[i] = pow(10.0, logp0 + (i+0.5) * dlogp);
      dp_arr[i] = pow(10.0, logp0 + (i+1) * dlogp) - pow(10.0, logp0 + i * dlogp);
   };
// Radius
   for (i = 0; i < Nr; i++) {
      if (log_rbins) {
         r_arr[i] = pow(10.0, log10(r0) + (i+0.5) * dr);
         dr_arr[i] = pow(10.0, log10(r0) + (i+1) * dr) - pow(10.0, log10(r0) + i * dr);
      }
      else {
         r_arr[i] = r0 + (i+0.5) * dr;
         dr_arr[i] = dr;
      };
   };
// Time
   for (i = 0; i < Nt; i++) t_arr[i] = pow(10.0, logt0 + i * dlogt);
};

#endif
