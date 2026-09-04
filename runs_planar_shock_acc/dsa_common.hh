#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#ifndef DSA_COMMON_HH
#define DSA_COMMON_HH

using namespace Spectrum;

// AT MOST 1 of the following 3 flags can be defined for any run
// Flag to enable particle splitting
// #define ENABLE_SPLITTING 0
// Flag to enable importance sampling
// #define ENABLE_IMP_SAMP 0
// Flag to enable likelihood test
// #define LIKELIHOOD_TEST

// Importance Sampling Options:
// 0: time-dependent, numerical sampling intended to match particle-splitting
// 1: time-independent, analytic expression from Prinsloo's thesis

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
double Qtf;                   // Total number of particles injected by tf
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
double p_arr_edges[Np+1];     // Momentum bin edges
double dp_arr[Np];            // Momentum bin sizes
double z_arr[Nz];             // Spatial bin centers
double z_arr_edges[Nz+1];     // Spatial bin edges
double z_spectrum;            // Spatial location where to plot spectrum
double t_arr[Nt];             // Time bin centers

const int N_params = 20;      // Number of parameters to read from file
double params[N_params];      // Array of parameters
std::ifstream params_file;    // Parameter file

double alpha;                 // Splitting strength (float)
std::string alpha_string;     // Splitting strength (string)
int n_thrs;                   // Number of momentum thresholds to split in discrete scheme
double A0;                    // Importance sampling strength
const int Nt_IS = 4001;       // Number of time slices in likelihood solution
const int Nx_IS = 10000;      // Number of spatial point in likelihood solution
const int Nx0_IS = 5000;      // Index of x = 0 in likelihood solution
double dt_IS;                 // Time step for importance sampling
double *t_arr_IS;             // Time slices array for importance sampling
double *x_arr_IS;             // Space points array for importance sampling
double *S_CDF;                // CDF of source function
double **Atx;                 // Separable factor of the auxiliary function that depends on t and x
double **dlnAtxdt;            // First derivative of log Atx with respect to t
double **dlnAtxdx;            // First derivative of log Atx with respect to x
double **d2lnAtxdx2;          // Second derivative of log Atx with respect to x

// Declare function here to use before definition
double A(double t, double x, double p);

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
   alpha = params[17];
   n_thrs = params[18];
   A0 = params[19];
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2) << alpha;
   alpha_string = stream.str();
};

// Define initialize momentum, position, and time arrays
void DefineArrays(void)
{
   int i;
// Momentum
   for (i = 0; i < Np; i++) {
      p_arr_edges[i] = pow(10.0, logp0 + i * dlogp);
      p_arr[i] = pow(10.0, logp0 + (i + 0.5) * dlogp);
      dp_arr[i] = pow(10.0, logp0 + (i + 1) * dlogp) - pow(10.0, logp0 + i * dlogp);
   };
   p_arr_edges[Np] = pow(10.0, logpf);
// Position
   for (i = 0; i < Nz; i++) {
      z_arr[i] = z0 + (i + 0.5) * dz;
      z_arr_edges[i] = z0 + i * dz;
   };
   z_arr_edges[Nz] = zf;
// Time
   t_arr[0] = 1.0 * one_day;
   t_arr[1] = 3.0 * one_day;
   t_arr[2] = 10.0 * one_day;
   t_arr[3] = 30.0 * one_day;
   t_arr[4] = 100.0 * one_day;
   tf = t_arr[3];
// Baseline value of Qtf
   Qtf = Q * tf;
};

// Phase-space density upstream (analytic)
inline double N1(double z, double p, double t)
{
   double a = pow(p0 / p, beta);
   double b = 2.0 * z / (tau * U_up);
   double c = sqrt(t / tau);
   double d = 0.5 * beta / c * log(p / p0) - z / (sqrt(t * tau) * U_up);
   return amp * sqrt(Cube(p0 / p)) * exp(U_up * z / (2.0 * kappa_up)) * (exp(b) * erfc(d - c) * a + exp(-b) * erfc(d + c) / a);
};

// Phase-space density downstream (analytic)
inline double N2(double z, double p, double t)
{
   double a = pow(p0 / p, beta);
   double b = 2.0 * z / (tau * U_dn);
   double c = sqrt(t / tau);
   double d = 0.5 * beta / c * log(p / p0) + z / (sqrt(t * tau) * U_dn);
   return amp * sqrt(Cube(p0 / p)) * exp(U_dn * z / (2.0 * kappa_dn)) * (exp(-b) * erfc(d - c) * a + exp(b) * erfc(d + c) / a);
};

// Phase-space density upstream AND downstream (analytic)
inline double N12(double z, double p, double t)
{
   if (z < 0.0) return N1(z, p, t);
   else return N2(z, p, t);
};

// Integral in momentum (analytic + Riemann sum)
double MomentumIntegral(double p1_in, double p2_in, int Np_in,
                        double z, double t)
{
   int i;
   double S = 0.0;

// Initialize arrays to integrate
   double *_p_arr = new double[Np_in];
   double *_dp_arr = new double[Np_in];
   double _logp1 = log10(p1_in);
   double _logp2 = log10(p2_in);
   double _dlogp = (_logp2 - _logp1) / Np_in;
   for (i = 0; i < Np_in; i++) {
      _p_arr[i] = pow(10.0, _logp1 + (i + 0.5) * _dlogp);
      _dp_arr[i] = pow(10.0, _logp1 + (i + 1) * _dlogp) - pow(10.0, _logp1 + i * _dlogp);
   };

// Integrate
   for (i = 0; i < Np_in; i++) S += N12(z, _p_arr[i], t) * Sqr(_p_arr[i]) * _dp_arr[i];

// Free memory
   delete[] _p_arr;
   delete[] _dp_arr;

   return S * M_4PI;
};

// Double integral in momentum and space (analytic + Riemann sum)
double SpaceMomentumIntegral(double z1_in, double z2_in, int Nz_in,
                             double p1_in, double p2_in, int Np_in,
                             double t)
{
   int i;
   double S = 0.0;

// Initialize arrays to integrate
   double *_z_arr = new double[Nz_in];
   double _dz = (z2_in - z1_in) / Nz_in;
   for (i = 0; i < Nz_in; i++) _z_arr[i] = z1_in + (i + 0.5) * _dz;

// Integrate
   for (i = 0; i < Nz_in; i++) S += MomentumIntegral(p1_in, p2_in, Np_in, _z_arr[i], t) * _dz;

// Free memory
   delete[] _z_arr;

   return S;
};

// Define importance sampling helper arrays
void ImportanceSamplingSetup(void)
{
#if defined(ENABLE_IMP_SAMP)
   int i, j, k;
#if ENABLE_IMP_SAMP == 1
   Qtf = Q * tf * A(0.0, 0.0, p0);
#else
   std::ifstream likelihood_file;

// Import time slice array
   t_arr_IS = new double[Nt_IS];
   likelihood_file.open("dsa_results/feynman_kac/dsa_likelihood_t_alpha=" + alpha_string + ".dat");
   for (i = 0; i < Nt_IS; i++) {
// Remember to swap rows (invert time) because h is solved for and saved from tf to 0
      k = Nt_IS - i - 1;
      likelihood_file >> t_arr_IS[k];
   };
   likelihood_file.close();

// Import space axis arrays
   x_arr_IS = new double[Nx_IS];
   likelihood_file.open("dsa_results/feynman_kac/dsa_likelihood_x_alpha=" + alpha_string + ".dat");
   for (i = 0; i < Nx_IS; i++) likelihood_file >> x_arr_IS[i];
   likelihood_file.close();

// Import likelihood
   Atx = Create2D<double>(Nt_IS, Nx_IS);
   likelihood_file.open("dsa_results/feynman_kac/dsa_likelihood_alpha=" + alpha_string + ".dat");
   for (i = 0; i < Nt_IS; i++) {
      k = Nt_IS - i - 1;
      for (j = 0; j < Nx_IS; j++) likelihood_file >> Atx[k][j];
   };
   likelihood_file.close();

// Define source CDF assuming t_arr_IS[Nt_IS-1] == tf AND x_arr_IS[Nx0_IS] == 0.
   S_CDF = new double[Nt_IS];
   S_CDF[0] = 0.0;
   dt_IS = t_arr_IS[1] - t_arr_IS[0];
// Integrate PDF to get CDF
   for (i = 1; i < Nt_IS; i++) {
      S_CDF[i] = S_CDF[i-1] + 0.5 * Q * dt_IS
               * (Atx[i-1][Nx0_IS] + Atx[i][Nx0_IS]);
   };
// Normalize and update Qtf
   Qtf = S_CDF[Nt_IS-1];
   for (i = 1; i < Nt_IS; i++) {
      S_CDF[i] /= Qtf;
   };

// Import derivatives of likelihood
   dlnAtxdt = Create2D<double>(Nt_IS, Nx_IS);
   likelihood_file.open("dsa_results/feynman_kac/dsa_likelihood_dt_alpha=" + alpha_string + ".dat");
   for (i = 0; i < Nt_IS; i++) {
      k = Nt_IS - i - 1;
      for (j = 0; j < Nx_IS; j++) likelihood_file >> dlnAtxdt[k][j];
   };
   likelihood_file.close();

   dlnAtxdx = Create2D<double>(Nt_IS, Nx_IS);
   likelihood_file.open("dsa_results/feynman_kac/dsa_likelihood_dx_alpha=" + alpha_string + ".dat");
   for (i = 0; i < Nt_IS; i++) {
      k = Nt_IS - i - 1;
      for (j = 0; j < Nx_IS; j++) likelihood_file >> dlnAtxdx[k][j];
   };
   likelihood_file.close();

   d2lnAtxdx2 = Create2D<double>(Nt_IS, Nx_IS);
   likelihood_file.open("dsa_results/feynman_kac/dsa_likelihood_dx2_alpha=" + alpha_string + ".dat");
   for (i = 0; i < Nt_IS; i++) {
      k = Nt_IS - i - 1;
      for (j = 0; j < Nx_IS; j++) likelihood_file >> d2lnAtxdx2[k][j];
   };
   likelihood_file.close();
#endif
#endif
};

// Free memory from importance sampling helper arrays
void ImportanceSamplingFree(void)
{
#if defined(ENABLE_IMP_SAMP)
#if ENABLE_IMP_SAMP == 0
   delete[] t_arr_IS;
   delete[] x_arr_IS;
   delete[] S_CDF;
   Delete2D(Atx);
   Delete2D(dlnAtxdt);
   Delete2D(dlnAtxdx);
   Delete2D(d2lnAtxdx2);
#endif
#endif
};

// Interpolation function for nonuniform distance between points
double InterpSCDFinv(double x)
{
   int i = LocateInArray(0, Nt_IS-1, S_CDF, x, false);
   if (i == -1) return 0.0;
   return ( t_arr_IS[i  ] * (S_CDF[i+1] - x)
          + t_arr_IS[i+1] * (x - S_CDF[i  ]) )
          / (S_CDF[i+1] - S_CDF[i]);
};

double InterpAtx_arr(double t, double x, double **Atx_arr)
{
   int i = t / dt_IS;
   if (i < 0 || Nt_IS-1 <= i) return 0.0;
   double dt_f = (t - t_arr_IS[i]) / dt_IS;
   int j = LocateInArray(0, Nx_IS-1, x_arr_IS, x, false);
   if (j == -1) return 0.0;
   double dx_f = (x - x_arr_IS[j]) / (x_arr_IS[j+1] - x_arr_IS[j]);
   return   Atx_arr[i][j  ] * (1.0 - dt_f) * (1.0 - dx_f) + Atx_arr[i+1][j  ] * dt_f * (1.0 - dx_f)
          + Atx_arr[i][j+1] * (1.0 - dt_f) * dx_f         + Atx_arr[i+1][j+1] * dt_f * dx_f;
};

// Initial time source function
inline double init_time(double s)
{
#if defined(ENABLE_IMP_SAMP)
#if ENABLE_IMP_SAMP == 1
   return s * tf;
#else
   return InterpSCDFinv(s);
#endif
#else
   return s * tf;
#endif
};

// Importance sampling auxiliary function and its derivatives
inline double A(double t, double x, double p)
{
#if defined(ENABLE_IMP_SAMP)
#if ENABLE_IMP_SAMP == 1
   return pow(cosh(x / w_sh), -A0 * w_sh);
#else
   double _A = 1.0;
   if (t < tf) _A = InterpAtx_arr(t, x, Atx);
   return _A * pow(p / p0, alpha);
#endif
#else
   return 1.0;
#endif
};
inline double dlnAdt(double t, double x)
{
#if ENABLE_IMP_SAMP == 1
   return 0.0;
#else
   return InterpAtx_arr(t, x, dlnAtxdt);
#endif
}
inline double dlnAdx(double t, double x)
{
#if ENABLE_IMP_SAMP == 1
   return -A0 * tanh(x / w_sh);
#else
   return InterpAtx_arr(t, x, dlnAtxdx);
#endif
};
inline double d2lnAdx2(double t, double x)
{
#if ENABLE_IMP_SAMP == 1
   return -(A0 / w_sh) * (1.0 - Sqr(tanh(x / w_sh)));
#else
   return InterpAtx_arr(t, x, d2lnAtxdx2);
#endif
};
inline double dlnAdp(double p)
{
#if ENABLE_IMP_SAMP == 1
   return 0.0;
#else
   return alpha / p;
#endif
};

#endif
