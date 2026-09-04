from scipy.interpolate import interp1d
import numpy as np

def InterpLikelihood(filename, x):
   """
   Load conditional path density from `filename` and interpolated it on new grid `x`.
   """
   h = np.loadtxt(filename)
   X_FD = h[0, :]
   ht = h[1, :]
   h_interp = interp1d(X_FD, ht, kind='linear')
   hint0 = h_interp(x)
   ht = h[2, :]
   h_interp = interp1d(X_FD, ht, kind='linear')
   hint1 = h_interp(x)
   ht = h[3, :]
   h_interp = interp1d(X_FD, ht, kind='linear')
   hint2 = h_interp(x)
   ht = h[4, :]
   h_interp = interp1d(X_FD, ht, kind='linear')
   hint3 = h_interp(x)
   return [hint0, hint1, hint2, hint3]

def CondPathDens(filename, factor_h, hint, x):
   """
   Compute conditional path density in spatial grid `x` at four different times by integrating over momentum.

   `filename` contains the post-processed path density.

   `factor_h` (True or False) controls whether or not the interpolated likelihood `hint` is factored into the integral.

   Outputs conditional path densities at four times, followed by the total number of particles at those same times.
   """
# Load parameters
   alpha = 1.0
   p = np.loadtxt("dsa_results/dsa_analytic_mom.dat")
   dp = np.loadtxt("dsa_results/dsa_analytic_dmom.dat")
   Nz = np.size(x)
   Np = np.size(p)
   dz = x[1] - x[0]

# Import conditional path densities
   CPD = np.loadtxt(filename)
   CPD0 = CPD[0:Nz,:]
   CPD1 = CPD[Nz:2*Nz,:]
   CPD2 = CPD[2*Nz:3*Nz,:]
   CPD3 = CPD[3*Nz:4*Nz,:]
# Integrate over momentum
   CPD0_x = np.zeros(Nz)
   CPD1_x = np.zeros(Nz)
   CPD2_x = np.zeros(Nz)
   CPD3_x = np.zeros(Nz)
   if factor_h:
      for j in range(Np):
         pj_alpha = p[j]**alpha
         CPD0_x = CPD0_x + CPD0[:,j] * hint[0] * pj_alpha * dp[j]
         CPD1_x = CPD1_x + CPD1[:,j] * hint[1] * pj_alpha * dp[j]
         CPD2_x = CPD2_x + CPD2[:,j] * hint[2] * pj_alpha * dp[j]
         CPD3_x = CPD3_x + CPD3[:,j] * hint[3] * pj_alpha * dp[j]
   else:
      for j in range(Np):
         CPD0_x = CPD0_x + CPD0[:,j] * dp[j]
         CPD1_x = CPD1_x + CPD1[:,j] * dp[j]
         CPD2_x = CPD2_x + CPD2[:,j] * dp[j]
         CPD3_x = CPD3_x + CPD3[:,j] * dp[j]
# Integrate number of particles as a sanity check
   CPD0_n = 0.0
   CPD1_n = 0.0
   CPD2_n = 0.0
   CPD3_n = 0.0
   for i in range(Nz):
      CPD0_n = CPD0_n + CPD0_x[i] * dz
      CPD1_n = CPD1_n + CPD1_x[i] * dz
      CPD2_n = CPD2_n + CPD2_x[i] * dz
      CPD3_n = CPD3_n + CPD3_x[i] * dz
# Return arrays
   return CPD0_x, CPD1_x, CPD2_x, CPD3_x, \
          CPD0_n, CPD1_n, CPD2_n, CPD3_n