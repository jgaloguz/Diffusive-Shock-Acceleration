# Import libraries
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import sys

# Load parameters
params = np.loadtxt("dsa_results/likelihood_params_FD.dat")
Nx_FD = params[0]
Nx_lft = params[1]
Nx_ctr = params[2]
Nx_rgt = params[3]
X1_FD = params[4]
X2_FD = params[5]
X3_FD = params[6]
X4_FD = params[7]
dt_FD = params[8]
Tf_FD = params[9]
Nt_out = params[10]
alpha = params[11]
t = np.loadtxt("dsa_results/dsa_analytic_time.dat")
x = np.loadtxt("dsa_results/dsa_analytic_pos.dat")
p = np.loadtxt("dsa_results/dsa_analytic_mom.dat")
dp = np.loadtxt("dsa_results/dsa_analytic_dmom.dat")
E = np.loadtxt("dsa_results/dsa_analytic_enr.dat")
Nz = np.size(x)
Np = np.size(p)

# Load coefficients
coeffs = np.loadtxt("dsa_results/likelihood_coeffs_FD.dat")
X_FD = coeffs[:,0]
dX_FD = X_FD[1:] - X_FD[:-1]
Cl = 1.0 / (dX_FD[1:] + dX_FD[:-1]) - 1.0 / dX_FD[:-1]
Cc = 1.0 / dX_FD[:-1] - 1.0 / dX_FD[1:]
Cr = 1.0 / dX_FD[1:] - 1.0 / (dX_FD[1:] + dX_FD[:-1])
U = coeffs[:,1]
K = coeffs[:,2]
dU = coeffs[:,3]

# Import data
t_idxs = [290, 270, 200, 0]
h = np.loadtxt("dsa_results/likelihood_solution_FD.dat")
ht = h[t_idxs[0],:]
h_interp = interp1d(X_FD, ht, kind='linear') 
hint0 = h_interp(x)
ht = h[t_idxs[1],:]
h_interp = interp1d(X_FD, ht, kind='linear') 
hint1 = h_interp(x)
ht = h[t_idxs[2],:]
h_interp = interp1d(X_FD, ht, kind='linear') 
hint2 = h_interp(x)
ht = h[t_idxs[3],:]
h_interp = interp1d(X_FD, ht, kind='linear') 
hint3 = h_interp(x)
base = np.loadtxt("dsa_results/dsa_forward_path_dens_pp_baseline.dat")
base0 = base[0:Nz,:]
base1 = base[Nz:2*Nz,:]
base2 = base[2*Nz:3*Nz,:]
base3 = base[3*Nz:4*Nz,:]
splt = np.loadtxt("dsa_results/dsa_forward_path_dens_pp_split.dat")
splt0 = splt[0:Nz,:]
splt1 = splt[Nz:2*Nz,:]
splt2 = splt[2*Nz:3*Nz,:]
splt3 = splt[3*Nz:4*Nz,:]
test = np.loadtxt("dsa_results/dsa_forward_path_dens_pp_ltest.dat")
test0 = test[0:Nz,:]
test1 = test[Nz:2*Nz,:]
test2 = test[2*Nz:3*Nz,:]
test3 = test[3*Nz:4*Nz,:]

# Integrate over momentum
base0_x = np.zeros(Nz)
base1_x = np.zeros(Nz)
base2_x = np.zeros(Nz)
base3_x = np.zeros(Nz)
baseh0_x = np.zeros(Nz)
baseh1_x = np.zeros(Nz)
baseh2_x = np.zeros(Nz)
baseh3_x = np.zeros(Nz)
splt0_x = np.zeros(Nz)
splt1_x = np.zeros(Nz)
splt2_x = np.zeros(Nz)
splt3_x = np.zeros(Nz)
test0_x = np.zeros(Nz)
test1_x = np.zeros(Nz)
test2_x = np.zeros(Nz)
test3_x = np.zeros(Nz)
for j in range(Np):
   pj_alpha = p[j]**alpha
   base0_x = base0_x + base0[:,j] * dp[j]
   base1_x = base1_x + base1[:,j] * dp[j]
   base2_x = base2_x + base2[:,j] * dp[j]
   base3_x = base3_x + base3[:,j] * dp[j]
   baseh0_x = baseh0_x + base0[:,j] * hint0 * pj_alpha * dp[j]
   baseh1_x = baseh1_x + base1[:,j] * hint1 * pj_alpha * dp[j]
   baseh2_x = baseh2_x + base2[:,j] * hint2 * pj_alpha * dp[j]
   baseh3_x = baseh3_x + base3[:,j] * hint3 * pj_alpha * dp[j]
   splt0_x = splt0_x + splt0[:,j] * dp[j]
   splt1_x = splt1_x + splt1[:,j] * dp[j]
   splt2_x = splt2_x + splt2[:,j] * dp[j]
   splt3_x = splt3_x + splt3[:,j] * dp[j]
   test0_x = test0_x + test0[:,j] * dp[j]
   test1_x = test1_x + test1[:,j] * dp[j]
   test2_x = test2_x + test2[:,j] * dp[j]
   test3_x = test3_x + test3[:,j] * dp[j]

# Integrate number of particles as a sanity check
base0_n = 0.0
base1_n = 0.0
base2_n = 0.0
base3_n = 0.0
splt0_n = 0.0
splt1_n = 0.0
splt2_n = 0.0
splt3_n = 0.0
dz = x[1] - x[0]
for i in range(Nz):
   base0_n = base0_n + base0_x[i] * dz
   base1_n = base1_n + base1_x[i] * dz
   base2_n = base2_n + base2_x[i] * dz
   base3_n = base3_n + base3_x[i] * dz
   splt0_n = splt0_n + splt0_x[i] * dz
   splt1_n = splt1_n + splt1_x[i] * dz
   splt2_n = splt2_n + splt2_x[i] * dz
   splt3_n = splt3_n + splt3_x[i] * dz

print("Baseline number of particles:")
print("{:10.3e}{:10.3e}{:10.3e}{:10.3e}".format(base0_n, base1_n, base2_n, base3_n))
print("Split number of particles:")
print("{:10.3e}{:10.3e}{:10.3e}{:10.3e}".format(splt0_n, splt1_n, splt2_n, splt3_n))

# Plot solution
p_idx = 1
fig = plt.figure(figsize=(15, 10), layout='tight')

ax1 = fig.add_subplot(221, projection='rectilinear')

ax1.semilogy(x, base0_x, linewidth=3, label="f")
ax1.semilogy(x, baseh0_x, linewidth=3, label="h * f")
ax1.semilogy(x, splt0_x, linewidth=3, linestyle="--", label="split")
# ax1.semilogy(x, test0_x, linewidth=3, label="ltest")
ax1.set_xlabel('$x$ (au)', fontsize=20)
ax1.set_ylabel('$f(x,t)$', fontsize=20)
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)
ax1.set_xlim(-3.0,3.0)
ax1.set_ylim(1.0e-7,1.0e-1)
ax1.legend(fontsize=20)
ax1.annotate("$t=${:.0f}".format(t[0]), (-2.8, 1.5e-7), fontsize=20)

ax2 = fig.add_subplot(222, projection='rectilinear')

ax2.semilogy(x, base1_x, linewidth=3, label="f")
ax2.semilogy(x, baseh1_x, linewidth=3, label="h * f")
ax2.semilogy(x, splt1_x, linewidth=3, linestyle="--", label="split")
# ax2.semilogy(x, test1_x, linewidth=3, label="ltest")
ax2.set_xlabel('$x$ (au)', fontsize=20)
ax2.set_ylabel('$f(x,t)$', fontsize=20)
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
ax2.set_xlim(-3.0,3.0)
ax2.set_ylim(1.0e-7,1.0e-1)
ax2.legend(fontsize=20)
ax2.annotate("$t=${:.0f}".format(t[1]), (-2.8, 1.5e-7), fontsize=20)

ax3 = fig.add_subplot(223, projection='rectilinear')

ax3.semilogy(x, base2_x, linewidth=3, label="f")
ax3.semilogy(x, baseh2_x, linewidth=3, label="h * f")
ax3.semilogy(x, splt2_x, linewidth=3, linestyle="--", label="split")
# ax3.semilogy(x, test2_x, linewidth=3, label="ltest")
ax3.set_xlabel('$x$ (au)', fontsize=20)
ax3.set_ylabel('$f(x,t)$', fontsize=20)
ax3.tick_params(axis='x', labelsize=20)
ax3.tick_params(axis='y', labelsize=20)
ax3.set_xlim(-3.0,3.0)
ax3.set_ylim(1.0e-7,1.0e-1)
ax3.legend(fontsize=20)
ax3.annotate("$t=${:.0f}".format(t[2]), (-2.8, 1.5e-7), fontsize=20)

ax4 = fig.add_subplot(224, projection='rectilinear')

ax4.semilogy(x, base3_x, linewidth=3, label="f")
ax4.semilogy(x, baseh3_x, linewidth=3, label="h * f")
ax4.semilogy(x, splt3_x, linewidth=3, linestyle="--", label="split")
# ax4.semilogy(x, test3_x, linewidth=3, label="ltest")
ax4.set_xlabel('$x$ (au)', fontsize=20)
ax4.set_ylabel('$f(x,t)$', fontsize=20)
ax4.tick_params(axis='x', labelsize=20)
ax4.tick_params(axis='y', labelsize=20)
ax4.set_xlim(-3.0,3.0)
ax4.set_ylim(1.0e-7,1.0e-1)
ax4.legend(fontsize=20)
ax4.annotate("$t=${:.0f}".format(t[3]), (-2.8, 1.5e-7), fontsize=20)

plt.savefig("dsa_results/likelihood_times_baseline_equals_splitting.png")
plt.show()
plt.close(fig)