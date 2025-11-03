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
t_idxs = [97, 90, 68, 0]
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
base = np.loadtxt("dsa_results/dsa_forward_path_dens_pp_analytic.dat")
base0 = base[0:Nz,:]
base1 = base[Nz:2*Nz,:]
base2 = base[2*Nz:3*Nz,:]
base3 = base[3*Nz:4*Nz,:]
splt = np.loadtxt("dsa_results/dsa_forward_path_dens_pp_split.dat")
splt0 = splt[0:Nz,:]
splt1 = splt[Nz:2*Nz,:]
splt2 = splt[2*Nz:3*Nz,:]
splt3 = splt[3*Nz:4*Nz,:]

# Plot solution
p_idx = 1
fig = plt.figure(figsize=(15, 10), layout='tight')

ax1 = fig.add_subplot(221, projection='rectilinear')

ax1.semilogy(x, splt0[:,p_idx], linewidth=3, label="split")
ax1.semilogy(x, base0[:,p_idx] * hint0 * p[p_idx]**alpha, linewidth=3, label="h * f")
ax1.set_xlabel('$x$ (au)', fontsize=20)
ax1.set_ylabel('$f(x,t)$', fontsize=20)
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)
ax1.set_xlim(-3.0,3.0)
ax1.set_ylim(1.0e-7,1.0e-2)
ax1.legend(fontsize=20)
ax1.annotate("$t=${:.1f}".format(t[0]), (1.8, 2.0e-3), fontsize=20)

ax2 = fig.add_subplot(222, projection='rectilinear')

ax2.semilogy(x, splt1[:,p_idx], linewidth=3, label="split")
ax2.semilogy(x, base1[:,p_idx] * hint1 * p[p_idx]**alpha, linewidth=3, label="h * f")
ax2.set_xlabel('$x$ (au)', fontsize=20)
ax2.set_ylabel('$f(x,t)$', fontsize=20)
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
ax2.set_xlim(-3.0,3.0)
ax2.set_ylim(1.0e-7,1.0e-2)
ax2.legend(fontsize=20)
ax2.annotate("$t=${:.1f}".format(t[1]), (1.8, 2.0e-3), fontsize=20)

ax3 = fig.add_subplot(223, projection='rectilinear')

ax3.semilogy(x, splt2[:,p_idx], linewidth=3, label="split")
ax3.semilogy(x, base2[:,p_idx] * hint2 * p[p_idx]**alpha, linewidth=3, label="h * f")
ax3.set_xlabel('$x$ (au)', fontsize=20)
ax3.set_ylabel('$f(x,t)$', fontsize=20)
ax3.tick_params(axis='x', labelsize=20)
ax3.tick_params(axis='y', labelsize=20)
ax3.set_xlim(-3.0,3.0)
ax3.set_ylim(1.0e-7,1.0e-2)
ax3.legend(fontsize=20)
ax3.annotate("$t=${:.1f}".format(t[2]), (1.8, 2.0e-3), fontsize=20)

ax4 = fig.add_subplot(224, projection='rectilinear')

ax4.semilogy(x, splt3[:,p_idx], linewidth=3, label="split")
ax4.semilogy(x, base3[:,p_idx] * hint3 * p[p_idx]**alpha, linewidth=3, label="h * f")
ax4.set_xlabel('$x$ (au)', fontsize=20)
ax4.set_ylabel('$f(x,t)$', fontsize=20)
ax4.tick_params(axis='x', labelsize=20)
ax4.tick_params(axis='y', labelsize=20)
ax4.set_xlim(-3.0,3.0)
ax4.set_ylim(1.0e-7,1.0e-2)
ax4.legend(fontsize=20)
ax4.annotate("$t=${:.1f}".format(t[3]), (1.8, 2.0e-3), fontsize=20)

plt.savefig("likelihood_times_baseline_equals_splitting.png")
plt.show()
plt.close(fig)