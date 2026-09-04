# Import libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys

# Import data
x = np.loadtxt("dsa_results/dsa_analytic_pos.dat")
E = np.loadtxt("dsa_results/dsa_analytic_enr.dat")
dmom = np.loadtxt("dsa_results/dsa_analytic_dmom.dat")
Nz = np.size(x)
Np = np.size(E)
analytic = np.loadtxt("dsa_results/dsa_forward_path_dens_pp_analytic.dat")
analytic_t3 = analytic[3*Nz+Nz//2,:]
baseline = np.loadtxt("dsa_results/dsa_forward_mom_pp_base.dat")
baseline_t3 = baseline[:,1]
dz = x[1] - x[0]
p0 = 14453.5
n_cpu = 128
Qtf = 1.73264
n_traj = 100
one_count_lim = np.ones(Np) * Qtf / (dz * p0 * dmom * n_traj * n_cpu)

# Plot data
fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='rectilinear')

ax.loglog(E, baseline_t3, linewidth=0, marker="o", markersize=7, label="numerical")
ax.loglog(E, analytic_t3, linewidth=4, label="analytic")
ax.loglog(E, one_count_lim, linewidth=2, linestyle="--", label="one-count-limit")
ax.loglog(E, 2.0*one_count_lim, linewidth=2, linestyle="--", label="two-count-limit")
ax.loglog(E, 3.0*one_count_lim, linewidth=2, linestyle="--", label="three-count-limit")
ax.set_xlabel('$E$ (MeV)', fontsize=24)
ax.set_ylabel('$F$ (particles cm$^{-3}$ MeV$^{-1}$)', fontsize=24)
ax.tick_params(axis='x', labelsize=24)
ax.tick_params(axis='y', labelsize=24)
ax.set_xlim(1.0e0, 1.0e6)
ax.legend(fontsize=20)

plt.savefig("dsa_results/dsa_one_count_limit.png")
plt.close(fig)