# Import libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys

# Import data
x = np.loadtxt("dsa_results/dsa_analytic_pos.dat")
E = np.loadtxt("dsa_results/dsa_analytic_enr.dat")
Nz = np.size(x)
Np = np.size(E)
analytic = np.loadtxt("dsa_results/dsa_forward_path_dens_pp_analytic.dat")
analytic_t3 = analytic[3*Nz:4*Nz,:]
baseline = np.loadtxt("dsa_results/dsa_forward_mom_4_pp_baseline.dat")
baseline_t3 = baseline[:,1]
split = np.loadtxt("dsa_results/dsa_forward_mom_4_pp_split.dat")
split_t3 = split[:,1]
imps = np.loadtxt("dsa_results/dsa_forward_mom_4_pp_imps.dat")
imps_t3 = imps[:,1]

# Plot data
fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='rectilinear')

ax.loglog(E, analytic_t3[Nz//2,:], linewidth=3, label="analytic")
ax.loglog(E, baseline_t3, linewidth=0, marker="o", markersize=5, label="base method")
ax.loglog(E, split_t3, linewidth=0, marker="s", markersize=5, label="splitting")
ax.loglog(E, imps_t3, linewidth=0, marker="^", markersize=5, label="imp. samp.")
ax.set_xlabel('$E$ (MeV)', fontsize=20)
ax.set_ylabel('$f$', fontsize=20)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.set_xlim(1.0, 100.0)
ax.legend(fontsize=20)

plt.savefig("dsa_results/dsa_analytic_numerical_match.png")
plt.show()
plt.close(fig)