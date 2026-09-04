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

filename = [
            "dsa_results/spectrum_comparison/dsa_forward_mom_pp_base.dat",
            "dsa_results/spectrum_comparison/dsa_forward_mom_pp_split_0.dat",
            "dsa_results/spectrum_comparison/dsa_forward_mom_pp_split_100.dat",
            "dsa_results/spectrum_comparison/dsa_forward_mom_pp_imps_0.dat",
            "dsa_results/spectrum_comparison/dsa_forward_mom_pp_imps_1.dat",
           ]

label = [
         "base method",
         "continuous splitting",
         "discrete splitting",
         "importance sampling splitting",
         "importance sampling Prinsloo",
        ]

exec_times = np.loadtxt("dsa_results/spectrum_comparison/time_stats.dat")

marker = ['o', 'x', '^', 's', '+', 'D']

analytic = np.loadtxt("dsa_results/dsa_forward_path_dens_pp_analytic.dat")
loc = 3*Nz + np.digitize([0.15], x)
analytic_t3 = analytic[loc[0],:]
numerical_t3 = []

n_cases = len(filename)
for run in range(n_cases):
   numerical = np.loadtxt(filename[run])
   numerical_t3.append(numerical[:,1])

# Plot data
fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='rectilinear')

for run in range(n_cases):
   ax.loglog(E, numerical_t3[run], linewidth=0, marker=marker[run], markersize=7, label=label[run]+" ({:.2f} $\\pm$ {:.2f})".format(exec_times[run,0], exec_times[run,1]))
ax.loglog(E, analytic_t3, linewidth=5, label="analytic")

ax.set_xlabel('$E$ (MeV)', fontsize=24)
ax.set_ylabel('$F$ (particles cm$^{-3}$ MeV$^{-1}$)', fontsize=24)
ax.tick_params(axis='x', labelsize=24)
ax.tick_params(axis='y', labelsize=24)
ax.legend(loc=3, fontsize=20)

plt.savefig("dsa_results/dsa_spectrum_comparison.png")
plt.close(fig)
