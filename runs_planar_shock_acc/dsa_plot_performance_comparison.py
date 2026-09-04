# Import libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import argparse
import sys

# Parse arguments
parser = argparse.ArgumentParser(description="Plot spectrum and e^2 metric (times unit cost) for the base method and a variety of enhanced sampling cases.")
parser.add_argument("case",
                    type=str,
                    choices=['split_0', 'split_1', 'imps_0', 'imps_1'],
                    help="This argument specifies the type of enhanced sampling method to use when plotting results.")
parser.add_argument("folder",
                    type=str,
                    help="This argument specifies the subfolder inside dsa_results where results are stored.")
args = parser.parse_args()

cost_per_particle = np.loadtxt("dsa_results/" + args.folder + "/steps_per_particle.dat")
n_cases = np.size(cost_per_particle, 0)
if args.case == "split_0" or args.case == "split_1" or args.case == "imps_0":
   label = ["$\\alpha$ = {:.2f}".format(cost_per_particle[i,0]) for i in range(n_cases)]
   filename = ["dsa_results/" + args.folder + "/dsa_forward_mom_pp_{:s}_{:.2f}.dat".format(args.case, cost_per_particle[i,0]) for i in range(n_cases)]
elif args.case == "imps_1":
   label = ["$A_0$ = {:.1f}".format(cost_per_particle[i,0]) for i in range(n_cases)]
   filename = ["dsa_results/" + args.folder + "/dsa_forward_mom_pp_{:s}_{:.1f}.dat".format(args.case, cost_per_particle[i,0]) for i in range(n_cases)]
label[0] = "base method"
filename[0] = "dsa_results/" + args.folder + "/dsa_forward_mom_pp_base.dat"

# Import data
x = np.loadtxt("dsa_results/dsa_analytic_pos.dat")
E = np.loadtxt("dsa_results/dsa_analytic_enr.dat")
Nz = np.size(x)
Np = np.size(E)

marker = ['o', 'x', '^', 's', '+', 'D']
Cpp = cost_per_particle[:,1] / cost_per_particle[0,1]

analytic = np.loadtxt("dsa_results/dsa_forward_path_dens_pp_analytic.dat")
loc = 3*Nz + np.digitize([0.15], x)
analytic_t3 = analytic[loc[0],:]
numerical_t3 = []
rmsrelerr_t3 = []
H_analytic_t3 = []
H_numerical_t3 = []
e2_metric_t3 = []
for run in range(n_cases):
   numerical = np.loadtxt(filename[run])
   numerical_t3.append(numerical[:,1])
   rmsrelerr_t3.append(numerical[:,3])
   H_analytic_t3.append(numerical[:,4])
   H_numerical_t3.append(numerical[:,5])
   e2_metric_t3.append(numerical[:,6])

# Plot data
fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='rectilinear')

for run in range(n_cases):
   ax.loglog(E, numerical_t3[run], linewidth=0, marker=marker[run], markersize=10, label=label[run])
ax.loglog(E, analytic_t3, linewidth=5, label="analytic")

ax.set_xlabel('$E$ (MeV)', fontsize=28)
ax.set_ylabel('$F$ (particles s MeV$^{-1}$)', fontsize=28)
ax.tick_params(axis='x', labelsize=28)
ax.tick_params(axis='y', labelsize=28)
ax.legend(fontsize=28)

plt.savefig("dsa_results/dsa_spectrum_comparison_{:s}.png".format(args.case))
plt.close(fig)

fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='rectilinear')

for run in range(n_cases):
   ax.loglog(E, rmsrelerr_t3[run], linewidth=0, marker=marker[run], markersize=10, label=label[run])

ax.set_xlabel('$E$ (MeV)', fontsize=28)
ax.set_ylabel('RMS rel. err.', fontsize=28)
ax.tick_params(axis='x', labelsize=28)
ax.tick_params(axis='y', labelsize=28)
ax.legend(fontsize=28)

plt.savefig("dsa_results/dsa_rms_rel_err_comparison_{:s}.png".format(args.case))
plt.close(fig)

fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='rectilinear')

for run in range(n_cases):
   ax.loglog(E, H_numerical_t3[run], linewidth=0, marker=marker[run], markersize=10, label=label[run])
ax.loglog(E, H_analytic_t3[0], linewidth=3, linestyle="-", label="analytic")

ax.set_xlabel('$E$ (MeV)', fontsize=28)
ax.set_ylabel('$H$', fontsize=28)
ax.tick_params(axis='x', labelsize=28)
ax.tick_params(axis='y', labelsize=28)
ax.legend(fontsize=28)

plt.savefig("dsa_results/dsa_H_analytic_numerical_comparison_{:s}.png".format(args.case))
plt.close(fig)

fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='rectilinear')

for run in range(n_cases):
   ax.loglog(E, e2_metric_t3[run], linewidth=0, marker=marker[run], markersize=10, label=label[run])
ax.set_xlabel('$E$ (MeV)', fontsize=28)
ax.set_ylabel('$\\hat{e}^2$', fontsize=28)
ax.tick_params(axis='x', labelsize=28)
ax.tick_params(axis='y', labelsize=28)
ax.set_ylim(bottom=1.0)
ax.legend(fontsize=28)

plt.savefig("dsa_results/dsa_e2_metric_comparison_{:s}.png".format(args.case))
plt.close(fig)

fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='rectilinear')

for run in range(n_cases):
   ax.loglog(E, Cpp[run] * e2_metric_t3[run], linewidth=0, marker=marker[run], markersize=10, label=label[run])
ax.set_xlabel('$E$ (MeV)', fontsize=28)
ax.set_ylabel('$\\hat{c}$', fontsize=28)
ax.tick_params(axis='x', labelsize=28)
ax.tick_params(axis='y', labelsize=28)
ax.set_ylim(bottom=1.0)
ax.legend(fontsize=28)

plt.savefig("dsa_results/dsa_e2_x_cost_metric_comparison_{:s}.png".format(args.case))
plt.close(fig)

# Compute maximum acceleration
print("Maximum Acceleration:")
for run in range(1,n_cases):
   max_val = np.max(e2_metric_t3[0] / (Cpp[run] * e2_metric_t3[run]))
   max_idx = np.argmax(e2_metric_t3[0] / (Cpp[run] * e2_metric_t3[run]))
   print("\t{:s}    E = {:.2f}    acc = {:.3f}x".format(label[run], E[max_idx], max_val))
