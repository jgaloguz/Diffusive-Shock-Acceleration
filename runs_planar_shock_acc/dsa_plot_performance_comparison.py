# Import libraries
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import argparse
import sys

# Parse arguments
parser = argparse.ArgumentParser(description="Plot spectrum and e^2 metric (times unit cost) for the base method and a variety of enhanced sampling cases.")
parser.add_argument("case",
                    type=str,
                    choices=['split_0', 'split_1', 'imps_0', 'imps_1', 'all'],
                    help="This argument specifies the type of enhanced sampling method to use when plotting results.")
parser.add_argument("--folder",
                    type=str,
                    default="e2metric",
                    help="This argument specifies the subfolder inside dsa_results where results are stored.")
args = parser.parse_args()

# Import data
x = np.loadtxt("dsa_results/dsa_analytic_pos.dat")
E = np.loadtxt("dsa_results/dsa_analytic_enr.dat")
E_l = np.loadtxt("dsa_results/dsa_analytic_enr_left.dat")
Nz = np.size(x)
Np = np.size(E)

analytic = np.loadtxt("dsa_results/dsa_forward_path_dens_pp_analytic.dat")
loc = 3*Nz + np.digitize([0.15], x)
analytic_t3 = analytic[loc[0],:]

marker = ['o', 'x', '^', 's', '+', 'D']

if args.case == 'all':
   case = ['split_0', 'split_1', 'imps_0', 'imps_1']
   row_title = ['$\\mathbf{continuous\\ splitting}$', '$\\mathbf{discrete\\ splitting}$', '$\\mathbf{splitting-matched}$\n$\\mathbf{importance\\ sampler}$', '$\\mathbf{shock-directed}$\n$\\mathbf{importance\\ sampler}$']
   folder = ['e2metric_split_0', 'e2metric_split_1', 'e2metric_imps_0', 'e2metric_imps_1']
else:
   case = [args.case]
   folder = [args.folder]


for j in range(len(case)):
   cost_per_particle = np.loadtxt("dsa_results/" + folder[j] + "/steps_per_particle.dat")
   n_cases = np.size(cost_per_particle, 0)
   if case[j] == "split_0" or case[j] == "split_1" or case[j] == "imps_0":
      label = ["$\\alpha$ = {:.2f}".format(cost_per_particle[i,0]) for i in range(n_cases)]
      filename = ["dsa_results/" + folder[j] + "/dsa_forward_mom_pp_{:s}_{:.2f}.dat".format(case[j], cost_per_particle[i,0]) for i in range(n_cases)]
   elif case[j] == "imps_1":
      label = ["$A_0$ = {:.1f}".format(cost_per_particle[i,0]) for i in range(n_cases)]
      filename = ["dsa_results/" + folder[j] + "/dsa_forward_mom_pp_{:s}_{:.1f}.dat".format(case[j], cost_per_particle[i,0]) for i in range(n_cases)]
   label[0] = "direct integration"
   filename[0] = "dsa_results/" + folder[j] + "/dsa_forward_mom_pp_base.dat"

   Cpp = cost_per_particle[:,1]

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
   fig = plt.figure(figsize=(18, 10), layout='tight')

   ax = fig.add_subplot(111, projection='rectilinear')

   for run in range(n_cases):
      ax.loglog(E, numerical_t3[run], linewidth=0, marker=marker[run], markersize=12, label=label[run])
   ax.loglog(E, analytic_t3, linewidth=5, label="analytic")

   if args.case == 'all' and j < len(case)-1:
      ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
   else:
      ax.set_xlabel('$E$ (MeV)', fontsize=36)
      ax.tick_params(axis='x', labelsize=36)
   ax.set_ylabel('$F$ (particles s MeV$^{-1}$)', fontsize=36)
   ax.tick_params(axis='y', labelsize=36)
   ax.legend(fontsize=36)

   plt.savefig("dsa_results/dsa_spectrum_comparison_{:s}.png".format(case[j]))
   plt.close(fig)

   fig = plt.figure(figsize=(18, 10), layout='tight')

   ax = fig.add_subplot(111, projection='rectilinear')

   for run in range(n_cases):
      ax.loglog(E, rmsrelerr_t3[run], linewidth=0, marker=marker[run], markersize=12, label=label[run])

   if args.case == 'all' and j < len(case)-1:
      ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
   else:
      ax.set_xlabel('$E$ (MeV)', fontsize=36)
      ax.tick_params(axis='x', labelsize=36)
   ax.set_ylabel('RMS rel. err.', fontsize=36)
   ax.tick_params(axis='y', labelsize=36)
   ax.legend(fontsize=36)

   plt.savefig("dsa_results/dsa_rms_rel_err_comparison_{:s}.png".format(case[j]))
   plt.close(fig)

   fig = plt.figure(figsize=(18, 10), layout='tight')

   ax = fig.add_subplot(111, projection='rectilinear')

   for run in range(n_cases):
      ax.loglog(E_l, H_numerical_t3[run], linewidth=0, marker=marker[run], markersize=12, label=label[run])
   ax.loglog(E_l, H_analytic_t3[0], linewidth=3, linestyle="-", label="analytic")

   if args.case == 'all' and j < len(case)-1:
      ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
   else:
      ax.set_xlabel('$E$ (MeV)', fontsize=36)
      ax.tick_params(axis='x', labelsize=36)
   ax.set_ylabel('$H$', fontsize=36)
   ax.tick_params(axis='y', labelsize=36)
   ax.legend(fontsize=36)

   plt.savefig("dsa_results/dsa_H_analytic_numerical_comparison_{:s}.png".format(case[j]))
   plt.close(fig)

   fig = plt.figure(figsize=(18, 10), layout='tight')

   ax = fig.add_subplot(111, projection='rectilinear')

   for run in range(n_cases):
      ax.loglog(E_l, e2_metric_t3[run], linewidth=0, marker=marker[run], markersize=12, label=label[run])

   if args.case == 'all':
      e2_string = '$\\hat{e}_{p,m}^2$'
      ax.set_ylabel('{:s}\n{:s}'.format(row_title[j],e2_string), fontsize=36)
      if j < len(case)-1:
         ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
      else:
         ax.set_xlabel('$E$ (MeV)', fontsize=36)
         ax.tick_params(axis='x', labelsize=36)
   else:
      ax.set_ylabel('$\\hat{e}_{p,m}^2$', fontsize=36)
      ax.set_xlabel('$E$ (MeV)', fontsize=36)
      ax.tick_params(axis='x', labelsize=36)
   ax.tick_params(axis='y', labelsize=36)
   ax.set_ylim(1.0,1.0e8)
   ax.legend(fontsize=36)

   plt.savefig("dsa_results/dsa_e2_metric_comparison_{:s}.png".format(case[j]))
   plt.close(fig)

   fig = plt.figure(figsize=(18, 10), layout='tight')

   ax = fig.add_subplot(111, projection='rectilinear')

   for run in range(n_cases):
      ax.loglog(E_l, Cpp[run] * e2_metric_t3[run], linewidth=0, marker=marker[run], markersize=12, label=label[run])

   if args.case == 'all' and j < len(case)-1:
      ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
   else:
      ax.set_xlabel('$E$ (MeV)', fontsize=36)
      ax.tick_params(axis='x', labelsize=36)
   ax.set_ylabel('$\\hat{c}_{c,m}$ (integration steps per seed)', fontsize=36)
   ax.tick_params(axis='y', labelsize=36)
   ax.set_ylim(1.0e5,1.0e13)
   if args.case != 'all':
      ax.legend(fontsize=36)

   plt.savefig("dsa_results/dsa_e2_x_cost_metric_comparison_{:s}.png".format(case[j]))
   plt.close(fig)

   # Compute maximum acceleration
   print("Maximum Acceleration:")
   for run in range(1,n_cases):
      interp = interp1d(E_l, Cpp[0] * e2_metric_t3[0] / (Cpp[run] * e2_metric_t3[run]), kind='linear')
      print("\t{:s}    acc (3e5 MeV) = {:.3f}x".format(label[run], interp(3.0e5)))
