# Import libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys

if len(sys.argv) > 1:
   which_acc = sys.argv[1]   
else:
   print("ERROR: No arguments were provided.")
   print("Valid arguments: analytic, baseline, split, imps")
   exit(1)

# Import data
maps = np.loadtxt("dsa_results/dsa_forward_path_dens_pp_" + which_acc + ".dat")
times = np.loadtxt("dsa_results/dsa_analytic_time.dat")
x = np.loadtxt("dsa_results/dsa_analytic_pos.dat")
E = np.loadtxt("dsa_results/dsa_analytic_enr.dat")
Nz = np.size(x)
Np = np.size(E)
map0 = maps[0:Nz,:]
map1 = maps[Nz:2*Nz,:]
map2 = maps[2*Nz:3*Nz,:]
map3 = maps[3*Nz:4*Nz,:]
mmax = np.max(map3)
mmin = mmax / 1.0e4

# Plot
fig = plt.figure(figsize=(15, 10), layout='tight')
ax1 = fig.add_subplot(221, projection='rectilinear')

pcm1 = ax1.pcolormesh(x, E, np.transpose(map0), shading='gouraud',
                      norm=LogNorm(vmin=mmin, vmax=mmax))
cb1 = fig.colorbar(pcm1, ax=ax1)
cb1.ax.tick_params(labelsize=20)
ax1.set_xlabel('$x$ (au)', fontsize=20)
ax1.set_ylabel('$E$ (MeV)', fontsize=20)
ax1.set_yscale('log')
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)
ax1.annotate("$t=${:.1f}".format(times[0]), (2.8, 60.0), fontsize=20)

ax2 = fig.add_subplot(222, projection='rectilinear')

pcm2 = ax2.pcolormesh(x, E, np.transpose(map1), shading='gouraud',
                      norm=LogNorm(vmin=mmin, vmax=mmax))
cb2 = fig.colorbar(pcm2, ax=ax2)
cb2.ax.tick_params(labelsize=20)
ax2.set_xlabel('$x$ (au)', fontsize=20)
ax2.set_ylabel('$E$ (MeV)', fontsize=20)
ax2.set_yscale('log')
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
ax2.annotate("$t=${:.1f}".format(times[1]), (2.8, 60.0), fontsize=20)

ax3 = fig.add_subplot(223, projection='rectilinear')

pcm3 = ax3.pcolormesh(x, E, np.transpose(map2), shading='gouraud',
                      norm=LogNorm(vmin=mmin, vmax=mmax))
cb3 = fig.colorbar(pcm3, ax=ax3)
cb3.ax.tick_params(labelsize=20)
ax3.set_xlabel('$x$ (au)', fontsize=20)
ax3.set_ylabel('$E$ (MeV)', fontsize=20)
ax3.set_yscale('log')
ax3.tick_params(axis='x', labelsize=20)
ax3.tick_params(axis='y', labelsize=20)
ax3.annotate("$t=${:.1f}".format(times[2]), (2.8, 60.0), fontsize=20)

ax4 = fig.add_subplot(224, projection='rectilinear')

pcm4 = ax4.pcolormesh(x, E, np.transpose(map3), shading='gouraud',
                      norm=LogNorm(vmin=mmin, vmax=mmax))
cb4 = fig.colorbar(pcm4, ax=ax4)
cb4.ax.tick_params(labelsize=20)
ax4.set_xlabel('$x$ (au)', fontsize=20)
ax4.set_ylabel('$E$ (MeV)', fontsize=20)
ax4.set_yscale('log')
ax4.tick_params(axis='x', labelsize=20)
ax4.tick_params(axis='y', labelsize=20)
ax4.annotate("$t=${:.1f}".format(times[3]), (2.8, 60.0), fontsize=20)

plt.savefig("dsa_path_density_forward_" + which_acc + ".png")
plt.show()
plt.close(fig)

# Test
analytic = np.loadtxt("dsa_results/dsa_forward_path_dens_pp_analytic.dat")
analytic_t3 = analytic[3*Nz:4*Nz,:]
baseline = np.loadtxt("dsa_results/dsa_forward_path_dens_pp_baseline.dat")
baseline_t3 = baseline[3*Nz:4*Nz,:]

fig = plt.figure(figsize=(15, 10), layout='tight')

ax = fig.add_subplot(111, projection='rectilinear')

ax.loglog(E, analytic_t3[Nz//2,:], linewidth=3, label="analytic")
ax.loglog(E, baseline_t3[Nz//2,:], linewidth=3, label="baseline")
ax.set_xlabel('$E$ (MeV)', fontsize=20)
ax.set_ylabel('$f$', fontsize=20)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.set_xlim(1.0, 100.0)
ax.legend(fontsize=20)

plt.savefig("baseline_analytic_match.png")
plt.show()
plt.close(fig)