# Import libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import argparse
import sys

# Parse arguments
parser = argparse.ArgumentParser(description="Plot trajectory heatmaps at 4 times from a chosen method")
parser.add_argument("--method",
                    type=str,
                    default='analytic',
                    choices=['analytic', 'base', 'split_0', 'split_10', 'split_100', 'split_1000', 'imps_0', 'imps_1'],
                    help="This argument specifies the simulation method, like analytic or base.")
args = parser.parse_args()

# Import data
maps = np.loadtxt("dsa_results/dsa_forward_path_dens_pp_" + args.method + ".dat")
times = np.loadtxt("dsa_results/dsa_analytic_time.dat")
x = np.loadtxt("dsa_results/dsa_analytic_pos.dat")
E = np.loadtxt("dsa_results/dsa_analytic_enr.dat")
Nz = np.size(x)
Np = np.size(E)
map0 = maps[0:Nz,:]
map1 = maps[Nz:2*Nz,:]
map2 = maps[2*Nz:3*Nz,:]
map3 = maps[3*Nz:4*Nz,:]
mmax = 6.0e-5
mmin = mmax / 1.0e3

# Plot data
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
ax1.annotate("$t=${:.0f}".format(times[0]), (-4.6, 2.0), fontsize=20)

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
ax2.annotate("$t=${:.0f}".format(times[1]), (-4.6, 2.0), fontsize=20)

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
ax3.annotate("$t=${:.0f}".format(times[2]), (-4.6, 2.0), fontsize=20)

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
ax4.annotate("$t=${:.0f}".format(times[3]), (-4.6, 2.0), fontsize=20)

plt.savefig("dsa_results/dsa_path_density_" + args.method + ".png")
plt.close(fig)