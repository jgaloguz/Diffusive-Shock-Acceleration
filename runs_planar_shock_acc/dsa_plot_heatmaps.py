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
                    choices=['analytic', 'base', 'split_0', 'split_10', 'split_100', 'split_1000', 'imps_0', 'imps_1', 'all'],
                    help="This argument specifies the simulation method, like analytic or base.")
args = parser.parse_args()

# Import data
if args.method == 'all':
   maps_base = np.loadtxt("dsa_results/dsa_forward_path_dens_pp_base.dat")
   maps_spl0 = np.loadtxt("dsa_results/dsa_forward_path_dens_pp_split_0.dat")
   maps_spl1 = np.loadtxt("dsa_results/dsa_forward_path_dens_pp_split_10.dat")
   maps_imps = np.loadtxt("dsa_results/dsa_forward_path_dens_pp_imps_1.dat")
   times = np.loadtxt("dsa_results/dsa_analytic_time.dat")
   x = np.loadtxt("dsa_results/dsa_analytic_pos.dat")
   E = np.loadtxt("dsa_results/dsa_analytic_enr.dat")
   Nz = np.size(x)
   Np = np.size(E)
   conv = 1.0 / 2.99792414e4 # conversion to physical units
   map_base_0 = conv * maps_base[0:Nz,:]
   map_base_1 = conv * maps_base[Nz:2*Nz,:]
   map_base_2 = conv * maps_base[2*Nz:3*Nz,:]
   map_base_3 = conv * maps_base[3*Nz:4*Nz,:]
   map_spl0_0 = conv * maps_spl0[0:Nz,:]
   map_spl0_1 = conv * maps_spl0[Nz:2*Nz,:]
   map_spl0_2 = conv * maps_spl0[2*Nz:3*Nz,:]
   map_spl0_3 = conv * maps_spl0[3*Nz:4*Nz,:]
   map_spl1_0 = conv * maps_spl1[0:Nz,:]
   map_spl1_1 = conv * maps_spl1[Nz:2*Nz,:]
   map_spl1_2 = conv * maps_spl1[2*Nz:3*Nz,:]
   map_spl1_3 = conv * maps_spl1[3*Nz:4*Nz,:]
   map_imps_0 = conv * maps_imps[0:Nz,:]
   map_imps_1 = conv * maps_imps[Nz:2*Nz,:]
   map_imps_2 = conv * maps_imps[2*Nz:3*Nz,:]
   map_imps_3 = conv * maps_imps[3*Nz:4*Nz,:]
   mmax = 1.0e-13
   mmin = mmax / 1.0e3

# Plot data
   fig = plt.figure(figsize=(15, 12), layout="constrained")
   pcms = []
   axs = []

   ax1 = fig.add_subplot(4, 4, 1, projection='rectilinear')

   pcm = ax1.pcolormesh(x, E, conv * np.transpose(map_base_0), shading='auto',
                         norm=LogNorm(vmin=mmin, vmax=mmax))
   pcms.append(pcm)
   axs.append(ax1)
   ax1.set_ylabel('$\\mathbf{direct}$\n$\\mathbf{integration}$\n$E$ (MeV)', fontsize=20)
   ax1.set_yscale('log')
   ax1.set_xticks([])
   ax1.tick_params(axis='y', labelsize=16)
   ax1.set_title("$t=${:.0f} day".format(times[0]), fontsize=20)

   ax2 = fig.add_subplot(4, 4, 2, projection='rectilinear')

   pcm = ax2.pcolormesh(x, E, conv * np.transpose(map_base_1), shading='auto',
                         norm=LogNorm(vmin=mmin, vmax=mmax))
   pcms.append(pcm)
   axs.append(ax2)
   ax2.set_yscale('log')
   ax2.set_xticks([])
   ax2.tick_params(axis='y', which='both', left=False, labelleft=False)
   ax2.set_title("$t=${:.0f} days".format(times[1]), fontsize=20)

   ax3 = fig.add_subplot(4, 4, 3, projection='rectilinear')

   pcm = ax3.pcolormesh(x, E, conv * np.transpose(map_base_2), shading='auto',
                         norm=LogNorm(vmin=mmin, vmax=mmax))
   pcms.append(pcm)
   axs.append(ax3)
   ax3.set_yscale('log')
   ax3.set_xticks([])
   ax3.tick_params(axis='y', which='both', left=False, labelleft=False)
   ax3.set_title("$t=${:.0f} days".format(times[2]), fontsize=20)

   ax4 = fig.add_subplot(4, 4, 4, projection='rectilinear')

   pcm = ax4.pcolormesh(x, E, conv * np.transpose(map_base_3), shading='auto',
                         norm=LogNorm(vmin=mmin, vmax=mmax))
   pcms.append(pcm)
   axs.append(ax4)
   ax4.set_yscale('log')
   ax4.set_xticks([])
   ax4.tick_params(axis='y', which='both', left=False, labelleft=False)
   ax4.set_title("$t=${:.0f} days".format(times[3]), fontsize=20)

   ax5 = fig.add_subplot(4, 4, 5, projection='rectilinear')

   pcm = ax5.pcolormesh(x, E, conv * np.transpose(map_spl0_0), shading='auto',
                         norm=LogNorm(vmin=mmin, vmax=mmax))
   pcms.append(pcm)
   axs.append(ax5)
   ax5.set_ylabel('$\\mathbf{continuous}$\n$\\mathbf{splitting}$\n$E$ (MeV)', fontsize=20)
   ax5.set_yscale('log')
   ax5.set_xticks([])
   ax5.tick_params(axis='y', labelsize=16)

   ax6 = fig.add_subplot(4, 4, 6, projection='rectilinear')

   pcm = ax6.pcolormesh(x, E, conv * np.transpose(map_spl0_1), shading='auto',
                         norm=LogNorm(vmin=mmin, vmax=mmax))
   pcms.append(pcm)
   axs.append(ax6)
   ax6.set_yscale('log')
   ax6.set_xticks([])
   ax6.tick_params(axis='y', which='both', left=False, labelleft=False)

   ax7 = fig.add_subplot(4, 4, 7, projection='rectilinear')

   pcm = ax7.pcolormesh(x, E, conv * np.transpose(map_spl0_2), shading='auto',
                         norm=LogNorm(vmin=mmin, vmax=mmax))
   pcms.append(pcm)
   axs.append(ax7)
   ax7.set_yscale('log')
   ax7.set_xticks([])
   ax7.tick_params(axis='y', which='both', left=False, labelleft=False)

   ax8 = fig.add_subplot(4, 4, 8, projection='rectilinear')

   pcm = ax8.pcolormesh(x, E, conv * np.transpose(map_spl0_3), shading='auto',
                         norm=LogNorm(vmin=mmin, vmax=mmax))
   pcms.append(pcm)
   axs.append(ax8)
   ax8.set_yscale('log')
   ax8.set_xticks([])
   ax8.tick_params(axis='y', which='both', left=False, labelleft=False)

   ax9 = fig.add_subplot(4, 4, 9, projection='rectilinear')

   pcm = ax9.pcolormesh(x, E, conv * np.transpose(map_spl1_0), shading='auto',
                         norm=LogNorm(vmin=mmin, vmax=mmax))
   pcms.append(pcm)
   axs.append(ax9)
   ax9.set_ylabel('$\\mathbf{discrete}$\n$\\mathbf{splitting}$\n$E$ (MeV)', fontsize=20)
   ax9.set_yscale('log')
   ax9.set_xticks([])
   ax9.tick_params(axis='y', labelsize=16)

   ax10 = fig.add_subplot(4, 4, 10, projection='rectilinear')

   pcm = ax10.pcolormesh(x, E, conv * np.transpose(map_spl1_1), shading='auto',
                         norm=LogNorm(vmin=mmin, vmax=mmax))
   pcms.append(pcm)
   axs.append(ax10)
   ax10.set_yscale('log')
   ax10.set_xticks([])
   ax10.tick_params(axis='y', which='both', left=False, labelleft=False)

   ax11 = fig.add_subplot(4, 4, 11, projection='rectilinear')

   pcm = ax11.pcolormesh(x, E, conv * np.transpose(map_spl1_2), shading='auto',
                         norm=LogNorm(vmin=mmin, vmax=mmax))
   pcms.append(pcm)
   axs.append(ax11)
   ax11.set_yscale('log')
   ax11.set_xticks([])
   ax11.tick_params(axis='y', which='both', left=False, labelleft=False)

   ax12 = fig.add_subplot(4, 4, 12, projection='rectilinear')

   pcm = ax12.pcolormesh(x, E, conv * np.transpose(map_spl1_3), shading='auto',
                         norm=LogNorm(vmin=mmin, vmax=mmax))
   pcms.append(pcm)
   axs.append(ax12)
   ax12.set_yscale('log')
   ax12.set_xticks([])
   ax12.tick_params(axis='y', which='both', left=False, labelleft=False)

   ax13 = fig.add_subplot(4, 4, 13, projection='rectilinear')

   pcm = ax13.pcolormesh(x, E, conv * np.transpose(map_imps_0), shading='auto',
                         norm=LogNorm(vmin=mmin, vmax=mmax))
   pcms.append(pcm)
   axs.append(ax13)
   ax13.set_xticks([-5,-2.5,0,2.5,5])
   ax13.set_xlabel('$x$ (au)', fontsize=20)
   ax13.set_ylabel('$\\mathbf{importance}$\n$\\mathbf{sampling}$\n$E$ (MeV)', fontsize=20)
   ax13.set_yscale('log')
   ax13.tick_params(axis='x', labelsize=16)
   ax13.tick_params(axis='y', labelsize=16)

   ax14 = fig.add_subplot(4, 4, 14, projection='rectilinear')

   pcm = ax14.pcolormesh(x, E, conv * np.transpose(map_imps_1), shading='auto',
                         norm=LogNorm(vmin=mmin, vmax=mmax))
   pcms.append(pcm)
   axs.append(ax14)
   ax14.set_xticks([-5,-2.5,0,2.5,5])
   ax14.set_xlabel('$x$ (au)', fontsize=20)
   ax14.set_yscale('log')
   ax14.tick_params(axis='x', labelsize=16)
   ax14.tick_params(axis='y', which='both', left=False, labelleft=False)

   ax15 = fig.add_subplot(4, 4, 15, projection='rectilinear')

   pcm = ax15.pcolormesh(x, E, conv * np.transpose(map_imps_2), shading='auto',
                         norm=LogNorm(vmin=mmin, vmax=mmax))
   pcms.append(pcm)
   axs.append(ax15)
   ax15.set_xticks([-5,-2.5,0,2.5,5])
   ax15.set_xlabel('$x$ (au)', fontsize=20)
   ax15.set_yscale('log')
   ax15.tick_params(axis='x', labelsize=16)
   ax15.tick_params(axis='y', which='both', left=False, labelleft=False)

   ax16 = fig.add_subplot(4, 4, 16, projection='rectilinear')

   pcm = ax16.pcolormesh(x, E, conv * np.transpose(map_imps_3), shading='auto',
                         norm=LogNorm(vmin=mmin, vmax=mmax))
   pcms.append(pcm)
   axs.append(ax16)
   ax16.set_xticks([-5,-2.5,0,2.5,5])
   ax16.set_xlabel('$x$ (au)', fontsize=20)
   ax16.set_yscale('log')
   ax16.tick_params(axis='x', labelsize=16)
   ax16.tick_params(axis='y', which='both', left=False, labelleft=False)

   cb = fig.colorbar(pcms[0], ax=axs, orientation='vertical')
   cb.set_label('Unweighted path intensity (particles cm$^{-3}$ (MeV/c)$^{-1}$)', fontsize=20)
   cb.ax.tick_params(labelsize=16)

   plt.savefig("dsa_results/dsa_path_density_" + args.method + ".png")
   plt.close(fig)

else:
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

   pcm1 = ax1.pcolormesh(x, E, conv * np.transpose(map0), shading='auto',
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

   pcm2 = ax2.pcolormesh(x, E, conv * np.transpose(map1), shading='auto',
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

   pcm3 = ax3.pcolormesh(x, E, conv * np.transpose(map2), shading='auto',
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

   pcm4 = ax4.pcolormesh(x, E, conv * np.transpose(map3), shading='auto',
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