# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import sys

# What to plot
# 0: spectral plots comparing baseline, splitting, and importance sampling
# 1: path density at several points in time for baseline and splitting
which_plot = 1

if which_plot == 0:
# Import data
   dsa_simulated_mom = [np.loadtxt("dsa_results/dsa_forward_mom_4_pp_baseline.dat"),
                        np.loadtxt("dsa_results/dsa_forward_mom_4_pp_split.dat"),
                        np.loadtxt("dsa_results/dsa_forward_mom_4_pp_imps.dat")]
   labels_sim = ["BS", "SP", "IS"]
   colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]

# Plot
   fig = plt.figure(figsize=(15, 10), layout='tight')
   ax = fig.add_subplot(111, projection='rectilinear')

   for t in range(3):
      ax.loglog(dsa_simulated_mom[t][:,0], dsa_simulated_mom[t][:,1],
                color=colors[t], label=labels_sim[t], linestyle="", marker="s")
   ax.set_xlabel('$E$ (MeV)', fontsize=20)
   ax.set_ylabel('$4\\pi fp^2$', fontsize=20)
   ax.tick_params(axis='x', labelsize=20)
   ax.tick_params(axis='y', labelsize=20)
   ax.legend(fontsize=20)

   plt.savefig("dsa_mom_forward_4_5_acc.png")
   plt.show()
   plt.close(fig)

elif which_plot == 1:
# Import data
   maps = np.loadtxt("dsa_results/dsa_forward_path_dens_pp_split.dat")
   Nz = np.size(maps,0) // 4
   Np = np.size(maps,1)
   map0 = maps[0:Nz,:]
   map1 = maps[Nz:2*Nz,:]
   map2 = maps[2*Nz:3*Nz,:]
   map3 = maps[3*Nz:4*Nz,:]

   x = np.linspace(-5.0, 5.0, num=Nz+1)
   x = 0.5 * (x[:-1] + x[1:])
   p = np.linspace(np.log10(1.0), np.log10(100.0), num=Np+1)
   p = np.power(10.0, 0.5 * (p[:-1] + p[1:]))

# Plot
   fig = plt.figure(figsize=(15, 10), layout='tight')
   ax1 = fig.add_subplot(221, projection='rectilinear')

   ax1.pcolormesh(x, p, np.log(np.transpose(map0)), shading='gouraud')
   ax1.set_xlabel('$x$ (au)', fontsize=20)
   ax1.set_ylabel('$E$ (MeV)', fontsize=20)
   ax1.tick_params(axis='x', labelsize=20)
   ax1.tick_params(axis='y', labelsize=20)
   ax1.annotate("$t=t_0$", (3.7, 90), fontsize=20)

   ax2 = fig.add_subplot(222, projection='rectilinear')

   ax2.pcolormesh(x, p, np.log(np.transpose(map1)), shading='gouraud')
   ax2.set_xlabel('$x$ (au)', fontsize=20)
   ax2.set_ylabel('$E$ (MeV)', fontsize=20)
   ax2.tick_params(axis='x', labelsize=20)
   ax2.tick_params(axis='y', labelsize=20)
   ax2.annotate("$t=t_1$", (3.7, 90), fontsize=20)

   ax3 = fig.add_subplot(223, projection='rectilinear')

   ax3.pcolormesh(x, p, np.log(np.transpose(map2)), shading='gouraud')
   ax3.set_xlabel('$x$ (au)', fontsize=20)
   ax3.set_ylabel('$E$ (MeV)', fontsize=20)
   ax3.tick_params(axis='x', labelsize=20)
   ax3.tick_params(axis='y', labelsize=20)
   ax3.annotate("$t=t_2$", (3.7, 90), fontsize=20)

   ax4 = fig.add_subplot(224, projection='rectilinear')

   ax4.pcolormesh(x, p, np.log(np.transpose(map3)), shading='gouraud')
   ax4.set_xlabel('$x$ (au)', fontsize=20)
   ax4.set_ylabel('$E$ (MeV)', fontsize=20)
   ax4.tick_params(axis='x', labelsize=20)
   ax4.tick_params(axis='y', labelsize=20)
   ax4.annotate("$t=t_3$", (3.7, 90), fontsize=20)

   plt.savefig("dsa_path_density_forward_split.png")
   plt.show()
   plt.close(fig)

else:
   print("Nothing to do...")