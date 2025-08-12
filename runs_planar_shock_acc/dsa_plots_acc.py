# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import sys

# Import data
which_time_flow = sys.argv[1]
dsa_simulated_mom = [np.loadtxt("dsa_results/dsa_" + which_time_flow + "_mom_4_pp_baseline.dat"),
                     np.loadtxt("dsa_results/dsa_" + which_time_flow + "_mom_4_pp_split.dat"),
                     np.loadtxt("dsa_results/dsa_" + which_time_flow + "_mom_4_pp_imps.dat")]
labels_sim = ["BS", "SP", "IS"]
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]

# Plot
fig = plt.figure(figsize=(15, 10), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')

for t in range(3):
   ax.loglog(dsa_simulated_mom[t][:,0], dsa_simulated_mom[t][:,1],
             color=colors[t], label=labels_sim[t], linestyle="", marker="s")
ax2.set_xlabel('$E$ (MeV)', fontsize=20)
ax2.set_ylabel('$4\\pi fp^2$', fontsize=20)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.legend(fontsize=20)

plt.savefig("dsa_mom_" + which_time_flow + "_4_5_acc.png")
plt.show()
plt.close(fig)